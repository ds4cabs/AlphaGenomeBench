"""
AlphaGenome Multi-Modal Variant Effect Prediction Pipeline
==========================================================
Author  : Zechuan Shi, Swarup Lab, UC Irvine
Contact : zechuas@uci.edu
Version : 1.1

Description:
    Batch variant effect prediction using the AlphaGenome model.
    Given a CSV of variants, this script runs multi-modal predictions
    (RNA-seq, ATAC, CTCF, H3K27ac) for each variant and saves a
    figure per variant to the specified output directory.

Usage:
    # Optional arguments shown with their default values
    python ./PATH_to_code/alphagenome_variant_pipeline.py \
        --csv variants_input.csv \
        --output_dir ./figures \
        --api_key YOUR_API_KEY \
        [--ontology UBERON:0001134] \
        [--zoom 32768] \
        [--interval_size 1048576] \
        [--aggregate_rna]

CSV Format (required columns):
    gene_name, rsID, chr, position, ref_allele, alt_allele
Optional columns:
    ontology_curie  (per-variant cell-type ontology; overrides --ontology if present)
Any extra columns are ignored safely.

Parameters:
    --interval_size : The genomic window (in bp) centered on the variant
                      that is fed to the model for prediction.
                      The model requires a large context window to make
                      accurate predictions — default is 1,048,576 bp (1 Mb).
                      Do not reduce this unless you have a specific reason.

    --zoom          : The sub-window (in bp) shown in the saved figure.
                      Default is 32,768 bp (32 kb) centered on the variant.
                      This only affects what region is displayed in the plot,
                      not the prediction itself.

    --aggregate_rna : Convert RNA-seq genomic curves to aggregated expression
                      (mean across interval). Plots show flat lines at total
                      expression level; track summary logs REF/ALT expression values.
"""

import os
import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for script usage
import matplotlib.pyplot as plt

from alphagenome.data import gene_annotation, genome, transcript as transcript_utils, track_data
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components


# Keywords that indicate a track is from a modified/non-wildtype condition
MODIFIED_KEYWORDS = ['crispr', 'modified', 'insertion', 'deletion', 'mutant', 'knockout']


# ---------------------------------------------------------------------------
# CLI Argument Parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run AlphaGenome variant effect predictions from a CSV file."
    )
    parser.add_argument(
        '--csv', required=True,
        help='Path to input CSV file with variant definitions.'
    )
    parser.add_argument(
        '--output_dir', required=True,
        help='Directory to save output figures.'
    )
    parser.add_argument(
        '--api_key', default=None,
        help=(
            'AlphaGenome API key. If omitted, the script will read it from the '
            'ALPHAGENOME_API_KEY environment variable.'
        )
    )
    parser.add_argument(
        '--ontology', default='UBERON:0001134',
        help='Ontology term for cell type (default: UBERON:0001134).'
    )
    parser.add_argument(
        '--zoom', type=int, default=32768,
        help=(
            'Sub-window (in bp) displayed in the saved figure (default: 32768 = 32 kb). '
            'Does not affect prediction accuracy, only what region is shown in the plot.'
        )
    )
    parser.add_argument(
        '--interval_size', type=int, default=1048576,
        help=(
            'Genomic window (in bp) centered on the variant, used as model input context '
            '(default: 1048576 = 1 Mb). The model needs a large window for accuracy -- '
            'only change if you have a specific reason.'
        )
    )
    parser.add_argument(
        '--aggregate_rna', action='store_true',
        help=(
            'Convert RNA-seq genomic curves to aggregated expression values (mean across '
            'interval). Plots show flat lines at total expression level instead of per-bp curves.'
        )
    )
    args = parser.parse_args()
    if not args.api_key:
        args.api_key = os.environ.get('ALPHAGENOME_API_KEY')
    if not args.api_key:
        raise SystemExit(
            "Missing API key. Provide --api_key or set environment variable ALPHAGENOME_API_KEY."
        )
    return args


# ---------------------------------------------------------------------------
# GTF / Annotation Setup (loaded once, shared across all variants)
# ---------------------------------------------------------------------------

def load_annotation():
    """Download and prepare the MANE Select transcript extractor."""
    print("Loading GTF annotation (this may take a moment on first run)...")
    gtf_url = (
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    gtf = pd.read_feather(gtf_url)
    gtf_mane = gene_annotation.filter_to_mane_select_transcript(gtf)
    extractor = transcript_utils.TranscriptExtractor(gtf_mane)
    print("Annotation loaded.")
    return extractor


# ---------------------------------------------------------------------------
# Track Helper Functions
# ---------------------------------------------------------------------------

def is_modified(name):
    """Return True if a track name contains any modified/non-wildtype keyword."""
    name_lower = name.lower()
    return any(kw in name_lower for kw in MODIFIED_KEYWORDS)


def select_best_track(tdata):
    """
    Smart track selection logic:
      - 1 track  → use as-is
      - Multiple tracks, all identical names → average them (true replicates)
      - Multiple tracks, different names     → prefer unmodified/wildtype tracks;
                                               if none found, fall back to the first track
    Returns (selected_tdata, log_lines) where log_lines is a list of strings
    describing what was selected and why.
    """
    log_lines = []
    n = tdata.values.shape[1]

    # No tracks at all — just return as-is and let caller decide whether to plot.
    if n == 0:
        log_lines.append("    [INFO] 0 tracks found for this modality — skipping.")
        return tdata, log_lines

    if n == 1:
        return tdata, log_lines

    names = tdata.metadata['name'].tolist() if 'name' in tdata.metadata.columns else []

    # Case: all names are identical → true replicates, average them
    if len(set(names)) == 1 and len(names) == n:
        log_lines.append(f"    → {n} identical tracks found — averaging as true replicates")
        log_lines.append(f"    → Track: {names[0]}")
        mean_values = tdata.values.mean(axis=1, keepdims=True)
        averaged = track_data.TrackData(
            values=mean_values,
            metadata=tdata.metadata.iloc[[0]].copy(),
            interval=tdata.interval,
            resolution=tdata.resolution,
        )
        return averaged, log_lines

    # Case: multiple tracks with different names → pick best unmodified one
    log_lines.append(f"    [WARNING] {n} different tracks found — selecting wildtype/unmodified:")
    for name in names:
        tag = '[MODIFIED - skipped]' if is_modified(name) else '[wildtype]'
        log_lines.append(f"      - {name}  {tag}")

    unmodified_idx = [i for i, name in enumerate(names) if not is_modified(name)]

    if unmodified_idx:
        chosen_idx = unmodified_idx[0]
        chosen_name = names[chosen_idx] if names else '<unknown>'
        log_lines.append(f"    → Selected: {chosen_name}")
    else:
        chosen_idx = 0
        fallback_name = names[0] if names else '<unknown>'
        log_lines.append(f"    [WARNING] No unmodified track found — falling back to: {fallback_name}")

    # Use explicit integer indices for robustness with underlying TrackData implementation.
    selected = tdata.filter_tracks([chosen_idx])
    return selected, log_lines


def filter_by_metadata(tracks, **filters):
    """Filter tracks by arbitrary metadata column key=value pairs."""
    mask = pd.Series([True] * len(tracks.metadata), index=tracks.metadata.index)
    for key, value in filters.items():
        if key in tracks.metadata.columns:
            mask &= (tracks.metadata[key] == value)
    # Convert boolean mask to explicit integer indices to avoid dtype issues.
    indices = tracks.metadata.index[mask].tolist()
    return tracks.filter_tracks(indices)


def filter_by_ontology(tracks, ontology_curie):
    """Filter tracks to a specific cell type ontology."""
    # Some Alphagenome outputs may not include an 'ontology_curie' column in
    # their metadata. In that case, fall back to returning all tracks.
    if 'ontology_curie' not in tracks.metadata.columns:
        return tracks
    mask = tracks.metadata['ontology_curie'] == ontology_curie
    indices = tracks.metadata.index[mask].tolist()
    return tracks.filter_tracks(indices)


def get_rna_track(tracks, rna_type='total'):
    """Filter RNA-seq tracks by type keyword (e.g. 'total', 'poly')."""
    mask = tracks.metadata['name'].str.contains(rna_type, case=False)
    return tracks.filter_tracks(mask.values)


def get_gtex_polya_track(tracks):
    """Filter RNA-seq tracks to GTEx polyA plus RNA-seq assays."""
    if 'name' not in tracks.metadata.columns:
        return tracks
    name_series = tracks.metadata['name'].astype(str).str.lower()
    mask = name_series.str.contains('gtex') & name_series.str.contains('polya plus rna-seq')
    return tracks.filter_tracks(mask.values)


def _exon_peaks_to_expression(tdata, exon_intervals, value_slice=None):
    """
    Compute mean of exon peaks: for each exon take max (peak) over that exon's
    bins, then average those peaks across exons. Returns array of shape (num_tracks,).
    If value_slice=(i_low, i_high) is given, only bins in that range are used
    (e.g. to restrict to the plotted window).
    """
    if not tdata.interval or not exon_intervals or tdata.values.size == 0:
        if value_slice is not None:
            i_l, i_h = value_slice
            return np.mean(tdata.values[i_l:i_h], axis=0)
        return np.mean(tdata.values, axis=0)
    inv = tdata.interval
    res = tdata.resolution
    n_bins = tdata.values.shape[0]
    slice_lo, slice_hi = value_slice if value_slice else (0, n_bins)
    peaks_per_track = []
    for exon in exon_intervals:
        if hasattr(exon, 'intersect') and hasattr(exon, 'overlaps') and exon.overlaps(inv):
            overlap = exon.intersect(inv)
        elif hasattr(exon, 'start'):
            ov_start = max(exon.start, inv.start)
            ov_end = min(exon.end, inv.end)
            overlap = genome.Interval(inv.chromosome, min(ov_start, ov_end), max(ov_start, ov_end))
        else:
            s, e = (exon[0], exon[1]) if isinstance(exon, (tuple, list)) else (exon['start'], exon['end'])
            ov_start = max(s, inv.start)
            ov_end = min(e, inv.end)
            overlap = genome.Interval(inv.chromosome, min(ov_start, ov_end), max(ov_start, ov_end))
        use_start = min(overlap.start, overlap.end)
        use_end = max(overlap.start, overlap.end)
        if use_start >= use_end:
            continue
        rel_start = use_start - inv.start
        rel_end = use_end - inv.start
        i_low = max(slice_lo, int(rel_start // res))
        i_high = min(slice_hi, int((rel_end + res - 1) // res))
        if i_low >= i_high:
            continue
        exon_slice = tdata.values[i_low:i_high]
        peak = np.max(exon_slice, axis=0)
        peaks_per_track.append(peak)
    if not peaks_per_track:
        return np.mean(tdata.values[slice_lo:slice_hi], axis=0)
    return np.mean(peaks_per_track, axis=0)


def _bin_range_for_interval(tdata, target_interval):
    """
    Return (i_low, i_high) such that tdata.values[i_low:i_high] covers the
    overlap of tdata.interval with target_interval. Returns (0, n_bins) if no overlap.
    """
    if not tdata.interval or not target_interval:
        return 0, tdata.values.shape[0]
    inv = tdata.interval
    if not hasattr(inv, 'overlaps') or not inv.overlaps(target_interval):
        return 0, tdata.values.shape[0]
    overlap = inv.intersect(target_interval)
    if overlap is None:
        return 0, tdata.values.shape[0]
    use_start = min(overlap.start, overlap.end)
    use_end = max(overlap.start, overlap.end)
    res = tdata.resolution
    n_bins = tdata.values.shape[0]
    rel_start = use_start - inv.start
    rel_end = use_end - inv.start
    i_low = max(0, int(rel_start // res))
    i_high = min(n_bins, int((rel_end + res - 1) // res))
    return i_low, max(i_low, i_high)


def aggregate_track_to_expression(tdata, exon_intervals=None, plot_interval=None):
    """
    Aggregate genomic track curves to a single expression value per track.
    If plot_interval is set, aggregation is restricted to that region (so reported
    values align with the plotted zoom window). If exon_intervals is provided,
    returns the mean of exon peaks (peak = max within each exon) in that region;
    otherwise the mean across the region.
    """
    if tdata.values.size == 0:
        return np.array([])
    i_low, i_high = _bin_range_for_interval(tdata, plot_interval)
    if plot_interval and exon_intervals and hasattr(plot_interval, 'overlaps'):
        exon_intervals = [e for e in exon_intervals if hasattr(e, 'overlaps') and e.overlaps(plot_interval)]
    if exon_intervals:
        return _exon_peaks_to_expression(tdata, exon_intervals, value_slice=(i_low, i_high))
    return np.mean(tdata.values[i_low:i_high], axis=0)


def track_to_aggregated_track(tdata, exon_intervals=None, plot_interval=None):
    """
    Convert genomic track curves to a constant track for plotting.
    If plot_interval is set, the constant is the aggregate over that region so
    the flat line aligns with the plotted window. If exon_intervals is provided,
    uses mean of exon peaks in that region; otherwise the mean across the region.
    """
    if tdata.values.size == 0:
        return tdata
    i_low, i_high = _bin_range_for_interval(tdata, plot_interval)
    if plot_interval and exon_intervals and hasattr(plot_interval, 'overlaps'):
        exon_intervals_in_plot = [e for e in exon_intervals if hasattr(e, 'overlaps') and e.overlaps(plot_interval)]
    else:
        exon_intervals_in_plot = exon_intervals
    if exon_intervals_in_plot:
        mean_per_track = _exon_peaks_to_expression(tdata, exon_intervals_in_plot, value_slice=(i_low, i_high))
    else:
        mean_per_track = np.mean(tdata.values[i_low:i_high], axis=0)
    mean_per_track = np.asarray(mean_per_track)
    if mean_per_track.ndim == 0:
        mean_per_track = mean_per_track[np.newaxis]
    mean_2d = np.broadcast_to(
        mean_per_track.reshape(1, -1), (tdata.values.shape[0], tdata.values.shape[1])
    )
    return track_data.TrackData(
        values=mean_2d.astype(tdata.values.dtype),
        metadata=tdata.metadata.copy(),
        interval=tdata.interval,
        resolution=tdata.resolution,
    )


# ---------------------------------------------------------------------------
# Track Summary Builder
# ---------------------------------------------------------------------------

def build_track_summary(outputs, ontology, gene_name, rsid, variant, exon_intervals=None, plot_interval=None):
    """
    Build a detailed track summary for a single variant as a list of strings,
    and apply select_best_track() to choose the appropriate track for each modality.

    exon_intervals: optional list of genome.Interval for exons. When provided,
        REF/ALT expression uses mean of exon peaks.
    plot_interval: optional genome.Interval for the plotted (zoom) window.
        When provided, expression is computed over this region only so reported
        values align with the expression plots.

    Returns
    -------
    summary_lines    : list of str
    selected         : dict of label -> (ref_track, alt_track)
    expression_report: list of dict with keys label, ref, alt (for REF vs ALT expression report)
    """
    ref = outputs.reference
    alt = outputs.alternate

    sep  = '=' * 60
    dash = chr(9472) * 60

    lines = [
        sep,
        f"  Variant : {gene_name} | {rsid} | {variant}",
        f"  Ontology: {ontology}",
        dash,
    ]

    expression_report = []

    def process_pair(label, key, ref_tdata, alt_tdata, log_expression=False):
        """Select best track for ref and alt, log results."""
        n_ref = ref_tdata.values.shape[1]
        n_alt = alt_tdata.values.shape[1]

        # If there are no tracks for this modality, log and skip adding it.
        if n_ref == 0 or n_alt == 0:
            lines.append(f"  {label:<30}   0 track(s)   [no matching tracks found — skipped]")
            return

        ref_selected, ref_log = select_best_track(ref_tdata)
        alt_selected, _       = select_best_track(alt_tdata)  # alt mirrors ref logic

        names = ref_tdata.metadata['name'].tolist() if 'name' in ref_tdata.metadata.columns else []
        name_str = ', '.join(names) if n_ref <= 4 else ', '.join(names[:4]) + f' ... (+{n_ref-4} more)'
        lines.append(f"  {label:<30} {n_ref:>3} track(s)   {name_str}")
        lines.extend(ref_log)
        if log_expression:
            ref_expr = aggregate_track_to_expression(
                ref_selected, exon_intervals=exon_intervals, plot_interval=plot_interval
            )
            alt_expr = aggregate_track_to_expression(
                alt_selected, exon_intervals=exon_intervals, plot_interval=plot_interval
            )
            ref_val = float(np.mean(ref_expr)) if ref_expr.size else np.nan
            alt_val = float(np.mean(alt_expr)) if alt_expr.size else np.nan
            region_note = " in plot window" if plot_interval else ""
            agg_note = f" (mean of exon peaks{region_note})" if exon_intervals else f" (mean{region_note})"
            lines.append(f"    → Aggregated expression{agg_note}: REF={ref_val:.4f}  ALT={alt_val:.4f}")
            expression_report.append({'label': label, 'ref': ref_val, 'alt': alt_val})
        selected[key] = (ref_selected, alt_selected)

    selected = {}

    process_pair(
        'RNA-seq (+) total',
        'rna_tot_pos',
        get_rna_track(ref.rna_seq.filter_to_positive_strand(), 'total'),
        get_rna_track(alt.rna_seq.filter_to_positive_strand(), 'total'),
        log_expression=True,
    )

    process_pair(
        'GTEx polyA (+)',
        'rna_poly_pos_gtex',
        get_gtex_polya_track(ref.rna_seq.filter_to_positive_strand()),
        get_gtex_polya_track(alt.rna_seq.filter_to_positive_strand()),
        log_expression=True,
    )

    process_pair(
        'RNA-seq (+) poly-A',
        'rna_poly_pos',
        get_rna_track(ref.rna_seq.filter_to_positive_strand(), 'poly'),
        get_rna_track(alt.rna_seq.filter_to_positive_strand(), 'poly'),
        log_expression=True,
    )

    process_pair(
        'RNA-seq (-) total',
        'rna_tot_neg',
        get_rna_track(ref.rna_seq.filter_to_negative_strand(), 'total'),
        get_rna_track(alt.rna_seq.filter_to_negative_strand(), 'total'),
        log_expression=True,
    )

    process_pair(
        'RNA-seq (-) poly-A',
        'rna_poly_neg',
        get_rna_track(ref.rna_seq.filter_to_negative_strand(), 'poly'),
        get_rna_track(alt.rna_seq.filter_to_negative_strand(), 'poly'),
        log_expression=True,
    )

    process_pair('ATAC', 'atac', ref.atac, alt.atac)

    process_pair(
        'CHIP-TF (CTCF)',
        'ctcf',
        filter_by_metadata(filter_by_ontology(ref.chip_tf, ontology), transcription_factor='CTCF'),
        filter_by_metadata(filter_by_ontology(alt.chip_tf, ontology), transcription_factor='CTCF'),
    )

    process_pair(
        'CHIP-Histone (H3K27ac)',
        'h3k27ac',
        filter_by_metadata(filter_by_ontology(ref.chip_histone, ontology), histone_mark='H3K27ac'),
        filter_by_metadata(filter_by_ontology(alt.chip_histone, ontology), histone_mark='H3K27ac'),
    )

    lines.append(dash)
    lines.append('')

    return lines, selected, expression_report


def format_expression_report(expression_report):
    """Format REF vs. ALT expression report as a list of strings (for printing and file output)."""
    if not expression_report:
        return []
    sep = '=' * 60
    lines = [
        sep,
        "  Track summary report: REF vs. ALT expression (aggregated mean)",
        sep,
        f"  {'Track':<28}   {'REF':>10}   {'ALT':>10}   {'ALT/REF':>10}",
        '  ' + '-' * 56,
    ]
    for row in expression_report:
        label, ref_val, alt_val = row['label'], row['ref'], row['alt']
        ratio = alt_val / ref_val if (ref_val != 0 and not np.isnan(ref_val)) else np.nan
        ratio_str = f"{ratio:.4f}" if not np.isnan(ratio) else "  — "
        ref_str = f"{ref_val:.4f}" if not np.isnan(ref_val) else "  — "
        alt_str = f"{alt_val:.4f}" if not np.isnan(alt_val) else "  — "
        lines.append(f"  {label:<28}   {ref_str:>10}   {alt_str:>10}   {ratio_str:>10}")
    lines.extend([sep, ''])
    return lines


def print_expression_report(expression_report):
    """Print the REF vs. ALT expression summary table."""
    lines = format_expression_report(expression_report)
    if lines:
        print('\n'.join(lines))


# ---------------------------------------------------------------------------
# Core Prediction & Plotting for a Single Variant
# ---------------------------------------------------------------------------

def run_variant(variant, gene_name, rsid, model, extractor, ontology,
                interval_size, zoom, output_dir, aggregate_rna=False):
    """
    Run the full prediction pipeline for a single variant and save the figure.

    Returns (figure_path, summary_lines)
    """
    print(f"\n{'='*60}")
    print(f"Processing: {gene_name} | {rsid} | {variant}")
    print(f"{'='*60}")

    interval = variant.reference_interval.resize(interval_size)

    # --- Run prediction ---
    outputs = model.predict_variant(
        interval=interval,
        variant=variant,
        ontology_terms=[ontology],
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.CHIP_TF,
            dna_client.OutputType.CHIP_HISTONE,
        },
    )

    # --- Transcript annotation and exon intervals (for exon-peak aggregation) ---
    transcripts = extractor.extract(interval)
    exon_intervals = []
    for t in transcripts:
        for exon in t.exons:
            exon_intervals.append(exon)

    # Plot window: expression values and flat-line aggregation use this so they align with the figure.
    plot_interval = variant.reference_interval.resize(zoom)

    # --- Build track summary and select best tracks ---
    summary_lines, selected, expression_report = build_track_summary(
        outputs, ontology, gene_name, rsid, variant,
        exon_intervals=exon_intervals,
        plot_interval=plot_interval,
    )
    print('\n'.join(summary_lines))
    print_expression_report(expression_report)

    # --- Assemble plot components (only for modalities that have data) ---
    def _rna_track(ref_t, alt_t):
        """Optionally convert RNA tracks to aggregated expression for plotting (over plot window)."""
        if aggregate_rna:
            return (
                track_to_aggregated_track(ref_t, exon_intervals=exon_intervals, plot_interval=plot_interval),
                track_to_aggregated_track(alt_t, exon_intervals=exon_intervals, plot_interval=plot_interval),
            )
        return ref_t, alt_t

    components = [plot_components.TranscriptAnnotation(transcripts)]

    if 'rna_tot_pos' in selected:
        r, a = _rna_track(*selected['rna_tot_pos'])
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': r, 'ALT': a},
                colors={'REF': 'dimgrey', 'ALT': 'red'},
                ylabel_template='Total RNA (+)',
            )
        )

    if 'rna_poly_pos_gtex' in selected:
        r, a = _rna_track(*selected['rna_poly_pos_gtex'])
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': r, 'ALT': a},
                colors={'REF': 'dimgrey', 'ALT': 'darkgreen'},
                ylabel_template='GTEx polyA RNA (+)',
            )
        )

    if 'rna_poly_pos' in selected:
        r, a = _rna_track(*selected['rna_poly_pos'])
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': r, 'ALT': a},
                colors={'REF': 'dimgrey', 'ALT': 'gold'},
                ylabel_template='Poly-A RNA (+)',
            )
        )

    if 'rna_tot_neg' in selected:
        r, a = _rna_track(*selected['rna_tot_neg'])
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': r, 'ALT': a},
                colors={'REF': 'dimgrey', 'ALT': 'red'},
                ylabel_template='Total RNA (-)',
            )
        )

    if 'rna_poly_neg' in selected:
        r, a = _rna_track(*selected['rna_poly_neg'])
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': r, 'ALT': a},
                colors={'REF': 'dimgrey', 'ALT': 'gold'},
                ylabel_template='Poly-A RNA (-)',
            )
        )

    if 'atac' in selected:
        atac_ref, atac_alt = selected['atac']
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': atac_ref, 'ALT': atac_alt},
                colors={'REF': 'dimgrey', 'ALT': 'blue'},
                ylabel_template='ATAC',
            )
        )

    if 'ctcf' in selected:
        ctcf_ref, ctcf_alt = selected['ctcf']
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': ctcf_ref, 'ALT': ctcf_alt},
                colors={'REF': 'dimgrey', 'ALT': 'orange'},
                ylabel_template='CTCF',
            )
        )

    if 'h3k27ac' in selected:
        h3k27ac_ref, h3k27ac_alt = selected['h3k27ac']
        components.append(
            plot_components.OverlaidTracks(
                tdata={'REF': h3k27ac_ref, 'ALT': h3k27ac_alt},
                colors={'REF': 'dimgrey', 'ALT': 'cyan'},
                ylabel_template='H3K27ac',
            )
        )

    cell_label = ontology.replace(':', '_')
    title = f"{gene_name} ({rsid}) — Multi-Modal Variant Effect [{cell_label}]"

    plot_components.plot(
        components,
        interval=variant.reference_interval.resize(zoom),
        annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
        title=title,
    )

    # --- Save figure ---
    safe_name = f"{gene_name}_{rsid}_{cell_label}_variant_effect.png"
    safe_name = safe_name.replace('/', '_').replace(' ', '_')
    save_path = os.path.join(output_dir, safe_name)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")

    return save_path, summary_lines, expression_report


# ---------------------------------------------------------------------------
# CSV Validation & Loading
# ---------------------------------------------------------------------------

REQUIRED_COLUMNS = {'gene_name', 'rsID', 'chr', 'position', 'ref_allele', 'alt_allele'}


def load_variants_csv(csv_path):
    """Load and validate the variants CSV. Returns a list of dicts."""
    try:
        df = pd.read_csv(csv_path)
    except UnicodeDecodeError:
        # Fallback for CSV files saved with non-UTF-8 encodings (common on Windows).
        df = pd.read_csv(csv_path, encoding='cp1252')
    df.columns = df.columns.str.strip()

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {missing}\n"
            f"Required: {REQUIRED_COLUMNS}\n"
            f"Found: {set(df.columns)}"
        )

    df['position']   = df['position'].astype(int)
    df['gene_name']  = df['gene_name'].astype(str).str.replace('\u00a0', ' ', regex=False).str.strip()
    df['rsID']       = df['rsID'].astype(str).str.replace('\u00a0', ' ', regex=False).str.strip()
    df['chr']        = df['chr'].astype(str).str.replace('\u00a0', ' ', regex=False).str.strip()
    df['ref_allele'] = df['ref_allele'].astype(str).str.replace('\u00a0', ' ', regex=False).str.strip()
    df['alt_allele'] = df['alt_allele'].astype(str).str.replace('\u00a0', ' ', regex=False).str.strip()
    if 'ontology_curie' in df.columns:
        df['ontology_curie'] = (
            df['ontology_curie']
            .astype(str)
            .str.replace('\u00a0', ' ', regex=False)
            .str.strip()
        )

    print(f"Loaded {len(df)} variant(s) from {csv_path}")
    return df.to_dict(orient='records')


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def main():
    args   = parse_args()
    run_ts = datetime.now().strftime('%Y%m%d_%H%M%S')   # e.g. 20250222_143501
    run_dt = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Validate output directory
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory : {os.path.abspath(args.output_dir)}")
    print(f"Run timestamp    : {run_dt}")

    # Load variants from CSV
    variants = load_variants_csv(args.csv)

    # Initialize model and annotation (once, shared)
    print("\nInitializing AlphaGenome model...")
    model     = dna_client.create(args.api_key)
    extractor = load_annotation()

    # Header block written once at the top of the combined summary file
    all_summary_lines = [
        '=' * 60,
        'AlphaGenome Variant Pipeline — Track Summary',
        '=' * 60,
        f"Run timestamp    : {run_dt}",
        f"Input CSV        : {os.path.abspath(args.csv)}",
        f"Output directory : {os.path.abspath(args.output_dir)}",
        f"Ontology         : {args.ontology}",
        f"Aggregate RNA    : {args.aggregate_rna}  (mean expression vs genomic curves)",
        f"Interval size    : {args.interval_size:,} bp  (model input context window)",
        f"Zoom             : {args.zoom:,} bp  (display window in figures)",
        f"Total variants   : {len(variants)}",
        '=' * 60,
        '',
    ]

    results = []
    errors  = []
    expression_report_rows = []  # Accumulate for REF vs. ALT expression CSV

    for i, row in enumerate(variants, start=1):
        gene_name = row['gene_name']
        rsid      = row['rsID']
        chrom     = row['chr']
        position  = row['position']
        ref       = row['ref_allele']
        alt       = row['alt_allele']
        # Per-variant ontology (if present) overrides global --ontology.
        variant_ontology = row.get('ontology_curie') or args.ontology

        try:
            variant = genome.Variant(
                chromosome=chrom,
                position=position,
                reference_bases=ref,
                alternate_bases=alt,
            )
            out_path, summary_lines, expression_report = run_variant(
                variant=variant,
                gene_name=gene_name,
                rsid=rsid,
                model=model,
                extractor=extractor,
                ontology=variant_ontology,
                interval_size=args.interval_size,
                zoom=args.zoom,
                output_dir=args.output_dir,
                aggregate_rna=args.aggregate_rna,
            )
            all_summary_lines.extend(summary_lines)
            results.append({
                'gene_name': gene_name, 'rsID': rsid,
                'status': 'success', 'output': out_path,
            })
            for er in expression_report:
                rv, av = er['ref'], er['alt']
                ratio = av / rv if (rv != 0 and not np.isnan(rv)) else np.nan
                expression_report_rows.append({
                    'gene_name': gene_name,
                    'rsID': rsid,
                    'ontology': variant_ontology,
                    'track': er['label'],
                    'REF': rv,
                    'ALT': av,
                    'ALT_over_REF': ratio if not np.isnan(ratio) else '',
                })

        except Exception as e:
            print(f"  ERROR processing {gene_name} ({rsid}): {e}")
            all_summary_lines.extend([
                '=' * 60,
                f"  Variant : {gene_name} | {rsid}",
                f"  STATUS  : ERROR — {e}",
                '=' * 60,
                '',
            ])
            errors.append({
                'gene_name': gene_name, 'rsID': rsid,
                'status': 'error', 'message': str(e),
            })

    # --- Write combined track summary file ---
    summary_txt_path = os.path.join(args.output_dir, f'track_summary_{run_ts}.txt')
    # Use UTF-8 so Unicode characters in log lines (e.g. em dashes) don't fail on Windows.
    with open(summary_txt_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(all_summary_lines))
    print(f"\nTrack summary saved : {summary_txt_path}")

    # --- Write REF vs. ALT expression report CSV ---
    expression_csv_path = os.path.join(args.output_dir, f'expression_report_{run_ts}.csv')
    if expression_report_rows:
        pd.DataFrame(expression_report_rows).to_csv(expression_csv_path, index=False)
        print(f"Expression report : {expression_csv_path}")

    # --- Write pipeline run summary CSV ---
    summary_csv_path = os.path.join(args.output_dir, f'pipeline_summary_{run_ts}.csv')
    pd.DataFrame(results + errors).to_csv(summary_csv_path, index=False)

    # --- Print final report ---
    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE  [{run_dt}]")
    print(f"  Success         : {len(results)}")
    print(f"  Errors          : {len(errors)}")
    print(f"  Figures         : {os.path.abspath(args.output_dir)}")
    print(f"  Track summary   : {summary_txt_path}")
    print(f"  Expression CSV  : {expression_csv_path if expression_report_rows else '(none)'}")
    print(f"  Run summary     : {summary_csv_path}")

    if errors:
        print("\nFailed variants:")
        for e in errors:
            print(f"  - {e['gene_name']} ({e['rsID']}): {e['message']}")


if __name__ == '__main__':
    main()
