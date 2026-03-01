"""
AlphaGenome Multi-Modal Variant Effect Prediction Pipeline
==========================================================
Author  : Zechuan Shi, Swarup Lab, UC Irvine
Contact : zechuas@uci.edu
Version : 1.1.5

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
        [--ontology EFO:0001187] \
        [--zoom 32768] \
        [--interval_size 1048576]

CSV Format (required columns):
    gene_name, rsID, chr, position, ref_allele, alt_allele
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
        '--api_key', required=True,
        help='AlphaGenome API key.'
    )
    parser.add_argument(
        '--ontology', default='EFO:0001187',
        help='Ontology term for cell type (default: EFO:0001187 = HepG2/Liver).'
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
    return parser.parse_args()


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
    extractor.gtf = gtf_mane  # Attach the dataframe for easy access later
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
      - 0 tracks → return safely (prevents index errors)
      - 1 track  → use as-is
      - Multiple tracks, all identical names → average them (replicates)
      - Multiple tracks, different names     → prefer unmodified/wildtype tracks
    """
    # SAFETY CHECK: If the track data is empty, return immediately
    if tdata.values.shape[1] == 0:
        return tdata, ["    → No tracks found for this modality."]

    log_lines = []
    n = tdata.values.shape[1]

    if n == 1:
        return tdata, log_lines

    names = tdata.metadata['name'].tolist() if 'name' in tdata.metadata.columns else []

    # Case: all names are identical → true replicates, average them
    if len(set(names)) == 1:
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
        log_lines.append(f"    → Selected: {names[chosen_idx]}")
    else:
        chosen_idx = 0
        log_lines.append(f"    [WARNING] No unmodified track found — falling back to: {names[0]}")

    selected = tdata.filter_tracks(
        [i == chosen_idx for i in range(n)]
    )
    return selected, log_lines

# def filter_by_metadata(tracks, **filters):
#     """Filter tracks by arbitrary metadata column key=value pairs."""
#     mask = pd.Series([True] * len(tracks.metadata), index=tracks.metadata.index)
#     for key, value in filters.items():
#         if key in tracks.metadata.columns:
#             mask &= (tracks.metadata[key] == value)
#     return tracks.filter_tracks(mask.values)
#
# #
def filter_by_metadata(tracks, **filters):
    """Filter tracks safely, even if metadata is empty or columns are missing."""
    if tracks.metadata.empty:
        return tracks

    mask = pd.Series([True] * len(tracks.metadata), index=tracks.metadata.index)
    for key, value in filters.items():
        if key in tracks.metadata.columns:
            mask &= (tracks.metadata[key] == value)
        else:
            # If the specific column is missing, we don't apply the filter
            print(f"  Note: Metadata column '{key}' not found. Skipping filter.")
    return tracks.filter_tracks(mask.values)
#

# def filter_by_ontology(tracks, ontology_curie):
#     """Filter tracks to a specific cell type ontology."""
#     mask = tracks.metadata['ontology_curie'] == ontology_curie
#     return tracks.filter_tracks(mask.values)
def filter_by_ontology(tracks, ontology_curie):
    """Specific check for ontology in both 'ontology_curie' and 'name' columns."""
    if tracks.metadata.empty:
        return tracks

    # Check common columns where ontology might be stored
    for col in ['ontology_curie', 'name']:
        if col in tracks.metadata.columns:
            mask = tracks.metadata[col].str.contains(ontology_curie, case=False)
            return tracks.filter_tracks(mask.values)

    print(f"  Warning: No ontology info found in metadata columns.")
    return tracks




# def get_rna_track(tracks, rna_type='total'):
#     """Filter RNA-seq tracks by type keyword (e.g. 'total', 'poly')."""
#     mask = tracks.metadata['name'].str.contains(rna_type, case=False)
#     return tracks.filter_tracks(mask.values)

# def get_rna_track(tracks, rna_type='total'):
#     """Filter RNA-seq tracks by type keyword safely."""
#     # If no tracks were returned by the API, return as-is
#     if tracks.metadata.empty:
#         return tracks
#     # Check if 'name' column exists before string filtering
#     if 'name' in tracks.metadata.columns:
#         mask = tracks.metadata['name'].str.contains(rna_type, case=False)
#         return tracks.filter_tracks(mask.values)
#     return tracks
def get_rna_track(tracks, rna_type='total'):
    """Filter RNA-seq tracks by type keyword safely."""
    # 1. Safety check: Handle None (if model output is missing) or empty metadata
    if tracks is None or tracks.metadata.empty:
        return tracks

    # 2. Check if 'name' column exists before string filtering
    if 'name' in tracks.metadata.columns:
        # 3. Use na=False to handle cases where names might be missing (NaN)
        # This prevents the filter from crashing on unexpected data gaps
        mask = tracks.metadata['name'].str.contains(rna_type, case=False, na=False)
        return tracks.filter_tracks(mask.values)

    return tracks

# ---------------------------------------------------------------------------
# Track Summary Builder
# ---------------------------------------------------------------------------

def build_track_summary(outputs, ontology, gene_name, rsid, variant):
    """
    Build a detailed track summary for a single variant as a list of strings,
    and apply select_best_track() to choose the appropriate track for each modality.

    Returns
    -------
    summary_lines : list of str
    selected      : dict of label -> (ref_track, alt_track)
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

    def process_pair(label, ref_tdata, alt_tdata):
        """Select best track for ref and alt, log results."""
        ref_selected, ref_log = select_best_track(ref_tdata)
        alt_selected, _       = select_best_track(alt_tdata)  # alt mirrors ref logic
        n = ref_tdata.values.shape[1]
        names = ref_tdata.metadata['name'].tolist() if 'name' in ref_tdata.metadata.columns else []
        name_str = ', '.join(names) if n <= 4 else ', '.join(names[:4]) + f' ... (+{n-4} more)'
        lines.append(f"  {label:<30} {n:>3} track(s)   {name_str}")
        lines.extend(ref_log)
        return ref_selected, alt_selected

    selected = {}

    r, a = process_pair(
        'RNA-seq (+) total',
        get_rna_track(ref.rna_seq.filter_to_positive_strand(), 'total'),
        get_rna_track(alt.rna_seq.filter_to_positive_strand(), 'total'),
    )
    selected['rna_tot_pos'] = (r, a)

    r, a = process_pair(
        'RNA-seq (-) total',
        get_rna_track(ref.rna_seq.filter_to_negative_strand(), 'total'),
        get_rna_track(alt.rna_seq.filter_to_negative_strand(), 'total'),
    )
    selected['rna_tot_neg'] = (r, a)

    # r, a = process_pair(
    #     'RNA-seq (-) poly-A',
    #     get_rna_track(ref.rna_seq.filter_to_negative_strand(), 'poly'),
    #     get_rna_track(alt.rna_seq.filter_to_negative_strand(), 'poly'),
    # )
    # selected['rna_poly_neg'] = (r, a)
    # 3. Poly-A RNA (Flexible Strand for Brain/GTEx compatibility)
    r_poly = get_rna_track(ref.rna_seq, 'poly')
    a_poly = get_rna_track(alt.rna_seq, 'poly')

    r, a = process_pair(
        'RNA-seq poly-A',
        r_poly,
        a_poly
    )
    selected['rna_poly_neg'] = (r, a)  # Keep the key 'rna_poly_neg' so the plotting code doesn't need to change


    r, a = process_pair('ATAC', ref.atac, alt.atac)
    selected['atac'] = (r, a)

    r, a = process_pair(
        'CHIP-TF (CTCF)',
        filter_by_metadata(filter_by_ontology(ref.chip_tf, ontology), transcription_factor='CTCF'),
        filter_by_metadata(filter_by_ontology(alt.chip_tf, ontology), transcription_factor='CTCF'),
    )
    selected['ctcf'] = (r, a)

    r, a = process_pair(
        'CHIP-Histone (H3K27ac)',
        filter_by_metadata(filter_by_ontology(ref.chip_histone, ontology), histone_mark='H3K27ac'),
        filter_by_metadata(filter_by_ontology(alt.chip_histone, ontology), histone_mark='H3K27ac'),
    )
    selected['h3k27ac'] = (r, a)

    lines.append(dash)
    lines.append('')

    return lines, selected

# ---------------------------------------------------------------------------
# Quantification Function - Gene expression only
# ---------------------------------------------------------------------------

def quantify_gene_expression(outputs, interval, gtf_mane):
    """
    Calculates signal for every transcript in the 1Mb window
    and returns both the transcript/gene-level summaries.
    """
    interval_start = interval.start
    interval_end = interval.end

    # 1. Filter GTF to find all features within this window
    genes_in_window = gtf_mane[
        (gtf_mane['Chromosome'] == interval.chromosome) &
        (gtf_mane['Start'] < interval_end) &
        (gtf_mane['End'] > interval_start)
    ].copy()

    # This ensures you only have one row per unique transcript_id
    if 'Feature' in genes_in_window.columns:
        genes_in_window = genes_in_window[genes_in_window['Feature'] == 'transcript']

    gene_results = []

    # Use Column 2 (Poly-A RNA-seq) for mature mRNA
    ref_values = outputs.reference.rna_seq.values[:, 2]
    alt_values = outputs.alternate.rna_seq.values[:, 2]

    for _, entry in genes_in_window.iterrows():
        # 2. Map genomic coordinates to array indices
        rel_start = int(max(0, entry['Start'] - interval_start))
        rel_end = int(min(len(ref_values), entry['End'] - interval_start))

        # 3. Sum signal across the transcript body
        ref_sum = ref_values[rel_start:rel_end].sum()
        alt_sum = alt_values[rel_start:rel_end].sum()

        # 4. Calculate Log2 Fold Change
        if ref_sum > 0 and alt_sum > 0:
            log2fc = np.log2(alt_sum / ref_sum)
        else:
            log2fc = 0.0

        # Capture transcript-level detail (ISOFORMS)
        gene_results.append({
            'gene_name': entry['gene_name'],
            'transcript_id': entry.get('transcript_id', 'N/A'), # Keep transcript_id info
            'gene_type': entry['gene_type'],
            'Start': entry['Start'],
            'End': entry['End'],
            'REF_sum': round(float(ref_sum), 2),
            'ALT_sum': round(float(alt_sum), 2),
            'log2FC': round(float(log2fc), 6)
        })

    # VERSION 1: The Detailed Table (Isoforms)
    isoform_df = pd.DataFrame(gene_results) ### NOTE let's just keep this for now, in case we will need to deal with isoform later, but now isoform in this RNA-seq

    # # VERSION 2: The Grouped Table (One row per Gene)
    # gene_summary = isoform_df.groupby(['gene_name', 'gene_type']).agg({
    #     'REF_sum': 'sum',
    #     'ALT_sum': 'sum'
    # }).reset_index()

    # --- GROUPING TO GENE LEVEL ---
    # We aggregate to keep all your requested columns
    gene_summary = isoform_df.groupby(['gene_name', 'gene_type']).agg({
        'transcript_id': 'first', # Keeps the primary transcript ID
        'Start': 'min',           # Smallest start site
        'End': 'max',             # Largest end site
        'REF_sum': 'sum',         # Total signal
        'ALT_sum': 'sum'          # Total signal
    }).reset_index()

    # Final recalculation of log2FC for the whole gene
    gene_summary['log2FC'] = gene_summary.apply(
        lambda x: round(np.log2(x['ALT_sum']/x['REF_sum']), 6) if (x['REF_sum'] > 0 and x['ALT_sum'] > 0) else 0.0,
        axis=1
    )

    return gene_summary

    # # Re-calculate Log2FC for the whole gene
    # def calc_log2fc(row):
    #     if row['REF_sum'] > 0 and row['ALT_sum'] > 0:
    #         return np.log2(row['ALT_sum'] / row['REF_sum'])
    #     return 0.0
    #
    # gene_summary['log2FC'] = gene_summary.apply(calc_log2fc, axis=1)

    # # Return BOTH versions to the main script
    # return gene_summary, isoform_df


# ---------------------------------------------------------------------------
# Core Prediction & Plotting for a Single Variant
# ---------------------------------------------------------------------------

def run_variant(variant, gene_name, rsid, model, extractor, ontology,
                interval_size, zoom, output_dir):
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

    # --- Build track summary and select best tracks ---
    summary_lines, selected = build_track_summary(
        outputs, ontology, gene_name, rsid, variant
    )
    print('\n'.join(summary_lines))

    # --- Unpack selected tracks ---
    rna_tot_pos_ref,  rna_tot_pos_alt  = selected['rna_tot_pos']
    rna_tot_neg_ref,  rna_tot_neg_alt  = selected['rna_tot_neg']
    rna_poly_neg_ref, rna_poly_neg_alt = selected['rna_poly_neg']
    atac_ref,         atac_alt         = selected['atac']
    ctcf_ref,         ctcf_alt         = selected['ctcf']
    h3k27ac_ref,      h3k27ac_alt      = selected['h3k27ac']

    # --- Transcript annotation ---
    transcripts = extractor.extract(interval)

    # # --- NEW v1.1.2: Gene Expression Quantification ---
    # # We use extractor.gtf which is the filtered MANE gtf loaded at start
    # quant_df = quantify_gene_expression(outputs, interval, extractor.gtf)
    #
    # # Save the quantification table
    # cell_label = ontology.replace(':', '_')
    # quant_csv_name = f"{gene_name}_{rsid}_{cell_label}_quantification.csv"
    # quant_path = os.path.join(output_dir, quant_csv_name)
    # quant_df.sort_values('REF_sum', ascending=False).to_csv(quant_path, index=False)
    # print(f"Quantification saved: {quant_path}")

    # # --- NEW v1.1.3: Gene & Isoform Expression Quantification ---
    # # Unpack both the Gene-level and Transcript-level DataFrames
    # gene_df, isoform_df = quantify_gene_expression(outputs, interval, extractor.gtf)
    #
    # # Save the Gene Summary (Grouped)
    # cell_label = ontology.replace(':', '_')
    # gene_csv_name = f"{gene_name}_{rsid}_{cell_label}_GENE_summary.csv"
    # gene_path = os.path.join(output_dir, gene_csv_name)
    # gene_df.sort_values('REF_sum', ascending=False).to_csv(gene_path, index=False)
    #
    # # Save the Isoform Details (Individual transcripts)
    # iso_csv_name = f"{gene_name}_{rsid}_{cell_label}_ISOFORM_details.csv"
    # iso_path = os.path.join(output_dir, iso_csv_name)
    # isoform_df.sort_values('REF_sum', ascending=False).to_csv(iso_path, index=False)
    #
    # print(f"Quantification saved:")
    # print(f"  → Gene level: {gene_path}")
    # print(f"  → Isoform level: {iso_path}")

    # --- NEW v1.1.5: Gene & Isoform Expression Quantification ---
    # Unpack both the Gene-level and Transcript-level DataFrames
    gene_df = quantify_gene_expression(outputs, interval, extractor.gtf)

    # # Save the Gene Summary (Grouped)
    # cell_label = ontology.replace(':', '_')
    # gene_csv_name = f"{gene_name}_{rsid}_{cell_label}_GENE_summary.csv"
    # gene_path = os.path.join(output_dir, gene_csv_name)
    # gene_df.sort_values('REF_sum', ascending=False).to_csv(gene_path, index=False)

    # Save only the single summary CSV
    cell_label = ontology.replace(':', '_')
    csv_path = os.path.join(output_dir, f"{gene_name}_{rsid}_{cell_label}_summary.csv")
    gene_df.to_csv(csv_path, index=False)

    print(f"Quantification saved:")
    print(f"  → Gene level: {csv_path}")

    # # NOTE : modify the run_variant function.
    # # Currently, your script assumes every data type (RNA, ATAC, ChIP) is always available,
    # # which causes it to crash when a tissue like Brain is missing those specific tracks.
    # # --- Assemble plot components ---
    # components = [
    #     plot_components.TranscriptAnnotation(transcripts),
    #     plot_components.OverlaidTracks(
    #         tdata={'REF': rna_tot_pos_ref, 'ALT': rna_tot_pos_alt},
    #         colors={'REF': 'dimgrey', 'ALT': 'red'},
    #         ylabel_template='Total RNA (+)',
    #     ),
    #     plot_components.OverlaidTracks(
    #         tdata={'REF': rna_tot_neg_ref, 'ALT': rna_tot_neg_alt},
    #         colors={'REF': 'dimgrey', 'ALT': 'red'},
    #         ylabel_template='Total RNA (-)',
    #     ),
    #     plot_components.OverlaidTracks(
    #         tdata={'REF': rna_poly_neg_ref, 'ALT': rna_poly_neg_alt},
    #         colors={'REF': 'dimgrey', 'ALT': 'gold'},
    #         ylabel_template='Poly-A RNA (-)',
    #     ),
    #     plot_components.OverlaidTracks(
    #         tdata={'REF': atac_ref, 'ALT': atac_alt},
    #         colors={'REF': 'dimgrey', 'ALT': 'blue'},
    #         ylabel_template='ATAC',
    #     ),
    #     plot_components.OverlaidTracks(
    #         tdata={'REF': ctcf_ref, 'ALT': ctcf_alt},
    #         colors={'REF': 'dimgrey', 'ALT': 'orange'},
    #         ylabel_template='CTCF',
    #     ),
    #     plot_components.OverlaidTracks(
    #         tdata={'REF': h3k27ac_ref, 'ALT': h3k27ac_alt},
    #         colors={'REF': 'dimgrey', 'ALT': 'cyan'},
    #         ylabel_template='H3K27ac',
    #     ),
    # ]
    # --- Assemble plot components ---
    # Start with the transcript annotation (always present)
    components = [plot_components.TranscriptAnnotation(transcripts)]

    # Add RNA Total (+) track if data exists
    if not rna_tot_pos_ref.metadata.empty:
        components.append(plot_components.OverlaidTracks(
            tdata={'REF': rna_tot_pos_ref, 'ALT': rna_tot_pos_alt},
            colors={'REF': 'dimgrey', 'ALT': 'red'},
            ylabel_template='Total RNA (+)',
        ))

    # Add RNA Total (-) track if data exists
    if not rna_tot_neg_ref.metadata.empty:
        components.append(plot_components.OverlaidTracks(
            tdata={'REF': rna_tot_neg_ref, 'ALT': rna_tot_neg_alt},
            colors={'REF': 'dimgrey', 'ALT': 'red'},
            ylabel_template='Total RNA (-)',
        ))

    # Add RNA Poly-A track if data exists
    if not rna_poly_neg_ref.metadata.empty:
        components.append(plot_components.OverlaidTracks(
            tdata={'REF': rna_poly_neg_ref, 'ALT': rna_poly_neg_alt},
            colors={'REF': 'dimgrey', 'ALT': 'gold'},
            ylabel_template='Poly-A RNA (-)',
        ))

    # Add ATAC track if data exists
    if not atac_ref.metadata.empty:
        components.append(plot_components.OverlaidTracks(
            tdata={'REF': atac_ref, 'ALT': atac_alt},
            colors={'REF': 'dimgrey', 'ALT': 'blue'},
            ylabel_template='ATAC',
        ))

    # Add CTCF track if data exists
    if not ctcf_ref.metadata.empty:
        components.append(plot_components.OverlaidTracks(
            tdata={'REF': ctcf_ref, 'ALT': ctcf_alt},
            colors={'REF': 'dimgrey', 'ALT': 'orange'},
            ylabel_template='CTCF',
        ))

    # Add H3K27ac track if data exists
    if not h3k27ac_ref.metadata.empty:
        components.append(plot_components.OverlaidTracks(
            tdata={'REF': h3k27ac_ref, 'ALT': h3k27ac_alt},
            colors={'REF': 'dimgrey', 'ALT': 'cyan'},
            ylabel_template='H3K27ac',
        ))

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

    # return save_path, summary_lines, gene_df, isoform_df  # v1.1.3 -- Add gene_df and isoform_df here
    return save_path, summary_lines, gene_df  # v1.1.5


# ---------------------------------------------------------------------------
# CSV Validation & Loading
# ---------------------------------------------------------------------------

REQUIRED_COLUMNS = {'gene_name', 'rsID', 'chr', 'position', 'ref_allele', 'alt_allele'}


def load_variants_csv(csv_path):
    """Load and validate the variants CSV. Returns a list of dicts."""
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {missing}\n"
            f"Required: {REQUIRED_COLUMNS}\n"
            f"Found: {set(df.columns)}"
        )

    df['position']   = df['position'].astype(int)
    df['chr']        = df['chr'].astype(str).str.strip()
    df['ref_allele'] = df['ref_allele'].astype(str).str.strip()
    df['alt_allele'] = df['alt_allele'].astype(str).str.strip()

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
        f"Interval size    : {args.interval_size:,} bp  (model input context window)",
        f"Zoom             : {args.zoom:,} bp  (display window in figures)",
        f"Total variants   : {len(variants)}",
        '=' * 60,
        '',
    ]

    results = []
    errors  = []

    for i, row in enumerate(variants, start=1):
        gene_name = row['gene_name']
        rsid      = row['rsID']
        chrom     = row['chr']
        position  = row['position']
        ref       = row['ref_allele']
        alt       = row['alt_allele']

        try:
            variant = genome.Variant(
                chromosome=chrom,
                position=position,
                reference_bases=ref,
                alternate_bases=alt,
            )
            # # v1.1.2
            # out_path, summary_lines = run_variant(
            #     variant=variant,
            #     gene_name=gene_name,
            #     rsid=rsid,
            #     model=model,
            #     extractor=extractor,
            #     ontology=args.ontology,
            #     interval_size=args.interval_size,
            #     zoom=args.zoom,
            #     output_dir=args.output_dir,
            # )
            # all_summary_lines.extend(summary_lines)
            # results.append({
            #     'gene_name': gene_name, 'rsID': rsid,
            #     'status': 'success', 'output': out_path,
            # })

            # # v1.1.3
            # Unpack all 4 returned objects
            # out_path, summary_lines, gene_df, isoform_df = run_variant(
            # # v1.1.5
            out_path, summary_lines, gene_df = run_variant(
                variant=variant,
                gene_name=gene_name,
                rsid=rsid,
                model=model,
                extractor=extractor,
                ontology=args.ontology,
                interval_size=args.interval_size,
                zoom=args.zoom,
                output_dir=args.output_dir,
            )
            all_summary_lines.extend(summary_lines)

            # Identify the "Top Hit" from the gene-level results
            top_hit = gene_df.iloc[gene_df['log2FC'].abs().argmax()]

            results.append({
                'gene_name': gene_name,
                'rsID': rsid,
                'chr': chrom,
                'pos': position,
                'top_affected_gene': top_hit['gene_name'],
                'max_log2FC': round(float(top_hit['log2FC']), 4),
                'ref_signal': round(float(top_hit['REF_sum']), 2),
                'status': 'success',
                'output_figure': out_path,
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
    with open(summary_txt_path, 'w') as f:
        f.write('\n'.join(all_summary_lines))
    print(f"\nTrack summary saved : {summary_txt_path}")

    # --- Write pipeline run summary CSV ---
    summary_csv_path = os.path.join(args.output_dir, f'pipeline_summary_{run_ts}.csv')
    pd.DataFrame(results + errors).to_csv(summary_csv_path, index=False)

    # --- Print final report ---
    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE  [{run_dt}]")
    print(f"  Success       : {len(results)}")
    print(f"  Errors        : {len(errors)}")
    print(f"  Figures       : {os.path.abspath(args.output_dir)}")
    print(f"  Track summary : {summary_txt_path}")
    print(f"  Run summary   : {summary_csv_path}")

    if errors:
        print("\nFailed variants:")
        for e in errors:
            print(f"  - {e['gene_name']} ({e['rsID']}): {e['message']}")


if __name__ == '__main__':
    main()
