"""
AlphaGenome Multi-Modal Variant Effect Prediction Pipeline
==========================================================
Author  : Zechuan Shi, UC Irvine
Contact : zechuas@uci.edu
Version : 1.0

Description:
    Batch variant effect prediction using the AlphaGenome model.
    Given a CSV of variants, this script runs multi-modal predictions
    (RNA-seq, ATAC, CTCF, H3K27ac) for each variant and saves a
    figure per variant to the specified output directory.

Usage:
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
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for script usage
import matplotlib.pyplot as plt

from alphagenome.data import gene_annotation, genome, transcript as transcript_utils, track_data
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components


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
        help='Zoom window size in bp for the final plot (default: 32768 = 32kb).'
    )
    parser.add_argument(
        '--interval_size', type=int, default=1048576,
        help='Full prediction interval size in bp (default: 1048576 = 1Mb).'
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
    print("Annotation loaded.")
    return extractor


# ---------------------------------------------------------------------------
# Track Helper Functions
# ---------------------------------------------------------------------------

def get_average_track(tdata):
    """Average multiple replicates into a single track."""
    if tdata.values.shape[1] <= 1:
        return tdata
    mean_values = tdata.values.mean(axis=1, keepdims=True)
    return track_data.TrackData(
        values=mean_values,
        metadata=tdata.metadata.iloc[[0]].copy(),
        interval=tdata.interval,
        resolution=tdata.resolution,
    )


def filter_by_metadata(tracks, **filters):
    """Filter tracks by arbitrary metadata column key=value pairs."""
    mask = pd.Series([True] * len(tracks.metadata), index=tracks.metadata.index)
    for key, value in filters.items():
        if key in tracks.metadata.columns:
            mask &= (tracks.metadata[key] == value)
    return tracks.filter_tracks(mask.values)


def filter_by_ontology(tracks, ontology_curie):
    """Filter tracks to a specific cell type ontology."""
    mask = tracks.metadata['ontology_curie'] == ontology_curie
    return tracks.filter_tracks(mask.values)


def get_rna_track(tracks, rna_type='total'):
    """Filter RNA-seq tracks by type keyword (e.g. 'total', 'poly')."""
    mask = tracks.metadata['name'].str.contains(rna_type, case=False)
    return tracks.filter_tracks(mask.values)


# ---------------------------------------------------------------------------
# Core Prediction & Plotting for a Single Variant
# ---------------------------------------------------------------------------

def run_variant(variant, gene_name, rsid, model, extractor, ontology, interval_size, zoom, output_dir):
    """
    Run the full prediction pipeline for a single variant and save the figure.

    Parameters
    ----------
    variant     : genome.Variant
    gene_name   : str
    rsid        : str
    model       : dna_client model handle
    extractor   : TranscriptExtractor
    ontology    : str  (e.g. 'EFO:0001187')
    interval_size : int
    zoom        : int
    output_dir  : str
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

    ref = outputs.reference
    alt = outputs.alternate

    # --- RNA-seq tracks ---
    rna_tot_pos_ref = get_average_track(get_rna_track(ref.rna_seq.filter_to_positive_strand(), 'total'))
    rna_tot_pos_alt = get_average_track(get_rna_track(alt.rna_seq.filter_to_positive_strand(), 'total'))
    rna_tot_neg_ref = get_average_track(get_rna_track(ref.rna_seq.filter_to_negative_strand(), 'total'))
    rna_tot_neg_alt = get_average_track(get_rna_track(alt.rna_seq.filter_to_negative_strand(), 'total'))
    rna_poly_neg_ref = get_average_track(get_rna_track(ref.rna_seq.filter_to_negative_strand(), 'poly'))
    rna_poly_neg_alt = get_average_track(get_rna_track(alt.rna_seq.filter_to_negative_strand(), 'poly'))

    # --- Epigenetic tracks ---
    # Filter to the requested cell type ontology first, then by mark/factor
    ref_chip_tf_ont  = filter_by_ontology(ref.chip_tf,      ontology)
    alt_chip_tf_ont  = filter_by_ontology(alt.chip_tf,      ontology)
    ref_chip_his_ont = filter_by_ontology(ref.chip_histone, ontology)
    alt_chip_his_ont = filter_by_ontology(alt.chip_histone, ontology)

    atac_ref    = get_average_track(ref.atac)
    atac_alt    = get_average_track(alt.atac)
    ctcf_ref    = get_average_track(filter_by_metadata(ref_chip_tf_ont,  transcription_factor='CTCF'))
    ctcf_alt    = get_average_track(filter_by_metadata(alt_chip_tf_ont,  transcription_factor='CTCF'))
    h3k27ac_ref = get_average_track(filter_by_metadata(ref_chip_his_ont, histone_mark='H3K27ac'))
    h3k27ac_alt = get_average_track(filter_by_metadata(alt_chip_his_ont, histone_mark='H3K27ac'))

    # --- Transcript annotation ---
    transcripts = extractor.extract(interval)

    # --- Assemble plot components ---
    components = [
        plot_components.TranscriptAnnotation(transcripts),
        plot_components.OverlaidTracks(
            tdata={'REF': rna_tot_pos_ref, 'ALT': rna_tot_pos_alt},
            colors={'REF': 'dimgrey', 'ALT': 'red'},
            ylabel_template='Total RNA (+)',
        ),
        plot_components.OverlaidTracks(
            tdata={'REF': rna_tot_neg_ref, 'ALT': rna_tot_neg_alt},
            colors={'REF': 'dimgrey', 'ALT': 'red'},
            ylabel_template='Total RNA (-)',
        ),
        plot_components.OverlaidTracks(
            tdata={'REF': rna_poly_neg_ref, 'ALT': rna_poly_neg_alt},
            colors={'REF': 'dimgrey', 'ALT': 'gold'},
            ylabel_template='Poly-A RNA (-)',
        ),
        plot_components.OverlaidTracks(
            tdata={'REF': atac_ref, 'ALT': atac_alt},
            colors={'REF': 'dimgrey', 'ALT': 'blue'},
            ylabel_template='ATAC',
        ),
        plot_components.OverlaidTracks(
            tdata={'REF': ctcf_ref, 'ALT': ctcf_alt},
            colors={'REF': 'dimgrey', 'ALT': 'orange'},
            ylabel_template='CTCF',
        ),
        plot_components.OverlaidTracks(
            tdata={'REF': h3k27ac_ref, 'ALT': h3k27ac_alt},
            colors={'REF': 'dimgrey', 'ALT': 'cyan'},
            ylabel_template='H3K27ac',
        ),
    ]

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
    return save_path


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

    # Basic type coercion
    df['position'] = df['position'].astype(int)
    df['chr'] = df['chr'].astype(str).str.strip()
    df['ref_allele'] = df['ref_allele'].astype(str).str.strip()
    df['alt_allele'] = df['alt_allele'].astype(str).str.strip()

    print(f"Loaded {len(df)} variant(s) from {csv_path}")
    return df.to_dict(orient='records')


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Validate output directory
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory: {os.path.abspath(args.output_dir)}")

    # Load variants from CSV
    variants = load_variants_csv(args.csv)

    # Initialize model and annotation (once, shared)
    print("\nInitializing AlphaGenome model...")
    model = dna_client.create(args.api_key)

    extractor = load_annotation()

    # Process each variant
    results = []
    errors = []

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
            out_path = run_variant(
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
            results.append({'gene_name': gene_name, 'rsID': rsid, 'status': 'success', 'output': out_path})

        except Exception as e:
            print(f"  ERROR processing {gene_name} ({rsid}): {e}")
            errors.append({'gene_name': gene_name, 'rsID': rsid, 'status': 'error', 'message': str(e)})

    # --- Summary ---
    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE")
    print(f"  Success : {len(results)}")
    print(f"  Errors  : {len(errors)}")
    print(f"  Figures saved to: {os.path.abspath(args.output_dir)}")

    if errors:
        print("\nFailed variants:")
        for e in errors:
            print(f"  - {e['gene_name']} ({e['rsID']}): {e['message']}")

    # Optionally save a run summary CSV
    summary_path = os.path.join(args.output_dir, 'pipeline_summary.csv')
    pd.DataFrame(results + errors).to_csv(summary_path, index=False)
    print(f"\nRun summary: {summary_path}")


if __name__ == '__main__':
    main()
