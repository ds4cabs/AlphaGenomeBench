

# AlphaGenome Variant Effect Pipeline

A batch pipeline for multi-modal variant effect prediction using [AlphaGenome](https://deepmind.google/blog/alphagenome-ai-for-better-understanding-the-genome/). Given a CSV of variants, the script predicts the effect of each variant on RNA-seq, chromatin accessibility (ATAC), CTCF binding, and histone modification (H3K27ac) in a specified cell type, and saves a figure for each variant.

---

## Requirements

- Python 3.10+
- AlphaGenome package and a valid API key ([request access here](https://deepmind.google.com/science/alphagenome/))
- Dependencies: `alphagenome`

---

## Input CSV Format

Prepare a CSV file with the following **required columns**:

**Note**: The **gene_name and rsID values** in the example CSV **may not be accurate**, they are placeholder examples for testing the code only. Please replace them with your actual variants before running.

| gene_name | rsID | chr | position | ref_allele | alt_allele |
|-----------|------|-----|----------|------------|------------|
| TBK1 | rs149000064 | chr12 | 64488488 | GAATT | G |
| TP53 | rs28934578 | chr17 | 7674220 | G | A |

- `position`: 1-based genomic coordinate (hg38)
- `ref_allele` / `alt_allele`: can be SNVs or indels
- Any extra columns in the CSV are safely ignored

A template CSV (`variants_input.csv`) is included in this repository.


---

## Usage

```bash
python ./PATH_to_code/alphagenome_variant_pipeline.py \
    --csv variants_input.csv \
    --output_dir ./alphagenome/figures \
    --ontology EFO:0001187 \
    --api_key YOUR_API_KEY
```

Code Running:
```
Output directory : ./alphagenome/figures
Run timestamp    : 2026-02-22 17:19:45
Loaded 4 variant(s) from variants_input.csv

Initializing AlphaGenome model...
Loading GTF annotation (this may take a moment on first run)...
Annotation loaded.

============================================================
Processing: TBK1 | rs149000064 | chr12:64488488:GAATT>G
============================================================
============================================================
  Variant : TBK1 | rs149000064 | chr12:64488488:GAATT>G
  Ontology: EFO:0001187
────────────────────────────────────────────────────────────
  RNA-seq (+) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) poly-A               1 track(s)   EFO:0001187 polyA plus RNA-seq
  ATAC                             1 track(s)   EFO:0001187 ATAC-seq
  CHIP-TF (CTCF)                   2 track(s)   EFO:0001187 TF ChIP-seq CTCF, EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF
    [WARNING] 2 different tracks found — selecting wildtype/unmodified:
      - EFO:0001187 TF ChIP-seq CTCF  [wildtype]
      - EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF  [MODIFIED - skipped]
    → Selected: EFO:0001187 TF ChIP-seq CTCF
  CHIP-Histone (H3K27ac)           1 track(s)   EFO:0001187 Histone ChIP-seq H3K27ac
────────────────────────────────────────────────────────────

Saved: ./alphagenome/figures/TBK1_rs149000064_EFO_0001187_variant_effect.png

============================================================
Processing: TP53 | rs28934578 | chr17:7674220:G>A
============================================================
============================================================
  Variant : TP53 | rs28934578 | chr17:7674220:G>A
  Ontology: EFO:0001187
────────────────────────────────────────────────────────────
  RNA-seq (+) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) poly-A               1 track(s)   EFO:0001187 polyA plus RNA-seq
  ATAC                             1 track(s)   EFO:0001187 ATAC-seq
  CHIP-TF (CTCF)                   2 track(s)   EFO:0001187 TF ChIP-seq CTCF, EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF
    [WARNING] 2 different tracks found — selecting wildtype/unmodified:
      - EFO:0001187 TF ChIP-seq CTCF  [wildtype]
      - EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF  [MODIFIED - skipped]
    → Selected: EFO:0001187 TF ChIP-seq CTCF
  CHIP-Histone (H3K27ac)           1 track(s)   EFO:0001187 Histone ChIP-seq H3K27ac
────────────────────────────────────────────────────────────

Saved: ./alphagenome/figures/TP53_rs28934578_EFO_0001187_variant_effect.png

============================================================
Processing: BRCA1 | rs80357713 | chr17:43071077:A>T
============================================================
============================================================
  Variant : BRCA1 | rs80357713 | chr17:43071077:A>T
  Ontology: EFO:0001187
────────────────────────────────────────────────────────────
  RNA-seq (+) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) poly-A               1 track(s)   EFO:0001187 polyA plus RNA-seq
  ATAC                             1 track(s)   EFO:0001187 ATAC-seq
  CHIP-TF (CTCF)                   2 track(s)   EFO:0001187 TF ChIP-seq CTCF, EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF
    [WARNING] 2 different tracks found — selecting wildtype/unmodified:
      - EFO:0001187 TF ChIP-seq CTCF  [wildtype]
      - EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF  [MODIFIED - skipped]
    → Selected: EFO:0001187 TF ChIP-seq CTCF
  CHIP-Histone (H3K27ac)           1 track(s)   EFO:0001187 Histone ChIP-seq H3K27ac
────────────────────────────────────────────────────────────

Saved: ./alphagenome/figures/BRCA1_rs80357713_EFO_0001187_variant_effect.png

============================================================
Processing: APOE | rs429358 | chr19:44908684:T>C
============================================================
============================================================
  Variant : APOE | rs429358 | chr19:44908684:T>C
  Ontology: EFO:0001187
────────────────────────────────────────────────────────────
  RNA-seq (+) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) total                1 track(s)   EFO:0001187 total RNA-seq
  RNA-seq (-) poly-A               1 track(s)   EFO:0001187 polyA plus RNA-seq
  ATAC                             1 track(s)   EFO:0001187 ATAC-seq
  CHIP-TF (CTCF)                   2 track(s)   EFO:0001187 TF ChIP-seq CTCF, EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF
    [WARNING] 2 different tracks found — selecting wildtype/unmodified:
      - EFO:0001187 TF ChIP-seq CTCF  [wildtype]
      - EFO:0001187 TF ChIP-seq CTCF genetically modified (insertion) using CRISPR targeting H. sapiens CTCF  [MODIFIED - skipped]
    → Selected: EFO:0001187 TF ChIP-seq CTCF
  CHIP-Histone (H3K27ac)           1 track(s)   EFO:0001187 Histone ChIP-seq H3K27ac
────────────────────────────────────────────────────────────

Saved: ./alphagenome/figures/APOE_rs429358_EFO_0001187_variant_effect.png

Track summary saved : ./alphagenome/figures/track_summary_20260222_171945.txt

============================================================
PIPELINE COMPLETE  [2026-02-22 17:19:45]
  Success       : 4
  Errors        : 0
  Figures       : ./alphagenome/figures
  Track summary : ./alphagenome/figures/track_summary_20260222_171945.txt
  Run summary   : ./alphagenome/figures/pipeline_summary_20260222_171945.csv

```

### All Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--csv` | ✅ | — | Path to input CSV file |
| `--output_dir` | ✅ | — | Directory to save output figures |
| `--api_key` | ✅ | — | Your AlphaGenome API key |
| `--ontology` | ❌ | `EFO:0001187` | Cell type ontology term (default: HepG2/Liver) |
| `--interval_size` | ❌ | `1048576` | Model input context window in bp (1 Mb). The model needs this large window for accurate predictions — only change if you have a specific reason |
| `--zoom` | ❌ | `32768` | Region displayed in the output figure in bp (32 kb). Does not affect prediction accuracy |

---

## Output

For each variant, the pipeline saves:
- A `.png` figure named `{gene_name}_{rsID}_{ontology}_variant_effect.png` showing REF vs ALT tracks for all modalities
- A `track_summary_{timestamp}.txt` logging the tracks available for each variant, which track was selected, and any warnings (e.g. CRISPR-modified tracks that were skipped)
- A `pipeline_summary_{timestamp}.csv` listing the status (success/error) of each variant

All output files are timestamped so multiple runs never overwrite each other.

Example figure for TBK1:

![TBK1 example](Plotting_RNA_ATAC_CTCF_H3K27ac/Figures/TBK1_rs149000064_EFO_0001187_variant_effect.png)

> *各位老师欢迎提出意见 More than happy to make changes*
