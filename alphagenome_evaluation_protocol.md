# Independent Evaluation of AlphaGenome: A Comprehensive Research Protocol

## Consortium Research Protocol v1.0
**Prepared for:** International Statistical and Population Genetics Consortium  
**Date:** February 2026  
**Target Publication:** Nature Methods / Nature Genetics / Genome Biology

---

## Executive Summary

This protocol outlines a comprehensive, systematic independent evaluation of AlphaGenome, the unified deep learning model for regulatory variant effect prediction published by Google DeepMind (Nature, January 2026). Our consortium of 20 researchers will conduct evaluations across multiple dimensions not fully explored in the original publication, with emphasis on:

1. Cross-ancestry generalizability
2. Novel disease cohort applications
3. Comparison with emerging models
4. Clinical utility assessment
5. Computational reproducibility

---

## Part I: Study Design Framework

### 1.1 Evaluation Objectives

| Objective | Priority | Novel Contribution |
|-----------|----------|-------------------|
| Cross-ancestry variant effect prediction | High | Original paper primarily used European ancestry data |
| Rare disease diagnostic utility | High | Clinical application not systematically evaluated |
| Comparison with post-publication models | High | Evo 2, Nucleotide Transformer v2, etc. |
| Long-read sequencing integration | Medium | ONT/PacBio structural variant prediction |
| Pharmacogenomic variant interpretation | Medium | Drug response prediction |
| Non-human primate transferability | Medium | Cross-species validation for preclinical research |
| Single-cell context predictions | Medium | Cell-type specificity beyond bulk data |
| Computational benchmark | High | Runtime, memory, reproducibility |

### 1.2 Consortium Structure

#### Work Package Organization (20 members)

```
WP1: Track Prediction Evaluation (4 members)
├── WP1.1: Expression track validation (2)
└── WP1.2: Chromatin/epigenomic tracks (2)

WP2: Variant Effect Prediction (6 members)
├── WP2.1: eQTL/sQTL evaluation (2)
├── WP2.2: Chromatin QTL evaluation (2)
└── WP2.3: Clinical variant interpretation (2)

WP3: Cross-Ancestry & Population Genetics (4 members)
├── WP3.1: Multi-ancestry QTL analysis (2)
└── WP3.2: Population-specific regulatory variation (2)

WP4: Disease Applications (3 members)
├── WP4.1: Mendelian disease variants (1.5)
└── WP4.2: Complex trait fine-mapping (1.5)

WP5: Benchmarking & Methods (3 members)
├── WP5.1: Model comparison framework (1.5)
└── WP5.2: Computational reproducibility (1.5)
```

### 1.3 Timeline (12 months)

| Phase | Duration | Activities |
|-------|----------|------------|
| Phase 1 | Months 1-2 | Data acquisition, pipeline setup, API access |
| Phase 2 | Months 3-6 | Primary evaluations across all WPs |
| Phase 3 | Months 7-9 | Secondary analyses, cross-WP integration |
| Phase 4 | Months 10-11 | Manuscript preparation |
| Phase 5 | Month 12 | Internal review, submission |

---

## Part II: Data Resources

### 2.1 Independent Datasets (Not Used in Original Training)

#### Expression QTLs

| Dataset | Ancestry | Tissues | Samples | Access |
|---------|----------|---------|---------|--------|
| eQTLGen Phase 2 | European | Blood | >30,000 | Public |
| GENOA | African American | Blood, lymphocytes | 1,263 | dbGaP |
| MESA | Multi-ethnic | Monocytes | 1,264 | dbGaP |
| GTEx v10 (if available) | Multi-ethnic | 54 tissues | ~1,000 | dbGaP |
| BioVU eQTL | Multi-ethnic | Multiple | >100,000 | Collaboration |
| Japanese eQTL | East Asian | Blood | 10,000+ | Collaboration |
| INTERVAL | European | Blood | 3,301 | EGA |
| BLUEPRINT | European | Immune cells | 197 | EGA |

### 2.2 Biobank GWAS & PheWAS Resources

| Resource | Population | Phenotypes | Samples | Key Features |
|----------|------------|------------|---------|--------------|
| **FinnGen R12** | Finnish | 2,400+ | 500,000+ | Summary stats, fine-mapping, founder enrichment |
| **UK Biobank** | European | 7,000+ | 500,000 | Exome + imputed, deep phenotyping |
| **AZ PheWAS** | Multi-ethnic | 1,500+ | 450,000+ | AstraZeneca internal + UKB, rare variants |
| **Genebass** | European | 4,500+ | 400,000 | Gene-burden tests, exome-wide associations |
| Biobank Japan | Japanese | 200+ | 200,000+ | East Asian GWAS, fine-mapping |
| All of Us | Multi-ethnic | 2,000+ | 300,000+ | Diverse ancestry, WGS available |

### 2.3 Variant Annotation & Population Databases

| Resource | Content | Size | Key Use Cases |
|----------|---------|------|---------------|
| **gnomAD v4** | Population frequencies | 800K exomes, 76K genomes | AF filtering, constraint (pLI, LOEUF) |
| **ClinVar (2026)** | Clinical interpretations | >2.5M submissions | Pathogenic/benign classification |
| **OMIM** | Gene-disease relationships | >16,000 genes | Mendelian disease annotation |
| HGMD Professional | Disease mutations | >350K variants | Literature-curated pathogenic |
| dbSNP b156 | Variant catalog | >1B variants | rsID mapping, validation |
| ClinGen | Gene/variant curation | >1,500 genes | Expert-curated validity |

### 2.4 Integrated Genetics Platforms

| Platform | Data Types | Key Features | Evaluation Use |
|----------|------------|--------------|----------------|
| **Open Targets Genetics** | GWAS, eQTL, chromatin | L2G scores, coloc, fine-mapping | Variant-to-gene validation |
| **Open Targets Platform** | Drug targets, diseases | Target tractability, safety | Clinical relevance scoring |
| GWAS Catalog | Published associations | >500K associations | Benchmark variant sets |
| PhenoScanner v2 | Cross-phenotype lookup | GWAS + eQTL + pQTL | Pleiotropy analysis |
| Ensembl VEP | Variant annotation | Consequence prediction | Functional annotation |

### 2.5 Splicing QTLs

| Dataset | Description | Variants |
|---------|-------------|----------|
| sQTLseeker2 results | Pan-tissue sQTL catalog | ~2M |
| Leafcutter GTEx reanalysis | Intron excision ratios | ~1.5M |
| ROSMAP brain sQTLs | Alzheimer's cohort | ~50K |
| CommonMind Consortium | Brain psychiatric disorders | ~100K |

### 2.6 Chromatin QTLs

| Dataset | Assay Type | Cell Types |
|---------|------------|------------|
| ENCODE4 | DNase, ATAC, ChIP | 200+ |
| Roadmap Epigenomics (new) | Histone ChIP | 127 |
| HaploReg v5 | Integrated chromatin | Multi-tissue |
| EpiMap | Chromatin states | 833 samples |
| CATlas | Single-cell ATAC | 30+ cell types |

### 2.7 Clinical Variants

| Resource | Variant Type | Count | Access |
|----------|--------------|-------|--------|
| ClinVar (2026 release) | Pathogenic/Benign/VUS | >2.5M | Public FTP |
| OMIM | Gene-phenotype | >16,000 genes | API/License |
| HGMD Professional | Disease mutations | >350K | License |
| Deciphering Developmental Disorders | De novo | ~10K | EGA |
| SPARK autism | De novo + inherited | ~15K families | SFARI Base |
| ClinGen Expert Panels | Curated assertions | >20K variants | Public |

### 2.8 GWAS Fine-Mapping

| Resource | Traits | Credible Sets | Method |
|----------|--------|---------------|--------|
| Open Targets Genetics v8 | >5,000 | >50,000 | SuSiE + FINEMAP |
| FinnGen R12 Fine-mapping | 2,400+ | >15,000 | SuSiE, founder LD |
| GWAS Catalog fine-mapping | >1,000 | >20,000 | Multiple methods |
| Biobank Japan | 200+ | >5,000 | FINEMAP |
| UKBB (Weissbrod) | 94 | ~2,000 | PolyFun + SuSiE |
| Genebass credible sets | 4,500+ | Gene-level | Burden test |

### 2.9 Novel Evaluation Datasets

#### Functional Assays Not Used in Original Paper

| Assay | Description | Variants Tested |
|-------|-------------|-----------------|
| MPRA (Tewhey) | Liver enhancer variants | 30,000 |
| MPRA (Ulirsch) | Erythroid variants | 2,756 |
| STARR-seq (Arnold) | Enhancer activity | 100,000+ |
| CRISPRi screens (Gasperini) | Enhancer-gene pairs | 5,920 |
| Perturb-seq | Single-cell perturbations | 10,000+ |
| BASE editing saturation | Coding + regulatory | Variable |

#### Structural Variant Resources

| Dataset | SV Types | Source |
|---------|----------|--------|
| gnomAD-SV v4 | Deletions, duplications, inversions | 76,156 genomes |
| 1000 Genomes SV | Multi-ancestry SVs | 2,504 individuals |
| HGSVC2 | Long-read resolved | 35 individuals |

---

## Part III: Database-Specific Evaluation Strategies

### 3.1 FinnGen-Specific Analyses

FinnGen provides unique opportunities due to Finnish founder population enrichment for rare variants:

**Evaluation Approaches:**
- **Founder variant effect prediction**: Test on variants enriched 10-100x in Finland vs. other populations
- **Fine-mapping validation**: Compare AlphaGenome scores to SuSiE posterior inclusion probabilities (PIPs)
- **Loss-of-function variant prioritization**: Predict regulatory LoF effects for enriched null alleles
- **Phenome-wide association**: Correlate predicted effects with PheWAS results across 2,400+ phenotypes

**Specific Tests:**
```
Analysis: Finnish Founder Variant Evaluation
1. Identify variants with AF(FinnGen) / AF(gnomAD-NFE) > 10
2. Score with AlphaGenome
3. Compare predicted effects to observed PheWAS associations
4. Evaluate enrichment of high-score variants in disease associations
```

### 3.2 Open Targets Integration

Open Targets Genetics provides variant-to-gene (V2G) assignments that can benchmark AlphaGenome:

**L2G Score Comparison:**
- Does AlphaGenome improve gene assignment beyond Open Targets L2G?
- Correlation between AlphaGenome expression scores and L2G scores
- Performance on credible sets where L2G assignment is uncertain

**Colocalization Validation:**
- Compare eQTL-GWAS colocalization results with AlphaGenome predictions
- Test if AlphaGenome can identify causal variants in colocalized loci
- Evaluate enhancer-gene predictions vs. Open Targets enhancer annotations

### 3.3 gnomAD Constraint Integration

gnomAD constraint metrics provide complementary information to sequence-based predictions:

**Evaluation Approaches:**
- **Constraint-stratified evaluation**: Separate performance analysis for high-pLI (>0.9) vs low-pLI genes
- **LOEUF correlation**: Test whether AlphaGenome regulatory scores correlate with loss-of-function intolerance
- **Rare variant AF filtering**: Evaluate predictions on gnomAD rare variants (AF < 0.001)
- **Regional constraint (RMC)**: Test within-gene regional differences using gnomAD regional missense constraint

### 3.4 ClinVar Clinical Validation

Systematic evaluation on ClinVar variants for clinical utility assessment:

**Classification Tasks:**
| Task | Gold Standard | Target Performance |
|------|---------------|-------------------|
| Pathogenic vs Benign | Expert-reviewed submissions | auROC > 0.85 |
| VUS prioritization | Later-resolved VUS | Top-quartile enrichment > 2x |
| Non-coding pathogenic | Regulatory region variants | Sensitivity > 70% |

**Temporal Validation:**
- Train/test split by submission date (pre-2024 vs 2024-2026)
- Evaluate on variants submitted after model training

### 3.5 OMIM Mendelian Disease Analysis

OMIM provides curated gene-disease relationships for Mendelian conditions:

**Evaluation Approaches:**
- **Regulatory variant discovery**: Identify non-coding pathogenic variants in OMIM disease genes
- **Promoter mutation prediction**: Benchmark on known promoter pathogenic variants from literature
- **Enhancer mutation detection**: Test long-range regulatory variant effects in OMIM genes
- **Tissue-specificity validation**: Match predicted tissue expression to disease manifestation patterns

### 3.6 Genebass Rare Variant Analysis

Genebass exome-wide associations provide gene-level validation:

**Evaluation Approaches:**
- **Burden test weighting**: Use AlphaGenome scores to weight rare variants in burden tests
- **Regulatory burden extension**: Extend burden tests to non-coding regulatory variants
- **Splice variant contribution**: Validate AlphaGenome splicing predictions with Genebass exome data
- **Gene-level effect correlation**: Compare predicted gene-level effects with observed associations

### 3.7 AZ PheWAS Integration

AstraZeneca PheWAS data provides industry-scale rare variant associations:

**Unique Opportunities:**
- Large-scale rare variant effects across diverse phenotypes
- Integration of AstraZeneca internal data with UK Biobank
- Cross-validation with pharmaceutical target knowledge

---

## Part IV: Primary Evaluation Strategies

### 4.1 Track Prediction Evaluation (WP1)

#### 3.1.1 Expression Prediction

**Objective:** Assess AlphaGenome's ability to predict gene expression across independent datasets.

**Metrics:**
- Pearson correlation (per-gene, per-track)
- Spearman correlation
- Mean squared error
- Tissue-specificity accuracy

**Analyses:**

```
Analysis 1.1: ENCODE4 RNA-seq Validation
- Compare predictions to new ENCODE4 RNA-seq data
- Focus on cell types not in training (e.g., new ENCODE cell lines)
- Stratify by expression level, gene length, GC content

Analysis 1.2: Single-Cell Expression Prediction
- Use pseudo-bulk from single-cell atlases (Human Cell Atlas)
- Evaluate cell-type specificity predictions
- Test on rare cell populations

Analysis 1.3: Condition-Specific Expression
- Stimulus-response data (e.g., immune activation)
- Disease vs. control comparisons
- Drug treatment effects
```

**Key Questions:**
1. Does AlphaGenome generalize to cell types not in training?
2. How well does it predict expression changes upon perturbation?
3. What is the performance on lowly-expressed genes?

#### 3.1.2 Chromatin Track Prediction

**Analyses:**

```
Analysis 1.4: ENCODE4 Chromatin Accessibility
- New ATAC-seq and DNase-seq datasets
- Single-cell ATAC (scATAC) pseudo-bulk validation
- Peak-level vs. genome-wide evaluation

Analysis 1.5: Histone Modification Prediction
- CUT&Tag data (higher resolution than ChIP)
- Novel histone marks (H3K27me2, H3K36me2)
- Bivalent domain prediction

Analysis 1.6: 3D Chromatin Architecture
- Micro-C data (higher resolution than Hi-C)
- HiChIP enhancer-promoter loops
- CTCF loop anchor prediction
```

### 4.2 Variant Effect Prediction (WP2)

#### 3.2.1 Expression QTL Evaluation

**Objective:** Systematically benchmark eQTL effect prediction across ancestries and tissues.

**Primary Metrics:**
- Sign prediction accuracy (auROC)
- Effect size correlation (Spearman ρ)
- Causality discrimination (auPRC)

**Analyses:**

```
Analysis 2.1: Multi-Ancestry eQTL Prediction
- Compare performance: EUR vs AFR vs EAS vs SAS vs AMR
- Test on ancestry-specific eQTLs
- Evaluate transferability of predictions

Analysis 2.2: Tissue-Specific eQTL Prediction
- Stratify by tissue specificity (τ index)
- Focus on hard-to-predict tissues (brain, testis)
- Cell-type-specific eQTLs from sorted populations

Analysis 2.3: Context-Dependent eQTLs
- Response eQTLs (reQTLs) - stimulus-dependent
- Age-dependent eQTLs
- Sex-specific eQTLs

Analysis 2.4: Distance-to-TSS Analysis
- Fine-grained distance bins (0-1kb, 1-5kb, 5-10kb, 10-50kb, 50-100kb, >100kb)
- Enhancer eQTLs specifically
- Comparison with Activity-by-Contact predictions
```

**Novel Evaluation: eQTL Fine-Mapping Integration**

```python
# Pseudo-code for fine-mapping evaluation
def evaluate_fine_mapping(credible_sets, alphagenome_scores):
    """
    Assess whether AlphaGenome prioritizes causal variants
    within statistically-defined credible sets
    """
    metrics = {}
    
    for cs in credible_sets:
        # Rank variants by AlphaGenome score
        ag_ranks = rank_by_score(cs.variants, alphagenome_scores)
        
        # Compare to PIP (posterior inclusion probability)
        pip_correlation = spearman(ag_ranks, cs.pip_values)
        
        # Check if top AG variant matches top PIP variant
        top_match = (ag_ranks[0] == cs.top_pip_variant)
        
        # Enrichment of high-PIP variants in top AG decile
        enrichment = calculate_enrichment(ag_ranks, cs.pip_values, 
                                          threshold=0.1)
    
    return metrics
```

#### 3.2.2 Splicing QTL Evaluation

**Analyses:**

```
Analysis 2.5: sQTL Effect Prediction
- Leafcutter intron excision ratios
- PSI (percent spliced in) prediction
- Splice site usage vs. junction count predictions

Analysis 2.6: Rare Splicing Variants
- MFASS functional validation data
- SpliceAI benchmarks (independent test set)
- Cryptic splice site activation

Analysis 2.7: Tissue-Specific Splicing
- Brain-specific isoforms
- Muscle-specific splicing
- Disease-relevant alternative splicing
```

#### 3.2.3 Chromatin QTL Evaluation

**Analyses:**

```
Analysis 2.8: Chromatin Accessibility QTLs
- ATAC-QTLs from multiple ancestries
- DNase-QTLs with matched cell types
- Single-cell caQTL validation

Analysis 2.9: Histone QTLs
- H3K27ac QTLs (enhancer activity)
- H3K4me3 QTLs (promoter activity)
- H3K4me1 QTLs (enhancer marking)

Analysis 2.10: 3D Contact QTLs
- Hi-C QTLs (loop strength)
- Effect on TAD boundaries
- CTCF binding variants
```

### 4.3 Cross-Ancestry Evaluation (WP3)

#### 3.3.1 Multi-Ancestry QTL Analysis

**Objective:** Evaluate whether AlphaGenome predictions generalize across human populations.

**Key Hypothesis:** Sequence-based models should perform equally across ancestries if trained adequately.

**Analyses:**

```
Analysis 3.1: Ancestry-Stratified Performance
Datasets:
- GTEx (primarily European)
- GENOA (African American)
- MESA (Hispanic, African American, European, Chinese)
- Biobank Japan (East Asian)
- FinnGen (Finnish)

Metrics per ancestry:
- eQTL sign prediction accuracy
- Effect size correlation
- Rare variant prediction

Analysis 3.2: Ancestry-Specific Regulatory Variants
- Variants at high frequency in one ancestry but rare/absent in others
- Population-specific enhancers
- Adaptive regulatory evolution

Analysis 3.3: Fine-Mapping Transferability
- Compare credible sets across ancestries
- Trans-ancestry fine-mapping validation
- GWAS loci with different causal variants by ancestry
```

**Critical Test: African Ancestry Performance**

```
Rationale: African populations have:
1. Greater genetic diversity
2. Shorter LD blocks (better fine-mapping)
3. Under-represented in training data

Specific Tests:
- AFR-specific eQTLs not found in EUR
- Variants in AFR-specific regulatory elements
- Performance on AFR GWAS fine-mapping
```

#### 3.3.2 Population-Specific Regulatory Variation

**Analyses:**

```
Analysis 3.4: Regulatory Element Evolution
- Compare predictions at evolutionarily divergent regulatory elements
- Test on regions with accelerated evolution (HARs)
- Archaic introgression variants (Neanderthal/Denisovan)

Analysis 3.5: Population Genetics Metrics Integration
- Correlation with selection signatures (iHS, XP-EHH)
- Prediction of deleterious regulatory variants (LINSIGHT, Eigen)
- Constraint metrics for regulatory elements (CDTS)
```

### 4.4 Disease Applications (WP4)

#### 3.4.1 Mendelian Disease Variants

**Objective:** Evaluate AlphaGenome for clinical variant interpretation.

**Analyses:**

```
Analysis 4.1: ClinVar Pathogenic vs. Benign Classification
- Non-coding pathogenic variants
- VUS (Variants of Uncertain Significance) prioritization
- Comparison with existing tools (CADD, DANN, Eigen)

Analysis 4.2: De Novo Mutation Prioritization
- DDD study variants
- SPARK autism cohort
- Congenital heart disease cohorts

Analysis 4.3: Specific Disease Categories
- Hemoglobinopathies (extensive regulatory variant data)
- Cancer predisposition (BRCA1/2 regulatory variants)
- Neurological disorders (brain-specific regulation)
- Immunodeficiencies (immune-specific enhancers)

Analysis 4.4: Diagnostic Yield Improvement
- Unsolved rare disease cases
- Reanalysis of negative exomes
- Integration with phenotype-driven prioritization
```

**Validation Framework:**

```
Gold Standard Datasets:
1. Literature-curated regulatory pathogenic variants (n≈500)
2. HGMD regulatory mutation class
3. Functional validation studies

Classification Task:
- Pathogenic regulatory vs. benign polymorphisms
- Metrics: auROC, auPRC, sensitivity at 90% specificity

Comparison Models:
- CADD v1.7
- DANN
- Eigen
- ReMM
- NCBoost
- regBase
```

#### 3.4.2 Complex Trait Fine-Mapping

**Objective:** Assess utility for GWAS variant prioritization.

**Analyses:**

```
Analysis 4.5: GWAS Credible Set Prioritization
- Open Targets Genetics credible sets
- Compare AlphaGenome ranking to:
  - Statistical fine-mapping (SuSiE, FINEMAP)
  - Functional annotations (PolyFun)
  - Existing variant effect predictors

Analysis 4.6: Trait-Specific Evaluation
- Metabolic traits (extensive fine-mapping)
- Autoimmune diseases (cell-type specific)
- Neuropsychiatric traits (brain enhancers)
- Cancer susceptibility loci

Analysis 4.7: Colocalization Enhancement
- eQTL-GWAS colocalization (coloc)
- Use AlphaGenome to prioritize colocalizing variants
- Multi-tissue colocalization
```

**Novel Analysis: Polygenic Score Enhancement**

```
Hypothesis: AlphaGenome scores can improve PRS by:
1. Weighting variants by predicted functional impact
2. Identifying causal variants in LD blocks
3. Improving cross-ancestry PRS transferability

Evaluation:
- UK Biobank held-out samples
- Comparison: standard PRS vs. AlphaGenome-weighted PRS
- Traits: height, BMI, lipids, blood pressure, T2D, CAD
```

### 4.5 Model Comparison & Benchmarking (WP5)

#### 3.5.1 Systematic Model Comparison

**Models to Compare:**

| Model | Type | Input Length | Outputs |
|-------|------|--------------|---------|
| AlphaGenome | Multi-modal | 1 Mb | 7,000+ tracks |
| Borzoi | Multi-modal | 500 kb | 7,000+ tracks |
| Enformer | Multi-modal | 200 kb | 5,000+ tracks |
| Evo 2 | Foundation | 1 Mb | Embeddings |
| Nucleotide Transformer v2 | Foundation | 12 kb | Embeddings |
| HyenaDNA | Foundation | 1 Mb | Embeddings |
| DNABERT-2 | Foundation | 512 bp | Embeddings |
| SpliceAI | Splicing | 10 kb | Splice scores |
| Pangolin | Splicing | 5 kb | Splice scores |
| ChromBPNet | Accessibility | 2 kb | Chromatin |
| ProCapNet | TSS | 2 kb | PRO-cap |

**Comparison Framework:**

```
Analysis 5.1: Head-to-Head Track Prediction
- Standardized test intervals
- Matched evaluation metrics
- Statistical significance testing (paired tests)

Analysis 5.2: Variant Effect Prediction Comparison
- Standardized variant sets
- Multiple scoring strategies per model
- Ensemble approaches

Analysis 5.3: Computational Efficiency
- Inference time per variant
- Memory requirements
- GPU vs. CPU performance
- Batch processing efficiency

Analysis 5.4: Ensemble and Integration
- Model ensembling strategies
- Feature combination approaches
- When does adding models help vs. hurt?
```

#### 3.5.2 Computational Reproducibility

**Analyses:**

```
Analysis 5.5: Reproducibility Assessment
- Run predictions multiple times
- Test across different hardware (GPU types)
- Version stability (API vs. local)

Analysis 5.6: Edge Case Evaluation
- Repetitive sequences
- GC-extreme regions
- Centromeric/telomeric regions
- Segmental duplications

Analysis 5.7: Input Sensitivity
- Effect of reference genome version (GRCh38 vs. T2T)
- N-masking effects
- Soft-masking effects
```

---

## Part V: Novel Evaluation Strategies

### 4.1 Cross-Species Validation

**Objective:** Evaluate transferability to non-human species relevant to preclinical research.

**Species:**
- Cynomolgus macaque (Macaca fascicularis) - key preclinical model
- Rhesus macaque (Macaca mulatta)
- Rat (Rattus norvegicus)
- Pig (Sus scrofa)

**Approach:**
```
Strategy 4.1.1: Orthologous Regulatory Element Prediction
- Identify conserved regulatory elements across species
- Compare predictions at orthologous positions
- Evaluate human-trained model on primate sequences

Strategy 4.1.2: Cross-Species QTL Validation
- Mouse eQTLs (BXD panel, Diversity Outbred)
- Rat eQTLs (HXB/BXH panel)
- Compare to human predictions at orthologous loci

Strategy 4.1.3: Evolutionary Constraint Validation
- PhyloP/phastCons scores vs. AlphaGenome predictions
- Species-specific regulatory elements
- Recently evolved human-specific enhancers
```

### 4.2 Structural Variant Prediction

**Objective:** Extend evaluation to larger genomic alterations.

```
Strategy 4.2.1: Deletion Effect Prediction
- gnomAD-SV deletions with expression effects
- Predict effect by "deleting" sequence (replacing with N's or flanking)
- Compare to measured expression changes

Strategy 4.2.2: Duplication Effect Prediction
- Copy number variant effects on expression
- Tandem vs. dispersed duplications

Strategy 4.2.3: Inversion Breakpoint Effects
- Effect on TAD boundaries
- Enhancer hijacking predictions
- Comparison with actual Hi-C data
```

### 4.3 Pharmacogenomic Applications

**Objective:** Evaluate utility for drug response prediction.

```
Strategy 4.3.1: Pharmacogene Regulatory Variants
- CYP450 enzyme regulation (CYP2D6, CYP3A4, CYP2C19)
- Transporter gene regulation (ABCB1, SLCO1B1)
- PharmGKB annotated variants

Strategy 4.3.2: Drug Response QTLs
- Statin response variants
- Warfarin sensitivity variants
- Chemotherapy response variants

Strategy 4.3.3: Drug-Induced Gene Expression
- L1000 connectivity map data
- Predict drug-induced expression changes
- Identify regulatory variants affecting drug response
```

### 4.4 Single-Cell Integration

**Objective:** Evaluate cell-type-specific predictions at single-cell resolution.

```
Strategy 4.4.1: Pseudo-Bulk Validation
- Human Cell Atlas datasets
- Cell-type-specific expression vs. predictions
- Rare cell population predictions

Strategy 4.4.2: Single-Cell ATAC Validation
- scATAC-seq from multiple tissues
- Cell-type-specific accessibility prediction
- Peak-level evaluation

Strategy 4.4.3: Single-Cell eQTL Validation
- Cell-type-specific eQTLs from sc-eQTL studies
- Dynamic eQTLs during differentiation
- Context-dependent genetic effects
```

### 4.5 In Silico Saturation Mutagenesis Validation

**Objective:** Systematically evaluate ISM predictions.

```
Strategy 4.5.1: MPRA Correlation
- Compare ISM scores to MPRA effect sizes
- Genome-wide MPRA datasets
- Cell-type-specific MPRAs

Strategy 4.5.2: DMS (Deep Mutational Scanning) Correlation
- Regulatory element DMS data
- Promoter DMS studies
- Enhancer saturation mutagenesis

Strategy 4.5.3: CRISPR Screen Validation
- CRISPRi tiling screens
- Base editing saturation screens
- Enhancer deletion screens
```

### 4.6 Temporal and Developmental Predictions

**Objective:** Evaluate predictions across developmental time.

```
Strategy 4.6.1: Developmental Stage Prediction
- Fetal vs. adult tissue expression
- Developmental enhancer activity
- Stage-specific regulatory variants

Strategy 4.6.2: Aging Effects
- Age-dependent gene expression changes
- Epigenetic clock-related variants
- Age-dependent eQTLs

Strategy 4.6.3: Differentiation Trajectories
- iPSC differentiation time courses
- Lineage-specific regulatory predictions
- Dynamic chromatin accessibility
```

---

## Part VI: Statistical Analysis Plan

### 5.1 Primary Metrics

| Metric | Use Case | Interpretation |
|--------|----------|----------------|
| Pearson r | Track prediction | Linear correlation with observed |
| Spearman ρ | Variant effect sizes | Rank correlation with QTL effects |
| auROC | Binary classification | Discrimination ability |
| auPRC | Imbalanced classification | Precision-recall for rare positives |
| MSE/RMSE | Continuous prediction | Prediction error magnitude |
| R² | Variance explained | Proportion of variance captured |

### 5.2 Statistical Testing

```
Comparison Tests:
- Paired t-test / Wilcoxon signed-rank: Compare models on same variants
- DeLong test: Compare auROC between models
- Bootstrap confidence intervals: Uncertainty quantification
- Permutation tests: Null distribution for enrichments

Multiple Testing Correction:
- Bonferroni for primary analyses
- FDR (Benjamini-Hochberg) for exploratory analyses
- Family-wise error rate control within evaluation categories
```

### 5.3 Power Analysis

```python
# Sample size requirements for detecting performance differences

def power_analysis_auroc(auroc1, auroc2, power=0.8, alpha=0.05):
    """
    Calculate required sample size to detect auROC difference
    """
    from scipy import stats
    import numpy as np
    
    # Effect size (Cohen's d approximation for auROC)
    effect_size = (auroc1 - auroc2) / 0.1  # Approximate SD
    
    # Two-sample comparison
    n = stats.power.TTestIndPower().solve_power(
        effect_size=effect_size,
        power=power,
        alpha=alpha
    )
    
    return int(np.ceil(n))

# Example: Detect 0.02 auROC difference
n_required = power_analysis_auroc(0.85, 0.83)
# Approximately 400 variants per group needed
```

### 5.4 Stratified Analysis Plan

All primary analyses will be stratified by:

1. **Variant characteristics:**
   - SNV vs. indel
   - Transition vs. transversion
   - Distance to TSS (bins)
   - MAF bins
   
2. **Genomic context:**
   - Promoter vs. enhancer vs. intronic vs. intergenic
   - GC content
   - Chromatin state
   - Conservation level

3. **Gene characteristics:**
   - Expression level
   - Tissue specificity
   - Gene length
   - Constraint (pLI, LOEUF)

4. **Population genetics:**
   - Ancestry
   - Allele frequency spectrum
   - Selection signals

---

## Part VII: Deliverables

### 6.1 Primary Manuscript

**Target Journal:** Nature Methods or Genome Biology

**Manuscript Structure:**
1. Introduction (AlphaGenome significance, need for independent evaluation)
2. Results
   - Track prediction performance on independent data
   - Cross-ancestry variant effect prediction
   - Clinical variant interpretation
   - Comparison with other models
   - Novel applications
3. Discussion (strengths, limitations, recommendations)
4. Methods (detailed protocols)

**Figures (planned):**
1. Study design and data overview
2. Track prediction benchmarks
3. eQTL/sQTL prediction performance by ancestry
4. Clinical variant classification
5. Model comparison across tasks
6. Novel application results
7. Computational benchmarks

### 6.2 Supplementary Resources

1. **Supplementary Data:**
   - Complete benchmark results tables
   - Per-variant predictions for all test sets
   - Model comparison matrices

2. **Code Repository:**
   - All evaluation scripts (GitHub)
   - Standardized benchmark datasets
   - Reproducible analysis pipelines

3. **Interactive Resource:**
   - Web portal for exploring results
   - Variant lookup tool
   - Model comparison dashboard

### 6.3 Preprint and Data Sharing

- Preprint on bioRxiv upon completion
- All evaluation datasets deposited (where permitted)
- Code available under open-source license
- Model predictions for community use

---

## Part VIII: Risk Assessment and Mitigation

### 7.1 Potential Challenges

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| API access limitations | Medium | High | Request research access; prepare local deployment |
| Computational costs | Medium | Medium | Cloud credits; optimize batch processing |
| Dataset access delays | Medium | Medium | Start applications early; have backup datasets |
| Model updates during study | Low | High | Version lock; document any changes |
| Negative results | Medium | Low | Pre-register; frame as informative either way |

### 7.2 Quality Control

```
QC Checkpoints:
1. Week 4: Data acquisition complete, pipelines tested
2. Week 12: Primary analyses 50% complete, interim review
3. Week 20: All analyses complete, cross-WP integration
4. Week 28: Manuscript draft for internal review
5. Week 36: Final submission preparation

QC Measures:
- Duplicate analysis by independent team members
- Code review for all analysis scripts
- Blinded evaluation where possible
- External advisory review at key milestones
```

---

## Part IX: Budget Estimate

### 8.1 Computational Costs

| Resource | Estimated Cost | Duration |
|----------|---------------|----------|
| Cloud GPU (inference) | $15,000 | 6 months |
| Cloud storage | $2,000 | 12 months |
| API costs (if applicable) | $5,000 | 6 months |
| **Subtotal** | **$22,000** | |

### 8.2 Personnel (In-Kind)

| Role | FTE | Months | Value* |
|------|-----|--------|--------|
| Senior investigators (5) | 0.1 each | 12 | $50,000 |
| Postdocs/Scientists (10) | 0.2 each | 12 | $200,000 |
| Students/Trainees (5) | 0.3 each | 12 | $75,000 |
| **Subtotal** | | | **$325,000** |

*Estimated salary + benefits contribution

### 8.3 Other Costs

| Item | Cost |
|------|------|
| Data access fees | $3,000 |
| Publication costs | $5,000 |
| Conference presentation | $5,000 |
| **Subtotal** | **$13,000** |

### 8.4 Total Estimated Budget

| Category | Amount |
|----------|--------|
| Computational | $22,000 |
| Personnel (in-kind) | $325,000 |
| Other | $13,000 |
| **Total** | **$360,000** |

---

## Appendices

### Appendix A: Detailed Evaluation Protocols

#### A.1 eQTL Evaluation Protocol

```python
"""
Standardized eQTL evaluation protocol
"""

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr, pearsonr

class eQTLEvaluator:
    def __init__(self, eqtl_data, alphagenome_scores):
        """
        eqtl_data: DataFrame with columns:
            - variant_id
            - gene_id
            - effect_size (beta)
            - effect_sign (+1 or -1)
            - p_value
            - tissue
            - ancestry
        
        alphagenome_scores: DataFrame with columns:
            - variant_id
            - gene_id
            - ag_score (AlphaGenome variant effect score)
        """
        self.eqtl = eqtl_data
        self.scores = alphagenome_scores
        self.merged = self._merge_data()
    
    def _merge_data(self):
        return pd.merge(
            self.eqtl, 
            self.scores,
            on=['variant_id', 'gene_id'],
            how='inner'
        )
    
    def evaluate_sign_prediction(self, stratify_by=None):
        """
        Evaluate sign (direction) prediction
        """
        results = {}
        
        if stratify_by is None:
            # Overall evaluation
            y_true = (self.merged['effect_sign'] > 0).astype(int)
            y_score = self.merged['ag_score']
            
            results['overall'] = {
                'auroc': roc_auc_score(y_true, y_score),
                'n_variants': len(y_true)
            }
        else:
            # Stratified evaluation
            for group, data in self.merged.groupby(stratify_by):
                y_true = (data['effect_sign'] > 0).astype(int)
                y_score = data['ag_score']
                
                if len(y_true) >= 50:  # Minimum sample size
                    results[group] = {
                        'auroc': roc_auc_score(y_true, y_score),
                        'n_variants': len(y_true)
                    }
        
        return results
    
    def evaluate_effect_size(self, stratify_by=None):
        """
        Evaluate effect size correlation
        """
        results = {}
        
        if stratify_by is None:
            rho, pval = spearmanr(
                self.merged['effect_size'].abs(),
                self.merged['ag_score'].abs()
            )
            results['overall'] = {
                'spearman_rho': rho,
                'p_value': pval,
                'n_variants': len(self.merged)
            }
        else:
            for group, data in self.merged.groupby(stratify_by):
                if len(data) >= 50:
                    rho, pval = spearmanr(
                        data['effect_size'].abs(),
                        data['ag_score'].abs()
                    )
                    results[group] = {
                        'spearman_rho': rho,
                        'p_value': pval,
                        'n_variants': len(data)
                    }
        
        return results
    
    def evaluate_causality(self, pip_threshold=0.5):
        """
        Evaluate causal variant prioritization
        Requires 'pip' column in eqtl_data
        """
        if 'pip' not in self.merged.columns:
            raise ValueError("PIP values required for causality evaluation")
        
        # Binary classification: causal (PIP > threshold) vs. non-causal
        y_true = (self.merged['pip'] > pip_threshold).astype(int)
        y_score = self.merged['ag_score'].abs()
        
        return {
            'auroc': roc_auc_score(y_true, y_score),
            'auprc': average_precision_score(y_true, y_score),
            'n_causal': y_true.sum(),
            'n_total': len(y_true)
        }
    
    def full_evaluation(self):
        """
        Run complete evaluation battery
        """
        return {
            'sign_overall': self.evaluate_sign_prediction(),
            'sign_by_ancestry': self.evaluate_sign_prediction('ancestry'),
            'sign_by_tissue': self.evaluate_sign_prediction('tissue'),
            'effect_overall': self.evaluate_effect_size(),
            'effect_by_ancestry': self.evaluate_effect_size('ancestry'),
        }
```

#### A.2 Model Comparison Protocol

```python
"""
Standardized model comparison framework
"""

class ModelComparator:
    def __init__(self, models_dict, test_data):
        """
        models_dict: {model_name: scorer_function}
        test_data: standardized test set
        """
        self.models = models_dict
        self.test_data = test_data
        self.predictions = {}
    
    def run_all_models(self):
        """Generate predictions from all models"""
        for name, scorer in self.models.items():
            print(f"Running {name}...")
            self.predictions[name] = scorer(self.test_data)
    
    def compare_auroc(self, y_true):
        """Compare auROC across models"""
        results = {}
        for name, preds in self.predictions.items():
            results[name] = roc_auc_score(y_true, preds)
        return pd.Series(results).sort_values(ascending=False)
    
    def statistical_comparison(self, y_true, model1, model2):
        """
        DeLong test for comparing two auROCs
        """
        from scipy import stats
        
        auc1 = roc_auc_score(y_true, self.predictions[model1])
        auc2 = roc_auc_score(y_true, self.predictions[model2])
        
        # Bootstrap confidence interval for difference
        n_bootstrap = 1000
        diffs = []
        n = len(y_true)
        
        for _ in range(n_bootstrap):
            idx = np.random.choice(n, n, replace=True)
            y_boot = y_true[idx]
            p1_boot = self.predictions[model1][idx]
            p2_boot = self.predictions[model2][idx]
            
            try:
                auc1_boot = roc_auc_score(y_boot, p1_boot)
                auc2_boot = roc_auc_score(y_boot, p2_boot)
                diffs.append(auc1_boot - auc2_boot)
            except:
                continue
        
        ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
        p_value = np.mean(np.array(diffs) < 0) * 2  # Two-sided
        
        return {
            'auc_diff': auc1 - auc2,
            'ci_95': (ci_low, ci_high),
            'p_value': min(p_value, 1 - p_value) * 2
        }
```

### Appendix B: Dataset Access Procedures

| Dataset | Application Process | Timeline |
|---------|-------------------|----------|
| GTEx v10 | dbGaP application | 2-4 weeks |
| UK Biobank | Research application | 4-8 weeks |
| ENCODE4 | Direct download | Immediate |
| ClinVar | Direct download | Immediate |
| FinnGen | Research agreement | 2-4 weeks |
| MESA | dbGaP application | 2-4 weeks |
| Biobank Japan | Collaboration request | 4-8 weeks |

### Appendix C: Computational Environment Specification

```yaml
# environment.yml
name: alphagenome_eval
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.11
  - numpy>=2.0
  - pandas>=2.0
  - scipy>=1.10
  - scikit-learn>=1.3
  - pytorch>=2.0
  - tensorflow>=2.15
  - jax>=0.4
  - flax>=0.8
  
  # Bioinformatics
  - pysam>=0.22
  - biopython>=1.82
  - pybedtools>=0.9
  - pyranges>=0.0.120
  
  # Genomics ML
  - kipoiseq>=0.7
  - selene-sdk>=0.5
  
  # Visualization
  - matplotlib>=3.8
  - seaborn>=0.13
  - plotly>=5.18
  
  # Statistics
  - statsmodels>=0.14
  - lifelines>=0.27
  
  # Utilities
  - tqdm>=4.66
  - joblib>=1.3
  - h5py>=3.10
  - zarr>=2.16
```

### Appendix D: Pre-Registration Template

```markdown
# Pre-Registration: Independent Evaluation of AlphaGenome

## Study Information
- Title: Comprehensive Independent Evaluation of AlphaGenome for 
         Regulatory Variant Effect Prediction
- Authors: [Consortium members]
- Date: [Registration date]
- Registry: OSF Registries / AsPredicted

## Hypotheses
1. AlphaGenome will achieve comparable performance to published results 
   on independent test sets
2. Performance will vary across ancestries, with potentially lower 
   accuracy in non-European populations
3. AlphaGenome will outperform single-modality models on tasks 
   requiring integration of multiple genomic features

## Primary Outcomes
1. auROC for eQTL sign prediction
2. Spearman correlation for eQTL effect size
3. auPRC for clinical variant classification

## Analysis Plan
[Detailed as specified in protocol]

## Data Collection
- All datasets specified in Section II
- No additional data collection after registration

## Sample Size
- Minimum 10,000 variants per evaluation category
- Power analysis as specified in Section V.3
```

---

## Document Control

| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 0.1 | Feb 2026 | Initial draft | [Lead author] |
| 0.2 | - | WP lead review | [WP leads] |
| 1.0 | - | Final protocol | Consortium |

---

*This protocol is intended for internal consortium use. Please contact the coordinating center before external distribution.*
