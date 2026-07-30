# uc-regulatory-variants
# Dark Matter Regulatory Variants in Ulcerative Colitis

**Vaishnavi Nagesh** · ML Engineer, Bioinformatics · [vnageshbio](https://github.com/vnageshbio)

> Applying deep learning sequence models (ChromBPNet) to identify which genome-wide significant UC GWAS variants are likely disrupting regulatory elements in disease-relevant cell types — and whether those effects are cell-type specific.

---

## Table of Contents

1. [Biological Question](#1-biological-question)
2. [Background](#2-background)
3. [Data Sources](#3-data-sources)
4. [Methods](#4-methods)
   - [M1 — Literature & Foundations](#m1--literature--foundations)
   - [M2 — Variant Filtering & Intersection](#m2--variant-filtering--intersection)
   - [M3 — ChromBPNet Training](#m3--chrombpnet-training)
5. [Results](#5-results)
6. [Computational Environment](#6-computational-environment)
7. [Limitations](#7-limitations)
8. [Next Steps](#8-next-steps)
9. [Repository Structure](#9-repository-structure)
10. [References](#10-references)

---

## 1. Biological Question

Genome-wide association studies (GWAS) for ulcerative colitis (UC) have identified hundreds of significant loci, yet the vast majority of these variants fall outside protein-coding sequences — in the non-coding "dark matter" of the genome. This raises a fundamental question:

**Why do most UC GWAS variants fall in non-coding DNA, and which ones are likely disrupting regulatory elements in disease-relevant cell types?**

The leading hypothesis is that these non-coding variants modulate gene expression by altering transcription factor (TF) binding sites or chromatin accessibility within regulatory elements such as enhancers and promoters. If true, the functional effect should be observable in the cell types most relevant to UC pathophysiology — specifically, immune cells mediating the inflammatory response (CD4+ T cells) and the intestinal epithelium (sigmoid colon).

This project tests that hypothesis computationally: using chromatin accessibility profiles and a trained sequence model (ChromBPNet) to predict which UC variants disrupt open chromatin, and whether those disruptions are cell-type specific.

---

## 2. Background

### Why non-coding variants are hard to interpret

Unlike coding variants, which can be interpreted through known gene function, non-coding variants have no simple functional annotation. A SNP 40kb upstream of a gene may be in an active enhancer in one cell type and entirely inert in another. Standard annotation tools (e.g., VEP) flag proximity to regulatory features but cannot model the sequence-level effect on TF binding or chromatin state.

### ChromBPNet as the modeling approach

ChromBPNet (Pampari et al., 2025) is a deep learning model trained on ATAC-seq data to predict chromatin accessibility from DNA sequence. Given a genomic sequence as input, it outputs a predicted accessibility profile. Critically, it can be used for *in silico* mutagenesis: by scoring the reference and alternate alleles of a variant, the model predicts the allelic effect on chromatin accessibility. This provides a sequence-resolved, cell-type-specific functional score for each variant.

### Cell type selection rationale

UC is an inflammatory bowel disease driven by dysregulation at the intersection of the immune system and intestinal mucosa. Two cell types are central:

- **CD4+ T cells** — orchestrate the adaptive immune response in UC; enriched for UC GWAS heritability in prior studies (Calderon et al., 2019)
- **Sigmoid colon** — the primary site of UC inflammation; open chromatin in colon epithelium is the tissue-proximal regulatory context

Testing both cell types allows a direct comparison: variants that affect CD4+ T cell chromatin but not colon (or vice versa) are strong candidates for cell-type-specific regulatory mechanisms. Variants affecting both are candidates for shared regulatory logic across the immune-epithelial axis.

---

## 3. Data Sources

| Resource | Accession | Description |
|---|---|---|
| GWAS summary statistics | GCST90473823 | UK Biobank UC GWAS; 458,440 European-ancestry participants; ICD10 K51 |
| CD4+ T cell ATAC-seq | ENCSR841LHT | ENCODE; paired-end; hg38; peaks (IDR-filtered) + BAM |
| Sigmoid colon ATAC-seq | ENCSR355SGJ | ENCODE; paired-end; hg38; peaks (IDR-filtered) + BAM |
| Reference genome | hg38 | GENCODE v44 exon annotations; ENCODE blacklist v2 |

All data are publicly available. GWAS summary statistics were downloaded from the NHGRI-EBI GWAS Catalog. ATAC-seq data were downloaded from the ENCODE portal (encodeproject.org).

---

## 4. Methods

### M1 — Literature & Foundations

Prior to analysis, the following foundational literature was reviewed:

- **Thurman et al. (2012)** — *Nature* — The accessible chromatin landscape of the human genome. Establishes DNase hypersensitivity as a marker of regulatory activity across 125 cell types.
- **Buenrostro et al. (2013)** — *Nature Methods* — ATAC-seq: Assay for Transposase-Accessible Chromatin using sequencing. The original ATAC-seq paper; defines the methodology underlying all chromatin accessibility data used here.
- **Kundaje et al. (2015)** — *Nature* — Integrative analysis of 111 reference human epigenomes (Roadmap Epigenomics). Provides the chromatin state framework and establishes cell-type specificity of regulatory elements.
- **Pampari et al. (2025)** — *bioRxiv* — ChromBPNet: bias factorized, baseresolution deep learning models of chromatin accessibility reveal cis-regulatory sequences in human cells. The primary modeling approach used in this project.

Conceptual grounding established in: ATAC-seq data processing, chromatin accessibility as a proxy for regulatory activity, TF motif analysis (JASPAR/HOCOMOCO), deep learning for sequence-to-function prediction, and GWAS variant interpretation.

---

### M2 — Variant Filtering & Intersection

**Goal:** Identify genome-wide significant UC GWAS variants that fall within open chromatin in CD4+ T cells and/or sigmoid colon.

#### 2.1 GWAS Variant Filtering

Starting from GCST90473823 (harmonized summary statistics, hg38):

1. **Genome-wide significance filter:** Retained variants with p < 5×10⁻⁸
   - Input: full summary statistics
   - Output: 4,081 significant variants

2. **Standard chromosome filter:** Retained chr1–chr22, chrX only (removed chrY, chrM, unplaced scaffolds)

3. **Non-coding filter:** Removed variants overlapping GENCODE v44 protein-coding exons (bedtools intersect -v)
   - Output: **3,797 non-coding variants**

4. **Blacklist filter:** Removed variants overlapping ENCODE hg38 blacklist v2 (regions of anomalous signal in ATAC-seq/ChIP-seq)
   - Final filtered set: 3,797 variants (minimal blacklist overlap in this dataset)

#### 2.2 ATAC-seq Peak Processing

For each cell type (CD4+ T cells, sigmoid colon):

1. Downloaded IDR-filtered consensus peak BED files from ENCODE
2. No additional merging required — ENCODE IDR peaks are already consensus-called across replicates
3. Peaks represent open chromatin regions; width varies (typically 150–500bp, summit-centered)

#### 2.3 Variant–Peak Intersection

Intersected filtered GWAS variants with ATAC-seq peaks per cell type using `bedtools intersect`:

```bash
bedtools intersect -a uc_variants_noncoding.bed -b cd4_peaks.bed -wa -u > cd4_variants.bed
bedtools intersect -a uc_variants_noncoding.bed -b colon_peaks.bed -wa -u > colon_variants.bed
```

#### 2.4 Results Summary

| Filter step | Variants |
|---|---|
| Genome-wide significant (p < 5×10⁻⁸) | 4,081 |
| Non-coding (exon-filtered) | 3,797 |
| In CD4+ T cell open chromatin | **174** |
| In sigmoid colon open chromatin | **77** |
| In both cell types | **50** |

The 50 variants overlapping open chromatin in both CD4+ T cells and sigmoid colon represent the highest-priority candidates for shared regulatory disruption across the immune-epithelial axis. The 124 CD4-specific and 27 colon-specific variants are candidates for cell-type-restricted regulatory effects.

All intermediate and final BED files are in `results/`.

---

### M3 — ChromBPNet Training

**Goal:** Train cell-type-specific ChromBPNet models on CD4+ T cell and sigmoid colon ATAC-seq data, enabling in silico prediction of chromatin accessibility from sequence.

#### 3.1 Overview of ChromBPNet Architecture

ChromBPNet is a bias-factorized deep learning model. It separates the sequence-specific accessibility signal from Tn5 transposase insertion bias (which confounds ATAC-seq signal at certain sequence motifs). The model consists of:

- **Bias model:** Trained first on GC-matched genomic negatives to learn Tn5 sequence bias
- **ChromBPNet main model:** Trained on ATAC-seq signal with the bias model held fixed, learning the residual sequence-specific accessibility signal
- **Architecture:** Dilated convolutional network (BPNet-style); predicts both profile shape and total count

#### 3.2 Input Preparation

For each cell type:

- **Positive regions:** ATAC-seq peaks (IDR-filtered), extended to 2,114bp centered on summit
- **Negative regions:** GC-matched genomic regions with no ATAC-seq signal (used to train bias model)
- **Sequence input:** hg38 FASTA, one-hot encoded (A/C/G/T)
- **Signal input:** ATAC-seq fragment coverage (from ENCODE BAM files, converted to bigWig)

#### 3.3 Training Procedure

Both models (CD4+ T cell, sigmoid colon) were trained using the ChromBPNet CLI via micromamba Python 3.8 environment on Google Colab Pro (NVIDIA A100, 40GB).

**Bias model training:**

```bash
LD_LIBRARY_PATH=/root/.local/share/mamba/envs/chrombpnet_v2/lib:$LD_LIBRARY_PATH \
./bin/micromamba run -n chrombpnet_v2 \
chrombpnet bias pipeline \
  -ibam $BAM \
  -d ATAC \
  -g hg38.fa \
  -c hg38.chrom.sizes \
  -p $PEAKS \
  -n $NEGATIVES \
  -fl $FOLD \
  -o $OUTPUT_DIR/bias
```

Early stopping applied; bias models converged in ~14 epochs.

**ChromBPNet main model training:**

```bash
LD_LIBRARY_PATH=/root/.local/share/mamba/envs/chrombpnet_v2/lib:$LD_LIBRARY_PATH \
./bin/micromamba run -n chrombpnet_v2 \
chrombpnet pipeline \
  -ibam $BAM \
  -d ATAC \
  -g hg38.fa \
  -c hg38.chrom.sizes \
  -p $PEAKS \
  -n $NEGATIVES \
  -fl $FOLD \
  -bl $BIAS_MODEL \
  -o $OUTPUT_DIR/chrombpnet
```

#### 3.4 Trained Model Outputs

For each cell type, the following model files were produced and saved to Google Drive:

| File | Description |
|---|---|
| `chrombpnet.h5` | Full ChromBPNet model (with bias) |
| `chrombpnet_nobias.h5` | Sequence-only model (bias factored out) — used for variant scoring |
| `bias_model_scaled.h5` | Trained Tn5 bias model |

Training loss curves and GC-matched negative plots are in `figures/`.

#### 3.5 Environment

ChromBPNet requires TensorFlow 2.8.0 and Python 3.8. On Google Colab Pro (which ships with newer CUDA/cuDNN libraries), the environment requires an explicit `LD_LIBRARY_PATH` override to resolve shared library conflicts:

```bash
# Full micromamba invocation required for every ChromBPNet command on Colab
LD_LIBRARY_PATH=/root/.local/share/mamba/envs/chrombpnet_v2/lib:$LD_LIBRARY_PATH \
./bin/micromamba run -n chrombpnet_v2 chrombpnet [command]
```

See `environment/` for full micromamba setup instructions and the `chrombpnet_v2` environment YAML.

---

## 5. Results

### Variant landscape

Of 4,081 genome-wide significant UC variants, 3,797 (93%) are non-coding — consistent with the established pattern that most GWAS signal lies outside protein-coding sequence. The small fraction in open chromatin (174 CD4+, 77 colon) reflects the cell-type specificity of ATAC-seq peaks: most regulatory elements are not accessible in all cell types, and only a subset of GWAS loci fall within accessible chromatin in any given cell type.

The 50 variants shared across both CD4+ T cells and sigmoid colon are particularly noteworthy: these loci show accessible chromatin in both an immune and an epithelial context, suggesting regulatory activity relevant to both arms of UC pathogenesis.

### ChromBPNet training

Both cell-type models trained successfully. Bias models converged in approximately 14 epochs with early stopping. Training loss curves show expected convergence behavior (available in `figures/training_loss/`). GC-matched negative plots confirm that the negative set is appropriately matched to the GC content distribution of ATAC-seq peaks, a prerequisite for unbiased bias model training.

### Pending: variant effect scores

Variant effect scoring (snp_score) and TF motif analysis (TF-MoDISco) are pending and will complete the functional interpretation of the 174 + 77 variant sets. Results will be added here upon completion.

---

## 6. Computational Environment

| Component | Specification |
|---|---|
| Hardware | Google Colab Pro, NVIDIA A100 40GB |
| OS | Ubuntu 20.04 (Colab) |
| Python | 3.8.x (micromamba environment: `chrombpnet_v2`) |
| ChromBPNet | v0.1.x |
| TensorFlow | 2.8.0 |
| Package manager | micromamba (via `./bin/micromamba`) |
| Genome reference | hg38 (GRCh38) |
| Genome tools | bedtools 2.30, samtools 1.17, pyBigWig |

**Known issue — Colab GPU compatibility:** Google Colab Pro ships with CUDA/cuDNN versions that conflict with the TF 2.8.0 environment required by ChromBPNet. The following prefix must be prepended to every `chrombpnet` command:

```bash
LD_LIBRARY_PATH=/root/.local/share/mamba/envs/chrombpnet_v2/lib:$LD_LIBRARY_PATH \
./bin/micromamba run -n chrombpnet_v2 chrombpnet [command]
```

Full setup instructions are in `environment/README.md`.

---

## 7. Limitations

**Bulk ATAC-seq resolution.** ENCODE ATAC-seq data represent bulk populations of CD4+ T cells and sigmoid colon cells, averaging across heterogeneous cell states (e.g., naive vs. activated T cells, crypt vs. villus colonocytes). Single-cell ATAC-seq (scATAC-seq) would provide finer resolution, but is not used here.

**European-ancestry GWAS.** GCST90473823 is derived from UK Biobank participants of European ancestry. GWAS effect sizes and significant loci may not generalize to other ancestries. The set of variants identified here may miss UC-associated regulatory variants that are only detected in multi-ancestry analyses.

**In silico variant scoring.** ChromBPNet predictions are sequence-based and have not been validated experimentally (e.g., by MPRA or CRISPRi) in this project. Predicted allelic effects should be interpreted as prioritization scores, not functional proofs.

**Linkage disequilibrium.** GWAS variants are association signals; the causal variant at each locus may be in LD with the reported lead SNP. This analysis uses lead SNPs; a LD-expansion step would increase sensitivity.

**Pending analyses.** Variant effect scoring (snp_score) and TF-MoDISco motif analysis are not yet complete. The biological interpretation of which variants disrupt which TF binding sites remains to be determined.

---

## 8. Next Steps

### Immediate (M3 completion)

- [ ] Run `chrombpnet snp_score` on all 174 CD4+ variants using the trained CD4+ ChromBPNet model
- [ ] Run `chrombpnet snp_score` on all 77 colon variants using the trained colon ChromBPNet model
- [ ] Compute log₂(ALT/REF) allelic imbalance scores for each variant in each cell type
- [ ] Rank variants by predicted effect size; identify top candidates
- [ ] Cross-cell-type comparison: which variants show discordant effects between CD4+ and colon?

### Downstream (M4)

- [ ] Run DeepLIFT contribution score analysis on top-ranked variants
- [ ] Run TF-MoDISco on contribution scores to identify disrupted TF binding motifs
- [ ] Match discovered motifs against JASPAR/HOCOMOCO
- [ ] Annotate each high-impact variant: which TF motif is disrupted, in which cell type

### Stretch

- [ ] Compare ChromBPNet variant scores to eQTL data from GTEx (colon, whole blood) as external validation
- [ ] Explore JEPA-DNA (Larey et al., 2026) as an alternative embedding approach for variant scoring
- [ ] Explore ProteinJEPA (Ofer et al., 2026) representations for TF proteins whose motifs are disrupted — a cross-modal regulatory mechanism question

---

## 9. Repository Structure

```
uc-regulatory-variants/
├── README.md                    # This document
├── data/
│   └── raw/                     # GWAS summary stats (not committed; see accessions above)
├── results/
│   ├── uc_variants_filtered.bed # 3,797 non-coding variants
│   ├── cd4_variants.bed         # 174 variants in CD4+ open chromatin
│   ├── colon_variants.bed       # 77 variants in colon open chromatin
│   └── shared_variants.bed      # 50 variants in both
├── models/
│   └── [stored on Google Drive — see environment/README.md for access]
├── figures/
│   ├── training_loss/           # ChromBPNet training curves
│   └── gc_matched_negatives/    # GC content distribution plots
├── notebooks/
│   ├── 01_variant_filtering.ipynb
│   ├── 02_atac_intersection.ipynb
│   └── 03_chrombpnet_training.ipynb
├── environment/
│   ├── README.md                # Setup instructions
│   └── chrombpnet_v2.yml        # micromamba environment spec
└── scripts/
    ├── filter_variants.sh
    ├── intersect_peaks.sh
    └── run_chrombpnet.sh
```

---

## 10. References

Buenrostro JD, et al. (2013). Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. *Nature Methods*, 10(12), 1213–1218.

Calderon D, et al. (2019). Landscape of stimulation-responsive chromatin across diverse human immune cells. *Nature Genetics*, 51, 1494–1505.

Larey A, et al. (2026). JEPA-DNA: Grounding Genomic Foundation Models through Joint-Embedding Predictive Architectures. *arXiv:2602.17162*.

Ofer D, Linial M, Shahaf D. (2026). ProteinJEPA: Latent prediction complements protein language models. *arXiv:2605.07554*.

Pampari A, et al. (2025). ChromBPNet: bias factorized, baseresolution deep learning models of chromatin accessibility reveal cis-regulatory sequences in human cells. *bioRxiv*.

Roadmap Epigenomics Consortium, Kundaje A, et al. (2015). Integrative analysis of 111 reference human epigenomes. *Nature*, 518, 317–330.

Thurman RE, et al. (2012). The accessible chromatin landscape of the human genome. *Nature*, 489, 75–82.

---

*Last updated: July 2026. Models trained on Google Colab Pro A100. Variant scoring and motif analysis in progress.*
