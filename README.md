# uc-regulatory-variants
dark-matter-uc
# Dark Matter Regulatory Variants in Ulcerative Colitis

## Biological Question
Why do most UC GWAS variants fall in non-coding DNA, and which ones 
are likely disrupting regulatory elements in disease-relevant cell types?

## Summary of Findings
- 4,081 significant UC variants (p<5e-8) from UK Biobank GWAS
- 3,797 non-coding (exon-filtered)
- 174 overlapping CD4+ T cell open chromatin
- 77 overlapping sigmoid colon open chromatin  
- 50 overlapping both — highest-priority regulatory candidates

## Data
| Resource | Accession | Description |
|---|---|---|
| GWAS | GCST90473823 | UK Biobank UC, 458,440 European, ICD10 K51 |
| CD4 ATAC-seq | ENCSR841LHT | ENCODE CD4+ T cell peaks + BAM |
| Colon ATAC-seq | ENCSR355SGJ | ENCODE sigmoid colon peaks + BAM |
| Reference | hg38 | GENCODE v44 exons, blacklist v2 |

## Methods
### M1 — Foundations
Literature reviewed: Thurman 2012, Buenrostro 2013, Kundaje 2015.
Conceptual grounding in ATAC-seq, chromatin accessibility, and 
deep learning for regulatory genomics.

### M2 — Variant Intersection
1. Downloaded harmonised GWAS summary statistics
2. Filtered: p<5e-8, standard chromosomes, non-coding (exon intersection)
3. Applied hg38 blacklist filter
4. Intersected with ENCODE ATAC-seq peaks per cell type
5. Results in results/ directory

### M3 — ChromBPNet Training
- Bias model trained per cell type (early stopping ~14 epochs)
- ChromBPNet main model trained per cell type
- Models saved: chrombpnet.h5, chrombpnet_nobias.h5, bias_model_scaled.h5
- Environment: Google Colab Pro A100, micromamba Python 3.8, TF 2.8.0

## Results
[Figures: GC-matched negatives plot, training loss curves, variant counts]

## Limitations
- Bulk ATAC-seq — averages across mixed cell states
- European-only GWAS — population generalizability limited
- Variant scoring (snp_score) and TF-MoDISco pending

## Next Steps
- Run chrombpnet snp_score on 174 CD4 + 77 colon variants
- Compute log2(ALT/REF) allelic imbalance scores
- TF-MoDISco contribution scores for motif discovery

## Environment
See environment/ for micromamba setup instructions including
LD_LIBRARY_PATH fix for Colab GPU compatibility.
