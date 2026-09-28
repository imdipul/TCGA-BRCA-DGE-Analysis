# Convergent and Tissue-Specific Transcriptomic Programs in Breast and Lung Adenocarcinoma

Reproducible R pipeline for a comparative differential gene expression (DGE), functional enrichment, multi-variable survival, and protein-interaction network analysis of TCGA-BRCA (breast invasive carcinoma) and TCGA-LUAD (lung adenocarcinoma), using public TCGA RNA-sequencing data.

This repository accompanies the manuscript *"Cross-cancer transcriptomic comparison of breast and lung adenocarcinoma: a comparative RNA-sequencing, survival, and interaction-network analysis using TCGA data"*.

## Repository Contents

| File / Folder | Description |
|---|---|
| `analysis.R` / `analysis.Rmd` | Complete, end-to-end reproducible analysis pipeline from raw count acquisition to functional clustering. |
| `validation_full_cohort.R` | Cloud-optimized verification script executing the full, un-subsampled TCGA-BRCA pipeline on complete cohorts. |
| `TCGA_BRCA_DESeq2_results.csv` | Full processed BRCA DESeq2 output matrix (53,081 genes). |
| `TCGA_LUAD_DESeq2_results.csv` | Full processed LUAD DESeq2 output matrix (52,165 genes). |
| `TCGA_BRCA_GO_Enrichment_results.csv`| Full Gene Ontology Biological Process enrichment results for BRCA. |
| `GO_shared_pancancer.csv` | Baseline unified GO enrichment records for multi-cancer shared DEGs. |
| `GO_shared_up.csv` | GO Biological Process terms unique to shared *upregulated* genes. |
| `GO_shared_down.csv` | GO Biological Process terms unique to shared *downregulated* genes. |
| `GO_BRCA_specific.csv` | Mixed-direction GO enrichment terms specific to the BRCA cohort. |
| `GO_brca_only_up.csv` | GO Biological Process terms specific to BRCA *upregulated* genes. |
| `GO_brca_only_down.csv` | GO Biological Process terms specific to BRCA *downregulated* genes. |
| `GO_LUAD_specific.csv` | Mixed-direction GO enrichment terms specific to the LUAD cohort. |
| `GO_luad_only_up.csv` | GO Biological Process terms specific to LUAD *upregulated* genes. |
| `GO_luad_only_down.csv` | GO Biological Process terms specific to LUAD *downregulated* genes. |
| `results.rds` | Structured serialization R data object for prompt workspace restoration. |
| `figures/` | Visual asset folder containing volcano plots, heatmaps, and direction-split GO dot plots. |

## Methods Summary

Unstranded RNA-sequencing count data (STAR - Counts workflow) were downloaded dynamically via the Genomic Data Commons (GDC) API using `TCGAbiolinks`. Differential expression profiles were modeled using `DESeq2` under a simple design contrasting primary tumor against solid tissue normal samples. Transcripts matching adjusted p-value < 0.05 and absolute log2 fold change > 1 were designated as differentially expressed. 

Functional enrichment pipelines were executed via `clusterProfiler` targeting Gene Ontology (GO) Biological Processes. Survival analysis was performed using Kaplan-Meier overall survival testing and median-split two-sided log-rank analysis within the `survival` package framework to evaluate the prognostic relevance of candidate genes like MMP11 and COL11A1.

## Reproducing This Analysis

Ensure your local system runs R (>= 4.2) with the following baseline dependencies compiled:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "clusterProfiler", "org.Hs.eg.db",
                        "EnhancedVolcano", "ComplexHeatmap"))
install.packages(c("ggplot2", "dplyr", "pheatmap", "survival", "survminer", "UpSetR"))
```

## Full Cohort Cloud Validation
To run the full-scale, un-subsampled verification mapping without local desktop memory bounds, execute the standalone validation_full_cohort.R script within a cloud runtime (e.g., Kaggle kernel) with internet access enabled to pull down matching API manifests.


## Data Availability

All transcriptomic files are public, open-access assets managed by the National Cancer Institute Genomic Data Commons (https://portal.gdc.cancer.gov) under project identifiers **TCGA-BRCA** and **TCGA-LUAD**. No proprietary or protected health information (PHI) is included in this repository.

## License

Code assets are distributed under the open [MIT License](https://choosealicense.com/licenses/mit/). Primary genomic data usage follows the clinical policies outlined by the NIH Genomic Data Sharing initiatives.

## Contact

Dipul Poudel — ISMT College, Nepal  
University of Sunderland, UK  
dipulpoudel123@gmail.com · [GitHub Profile](https://github.com/imdipul)
