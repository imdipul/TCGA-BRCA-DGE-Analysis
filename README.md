# Convergent and Tissue-Specific Transcriptomic Programs in Breast and Lung Adenocarcinoma

Reproducible R pipeline for a comparative differential gene expression (DGE), functional enrichment, multi-variable survival, and protein-interaction network analysis of TCGA-BRCA (breast invasive carcinoma) and TCGA-LUAD (lung adenocarcinoma), using public TCGA RNA-sequencing data.

This repository accompanies the manuscript *"Convergent and tissue-specific transcriptomic programs in breast and lung adenocarcinoma: a comparative RNA-sequencing, survival, and protein-interaction network analysis using TCGA data"*.

## Key Findings

- **Differential Expression Landscape:** Of 43,813 genes evaluated in TCGA-BRCA, **11,226 were significantly differentially expressed** (\(p_{\text{adj}} < 0.05\), \(\vert{}\log_2\text{FC}\vert{} > 1\)), comprising 6,517 upregulated and 4,709 downregulated transcripts in tumor vs. normal tissue.
- **Top Stromal Drivers:** The two most significant dysregulated transcripts were **MMP11** (\(\log_2\text{FC} = 6.32\)) and **COL11A1** (\(\log_2\text{FC} = 6.30\)), identifying them as dominant stroma-remodeling biomarkers.
- **Cross-Cancer Overlaps:** Comparative filtering against TCGA-LUAD (14,774 significant DEGs) isolated **5,795 shared core DEGs** across both malignancies. Functional pathways diverge cleanly between shared, BRCA-specific (lipid/hormonal metabolism), and LUAD-specific (adaptive immunity and ciliary motility) modules.
- **Direction-Split Pathway Enrichment:** Strategic splitting of shared and tissue-specific gene sets by regulation direction revealed sharp functional partitioning. For example, LUAD-specific upregulated programs are strongly enriched for organelle fission, nuclear division, and meiotic cell cycle processes, whereas tissue-specific downregulated streams govern separate physiological profiles.
- **Confounded Prognostic Assessment:** Multivariate Cox proportional hazards modeling (adjusting for patient age and clinical TNM staging matrix criteria, \(N = 1082\)) demonstrated that high DGE significance does not uniformly imply survival liability. In adjusted models, high expression of **MMP11** (\(\text{HR} = 1.103\), \(95\%\text{ CI } [0.917 - 1.326]\), \(p = 0.298\)) and **COL11A1** (\(\text{HR} = 1.146\), \(95\%\text{ CI } [0.967 - 1.358]\), \(p = 0.115\)) were not independent prognostic predictors, confirming that crude transcriptional abundance is secondary to formal staging criteria.

## Repository Contents

| File / Folder | Description |
|---|---|
| `analysis.R` / `analysis.Rmd` | Complete, end-to-end reproducible analysis pipeline from raw count acquisition to functional clustering. |
| `TCGA_BRCA_DESeq2_results.csv` | Full processed BRCA DESeq2 output matrix (53,081 genes). |
| `TCGA_LUAD_DESeq2_results.csv` | Full processed LUAD DESeq2 output matrix (52,165 genes). |
| `GO_shared_pancancer.csv` | Baseline unified GO enrichment records for multi-cancer shared DEGs. |
| `GO_shared_up.csv` | **[NEW]** GO Biological Process terms unique to shared *upregulated* genes. |
| `GO_shared_down.csv` | **[NEW]** GO Biological Process terms unique to shared *downregulated* genes. |
| `GO_brca_only_up.csv` | **[NEW]** GO Biological Process terms specific to BRCA *upregulated* genes. |
| `GO_brca_only_down.csv` | **[NEW]** GO Biological Process terms specific to BRCA *downregulated* genes. |
| `GO_luad_only_up.csv` | **[NEW]** GO Biological Process terms specific to LUAD *upregulated* genes. |
| `GO_luad_only_down.csv` | **[NEW]** GO Biological Process terms specific to LUAD *downregulated* genes. |
| `df.rds` / `results.rds` | Structured serialization R data objects for prompt workspace restoration. |
| `figures/` | Visual asset folder containing volcano plots, heatmaps, and direction-split GO dot plots. |

## Methods Summary

Unstranded RNA-sequencing count data (STAR - Counts) were downloaded dynamically via the Genomic Data Commons (GDC) API using `TCGAbiolinks`. Differential expression profiles were modeled using `DESeq2` under a simple design contrasting primary tumor against solid tissue normal samples. Transcripts matching \(p_{\text{adj}} < 0.05\) and \(\vert{}\log_2\text{FC}\vert{} > 1\) were designated as differentially expressed. 

Functional enrichment pipelines were executed via `clusterProfiler` targeting Gene Ontology (GO) Biological Processes, with input lists stratified explicitly by regulatory direction (upregulated vs. downregulated) to prevent reciprocal signal cancellation. Survival analysis was performed using multivariate Cox proportional hazards modeling within the `survival` package framework, controlling for age at index and pathologic staging structures.

## Reproducing This Analysis

Ensure your local system runs R (\(\ge 4.2\)) with the following baseline dependencies compiled:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "clusterProfiler", "org.Hs.eg.db",
                        "EnhancedVolcano", "ComplexHeatmap"))
install.packages(c("ggplot2", "dplyr", "pheatmap", "survival", "survminer", "UpSetR"))
```

To execute, fetch the scripts locally, configure your RStudio active working directory to the repository home node, and step through `analysis.R` or compile the markdown script.

## Data Availability

All transcriptomic files are public, open-access assets managed by the National Cancer Institute Genomic Data Commons (https://portal.gdc.cancer.gov) under project identifiers **TCGA-BRCA** and **TCGA-LUAD**. No proprietary or protected health information (PHI) is included in this repository.

## License

Code assets are distributed under the open [MIT License](https://choosealicense.com/licenses/mit/). Primary genomic data usage follows the clinical policies outlined by the NIH Genomic Data Sharing initiatives.

## Contact

Dipul Poudel — ISMT College, Nepal  
University of Sunderland, UK  
dipulpoudel123@gmail.com · [GitHub Profile](https://github.com/imdipul)
