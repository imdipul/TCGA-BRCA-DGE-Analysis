# Convergent and Tissue-Specific Transcriptomic Programs in Breast and Lung Adenocarcinoma

Reproducible R pipeline for a comparative differential gene expression (DGE), functional enrichment, survival, and protein-interaction network analysis of TCGA-BRCA (breast invasive carcinoma) and TCGA-LUAD (lung adenocarcinoma), using public TCGA RNA-sequencing data.

This repository accompanies the manuscript *"Convergent and tissue-specific transcriptomic programs in breast and lung adenocarcinoma: a comparative RNA-sequencing, survival, and protein-interaction network analysis using TCGA data"* (Dipul Poudel, ISMT College, Chitwan, Nepal — submitted to PLOS ONE).

## Key findings

- Of 43,813 genes tested in TCGA-BRCA, **11,226 were significantly differentially expressed** (padj < 0.05, |log2FC| > 1): 6,517 upregulated, 4,709 downregulated in tumor vs. normal tissue.
- The two most significant genes were **MMP11** (log2FC = 6.32) and **COL11A1** (log2FC = 6.30), both established stromal remodeling markers in breast cancer.
- Comparison with TCGA-LUAD (14,774 significant DEGs) identified **5,795 DEGs shared between both cancers** — enriched for cell-cycle and extracellular-matrix processes — alongside 5,431 BRCA-specific DEGs (hormone/lipid-metabolism biology) and 8,979 LUAD-specific DEGs (adaptive-immune/ciliary-motility biology).
- Three MAGE-A cancer-testis antigens (**MAGEA3, MAGEA6, MAGEA12**) were independently among the most extreme shared genes in both cancers, consistent with epigenetic reactivation as a convergent pan-cancer mechanism.
- Kaplan-Meier analysis showed a non-significant survival trend for high **COL11A1** expression (p = 0.068) but no association for **MMP11** (p = 0.46) — despite MMP11 having the stronger differential expression statistics — demonstrating that DGE significance does not imply prognostic value.

## Repository contents

| File | Description |
|---|---|
| `analysis.Rmd` | Full, reproducible analysis pipeline: data download through comparative GO enrichment |
| `TCGA_BRCA_DESeq2_results.csv` | Full BRCA DESeq2 output (53,081 genes) |
| `TCGA_LUAD_DESeq2_results.csv` | Full LUAD DESeq2 output (52,165 genes) |
| `TCGA_BRCA_GO_Enrichment_results.csv` | GO Biological Process enrichment, all significantly upregulated BRCA genes |
| `GO_shared_pancancer.csv` | GO enrichment for the 5,795 DEGs shared between BRCA and LUAD |
| `GO_BRCA_specific.csv` | GO enrichment for the 5,431 DEGs unique to BRCA |
| `GO_LUAD_specific.csv` | GO enrichment for the 8,979 DEGs unique to LUAD |
| `results.rds` | Saved R results object |
| `figures/` | All generated figures: volcano plot, top-50 heatmap, GO dot plots, COL11A1/MMP11 Kaplan-Meier curves, STRING PPI network, BRCA-vs-LUAD DEG overlap (Venn + UpSet), supplementary STRING PubMed-enrichment figure |

## Methods summary

RNA-sequencing count data (STAR – Counts) for TCGA-BRCA (200 randomly subsampled tumors, seed = 42, plus all 113 available normals) and TCGA-LUAD (200 subsampled tumors plus all 59 normals) were obtained via the GDC using `TCGAbiolinks` and analyzed with `DESeq2` (padj < 0.05, |log2FC| > 1). Functional enrichment used `clusterProfiler` (GO Biological Process). Kaplan–Meier survival analysis used the full 1,218-patient BRCA cohort. A STRING protein–protein interaction network was built from the top 10 most significant BRCA DEGs. BRCA and LUAD DEG sets were intersected to identify shared versus tissue-specific genes and biological processes.

Full methodological detail, statistical results, and discussion are provided in the accompanying manuscript.

**Note on sample size:** the discovery differential expression analysis used a 200-tumor subsample per cancer (rather than the full 1,218/539-tumor cohorts) to keep DESeq2's dispersion estimation tractable on standard personal computing hardware. This is discussed explicitly as a limitation in the manuscript, together with the rationale and its expected impact on the results.

## Reproducing this analysis

Requirements: R (≥ 4.2) with the following packages:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "clusterProfiler", "org.Hs.eg.db",
                        "EnhancedVolcano", "ComplexHeatmap"))
install.packages(c("ggplot2", "dplyr", "pheatmap", "survival", "survminer", "UpSetR"))
```

Then open `analysis.Rmd` in RStudio and run all chunks in order. Data download (via `TCGAbiolinks`) requires an internet connection and will take approximately 20–40 minutes depending on connection speed; the full pipeline (download through comparative GO enrichment) takes roughly 1–2 hours on a standard laptop.

## Data availability

All raw data are public and unrestricted, available from the NCI Genomic Data Commons (https://portal.gdc.cancer.gov) under project identifiers **TCGA-BRCA** and **TCGA-LUAD**. No proprietary, restricted-access, or personally identifying data are used or redistributed in this repository — only derived summary statistics (DESeq2 and GO enrichment result tables).

## Citation

If you use this pipeline or its results, please cite:

> Poudel D. Convergent and tissue-specific transcriptomic programs in breast and lung adenocarcinoma: a comparative RNA-sequencing, survival, and protein-interaction network analysis using TCGA data. *Submitted to PLOS ONE, 2026.*

(Citation will be updated with full publication details upon acceptance.)

## License

Code in this repository is released under the [MIT License](https://choosealicense.com/licenses/mit/). TCGA data is subject to the [NIH Genomic Data Sharing Policy](https://gds.nih.gov/) and is publicly available without restriction for the project identifiers used here.

## Contact

Dipul Poudel — ISMT College, Chitwan, Nepal
dipulpoudel123@gmail.com · [GitHub](https://github.com/imdipul)
