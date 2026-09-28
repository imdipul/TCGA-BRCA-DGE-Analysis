library(TCGAbiolinks)
library(DESeq2)

# 1. Rebuild the query manifest so R knows where the files are on your disk
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

# 2. Load the downloaded data object into memory
data <- GDCprepare(query)

# 3. Separate tumor and normal sample names
all_tumor <- colnames(data)[colData(data)$sample_type == "Primary Tumor"]
all_normal <- colnames(data)[colData(data)$sample_type == "Solid Tissue Normal"]

# 4. Subset 200 tumor samples to keep your Mac's RAM completely safe
set.seed(42) 
selected_tumor <- sample(all_tumor, 200)

# 5. Combine the 113 normals and 200 selected tumors
keep_samples <- c(all_normal, selected_tumor)
data_subsampled <- data[, keep_samples]

# 6. Clean up factor levels
colData(data_subsampled)$sample_type <- factor(colData(data_subsampled)$sample_type, 
                                               levels = c("Solid Tissue Normal", "Primary Tumor"))

# 7. Build the dataset object
dds <- DESeqDataSet(data_subsampled, design = ~ sample_type)

# 8. Filter out low-count genes
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]

# 9. Run DESeq safely without parallel workers
dds <- DESeq(dds)

# 10. Extract and print results
res <- results(dds)
head(res)


# Save the results table as a CSV file in your project folder
write.csv(as.data.frame(res), file = "TCGA_BRCA_DESeq2_results.csv")


# Convert results to a proper dataframe
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Clean Ensembl IDs by removing version numbers (e.g., ENSG00000000003.15 -> ENSG00000000003)
res_df$ensembl_clean <- sub("\\..*", "", res_df$gene)

# Check your total significant, upregulated, and downregulated gene counts
print(paste("Total significant genes:", sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE)))
print(paste("Upregulated genes:", sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)))
print(paste("Downregulated genes:", sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)))





# 1. Clean the Ensembl IDs properly
res_df$ensembl_clean <- sub("\\..*", "", res_df$gene)

# 2. Map Ensembl IDs to official gene symbols
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = res_df$ensembl_clean,
                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res_df$symbol <- symbols[res_df$ensembl_clean]

# 3. Replace any missing symbols back with their original ID so there are no blanks
res_df$symbol <- ifelse(is.na(res_df$symbol), res_df$gene, res_df$symbol)



library(EnhancedVolcano)

EnhancedVolcano(res_df,
                lab = res_df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Tumour vs Normal - Differential Gene Expression',
                subtitle = 'TCGA-BRCA Subsampled',
                legendPosition = 'right'
)


library(pheatmap)
library(dplyr)

# 1. Isolate the top 50 most significant DEGs
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# 2. Normalize the counts using Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes$gene, ]

# 3. Rename matrix rows to friendly Gene Symbols
rownames(mat) <- sig_genes$symbol

# 4. Scale rows (Z-score normalization)
mat_scaled <- t(scale(t(mat)))

# 5. Build the sample annotation tracker column cleanly
annotation <- data.frame(
  Condition = colData(data_subsampled)$sample_type,
  row.names = colnames(mat)
)

# 6. Render the professional pheatmap
pheatmap(mat_scaled,
  annotation_col = annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = 'Top 50 DEGs: Tumour vs Normal'
)




library(pheatmap)
library(dplyr)

# 1. Isolate the top 50 most significant DEGs
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# 2. Normalize the counts using Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes$gene, ]

# 3. Rename matrix rows to friendly Gene Symbols
rownames(mat) <- sig_genes$symbol

# 4. Scale rows (Z-score normalization)
mat_scaled <- t(scale(t(mat)))

# 5. Build the sample annotation tracker column cleanly
annotation <- data.frame(
  Condition = colData(data_subsampled)$sample_type,
  row.names = colnames(mat)
)

# 6. Render the professional pheatmap
pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = 'Top 50 DEGs: Tumour vs Normal'
)



library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Filter out genes that are highly significantly upregulated in the tumor
upregulated_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.01, log2FoldChange > 1.5) %>%
  pull(ensembl_clean)

# 2. Run the Biological Process (BP) Gene Ontology enrichment analysis
go_results <- enrichGO(
  gene          = upregulated_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = 'ENSEMBL',
  ont           = "BP", # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# 3. View the top enriched biological pathways
head(go_results)

# 4. Generate a clean dotplot of the top 10 pathways
dotplot(go_results, showCategory = 10, title = "Enriched Biological Processes in BRCA Tumors")

# Save the full GO pathway analysis results
write.csv(as.data.frame(go_results), file = "TCGA_BRCA_GO_Enrichment_results.csv")



head(res_df[order(res_df$padj), c("symbol", "log2FoldChange", "padj")], 10)


library(clusterProfiler)
library(org.Hs.eg.db)

# Get significantly upregulated genes (clean Ensembl IDs, no version suffix)
up_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
  pull(ensembl_clean)

# Convert Ensembl IDs to Entrez IDs (required format for enrichGO)
gene_ids <- bitr(up_genes,
  fromType = 'ENSEMBL',
  toType = 'ENTREZID',
  OrgDb = org.Hs.eg.db)

# Run GO enrichment - Biological Process
ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego, showCategory = 20,
  title = 'GO Biological Process Enrichment')

library(clusterProfiler)
library(org.Hs.eg.db)

# Get significantly upregulated genes (clean Ensembl IDs, no version suffix)
up_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
  pull(ensembl_clean)

# Convert Ensembl IDs to Entrez IDs (required format for enrichGO)
gene_ids <- bitr(up_genes,
                 fromType = 'ENSEMBL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)

# Run GO enrichment - Biological Process
ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego, showCategory = 20,
        title = 'GO Biological Process Enrichment')



sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE)
sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)   # upregulated
sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)  # downregulated

library(pheatmap)
library(dplyr)

# 1. Isolate the top 50 most significant DEGs
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# 2. Normalize the counts using Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes$gene, ]

# 3. Rename matrix rows to friendly Gene Symbols
rownames(mat) <- sig_genes$symbol

# 4. Scale rows (Z-score normalization)
mat_scaled <- t(scale(t(mat)))

# 5. Build the sample annotation tracker column cleanly
annotation <- data.frame(
  Condition = colData(data_subsampled)$sample_type,
  row.names = colnames(mat)
)

# 6. Render the professional pheatmap
pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = 'Top 50 DEGs: Tumour vs Normal'
)

library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Filter out genes that are highly significantly upregulated in the tumor
upregulated_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.01, log2FoldChange > 1.5) %>%
  pull(ensembl_clean)

# 2. Run the Biological Process (BP) Gene Ontology enrichment analysis
go_results <- enrichGO(
  gene          = upregulated_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = 'ENSEMBL',
  ont           = "BP", # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# 3. View the top enriched biological pathways
head(go_results)

# 4. Generate a clean dotplot of the top 10 pathways
dotplot(go_results, showCategory = 10, title = "Enriched Biological Processes in BRCA Tumors")


install.packages(c("survival", "survminer"))

# Force a fresh download and overwrite the broken files
install.packages(c("survival", "survminer"), force = TRUE)

library(survival)
library(survminer)

# 1. Load all required biological and clinical libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survival)
library(survminer)

# 2. Reload the saved data object back into active RAM
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)
data <- GDCprepare(query)

# 3. Pull clinical/survival metadata
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# 4. Extract and normalize expressions for the full tumor cohort
all_tumor_data <- data[, colData(data)$sample_type == "Primary Tumor"]
dds_full <- DESeqDataSet(all_tumor_data, design = ~ 1)
dds_full <- dds_full[rowSums(counts(dds_full)) >= 10, ]
vsd_full <- vst(dds_full, blind = TRUE)

# 5. Map your lead genes
mmp11_id  <- "ENSG00000099953.10"
col11a1_id <- "ENSG00000060718.22"

expr_mmp11  <- assay(vsd_full)[mmp11_id, ]
expr_col11a1 <- assay(vsd_full)[col11a1_id, ]

# 6. Merge expression profiles with clinical survival tracks
expr_df <- data.frame(
  sample = colnames(vsd_full),
  patient = substr(colnames(vsd_full), 1, 12),
  MMP11 = expr_mmp11,
  COL11A1 = expr_col11a1
)
merged <- merge(expr_df, clinical, by.x = "patient", by.y = "submitter_id")

# 7. Construct clinical outcome survival data points
merged$time  <- ifelse(merged$vital_status == "Dead", merged$days_to_death, merged$days_to_last_follow_up)
merged$event <- ifelse(merged$vital_status == "Dead", 1, 0)
merged <- merged[!is.na(merged$time) & merged$time >= 0, ]

# 8. Render Kaplan-Meier Overall Survival curve for MMP11
merged$MMP11_group <- ifelse(merged$MMP11 > median(merged$MMP11), "High", "Low")
fit_mmp11 <- survfit(Surv(time, event) ~ MMP11_group, data = merged)

ggsurvplot(fit_mmp11, data = merged, pval = TRUE, risk.table = TRUE,
           title = "MMP11 Expression and Overall Survival - TCGA-BRCA",
           legend.title = "MMP11 Expression", legend.labs = c("High", "Low"))

# 9. Render Kaplan-Meier Overall Survival curve for COL11A1
merged$COL11A1_group <- ifelse(merged$COL11A1 > median(merged$COL11A1), "High", "Low")
fit_col11a1 <- survfit(Surv(time, event) ~ COL11A1_group, data = merged)

ggsurvplot(fit_col11a1, data = merged, pval = TRUE, risk.table = TRUE,
           title = "COL11A1 Expression and Overall Survival - TCGA-BRCA",
           legend.title = "COL11A1 Expression", legend.labs = c("High", "Low"))


# ==========================================
# EXTENSION 3: LUNG ADENOCARCINOMA (TCGA-LUAD)
# ==========================================
library(TCGAbiolinks)
library(DESeq2)

# 1. Query and download LUAD data
query_luad <- GDCquery(
  project = 'TCGA-LUAD',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)
GDCdownload(query_luad)
data_luad <- GDCprepare(query_luad)

# 2. Subsample tumor cohort for balance and RAM protection (59 Normal vs 200 Tumor)
all_tumor_luad  <- colnames(data_luad)[colData(data_luad)$sample_type == "Primary Tumor"]
all_normal_luad <- colnames(data_luad)[colData(data_luad)$sample_type == "Solid Tissue Normal"]

set.seed(42)
selected_tumor_luad <- sample(all_tumor_luad, 200)
keep_luad <- c(all_normal_luad, selected_tumor_luad)
data_luad_sub <- data_luad[, keep_luad]

colData(data_luad_sub)$sample_type <- factor(colData(data_luad_sub)$sample_type,
                                             levels = c("Solid Tissue Normal", "Primary Tumor"))

# 3. Run DESeq2 Pipeline safely
dds_luad <- DESeqDataSet(data_luad_sub, design = ~ sample_type)
dds_luad <- dds_luad[rowSums(counts(dds_luad)) >= 10, ]
dds_luad <- DESeq(dds_luad)

# 4. Extract and Save Results Table
res_luad <- results(dds_luad)
res_luad_df <- as.data.frame(res_luad)
res_luad_df$gene <- rownames(res_luad_df)
res_luad_df$ensembl_clean <- sub("\\..*", "", res_luad_df$gene)

write.csv(res_luad_df, file = "TCGA_LUAD_DESeq2_results.csv")

# 5. Print out total significant genes for your notes
print(paste("LUAD significant genes:", sum(res_luad_df$padj < 0.05 & abs(res_luad_df$log2FoldChange) > 1, na.rm = TRUE)))


# Clear background memory space
gc()

# Force a clean re-run of the LUAD query and download step
library(TCGAbiolinks)
library(DESeq2)

query_luad <- GDCquery(
  project = 'TCGA-LUAD',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

# This command will skip broken pieces and download cleanly
GDCdownload(query_luad, method = "client", files.per.chunk = 10)
data_luad <- GDCprepare(query_luad)


# 1. Subsample LUAD tumors for balance and RAM protection (all Normals + 200 Tumors)
all_tumor_luad  <- colnames(data_luad)[colData(data_luad)$sample_type == "Primary Tumor"]
all_normal_luad <- colnames(data_luad)[colData(data_luad)$sample_type == "Solid Tissue Normal"]

set.seed(42)
selected_tumor_luad <- sample(all_tumor_luad, 200)
keep_luad <- c(all_normal_luad, selected_tumor_luad)
data_luad_sub <- data_luad[, keep_luad]

colData(data_luad_sub)$sample_type <- factor(colData(data_luad_sub)$sample_type,
                                             levels = c("Solid Tissue Normal", "Primary Tumor"))

# 2. Run the DESeq2 pipeline safely on Lung Cancers
dds_luad <- DESeqDataSet(data_luad_sub, design = ~ sample_type)
dds_luad <- dds_luad[rowSums(counts(dds_luad)) >= 10, ]
dds_luad <- DESeq(dds_luad)

# 3. Extract and save results table
res_luad <- results(dds_luad)
res_luad_df <- as.data.frame(res_luad)
res_luad_df$gene <- rownames(res_luad_df)
res_luad_df$ensembl_clean <- sub("\\..*", "", res_luad_df$gene)

write.csv(res_luad_df, file = "TCGA_LUAD_DESeq2_results.csv")

# 4. Print total significant genes
print(paste("LUAD significant genes:", sum(res_luad_df$padj < 0.05 & abs(res_luad_df$log2FoldChange) > 1, na.rm = TRUE)))

# 1. Safely load your saved BRCA results from your disk
res_df <- read.csv("TCGA_BRCA_DESeq2_results.csv")

# 2. Extract significant DEGs for BRCA (using your clean ensembl column)
brca_sig_genes <- res_df$ensembl_clean[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1)]

# 3. Extract significant DEGs for LUAD
luad_sig_genes <- res_luad_df$ensembl_clean[which(res_luad_df$padj < 0.05 & abs(res_luad_df$log2FoldChange) > 1)]

# 4. Calculate the intersections and tissue-specific sets
shared_degs     <- intersect(brca_sig_genes, luad_sig_genes)
brca_only_degs  <- setdiff(brca_sig_genes, luad_sig_genes)
luad_only_degs  <- setdiff(luad_sig_genes, brca_sig_genes)

# 5. Print out the overlap breakdown metrics for your notes
print(paste("Shared core cancer genes across both tissues:", length(shared_degs)))
print(paste("Genes completely unique to Breast Cancer (BRCA):", length(brca_only_degs)))
print(paste("Genes completely unique to Lung Cancer (LUAD):", length(luad_only_degs)))

# 6. Check exactly how your top 5 BRCA hub genes are acting in Lung Cancer
hubs <- c("ENSG00000099953", "ENSG00000123508", "ENSG00000060718", "ENSG00000117650", "ENSG00000090889")
hub_names <- c("MMP11", "COL10A1", "COL11A1", "NEK2", "KIF4A")

print("--- Checking Top Hub Expression in Lung Adenocarcinoma ---")
for(i in 1:length(hubs)){
  match_idx <- which(res_luad_df$ensembl_clean == hubs[i])
  if(length(match_idx) > 0){
    cat(hub_names[i], "in LUAD -> log2FC:", round(res_luad_df$log2FoldChange[match_idx], 2), 
        " | padj:", format.pval(res_luad_df$padj[match_idx], eps=1e-10), "\n")
  }
}


cat("Shared genes:", length(shared_degs), "\n")
cat("BRCA unique:", length(brca_only_degs), "\n")
cat("LUAD unique:", length(luad_only_degs), "\n")


# Safely re-calculate the exact cross-cancer overlaps
res_brca <- read.csv("TCGA_BRCA_DESeq2_results.csv")
res_luad <- read.csv("TCGA_LUAD_DESeq2_results.csv")

brca_sig <- res_brca$ensembl_clean[which(res_brca$padj < 0.05 & abs(res_brca$log2FoldChange) > 1)]
luad_sig <- res_luad$ensembl_clean[which(res_luad$padj < 0.05 & abs(res_luad$log2FoldChange) > 1)]

shared_degs    <- intersect(brca_sig, luad_sig)
brca_only_degs <- setdiff(brca_sig, luad_sig)
luad_only_degs <- setdiff(luad_sig, brca_sig)

cat("Shared Core Pan-Cancer Genes:", length(shared_degs), "\n")
cat("Unique Breast Cancer (BRCA) Genes:", length(brca_only_degs), "\n")
cat("Unique Lung Cancer (LUAD) Genes:", length(luad_only_degs), "\n")


# 1. Reload the data tables fresh
res_brca <- read.csv("TCGA_BRCA_DESeq2_results.csv")
res_luad <- read.csv("TCGA_LUAD_DESeq2_results.csv")

# 2. Force-clean the Ensembl IDs for BOTH datasets from their raw row names
# (This ensures the columns look 100% identical)
res_brca$ensembl_clean <- sub("\\..*", "", res_brca$X)
if(!"X" %in% colnames(res_brca)) res_brca$ensembl_clean <- sub("\\..*", "", res_brca$gene)

res_luad$ensembl_clean <- sub("\\..*", "", res_luad$X)
if(!"X" %in% colnames(res_luad)) res_luad$ensembl_clean <- sub("\\..*", "", res_luad$gene)

# 3. Re-extract the true significant gene lists
brca_sig <- res_brca$ensembl_clean[which(res_brca$padj < 0.05 & abs(res_brca$log2FoldChange) > 1)]
luad_sig <- res_luad$ensembl_clean[which(res_luad$padj < 0.05 & abs(res_luad$log2FoldChange) > 1)]

# 4. Recalculate intersections
shared_degs    <- intersect(brca_sig, luad_sig)
brca_only_degs <- setdiff(brca_sig, luad_sig)
luad_only_degs <- setdiff(luad_sig, brca_sig)

# 5. Print out the real metrics
cat("Shared Core Pan-Cancer Genes:", length(shared_degs), "\n")
cat("Unique Breast Cancer (BRCA) Genes:", length(brca_only_degs), "\n")
cat("Unique Lung Cancer (LUAD) Genes:", length(luad_only_degs), "\n")



library(clusterProfiler)
library(org.Hs.eg.db)

# Convert each set to Entrez IDs
shared_entrez    <- bitr(shared_degs,    fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
brca_only_entrez <- bitr(brca_only_degs, fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
luad_only_entrez <- bitr(luad_only_degs, fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Hs.eg.db)

go_shared <- enrichGO(shared_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
go_brca   <- enrichGO(brca_only_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
go_luad   <- enrichGO(luad_only_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)

dotplot(go_shared, showCategory=15, title="Shared Pan-Cancer DEGs: GO Biological Process")
dotplot(go_brca,   showCategory=15, title="BRCA-Specific DEGs: GO Biological Process")
dotplot(go_luad,   showCategory=15, title="LUAD-Specific DEGs: GO Biological Process")

write.csv(as.data.frame(go_shared), "GO_shared_pancancer.csv")
write.csv(as.data.frame(go_brca), "GO_BRCA_specific.csv")
write.csv(as.data.frame(go_luad), "GO_LUAD_specific.csv")

library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Reload the files fresh from your disk
res_brca <- read.csv("TCGA_BRCA_DESeq2_results.csv")
res_luad <- read.csv("TCGA_LUAD_DESeq2_results.csv")

# 2. Clean up the columns uniformly
res_brca$ensembl_clean <- sub("\\..*", "", res_brca$X)
if(!"X" %in% colnames(res_brca)) res_brca$ensembl_clean <- sub("\\..*", "", res_brca$gene)

res_luad$ensembl_clean <- sub("\\..*", "", res_luad$X)
if(!"X" %in% colnames(res_luad)) res_luad$ensembl_clean <- sub("\\..*", "", res_luad$gene)

# 3. Safely extract your clean significance lists
brca_sig <- res_brca$ensembl_clean[which(res_brca$padj < 0.05 & abs(res_brca$log2FoldChange) > 1)]
luad_sig <- res_luad$ensembl_clean[which(res_luad$padj < 0.05 & abs(res_luad$log2FoldChange) > 1)]

# 4. Compute the intersections
shared_degs    <- intersect(brca_sig, luad_sig)
brca_only_degs <- setdiff(brca_sig, luad_sig)
luad_only_degs <- setdiff(luad_sig, brca_sig)

# 5. Translate them into Entrez IDs (Ignore any normal % mapping alerts)
shared_entrez    <- bitr(shared_degs,    fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
brca_only_entrez <- bitr(brca_only_degs, fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
luad_only_entrez <- bitr(luad_only_degs, fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Hs.eg.db)

# 6. Run the 3 separate GO Enrichment analyses
go_shared <- enrichGO(shared_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
go_brca   <- enrichGO(brca_only_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
go_luad   <- enrichGO(luad_only_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)

# 7. Write the CSV tables straight to your folder directory
write.csv(as.data.frame(go_shared), "GO_shared_pancancer.csv")
write.csv(as.data.frame(go_brca), "GO_BRCA_specific.csv")
write.csv(as.data.frame(go_luad), "GO_LUAD_specific.csv")

# 8. Render the Shared Pan-Cancer plot first
dotplot(go_shared, showCategory=15, title="Shared Pan-Cancer DEGs: GO Biological Process")


dotplot(go_brca, showCategory=15, title="BRCA-Specific DEGs: GO Biological Process")
dotplot(go_luad, showCategory=15, title="LUAD-Specific DEGs: GO Biological Process")

res_brca <- read.csv("TCGA_BRCA_DESeq2_results.csv")

# Filter rows matching our target symbols securely
brca_targets <- res_brca[res_brca$symbol %in% c("MMP11", "COL11A1") | res_brca$gene %in% c("MMP11", "COL11A1"), ]
print(brca_targets)

res_brca <- read.csv("TCGA_BRCA_DESeq2_results.csv")

# Clean the Ensembl IDs in column X by removing version numbers
res_brca$ensembl_clean <- sub("\\..*", "", res_brca$X)

# Define the exact Ensembl IDs for MMP11 and COL11A1
mmp11_ensembl  <- "ENSG00000099953"
col11a1_ensembl <- "ENSG00000060718"

# Print out your rows cleanly
brca_targets <- res_brca[res_brca$ensembl_clean %in% c(mmp11_ensembl, col11a1_ensembl), ]
print(brca_targets[, c("X", "baseMean", "log2FoldChange", "padj")])


res_luad <- read.csv("TCGA_LUAD_DESeq2_results.csv")

res_brca$ensembl_clean <- sub("\\..*", "", res_brca$X)
res_luad$ensembl_clean <- sub("\\..*", "", res_luad$X)

# Isolate significant upregulated genes with a strong signal (log2FC > 1.5)
brca_up <- res_brca$ensembl_clean[which(res_brca$padj < 0.05 & res_brca$log2FoldChange > 1.5)]
luad_up <- res_luad$ensembl_clean[which(res_luad$padj < 0.05 & res_luad$log2FoldChange > 1.5)]

shared_up <- intersect(brca_up, luad_up)
cat("Universally Upregulated Core Genes:", length(shared_up), "\n")


res_brca$combined_lfc <- abs(res_brca$log2FoldChange) + abs(res_luad$log2FoldChange[match(res_brca$ensembl_clean, res_luad$ensembl_clean)])
top_shared <- res_brca[res_brca$ensembl_clean %in% shared_up, ]
top_shared <- top_shared[order(-top_shared$combined_lfc), ]

# Print the top 5 rows
head(top_shared[, c("X", "baseMean", "log2FoldChange", "padj")], 5)



load(".RData")

library(survival); library(dplyr)

clinical$stage_simple <- case_when(
  grepl("Stage IV", clinical$ajcc_pathologic_stage) ~ "IV",
  grepl("Stage III", clinical$ajcc_pathologic_stage) ~ "III",
  grepl("Stage II", clinical$ajcc_pathologic_stage) ~ "II",
  grepl("Stage I$|Stage IA|Stage IB", clinical$ajcc_pathologic_stage) ~ "I",
  TRUE ~ NA_character_
)
clinical$stage_simple <- factor(clinical$stage_simple, levels = c("I","II","III","IV"))

build_cox_data <- function(gene_id) {
  expr <- assay(vsd_full)[gene_id, ]
  df <- data.frame(patient = substr(colnames(vsd_full), 1, 12), expr = as.numeric(expr))
  merged <- merge(df, clinical, by.x = "patient", by.y = "submitter_id")
  merged$time  <- ifelse(merged$vital_status == "Dead", merged$days_to_death, merged$days_to_last_follow_up)
  merged$event <- ifelse(merged$vital_status == "Dead", 1, 0)
  merged <- merged[!is.na(merged$time) & merged$time >= 0, ]
  merged$age <- merged$age_at_index
  merged <- merged[!is.na(merged$age) & !is.na(merged$stage_simple), ]
  merged$expr_z <- as.numeric(scale(merged$expr))
  merged
}

mmp11_df   <- build_cox_data("ENSG00000099953.10")
col11a1_df <- build_cox_data("ENSG00000060718.22")

cox_mmp11   <- coxph(Surv(time, event) ~ expr_z + age + stage_simple, data = mmp11_df)
cox_col11a1 <- coxph(Surv(time, event) ~ expr_z + age + stage_simple, data = col11a1_df)

extract_hr <- function(model, name, df) {
  s <- summary(model)
  hr <- s$coefficients["expr_z", "exp(coef)"]
  ci <- s$conf.int["expr_z", c("lower .95","upper .95")]
  p  <- s$coefficients["expr_z", "Pr(>|z|)"]
  cat(name, ": HR =", round(hr,3), " 95% CI [", round(ci[1],3), "-", round(ci[2],3), "]  p =", signif(p,4), "\n")
  cat(name, "model N =", nrow(df), "\n")
}
extract_hr(cox_mmp11, "MMP11", mmp11_df)
extract_hr(cox_col11a1, "COL11A1", col11a1_df)



library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)

# Re-create the query definition (this does NOT download files)
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

data <- GDCprepare(query)

tumor_full <- data[, colData(data)$sample_type == 'Primary Tumor'] 
dds_full <- DESeqDataSet(tumor_full, design = ~ 1) 
dds_full <- dds_full[rowSums(counts(dds_full)) >= 10, ] 
vsd_full <- vst(dds_full, blind = TRUE)

clinical <- GDCquery_clinic(project = 'TCGA-BRCA', type = 'clinical')


library(survival); library(dplyr)

clinical$stage_simple <- case_when(
  grepl("Stage IV", clinical$ajcc_pathologic_stage) ~ "IV",
  grepl("Stage III", clinical$ajcc_pathologic_stage) ~ "III",
  grepl("Stage II", clinical$ajcc_pathologic_stage) ~ "II",
  grepl("Stage I$|Stage IA|Stage IB", clinical$ajcc_pathologic_stage) ~ "I",
  TRUE ~ NA_character_
)
clinical$stage_simple <- factor(clinical$stage_simple, levels = c("I","II","III","IV"))

build_cox_data <- function(gene_id) {
  expr <- assay(vsd_full)[gene_id, ]
  df <- data.frame(patient = substr(colnames(vsd_full), 1, 12), expr = as.numeric(expr))
  merged <- merge(df, clinical, by.x = "patient", by.y = "submitter_id")
  merged$time  <- ifelse(merged$vital_status == "Dead", merged$days_to_death, merged$days_to_last_follow_up)
  merged$event <- ifelse(merged$vital_status == "Dead", 1, 0)
  merged <- merged[!is.na(merged$time) & merged$time >= 0, ]
  merged$age <- merged$age_at_index
  merged <- merged[!is.na(merged$age) & !is.na(merged$stage_simple), ]
  merged$expr_z <- as.numeric(scale(merged$expr))
  merged
}

mmp11_df   <- build_cox_data("ENSG00000099953.10")
col11a1_df <- build_cox_data("ENSG00000060718.22")

cox_mmp11   <- coxph(Surv(time, event) ~ expr_z + age + stage_simple, data = mmp11_df)
cox_col11a1 <- coxph(Surv(time, event) ~ expr_z + age + stage_simple, data = col11a1_df)

extract_hr <- function(model, name, df) {
  s <- summary(model)
  hr <- s$coefficients["expr_z", "exp(coef)"]
  ci <- s$conf.int["expr_z", c("lower .95","upper .95")]
  p  <- s$coefficients["expr_z", "Pr(>|z|)"]
  cat(name, ": HR =", round(hr,3), " 95% CI [", round(ci[1],3), "-", round(ci[2],3), "]  p =", signif(p,4), "\n")
  cat(name, "model N =", nrow(df), "\n")
}
extract_hr(cox_mmp11, "MMP11", mmp11_df)
extract_hr(cox_col11a1, "COL11A1", col11a1_df)


list.files()

temp_data <- readRDS("df.rds")
summary(temp_data)

# 1. Load the pre-calculated results matrices
res_brca <- read.csv("TCGA_BRCA_DESeq2_results.csv")
res_luad <- read.csv("TCGA_LUAD_DESeq2_results.csv")

# 2. Extract significant DEGs based on thresholds
brca_sig <- res_brca$ensembl_clean[which(res_brca$padj < 0.05 & abs(res_brca$log2FoldChange) > 1)]
luad_sig <- res_luad$ensembl_clean[which(res_luad$padj < 0.05 & abs(res_luad$log2FoldChange) > 1)]

# 3. Calculate intersections and specific groups
shared_degs     <- intersect(brca_sig, luad_sig)
brca_only_degs  <- setdiff(brca_sig, luad_sig)
luad_only_degs  <- setdiff(luad_sig, brca_sig)





library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Helper function to format tables safely using an explicit namespace mapping
prepare_lfc_table <- function(df) {
  if (!"ensembl_clean" %in% colnames(df)) {
    df$ensembl_clean <- rownames(df)
  }
  # Explicitly call dplyr::select to prevent conflicts with AnnotationDbi
  df %>% as.data.frame() %>% dplyr::select(ensembl_clean, log2FoldChange)
}

# 1. Format Log Fold Change tables (Fixed with namespace mapping)
b_lfc <- prepare_lfc_table(res_brca)
l_lfc <- prepare_lfc_table(res_luad)

# 2. Extract direction-split groups
shared_up      <- b_lfc %>% filter(ensembl_clean %in% shared_degs, log2FoldChange > 0) %>% pull(ensembl_clean)
shared_down    <- b_lfc %>% filter(ensembl_clean %in% shared_degs, log2FoldChange < 0) %>% pull(ensembl_clean)
brca_only_up   <- b_lfc %>% filter(ensembl_clean %in% brca_only_degs, log2FoldChange > 0) %>% pull(ensembl_clean)
brca_only_down <- b_lfc %>% filter(ensembl_clean %in% brca_only_degs, log2FoldChange < 0) %>% pull(ensembl_clean)
luad_only_up   <- l_lfc %>% filter(ensembl_clean %in% luad_only_degs, log2FoldChange > 0) %>% pull(ensembl_clean)
luad_only_down <- l_lfc %>% filter(ensembl_clean %in% luad_only_degs, log2FoldChange < 0) %>% pull(ensembl_clean)

# 3. Core GO Pathway Function
run_go_bp <- function(ids, title) {
  if(length(ids) == 0) {
    cat("Skipping", title, "- No genes in this group.\n")
    return(data.frame())
  }
  clean_ids <- gsub("\\..*$", "", ids)
  entrez <- bitr(clean_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  
  if(!is.null(ego) && nrow(ego) > 0) {
    print(dotplot(ego, showCategory = 12, title = title))
  } else {
    cat("No significant terms for:", title, "\n")
  }
  as.data.frame(ego)
}

# 4. Execute functional runs and generate plots
go_shared_up   <- run_go_bp(shared_up, "Shared - Upregulated")
go_shared_down <- run_go_bp(shared_down, "Shared - Downregulated")
go_brca_up     <- run_go_bp(brca_only_up, "BRCA-only - Upregulated")
go_brca_down   <- run_go_bp(brca_only_down, "BRCA-only - Downregulated")
go_luad_up     <- run_go_bp(luad_only_up, "LUAD-only - Upregulated")
go_luad_down   <- run_go_bp(luad_only_down, "LUAD-only - Downregulated")

# 5. Export results sheets to your directory
write.csv(go_shared_up, "GO_shared_up.csv"); write.csv(go_shared_down, "GO_shared_down.csv")
write.csv(go_brca_up, "GO_brca_only_up.csv"); write.csv(go_brca_down, "GO_brca_only_down.csv")
write.csv(go_luad_up, "GO_luad_only_up.csv"); write.csv(go_luad_down, "GO_luad_only_down.csv")
