###############################################################################
# BRB-seq K562 drug-treatment time course
# PCA, sample-distance heatmap, gene heatmap, DESeq2 and volcano plots
#
# Course: Omics in practice - Module 2b
# Author: Lindsey Bakx
#
# INPUT FILES (put these in the SAME folder as this script):
#   umi_counts.txt       - STARsolo UMI counts (1MM_Directional), genes x 24 barcodes
#   barcode_samples.tsv  - barcode -> sample number (1-24)
#   features.tsv         - Ensembl gene ID -> gene symbol (from Solo.out/Gene/raw)
#
# HOW TO RUN (RStudio):
#   1. Session > Set Working Directory > To Source File Location
#   2. Run one section at a time: put the cursor in a section and press
#      Ctrl+Alt+T, or select lines and press Ctrl+Enter.
#   All plots and tables are written to the folder "output".
#
# EXPERIMENTAL DESIGN (from the course slides, confirmed by the instructor):
#   8 timepoints x 3 replicates = 24 samples
#   Samples 1-8   = replicate 1 (0h, 6h, 1d, 2d, 3d, 4d, 7d, 8d)
#   Samples 9-16  = replicate 2 (same order)
#   Samples 17-24 = replicate 3 (same order)
###############################################################################


## 0. Settings ----------------------------------------------------------------
gene_of_interest <- "TBP"                 # your TF / gene of interest (symbol)
timepoints  <- c("0h", "6h", "1d", "2d", "3d", "4d", "7d", "8d")
reference   <- "0h"                       # start timepoint = control
padj_cutoff <- 0.05                       # adjusted p-value threshold
lfc_cutoff  <- 1                          # |log2FC| > 1  = more than 2-fold
min_count   <- 10                         # filtering rule: keep genes with
min_samples <- 3                          # >= 10 counts in >= 3 samples
n_top_var   <- 30                         # genes in the variable-gene heatmap

dir.create("output", showWarnings = FALSE)


## 1. Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
})


## 2. Read the input files ----------------------------------------------------
counts <- read.delim("umi_counts.txt", row.names = 1, check.names = FALSE)
features <- read.delim("features.tsv", header = FALSE, stringsAsFactors = FALSE)
barcode_samples <- read.table("barcode_samples.tsv", header = TRUE,
                              stringsAsFactors = FALSE)

cat("Count table:", nrow(counts), "genes x", ncol(counts), "barcodes\n")

# Ensembl ID -> gene symbol lookup
symbol_of <- setNames(features$V2, features$V1)


## 3. Sample metadata ---------------------------------------------------------
meta <- data.frame(barcode   = barcode_samples$Barcode,
                   sample_nr = barcode_samples$Sample,
                   stringsAsFactors = FALSE)
meta$replicate <- factor(paste0("rep", ceiling(meta$sample_nr / 8)))
meta$time      <- factor(timepoints[((meta$sample_nr - 1) %% 8) + 1],
                         levels = timepoints)
meta$time      <- relevel(meta$time, ref = reference)
meta$sample    <- sprintf("S%02d_%s_%s", meta$sample_nr, meta$time, meta$replicate)

# Fixed, transparent column order: by timepoint, then replicate
meta <- meta[order(meta$time, meta$replicate), ]
rownames(meta) <- meta$sample
print(meta)

write.table(meta, "output/sample_metadata.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)


## 4. Match counts to metadata + validation checks ----------------------------
stopifnot(setequal(colnames(counts), meta$barcode))   # same 24 barcodes
counts <- as.matrix(counts[, meta$barcode])           # put columns in meta order
colnames(counts) <- meta$sample                       # barcode -> sample name

stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))
stopifnot(identical(colnames(counts), rownames(meta)))
storage.mode(counts) <- "integer"

# Readable count table with sample names and gene symbols (for Excel / report)
named_counts <- data.frame(gene_id = rownames(counts),
                           symbol  = symbol_of[rownames(counts)],
                           counts, check.names = FALSE)
write.table(named_counts, "output/gene_counts_named.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("Total UMI counts per sample:\n")
print(colSums(counts))


## 5. Filter lowly expressed genes --------------------------------------------
keep <- rowSums(counts >= min_count) >= min_samples
counts_f <- counts[keep, ]
cat("Genes kept after filtering:", nrow(counts_f), "of", nrow(counts), "\n")


## 6. Build the DESeq2 object (raw counts!) -----------------------------------
# Design: expression depends on timepoint. Reference level = 0h, so every
# coefficient means "timepoint X versus 0h" (positive log2FC = higher at X).
# If the PCA shows that samples group by replicate, you can discuss with the
# instructor whether to use  design = ~ replicate + time  instead.
dds <- DESeqDataSetFromMatrix(countData = counts_f,
                              colData   = meta,
                              design    = ~ time)


## 7. Variance-stabilising transformation (for visualisation only) ------------
vsd <- vst(dds, blind = TRUE)
vmat <- assay(vsd)


## 8. PCA plot ----------------------------------------------------------------
pca <- plotPCA(vsd, intgroup = c("time", "replicate"), returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(PC1, PC2, color = time, shape = replicate)) +
  geom_point(size = 3.5) +
  geom_text_repel(aes(label = sub("_.*", "", name)), size = 3,
                  show.legend = FALSE) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  labs(title = "PCA of BRB-seq samples (VST, top 500 variable genes)",
       x = paste0("PC1: ", pct[1], "% variance"),
       y = paste0("PC2: ", pct[2], "% variance"),
       color = "Timepoint", shape = "Replicate") +
  theme_bw()
print(p_pca)
ggsave("output/PCA_plot.png", p_pca, width = 8, height = 6, dpi = 300)


## 9. Sample-to-sample distance heatmap ---------------------------------------
sample_dist <- dist(t(vmat))
pheatmap(as.matrix(sample_dist),
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         annotation_col = meta[, c("time", "replicate")],
         main = "Sample-to-sample distances (VST)",
         filename = "output/sample_distance_heatmap.png",
         width = 9, height = 8)


## 10. Heatmap of selected genes (top variable genes + gene of interest) ------
goi_id <- names(symbol_of)[symbol_of == gene_of_interest]
goi_id <- intersect(goi_id, rownames(vmat))
if (length(goi_id) == 0) {
  warning(gene_of_interest, " not found after filtering (too few counts?)")
}

gene_var <- apply(vmat, 1, var)
top_ids  <- names(sort(gene_var, decreasing = TRUE))[1:n_top_var]
sel_ids  <- unique(c(top_ids, goi_id))

hm <- vmat[sel_ids, ]
hm_labels <- ifelse(is.na(symbol_of[sel_ids]), sel_ids, symbol_of[sel_ids])
rownames(hm) <- make.unique(hm_labels)

pheatmap(hm,
         scale = "row",                      # row z-score: pattern per gene
         cluster_cols = FALSE,               # fixed order: 0h -> 8d
         cluster_rows = TRUE,
         gaps_col = seq(3, 21, by = 3),      # gap between timepoints
         annotation_col = meta[, c("time", "replicate")],
         main = paste0("Top ", n_top_var, " variable genes + ",
                       gene_of_interest, " (row z-score of VST)"),
         filename = "output/heatmap_selected_genes.png",
         width = 10, height = 9)


## 11. Run DESeq2 -------------------------------------------------------------
dds <- DESeq(dds)
print(resultsNames(dds))


## 12. Every timepoint versus 0h (includes 8d vs 0h) --------------------------
# This takes a few minutes: apeglm shrinkage is run for 7 comparisons.
all_results <- list()
summary_tab <- data.frame()

for (tp in setdiff(timepoints, reference)) {
  coef_name <- paste0("time_", tp, "_vs_", reference)
  stopifnot(coef_name %in% resultsNames(dds))
  message("Running ", coef_name, " ...")

  res    <- results(dds, name = coef_name, alpha = padj_cutoff)
  res_sh <- lfcShrink(dds, coef = coef_name, res = res,
                      type = "apeglm", quiet = TRUE)

  df <- as.data.frame(res_sh)
  df$gene_id    <- rownames(df)
  df$symbol     <- symbol_of[df$gene_id]
  df$comparison <- paste(tp, "vs", reference)
  df$category   <- "not significant"
  df$category[!is.na(df$padj) & df$padj < padj_cutoff &
                df$log2FoldChange >  lfc_cutoff] <- "up"
  df$category[!is.na(df$padj) & df$padj < padj_cutoff &
                df$log2FoldChange < -lfc_cutoff] <- "down"
  df <- df[order(df$padj), c("gene_id", "symbol", "baseMean",
                             "log2FoldChange", "lfcSE", "pvalue", "padj",
                             "category", "comparison")]

  write.csv(df, sprintf("output/DESeq2_%s_vs_%s.csv", tp, reference),
            row.names = FALSE)
  all_results[[tp]] <- df

  goi_row <- df[df$symbol %in% gene_of_interest, ]
  summary_tab <- rbind(summary_tab, data.frame(
    comparison   = paste(tp, "vs", reference),
    n_up         = sum(df$category == "up"),
    n_down       = sum(df$category == "down"),
    goi_log2FC   = if (nrow(goi_row)) round(goi_row$log2FoldChange[1], 3) else NA,
    goi_padj     = if (nrow(goi_row)) signif(goi_row$padj[1], 3) else NA))
}

names(summary_tab)[4:5] <- paste0(gene_of_interest, c("_log2FC", "_padj"))
print(summary_tab)
write.csv(summary_tab, "output/DE_summary_all_timepoints.csv", row.names = FALSE)


## 13. Volcano plot function --------------------------------------------------
make_volcano <- function(df, title, n_labels = 10) {
  df <- df[!is.na(df$padj), ]
  df$neg_log10_padj <- -log10(pmax(df$padj, 1e-300))
  top <- df[df$category != "not significant", ]
  top <- head(top[order(top$padj), ], n_labels)
  goi <- df[df$symbol %in% gene_of_interest, ]

  ggplot(df, aes(log2FoldChange, neg_log10_padj, color = category)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c(down = "#2166AC",
                                  "not significant" = "grey75",
                                  up = "#B2182B")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(data = top, aes(label = symbol), size = 3,
                    color = "black", max.overlaps = 30) +
    geom_point(data = goi, color = "darkorange", size = 3) +
    geom_label_repel(data = goi, aes(label = symbol), color = "darkorange",
                     size = 3.5, show.legend = FALSE) +
    labs(title = title,
         x = "Shrunken log2 fold change (apeglm)",
         y = "-log10(adjusted p-value)",
         color = NULL) +
    theme_bw()
}


## 14. Volcano plot: latest timepoint (8d) versus 0h --------------------------
p_volcano_8d <- make_volcano(all_results[["8d"]],
                             "Volcano plot: K562 8d vs 0h")
print(p_volcano_8d)
ggsave("output/volcano_8d_vs_0h.png", p_volcano_8d,
       width = 8, height = 6, dpi = 300)


## 15. Volcano plots: every timepoint versus 0h -------------------------------
for (tp in names(all_results)) {
  ggsave(sprintf("output/volcano_%s_vs_%s.png", tp, reference),
         make_volcano(all_results[[tp]],
                      paste0("Volcano plot: K562 ", tp, " vs ", reference)),
         width = 8, height = 6, dpi = 300)
}

# All seven in one figure (no gene labels, for a quick overview)
combined <- do.call(rbind, all_results)
combined <- combined[!is.na(combined$padj), ]
combined$comparison <- factor(combined$comparison,
                              levels = paste(setdiff(timepoints, reference),
                                             "vs", reference))
combined$neg_log10_padj <- -log10(pmax(combined$padj, 1e-300))

p_volcano_all <- ggplot(combined,
                        aes(log2FoldChange, neg_log10_padj, color = category)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c(down = "#2166AC",
                                "not significant" = "grey75",
                                up = "#B2182B")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  facet_wrap(~ comparison, nrow = 2) +
  labs(title = "Every timepoint versus 0h",
       x = "Shrunken log2 fold change", y = "-log10(adjusted p-value)",
       color = NULL) +
  theme_bw()
print(p_volcano_all)
ggsave("output/volcano_all_timepoints.png", p_volcano_all,
       width = 14, height = 7, dpi = 300)

# Number of DE genes over time (helps you choose which comparison to report)
de_long <- rbind(
  data.frame(comparison = summary_tab$comparison, direction = "up",
             n = summary_tab$n_up),
  data.frame(comparison = summary_tab$comparison, direction = "down",
             n = summary_tab$n_down))
de_long$comparison <- factor(de_long$comparison, levels = summary_tab$comparison)

p_de_counts <- ggplot(de_long, aes(comparison, n, fill = direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(down = "#2166AC", up = "#B2182B")) +
  labs(title = paste0("DE genes per timepoint (padj < ", padj_cutoff,
                      ", |log2FC| > ", lfc_cutoff, ")"),
       x = NULL, y = "Number of genes", fill = NULL) +
  theme_bw()
print(p_de_counts)
ggsave("output/DE_gene_counts_per_timepoint.png", p_de_counts,
       width = 8, height = 5, dpi = 300)


## 16. Gene of interest at replicate level ------------------------------------
if (length(goi_id) > 0) {
  goi_counts <- plotCounts(dds, gene = goi_id[1],
                           intgroup = c("time", "replicate"),
                           returnData = TRUE)
  p_goi <- ggplot(goi_counts, aes(time, count, color = replicate,
                                  group = replicate)) +
    geom_point(size = 3) +
    geom_line() +
    labs(title = paste0(gene_of_interest, " (", goi_id[1], ")"),
         x = "Timepoint", y = "Normalised count") +
    theme_bw()
  print(p_goi)
  ggsave(paste0("output/", gene_of_interest, "_counts_over_time.png"), p_goi,
         width = 7, height = 5, dpi = 300)
}


## 17. Save session info (reproducibility) ------------------------------------
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
message("Done. All results are in the folder: ", normalizePath("output"))
