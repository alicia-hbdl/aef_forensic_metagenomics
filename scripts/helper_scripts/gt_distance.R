# Goal: Evaluate and compare different metagenomic pipeline configurations by projecting their 
# species abundance profiles into reduced-dimensional space (PCA, t-SNE, UMAP), and identifying
# which configuration most closely matches a known ground truth profile.

# - Each column represents a pipeline configuration (identified by a unique run ID)
#   - One column is labeled "truth" and serves as the ground truth benchmark
#   - Run IDs are linked to their settings in a metadata table
# - Each row corresponds to a species
# - Matrix values are species read counts or relative abundances

# Currently running in "mock" mode on a combined Bracken report file 
library(tidyverse)
library(matrixStats)
library(Rtsne)
library(umap)
library(patchwork)
library(reshape2)  # for melt()
library(ggplotify)
library(pheatmap)
library(viridis)
library(ggrepel)
# Example manual color scale
sample_colors <- c(
  "ZC1_S4" = "red",
  "ZC2_S5" = "orange",
  "ZC3_S6" = "yellow",
  "ZC4_S7" = "green",
  "ZP1_S20" = "blue",
  "ZP2_S21" = "purple",
  "ZP3_S22" = "pink"
)

# Read CSV file with species as row names and samples as columns
file <- read.csv("/Users/aliciahobdell/Desktop/final_project/zymobiomics_folder/results/runs/run_0705_1620/combined_breports.csv", row.names = 1)

# TODO: Aggregate combined Bracken reports by run.
#       For each run ID, sum or average the 7 sample profiles to create a single vector.
#       Final format: rows = species, columns = run IDs (including "truth").

# Replace NAs with 0.
file[is.na(file)] <- 0

# Filter out low-abundance species
# Keep species (rows) that have at least 20 reads in at least 50% of the samples
to_keep <- rowSums(file >= 20) >= ncol(file) / 2
filtered <- file[to_keep,]

# Normalize
# Compute total classified reads per sample (column sums)
lib_sizes <- colSums(filtered)

# Visualize library sizes per sample
df_lib <- data.frame(Run = names(lib_sizes), Reads = lib_sizes)

# Ensure Run is a factor matching sample_colors
df_lib$Run <- factor(df_lib$Run, levels = names(sample_colors))

p1 <- ggplot(df_lib, aes(x = Run, y = Reads, fill = Run)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = sample_colors) +
  labs(x = "Pipeline Run", y = "Total Reads", title = "Library Size per Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Correlation: species abundance vs. library size (raw counts)
cor_vals_raw <- apply(filtered, 1, \(x) cor(x, lib_sizes))
df_cor_raw <- data.frame(Correlation = cor_vals_raw)
p2 <- ggplot(df_cor_raw, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_raw), sd = sd(cor_vals_raw)), color = "red", linewidth = 1) +
  labs(x = "Correlation Coefficient", y = "Density", title = "Raw Abundance vs. Library Size") +
  theme_minimal()

# Normalize to CPM
filtered_cpm <- t(t(filtered) / lib_sizes * 1e6)

# Correlation after CPM normalization
cor_vals_cpm <- apply(filtered_cpm, 1, \(x) cor(x, lib_sizes))
df_cor_cpm <- data.frame(Correlation = cor_vals_cpm)
p3 <- ggplot(df_cor_cpm, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(cor_vals_cpm), sd = sd(cor_vals_cpm)), color = "red", linewidth = 1) +
  labs(x = "Correlation (CPM)", y = "Number of Species", title = "CPM-Normalized Correlation") +
  theme_minimal()

# Log mean vs. log standard deviation (CPM scale)
df_var_cpm <- data.frame(Mean = log1p(rowMeans(filtered_cpm)), SD = log1p(rowSds(filtered_cpm)))
cor_cpm <- cor(df_var_cpm$Mean, df_var_cpm$SD)
p4 <- ggplot(df_var_cpm, aes(x = Mean, y = SD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  annotate("text", x = min(df_var_cpm$Mean), y = max(df_var_cpm$SD),
           label = paste0("r = ", round(cor_cpm, 3)), hjust = 0, size = 4, fontface = "italic") +
  labs(x = "log(Mean CPM + 1)", y = "log(SD CPM + 1)", title = "Species Variability (CPM)") +
  theme_minimal()

# Log2 transform for downstream analysis
filtered_logcpm <- log2(filtered_cpm+1)

# Log mean vs. log SD (logCPM scale)
df_var_logcpm <- data.frame(Mean = log1p(rowMeans(filtered_logcpm)), SD = log1p(rowSds(filtered_logcpm)))
cor_logcpm <- cor(df_var_logcpm$Mean, df_var_logcpm$SD)
p5 <- ggplot(df_var_logcpm, aes(x = Mean, y = SD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  annotate("text", x = max(df_var_logcpm$Mean)-0.2, y = max(df_var_logcpm$SD),
           label = paste0("r = ", round(cor_logcpm, 3)), hjust = 0, size = 4, fontface = "italic") +
  labs(x = "log(Mean logCPM + 1)", y = "log(SD logCPM + 1)", title = "Species Variability (logCPM)") +
  theme_minimal()

combined_transformations <- p1 / (p2 + p3) / (p4 + p5) 
ggsave("combined_transformations.png", combined_transformations, width = 10, height = 10) 

# Use PCA, t-SNE, and UMAP to reduce species-space vectors into 2D.
#  Plot the reduced 2D coordinates of each pipeline setting.
#  Annotate points with their setting name (including “truth”).

# PCA: visualize configurations in reduced space
pca <- prcomp(t(filtered_logcpm), scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>% rownames_to_column("sample")
pca_df$color <- sample_colors[pca_df$sample]
p6 <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample, color = sample)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = sample_colors) +
  labs(title = "PCA", x = "PC1", y = "PC2")

# t-SNE: visualize nonlinear separation
tsne <- Rtsne(t(filtered_logcpm), perplexity = 2)
tsne_df <- as.data.frame(tsne$Y)
colnames(tsne_df) <- c("tsne_1", "tsne_2")
tsne_df$sample <- colnames(filtered_logcpm)
tsne_df$color <- sample_colors[tsne_df$sample]
p7 <- ggplot(tsne_df, aes(x = tsne_1, y = tsne_2, label = sample, color = sample)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = sample_colors) +
  labs(title = "t-SNE", x = "t-SNE 1", y = "t-SNE 2")

# UMAP: another non-linear projection
umap_cfg <- umap.defaults; umap_cfg$n_neighbors <- 3
umap_res <- umap(t(filtered_logcpm), config = umap_cfg)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("umap_1", "umap_2")
umap_df$sample <- colnames(filtered_logcpm)
umap_df$color <- sample_colors[umap_df$sample]
p8 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, label = sample, color = sample)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = sample_colors) +
  labs(title = "UMAP", x = "UMAP 1", y = "UMAP 2") 

# Helper: wrap pheatmap as a ggplot-compatible grob
pheatmap_grob <- function(mat) {
  p <- pheatmap(
    mat,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color = viridis(100, option = "D", direction = -1),
    angle_col = 45,  # <<<<< Rotate bottom labels
    silent = TRUE
  )
  as.ggplot(p[[4]])  # Extract heatmap grob and wrap in ggplot-compatible form
}

# Distances in original logCPM space (species abundance profiles)
raw_dists   <- as.matrix(dist(t(filtered_logcpm)))
p9  <- pheatmap_grob(raw_dists)

# Distances in PCA space (first 2 principal components)
pca_coords <- pca_df[, c("PC1", "PC2")]
rownames(pca_coords) <- pca_df$sample
pca_dists <- as.matrix(dist(pca_coords))
pca_dists <- pca_dists / max(pca_dists)
p10 <- pheatmap_grob(pca_dists)

# Distances in UMAP space
rownames(umap_df) <- umap_df$sample
umap_dists  <- as.matrix(dist(umap_df[, c("umap_1", "umap_2")]))
umap_dists <- umap_dists / max(umap_dists)
p11 <- pheatmap_grob(umap_dists)

# Distances in t-SNE space
rownames(tsne_df) <- tsne_df$sample
tsne_dists  <- as.matrix(dist(tsne_df[, c("tsne_1", "tsne_2")]))
tsne_dists <- tsne_dists / max(tsne_dists)
p12 <- pheatmap_grob(tsne_dists)

combined_distances <- (p6 + p7 + p8) / (p10 + p11 + p12)
ggsave("combined_distances.png", combined_distances, width = 12, height = 8) 

# TODO: Use distances from "truth" column to rank best-matching configurations
