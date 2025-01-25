# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(scMetabolism)

# Load gene score data
load("Gene_scores_combine.Rdata")

# Extract high-scoring genes
gene_list <- Gene_scores[Gene_scores$Sum_Score >= 10, ]$Overlap_Genes
gene_list <- unique(all_data$Gene)

# Load single-cell RNA-seq data
load("pbmc_sub.RData")

# Check if genes exist in the Seurat object
valid_genes <- gene_list[gene_list %in% rownames(seurat_obj_sub)]
missing_genes <- gene_list[!gene_list %in% rownames(seurat_obj_sub)]

# Warn about missing genes
if (length(missing_genes) > 0) {
  warning("以下基因不在 Seurat 对象中：", paste(missing_genes, collapse = ", "))
}

# Create a subset of cells and genes
seurat_obj_sub <- subset(pbmc_sub, Diagnosis2 = "AD")
seurat_obj_sub2 <- subset(seurat_obj_sub, features = valid_genes)

# Set cell identities
Idents(seurat_obj_sub) <- "broad.cell.type"  # AD

# Set working directory
name <- "HD"

# Perform metabolic pathway analysis
seurat_obj_sub2 <- sc.metabolism.Seurat(
  obj = seurat_obj_sub2,
  method = "VISION",
  imputation = FALSE,
  ncores = 5,
  metabolism.type = "KEGG"
)

# Extract metabolic scores
metabolism.matrix_g <- seurat_obj_sub2@assays$METABOLISM$score

# Calculate absolute sum of metabolic scores
metabolism.matrix_g$abs_sum <- rowSums(abs(metabolism.matrix_g))

# Sort by absolute sum and select top 30 pathways
metabolism.matrix_g <- metabolism.matrix_g[order(metabolism.matrix_g$abs_sum, decreasing = TRUE), ]
input.pathway <- rownames(head(metabolism.matrix_g, 30))

# Plot metabolic pathways
pdf(paste0(name, "-disDotPlot.metabolism.gene.pdf"), height = 10, width = 8)
DotPlot.metabolism2(
  obj = seurat_obj_sub2,
  pathway = input.pathway,
  phenotype = cellgroup,
  norm = "y",
  gradient_colors = rev(c("#ffccd7", "#d7a5e9", "#9feae7", "#bfeeee", "#e5fcff", "#ffe0bd", "#ffe070"))
)
dev.off()

# Prepare data for visualization
seurat_obj_sub2$celltype <- seurat_obj_sub2[[cellgroup]][, 1]
mscore_data <- data.frame(
  t(seurat_obj_sub2@assays[["METABOLISM"]][["score"]]),
  seurat_obj_sub2$celltype
)

# Calculate average metabolic scores per cell type
avg_sM <- aggregate(
  mscore_data[, 1:(ncol(mscore_data) - 1)],
  list(mscore_data$seurat_obj_sub2.celltype),
  mean
)

# Format data for visualization
rownames(avg_sM) <- avg_sM$Group.1
avg_sM <- data.frame(t(avg_sM[, -1]))
avg_sM$KEGG <- rownames(seurat_obj_sub2@assays[["METABOLISM"]][["score"]])
rownames(avg_sM) <- avg_sM$KEGG

# Extract top 5 pathways per cell type
c_k_l <- c()
for (c in 1:ncol(avg_sM)) {
  c_k <- avg_sM[order(avg_sM[, c]), ]$KEGG[1:5]
  c_k_l <- c(c_k_l, c_k)
}
c_k_l <- unique(c_k_l)

# Prepare data for heatmap
c_k_d <- avg_sM[avg_sM$KEGG %in% c_k_l, ]
data_long <- melt(as.matrix(c_k_d[, -ncol(c_k_d)]), varnames = c("Pathway", "Cell_Type"), value.name = "value")

# Add significance labels
data_long <- data_long %>%
  mutate(Significance = case_when(
    abs(value) >= 0.7 ~ "***",
    abs(value) >= 0.4 ~ "**",
    abs(value) >= 0.2 ~ "*",
    TRUE ~ ""
  ))

# Plot heatmap
heatmap_dot_plot <- ggplot(data_long, aes(x = Cell_Type, y = Pathway, size = abs(value), fill = value)) +
  geom_point(shape = 21, color = "white", stroke = 0.5) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  geom_text(aes(label = Significance), color = "red", size = 4, fontface = "bold", vjust = -1) +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#c095e4", "#fcedf2", "#ffa0c5"))(111),
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  scale_size_continuous(
    range = c(5, 10),
    name = "Absolute Value"
  ) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  theme(
    axis.text.y = element_text(size = 10, face = "italic", color = c(
      rep("#ff595e", nrow(c_k_d) / 4),
      rep("#8ac926", nrow(c_k_d) / 4),
      rep("#1982c4", nrow(c_k_d) / 4),
      rep("#6a4c93", nrow(c_k_d) / 4)
    )),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    legend.position = "right"
  ) +
  guides(
    fill = guide_colorbar(barheight = 15, barwidth = 1, ticks.linewidth = 0.5),
    size = guide_legend(override.aes = list(fill = "grey50"))
  )

# Save heatmap
pdf(paste0(name, "-disDotplotsheatmap.metabolism.gene.pdf"), height = 8, width = 10)
print(heatmap_dot_plot)
dev.off()

# Save average metabolic scores
write.csv(avg_sM, file = paste0(name, '-ave_sMetaByCelltype.gene.csv'))

# Save metabolic data
save(metabolism.matrix_g, mscore_data, file = paste0(name, "_seurat_obj_metabolism.gene.Rdata"))