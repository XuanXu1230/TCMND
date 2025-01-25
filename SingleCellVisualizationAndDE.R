# Load required libraries
library(SCP)
library(Seurat)

# Load single-cell RNA-seq data
load("~/AD/PD/scRNA/seurat_obj_PD.Rdata")
name <- "PD"  # Dataset identifier (can be changed to "SZ" or other names)

# Extract counts matrix and metadata from SingleCellExperiment object
counts_matrix <- counts(notaras_brainOrganoids.sce)
print(dim(counts_matrix))  # Check dimensions of the counts matrix

# Extract normalized matrix
norm_matrix <- assay(notaras_brainOrganoids.sce, "log1p_sctransform")
print(dim(norm_matrix))  # Check dimensions of the normalized matrix

# Extract cell metadata
cell_metadata <- colData(notaras_brainOrganoids.sce)
head(cell_metadata)  # Preview cell metadata

# Extract gene metadata
gene_metadata <- rowData(notaras_brainOrganoids.sce)
head(gene_metadata)  # Preview gene metadata

# Standardize gene names (replace underscores with dashes)
rownames(counts_matrix) <- gsub("_", "-", rownames(counts_matrix))
rownames(gene_metadata) <- gsub("_", "-", rownames(gene_metadata))

# Align genes and cells between counts matrix and metadata
common_genes <- intersect(rownames(counts_matrix), rownames(gene_metadata))
counts_matrix <- counts_matrix[common_genes, , drop = FALSE]
gene_metadata <- gene_metadata[common_genes, , drop = FALSE]

common_cells <- intersect(colnames(counts_matrix), rownames(cell_metadata))
counts_matrix <- counts_matrix[, common_cells, drop = FALSE]
cell_metadata <- cell_metadata[common_cells, , drop = FALSE]

# Ensure counts_matrix is a matrix and metadata is a data frame
counts_matrix <- as.matrix(counts_matrix)
cell_metadata <- as.data.frame(cell_metadata)
gene_metadata <- as.data.frame(gene_metadata)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = cell_metadata
)

# Add gene metadata to Seurat object
seurat_obj[["RNA"]]@meta.features <- gene_metadata

# Check Seurat object
print(seurat_obj)

# Preprocessing steps
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Calculate mitochondrial percentage
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")

# Run SCTransform for normalization and variance stabilization
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimensionality reduction and clustering
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- RunTSNE(seurat_obj, check_duplicates = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# Visualize clustering results
DimPlot(seurat_obj, label = TRUE, reduction = "tsne", group.by = c("labels", "condition"))

# Standard SCP analysis
seurat_obj <- Standard_SCP(srt = seurat_obj)

# Save and visualize cell dimensionality reduction plots
pdf(paste0(name, "-CellDimPlot.pdf"), height = 5, width = 10)
CellDimPlot(
  srt = seurat_obj, group.by = c("labels", "condition"),
  reduction = "TSNE", theme_use = "theme_blank"
)
dev.off()

# Run PAGA for trajectory analysis
seurat_obj <- RunPAGA(
  srt = seurat_obj, group_by = "labels",
  linear_reduction = "pca", nonlinear_reduction = "umap"
)
pdf(paste0(name, "-PAGAPlot.pdf"))
PAGAPlot(srt = seurat_obj, reduction = "UMAP", label = TRUE, label_insitu = TRUE, label_repel = TRUE)
dev.off()

# Run RNA velocity analysis
sce.sub <- RunSCVELO(
  srt = seurat_obj, group_by = "labels",
  linear_reduction = "pca", nonlinear_reduction = "umap"
)
pdf(paste0(name, "-VelocityPlot.pdf"))
VelocityPlot(srt = pancreas_sub, reduction = "UMAP", group_by = "labels")
dev.off()

# Save Seurat object
save(seurat_obj, DEGs, file = "SZ_seurat_obj.Rdata")

# Differential expression analysis
seurat_obj <- RunDEtest(srt = seurat_obj, group_by = c("labels"), fc.threshold = 1, only.pos = FALSE)
pdf(paste0(name, "-VolcanoPlot.pdf"), height = 15, width = 15)
VolcanoPlot(srt = seurat_obj, group_by = "labels")
dev.off()

# Visualize marker genes
pdf(paste0(name, "-markerVocalno.pdf"), height = 8, width = 15)
markerVocalno(
  markers = DEGs,
  topn = 5,
  labelCol = colorRampPalette(ggsci::pal_npg("nrc")(10))(length(unique(seurat_obj$celltype)))
)
dev.off()

# Identify cluster markers
Idents(seurat_obj) <- "labels"
DEGs <- FindAllMarkers(seurat_obj, only.pos = FALSE,
                       min.pct = 0.25,
                       logfc.threshold = 0)

# Summarize cluster markers
table(DEGs$cluster)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(scRNAtoolVis)
library(ClusterGVis)
library(plot1cell)

# Define cell group and disease group columns
cellgroup <- "original_name"  # PD
disgroup <- "donor_status"

# Load sample table and define column names
sample_table <- read.csv("sample_table.csv")  # Replace with your sample table file
names(sample_table) <- c("Samples", "celltype", "CellNumber")

# Define sample order and colors
sample_table$Samples <- factor(sample_table$Samples, levels = c("CONTROL", "MCI", "AD"))  # AD
sample_table$Samples <- factor(sample_table$Samples, levels = c("Con", "HD"))  # HD
colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C")

# Plot cell type proportions
p1 <- ggplot(sample_table, aes(x = Samples, weight = CellNumber, fill = celltype)) +
  geom_bar(position = "fill", width = 0.7, size = 0.5, colour = '#222222') +
  scale_fill_manual(values = colors) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5, size = 16)) +
  labs(y = "Percentage") +
  RotatedAxis()

# Display the plot
p1

# Calculate cell type proportions
sample_table <- sample_table %>%
  group_by(Samples) %>%
  mutate(Proportion = CellNumber / sum(CellNumber))

# Plot alluvial diagram
p1 <- ggplot(sample_table, 
             aes(x = Samples, stratum = celltype, alluvium = celltype, y = Proportion, fill = celltype, label = celltype)) +  
  scale_y_continuous(label = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = colors) +  
  geom_flow() +  
  geom_stratum(width = 0.6, colour = "grey40") +  
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5, size = 16),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10)) +
  labs(y = "Percentage", fill = "Cell Type")

# Save the plot
pdf(paste0(name, "-Percentagegroup.pdf"), height = 8, width = 8)
p1
dev.off()

# Differential expression analysis
Idents(seurat_obj) <- disgroup  # HD
{
  seurat_obj$original_name <- as.factor(seurat_obj$original_name)
  seurat_obj <- subset(seurat_obj, subset = !is.na(original_name) & !is.na(donor_status))
  DefaultAssay(seurat_obj) <- "RNA"
  
  cell_types <- unique(seurat_obj$original_name)
  group_degs <- lapply(cell_types, function(cell_type) {
    degs <- FindMarkers(seurat_obj[, seurat_obj$original_name == cell_type], 
                        ident.1 = 'Healthy',
                        ident.2 = "Idiopathic Parkinson's disease")
    degs$cluster <- cell_type
    return(degs)
  })
}

# Combine differential expression results
degs_combined <- do.call(rbind, group_degs)
degs_combined$gene <- rownames(degs_combined)

# Volcano plot for differential expression results
pdf(paste0(name, "-jjVolcanogroup.pdf"), height = 8, width = 10)
jjVolcano(degs_combined, topGeneN = 0, tile.col = colors[1:length(cell_types)],
          pSize = 0.5, fontface = 'italic', flip = TRUE, aesCol = c('#c0c0fd', '#ffcfbc')) +  
  geom_text_repel(data = degs_combined[degs_combined$gene %in% mygene,],                  
                  aes(x = cluster, y = avg_log2FC, label = gene),                  
                  position = position_jitter(seed = 1),                   
                  show.legend = FALSE,                   
                  size = 2.5,                   
                  box.padding = unit(0.2, "lines"))
dev.off()

# Heatmap of top genes
pdf(paste0(name, "-FeatureDimPlottopgene.pdf"), height = 8, width = 10)
ht <- GroupHeatmap(
  srt = seurat_obj,
  features = topgene,
  group.by = cellgroup,
  heatmap_palette = "YlOrRd",
  cell_annotation = disgroup,
  cell_annotation_palette = c("Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
print(ht$plot)
dev.off()

# Feature annotation and heatmap
seurat_obj <- AnnotateFeatures(seurat_obj, species = "Homo_sapiens", db = c("Chromosome", "GeneType", "Enzyme", "TF", "CSPA", "VerSeDa"))
ht <- FeatureHeatmap(
  srt = seurat_obj, group.by = "broad.cell.type", features = DEGs$gene, feature_split = DEGs$cluster,
  species = "Homo_sapiens", db = c("GO_BP", "KEGG", "WikiPathway"), anno_terms = TRUE,
  label_color = "black",
  heatmap_palette = "PiYG",
  group_palette = "Set3",
  cell_annotation = NULL,
  cell_annotation_palette = "Set3",
  feature_annotation = NULL,
  feature_annotation_palette = "Set2",
  height = 5, width = 4
)

pdf(paste(name, "dis-FeatureHeatmap.pdf"), width = 20, height = 10)
print(ht$plot)
dev.off()

# Slingshot trajectory analysis
seurat_obj <- RunSlingshot(srt = seurat_obj, group.by = "cell_type", reduction = "UMAP")

# Group heatmap
ht <- GroupHeatmap(
  srt = seurat_obj,
  features = topgene,
  group.by = cellgroup,
  heatmap_palette = "YlOrRd",
  cell_annotation = disgroup,
  cell_annotation_palette = c("Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)

pdf(paste(name, "-GroupHeatmap.pdf"), width = 12, height = 10)
print(ht$plot)
dev.off()

# Feature dimension plot
pdf(paste0(name, "-GroupFeatureDimPlot.pdf"), height = 8, width = 8)
FeatureDimPlot(
  srt = seurat_obj_sub, features = topgene,
  compare_features = TRUE, label = TRUE, label_insitu = TRUE,
  reduction = "UMAP", theme_use = "theme_blank"
)
dev.off()

# Average heatmap
AverageHeatmap(object = seurat_obj, markerGene = topgene)

# Differential expression analysis for subsets
seurat_obj_sub <- subset(seurat_obj, Condition == "HD")
Idents(seurat_obj_sub) <- cellgroup
DEGs <- Seurat::FindAllMarkers(seurat_obj_sub,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# Save results
save(degs_combined, DEGs, file = paste0(name, "DEG_disgroup.Rdata"))

# Volcano plot for subset DEGs
pdf(paste0(name, "-disVolcanotopgene.pdf"), height = 8, width = 8)
jjVolcano(DEGs, topGeneN = 0, tile.col = colors[1:length(unique(seurat_obj_sub[[cellgroup]][,1]))],
          pSize = 0.5, fontface = 'italic', polar = TRUE, aesCol = c('#adb7dc', '#d7adbe')) +  
  geom_text_repel(data = DEGs[DEGs$gene %in% topgene,],                  
                  aes(x = cluster, y = avg_log2FC, label = gene),                  
                  position = position_jitter(seed = 1),                   
                  show.legend = FALSE,                   
                  size = 2.5,                   
                  box.padding = unit(0.2, "lines"))
dev.off()

# Enrichment analysis
seurat_obj_sub <- RunEnrichment(
  srt = seurat_obj_sub, group_by = cellgroup, db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)

pdf(paste0(name, "-disEnrichmentPlot.pdf"), height = 20, width = 40)
EnrichmentPlot(
  srt = seurat_obj_sub, group_by = cellgroup,
  plot_type = "network", palette = "Set2", padjustCutoff = 0.05
)
dev.off()

# GSEA analysis
seurat_obj_sub <- RunGSEA(
  srt = seurat_obj_sub, group_by = cellgroup, db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)

pdf(paste0(name, "-disGSEAPlot.pdf"), height = 20, width = 25)
GSEAPlot(
  srt = seurat_obj_sub, group_by = cellgroup, plot_type = "bar",
  direction = "both", topTerm = 20
)
dev.off()

# Feature statistics plot
pdf(paste0(name, "-disFeatureStatPlottopgene.pdf"), height = 10, width = 25)
FeatureStatPlot(
  srt = seurat_obj_sub, group.by = cellgroup, bg.by = cellgroup,
  stat.by = topgene, add_box = TRUE, nrow = 2, ncol = 4,
  comparisons = list(
    c("Astrocyte", "Microglia"),
    c("Microglia", "Neuron"),
    c("Neuron", "Oligodendrocyte"),
    c("Astrocyte", "Oligodendrocyte")
  )
)
dev.off()

# Circlize plot for cell type visualization
Idents(seurat_obj) <- "broad.cell.type"  # AD
circ_data <- prepare_circlize_data(seurat_obj, scale = 0.8)
cluster_colors <- brewer.pal(n = length(levels(seurat_obj)), "Set3")
group_colors <- brewer.pal(n = length(unique(seurat_obj$braaksc)), "Dark2")
rep_colors <- brewer.pal(n = length(unique(seurat_obj$Diagnosis2)), "Pastel1")

pdf(paste0(name, "-circlize_plot.pdf"), width = 6, height = 6)
plot_circlize(circ_data, do.label = TRUE, pt.size = 0.01, col.use = cluster_colors, bg.color = 'white', kde2d.n = 200, repel = TRUE, label.cex = 0.6)
add_track(circ_data, group = "braaksc", colors = group_colors, track_num = 2)
add_track(circ_data, group = "Diagnosis2", colors = rep_colors, track_num = 3)
dev.off()