# Load necessary R packages
library(linkET)
library(tidyverse)
library(WGCNA)

# Load preprocessed data
load("~/output.Rdata")
exp_sub <- exp_unique
num <- "GSE5747"  # Dataset identifier (can be changed to other datasets like GSE99039, GSE63063, etc.)

# Transpose the expression matrix (genes as columns, samples as rows)
mydata <- data.frame(t(exp_sub))
colnames(mydata) <- rownames(exp_sub)  # Set column names as gene names
rownames(mydata) <- colnames(exp_sub)  # Set row names as sample names

# Filter genes and samples for quality control
gsg <- goodSamplesGenes(mydata, verbose = 3)
if (!gsg$allOK) {
  # Print and remove problematic genes and samples
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing genes:", paste(names(mydata)[!gsg$goodGenes], collapse = ",")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing samples:", paste(rownames(mydata)[!gsg$goodSamples], collapse = ",")))
  }
  mydata <- mydata[gsg$goodSamples, gsg$goodGenes]  # Retain only good samples and genes
}

nGenes <- ncol(mydata)  # Number of genes
nSamples <- nrow(mydata)  # Number of samples

# Load phenotype data (sample traits)
traitData <- as.data.frame(cbind(colnames(exp_sub), Group))
colnames(traitData) <- c("Sample", "Group")
rownames(traitData) <- traitData$Sample

# Encode phenotype data into binary variables
traitData$AD <- ifelse(traitData$Group == "AD", 1, 0)
traitData$HC <- ifelse(traitData$Group == "AD", 0, 1)
traitData$PD <- ifelse(traitData$Group == 2, 1, 0)
traitData$HC <- ifelse(traitData$Group == 1, 1, 0)
traitData$SZ <- ifelse(traitData$Group == 2, 1, 0)
traitData$HC <- ifelse(traitData$Group == 1, 1, 0)

traitData <- traitData[, c(3, 4)]  # Retain only relevant columns

# Match phenotype data with expression data samples
fpkmSamples <- rownames(mydata)
traitSamples <- rownames(traitData)
traitRows <- match(fpkmSamples, traitSamples)
dataTraits <- traitData[traitRows, ]

# Plot sample clustering dendrogram with trait heatmap
sampleTree2 <- hclust(dist(mydata), method = "average")
traitColors <- numbers2colors(dataTraits, signed = FALSE)
pdf(paste0(num, "_sample-subtype-cluster.pdf"), width = 8, height = 6)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(dataTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2)
dev.off()

# Determine the optimal soft-thresholding power for network construction
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(mydata, powerVector = powers, verbose = 5)

# Plot scale-free topology and mean connectivity to select the best power value
pdf(paste0(num, "_beta-value.pdf"), width = 8, height = 6)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")  # Threshold for scale-free topology
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# Construct the co-expression network using the optimal power value
cor <- WGCNA::cor  # Use WGCNA's cor function to avoid conflicts
net <- blockwiseModules(mydata, power = sft$powerEstimate, maxBlockSize = nGenes,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE, saveTOMFileBase = paste(num, "-TOM"),
                        verbose = 3)
table(net$colors)  # View module assignments

# Convert module labels to colors for visualization
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs <- net$MEs  # Module eigengenes
geneTree <- net$dendrograms[[1]]  # Gene dendrogram

# Extract genes in each module
g <- merge(as.data.frame(table(net$colors)), as.data.frame(table(moduleColors)))
genes <- list()
for (i in 1:nrow(g)) {
  genes[[i]] <- names(net$colors)[net$colors == g$Var1[[i]]]
}
names(genes) <- g$moduleColors

# Plot the gene dendrogram with module colors
pdf(paste0(num, "_-genes-modules.pdf"), width = 8, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Calculate module-trait correlations
MEs0 <- moduleEigengenes(mydata, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
colnames(MEs) <- str_remove(colnames(MEs), "ME")  # Remove "ME" prefix
MEs <- dplyr::select(MEs, -grey)  # Exclude the grey module

moduleTraitCor <- cor(MEs, dataTraits, use = "p")  # Correlation matrix
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)  # P-values

# Plot module-trait correlation heatmap
pdf(paste0(num, "_wgcna_heatmap.pdf"), height = 10, width = 6)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dataTraits), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()

# Select significant modules based on correlation with traits
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
moduleTraitCor <- as.data.frame(moduleTraitCor)
selec_mod <- moduleTraitCor[order(abs(moduleTraitCor[, 2]), decreasing = TRUE), ]
selec_mod <- head(selec_mod, n = ceiling(0.5 * nrow(selec_mod)))  # Top 50% modules
print(selec_mod)  # View selected modules