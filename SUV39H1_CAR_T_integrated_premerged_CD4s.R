# SUV39H1 CAR-T CD4+ Cell Analysis Script
# This script processes and analyzes CD4+ CAR-T cell scRNA-seq data
# to investigate the effects of SUV39H1 knockout (gSUV) compared to Mock control

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)
library(paletteer)
library(RColorBrewer)
library(tidyverse)
library(EnhancedVolcano)
library(harmony)
library(Nebulosa)
library(UCell)
library(ggsignif)
library(reshape2)
library(viridis)
library(ggrepel)
library(colorspace)
library(Scillus) # For proportion plots

# Set up parallel processing for efficiency
library(future)
plan("multiprocess", workers = 8, strategy = "multicore")
options(future.globals.maxSize = 1500 * 2024^2)

# Color setup for visualization
mycolors <- pals::kelly(18)
mycolors.CD4 <- mycolors[4:16]
mycolors.CD4.KO <- darken(mycolors.CD4, 0.2) # Darker version for KO/gSUV

# Custom function for filtering non-TCR genes
"%!in%" <- function(x, y) !("%in%"(x, y))

# Load data (assuming the object is available)
# If not, you would need to load from a saved .RData file
# load("path/to/your/SUV39H1_CAR_T_integrated_premerged_CD4s_cleaned.Rdata")

#-----------------------------------------------------------------------------
# 1. Data Preprocessing (if needed) and Setting Up Cell Metadata
#-----------------------------------------------------------------------------

# Ensure genotype labels are consistent
genotype_mapping <- c(
  "Mock.1" = "Mock", "Mock.2" = "Mock", "Mock.3" = "Mock",
  "Suv.ko.1" = "gSUV", "Suv.ko.2" = "gSUV", "Suv.ko.3" = "gSUV"
)

# Apply genotype mapping to all cells
orig_idents <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$orig.ident
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype <- 
  plyr::mapvalues(orig_idents, names(genotype_mapping), genotype_mapping)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype <- 
  factor(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype, levels = c("Mock", "gSUV"))

# Create genotype_timepoint metadata column
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype_timepoint <- 
  paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype, 
        suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype_timepoint <- 
  factor(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype_timepoint, 
         levels = c("Mock D8", "gSUV D8", "Mock D28", "gSUV D28"))

# Extract cluster numbers by removing all descriptive suffixes
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- 
  gsub("-.*", "", Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))

#-----------------------------------------------------------------------------
# 2. UMAP Visualization of Cellular Clusters (Figure 1C)
#-----------------------------------------------------------------------------

# Rename clusters with more descriptive names
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) <- factor(x = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), levels = c(
  "C1-Active translation-IL7R", "C2-Stem/memory-KLF2", "C3-Stem/cytotoxic-GNLY", "C4-Resident-ITGA1", 
  "C5-Activated-CD69", "C6-Activated-JUN", "C7-Active translation-RPS27",
  "C8-Effector-KLRB1", "C9-Cytosolic stress-HSPA1B", "C10-Cytosolic stress-LncRNA", 
  "C11-Pro-apopototic-BAX", "C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"
))

# Plot UMAP with clusters
dimPlot_clusters <- DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
                          reduction = "umap", 
                          label = TRUE, 
                          label.size = 4, 
                          repel = TRUE,
                          pt.size = 0.8,
                          cols = mycolors.CD4) + 
  ggtitle("CD4+ CAR T cells") + 
  theme(plot.title = element_text(size = 14, face = "bold")) +
  NoLegend()

# Save the figure
pdf("Figure1C_CD4_UMAP_Clusters.pdf", width = 8, height = 7)
print(dimPlot_clusters)
dev.off()

#-----------------------------------------------------------------------------
# 3. Generate Signature Module Scores for T Cell States
#-----------------------------------------------------------------------------

# Define signature gene lists
# Load signature gene lists from files or define them directly
# Here we're defining a few as examples

# Stemness/Memory signature 
stem_memory_genes <- c("TCF7", "LEF1", "KLF2", "SELL", "IL7R", "CCR7", "CD27", "CD28", "FOXP1")

# Cytotoxic effector signature
cytotoxic_genes <- c("NKG7", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "IFNG", "PRF1")

# Tissue resident memory signature
resident_memory_genes <- c("ITGA1", "ITGAE", "CD69", "ZNF683", "CXCR6", "DUSP6", "RBPJ", "NOTCH1")

# Circulating signature
circulating_genes <- c("KLF2", "S1PR1", "SELL", "CCR7", "TCF7", "LEF1")

# Add module scores for each signature
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  features = list(stem_memory_genes), 
  name = "Stemness_Memory", 
  assay = "RNA"
)

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  features = list(cytotoxic_genes), 
  name = "Cytotoxic_Effector", 
  assay = "RNA"
)

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  features = list(resident_memory_genes), 
  name = "Resident_Memory", 
  assay = "RNA"
)

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  features = list(circulating_genes), 
  name = "Circulating", 
  assay = "RNA"
)

#-----------------------------------------------------------------------------
# 4. Generate Violin Plots for Signatures (Figure 1B)
#-----------------------------------------------------------------------------

# Create violin plots for each signature comparing genotypes and timepoints
vln_stemness <- VlnPlot(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  "Stemness_Memory1", 
  group.by = "genotype_timepoint", 
  pt.size = 0, 
  cols = c("darkgrey", "firebrick3", "black", "red")
) + 
  NoLegend() + 
  stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.6
  ) + 
  ggtitle("Stemness/memory") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

vln_circulating <- VlnPlot(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  "Circulating1", 
  group.by = "genotype_timepoint", 
  pt.size = 0, 
  cols = c("darkgrey", "firebrick3", "black", "red")
) + 
  NoLegend() + 
  stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.3
  ) + 
  ggtitle("Circulating") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

vln_cytotoxic <- VlnPlot(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  "Cytotoxic_Effector1", 
  group.by = "genotype_timepoint", 
  pt.size = 0, 
  cols = c("darkgrey", "firebrick3", "black", "red")
) + 
  NoLegend() + 
  stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 1.0
  ) + 
  ggtitle("Cytotoxic effectors") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

vln_resident <- VlnPlot(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  "Resident_Memory1", 
  group.by = "genotype_timepoint", 
  pt.size = 0, 
  cols = c("darkgrey", "firebrick3", "black", "red")
) + 
  NoLegend() + 
  stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.5
  ) + 
  ggtitle("Resident memory") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the violin plots into one figure
signatures_vlnplot <- vln_stemness + vln_circulating + vln_cytotoxic + vln_resident + 
  plot_layout(ncol = 4, widths = c(1, 1, 1, 1))

# Save the figure
pdf("Figure1B_CD4_Signatures_Violins.pdf", width = 12, height = 5)
print(signatures_vlnplot)
dev.off()

#-----------------------------------------------------------------------------
# 5. Differential Expression Analysis: Mock vs gSUV (Figure 1A)
#-----------------------------------------------------------------------------

# Find differentially expressed genes between Mock and gSUV
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype
degs_bulk <- FindMarkers(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  ident.1 = "gSUV", 
  ident.2 = "Mock", 
  logfc.threshold = 0.1, 
  test.use = "wilcox"
)

# Prepare volcano plot
colors <- rep("black", nrow(degs_bulk))
colors[which(degs_bulk$avg_log2FC > 0.25 & degs_bulk$p_val_adj < 0.001)] <- "salmon"
colors[which(degs_bulk$avg_log2FC < -0.25 & degs_bulk$p_val_adj < 0.001)] <- "grey50"
stem_genes <- c("KLF2", "S1PR1", "SELL", "KLF3", "LEF1", "TCF7", "IL7R", "CCR7")
colors[which(rownames(degs_bulk) %in% stem_genes & degs_bulk$avg_log2FC > 0)] <- "red"

# Create volcano plot
volcano_plot <- EnhancedVolcano(
  degs_bulk,
  lab = rownames(degs_bulk),
  labCol = colors,
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = TRUE,
  maxoverlapsConnectors = 18,
  pCutoff = 0.001,
  FCcutoff = 0.25,
  xlim = c(-1.2, 1.2),
  title = "CD4+ CAR T cells",
  subtitle = NULL,
  caption = NULL,
  legendLabSize = 10,
  legendPosition = "none"
)

# Save the volcano plot
pdf("Figure1A_CD4_Volcano_Mock_vs_gSUV.pdf", width = 8, height = 7)
print(volcano_plot)
dev.off()

#-----------------------------------------------------------------------------
# 6. Density Plots Showing Cell Distribution Changes (Figure 1E)
#-----------------------------------------------------------------------------

# Filter data by timepoint for density plots
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "umap"))
df <- cbind(df, 
            suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, 
            suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")

# Day 8 density plots
df.d8 <- df %>% filter(timepoint == "D8")
d8_cell_counts <- table(df.d8$genotype)
d8_labels <- paste0(names(d8_cell_counts), " - ", d8_cell_counts, " cells")
names(d8_labels) <- names(d8_cell_counts)

d8_density_plot <- ggplot(data = df.d8) +
  geom_point(aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), 
                 bins = 20, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  facet_grid(. ~ genotype, labeller = labeller(genotype = d8_labels)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black", face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  NoLegend() +
  annotate("text", x = 0, y = -5, label = "Day 8", size = 5)

# Day 28 density plots
df.d28 <- df %>% filter(timepoint == "D28")
d28_cell_counts <- table(df.d28$genotype)
d28_labels <- paste0(names(d28_cell_counts), " - ", d28_cell_counts, " cells")
names(d28_labels) <- names(d28_cell_counts)

d28_density_plot <- ggplot(data = df.d28) +
  geom_point(aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), 
                 bins = 30, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  facet_grid(. ~ genotype, labeller = labeller(genotype = d28_labels)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black", face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  NoLegend() +
  annotate("text", x = 0, y = -5, label = "Day 28", size = 5)

# Combine density plots
density_plots <- d8_density_plot / d28_density_plot

# Save the combined density plots
pdf("Figure1E_CD4_Density_Plots.pdf", width = 10, height = 10)
print(density_plots)
dev.off()

#-----------------------------------------------------------------------------
# 7. Fold Change in Cluster Proportions (Figure 1F)
#-----------------------------------------------------------------------------

# Calculate fold change for D8
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  subset = timepoint == "D8"
)
query.list.D8 <- SplitObject(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, split.by = "genotype")
norm.c.D8 <- table(Idents(query.list.D8[["Mock"]]))/sum(table(Idents(query.list.D8[["Mock"]])))
norm.q.D8 <- table(Idents(query.list.D8[["gSUV"]]))/sum(table(Idents(query.list.D8[["gSUV"]])))
foldchange.D8 <- norm.q.D8/norm.c.D8
foldchange.D8 <- sort(foldchange.D8, decreasing = TRUE)
tb.m.D8 <- melt(foldchange.D8)
colnames(tb.m.D8) <- c("Cluster", "Fold_change")

# Calculate fold change for D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  subset = timepoint == "D28"
)
query.list.D28 <- SplitObject(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, split.by = "genotype")
norm.c.D28 <- table(Idents(query.list.D28[["Mock"]]))/sum(table(Idents(query.list.D28[["Mock"]])))
norm.q.D28 <- table(Idents(query.list.D28[["gSUV"]]))/sum(table(Idents(query.list.D28[["gSUV"]])))
foldchange.D28 <- norm.q.D28/norm.c.D28
foldchange.D28 <- sort(foldchange.D28, decreasing = TRUE)
tb.m.D28 <- melt(foldchange.D28)
colnames(tb.m.D28) <- c("Cluster", "Fold_change")

# Create barplots for fold changes
fc_plot_D8 <- ggplot(tb.m.D8, aes(x = Cluster, y = Fold_change, fill = Cluster)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_fill_manual(values = mycolors.CD4) + 
  geom_hline(yintercept = 1) + 
  scale_y_continuous(trans = 'log2', limits = c(0.6, 4.3), breaks = seq(-20, 20, by = 1)) + 
  theme(axis.text.x = element_blank(), legend.position = "none") + 
  ggtitle("CD4 - gSUV vs. mock - D8") +
  labs(y = "Fold change")

fc_plot_D28 <- ggplot(tb.m.D28, aes(x = Cluster, y = Fold_change, fill = Cluster)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_fill_manual(values = mycolors.CD4) + 
  geom_hline(yintercept = 1) + 
  scale_y_continuous(trans = 'log2', limits = c(0.6, 4.3), breaks = seq(-20, 20, by = 1)) + 
  theme(axis.text.x = element_blank(), legend.position = "none") + 
  ggtitle("CD4 - gSUV vs. mock - D28") +
  labs(y = "Fold change")

# Combine fold change plots
fc_plots <- fc_plot_D8 / fc_plot_D28

# Save the fold change plots
pdf("Figure1F_CD4_Cluster_Fold_Change.pdf", width = 8, height = 10)
print(fc_plots)
dev.off()

#-----------------------------------------------------------------------------
# 8. Feature Plots for Signature Scores (Figure 1D)
#-----------------------------------------------------------------------------

# Create feature plots for each signature
fp_stemness <- FeaturePlot(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  "Stemness_Memory1", 
  pt.size = 0.5, 
  min.cutoff = "q3", 
  max.cutoff = "q97", 
  order = TRUE
) + 
  scale_color_paletteer_c("pals::coolwarm") + 
  NoAxes() +
  ggtitle("Stemness/memory")

fp_cytotoxic <- FeaturePlot(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  "Cytotoxic_Effector1", 
  pt.size = 0.5, 
  min.cutoff = "q3", 
  max.cutoff = "q97", 
  order = TRUE
) + 
  scale_color_paletteer_c("pals::coolwarm") + 
  NoAxes() +
  ggtitle("Cytotoxic effector")

fp_resident <- FeaturePlot(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, 
  "Resident_Memory1", 
  pt.size = 0.5, 
  min.cutoff = "q3", 
  max.cutoff = "q97", 
  order = TRUE
) + 
  scale_color_paletteer_c("pals::coolwarm") + 
  NoAxes() +
  ggtitle("Resident memory Milner et al. 2018")

fp_circulating <- FeaturePlot(
  suv.car.t.integrate
