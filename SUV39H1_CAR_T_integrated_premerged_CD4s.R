library(Seurat)
library(ggplot2)
library(sctransform)
library(patchwork)
library(pheatmap)
library(Scillus)
library(pals)
library(RColorBrewer)
library(tidyverse)
library(paletteer)
library(EnhancedVolcano)
library(readxl)
library(ggsignif)
library(harmony)
library(Nebulosa)
library(data.table)
library(gridExtra)
library(CytoTRACE)

library(future)
plan("multiprocess", workers = 16, strategy = "multicore")
options(future.globals.maxSize = 1500 * 2024^2)

# Signatures
sign_term_trem2 <- toupper(readLines("E:/Work/Signatures/Resident Memory/Synder_et_al_2019/sign_Term_Trem_Synder_bis.txt"))
sign_memory <- toupper(read.csv("E:/Work/Signatures/Stemness/signature_memory_stemness.csv", h = F)[,1])
sign_senescence <- toupper(read.csv("E:/Work/Signatures/Senescence/sign_sag.txt", h = F)[,1])
sign_apoptosis <- toupper(read.csv("E:/Work/Signatures/Signatures Papier Leticia/Apoptosis.txt", h = F)[,1])
sign_DNA_repair <- toupper(read.csv("E:/Work/Signatures/DNA repair/sign_DNA_repair.txt", h = F)[,1])
sign_core_DNA_repair <- toupper(read.csv("E:/Work/Signatures/DNA repair/Core_DNA_damage_repair.txt", h = F)[,1])
sign_DNA_damage_repair <- toupper(read.csv("E:/Work/Signatures/DNA repair/DNA_damage_repair_Knijnenburg.txt", h = F)[,1])

# Colors setup
mycolors <- pals::kelly(18)
mycolors.CD4 <- mycolors[5:18]

# Saving
save(suv.car.t.integrated.merged.timepoints.premerged.cd4s, file = "E:/Work/A_R_files/SUV39H1_CAR_T_integrated_premerged_CD4s.Rdata")
save(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, file = "E:/Work/A_R_files/SUV39H1_CAR_T_integrated_premerged_CD4s_cleaned.Rdata")

# Merging CD4s only
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- subset(suv.car.t.integrated.merged.timepoints.premerged, subset = cd4A > 0 & cd4B > 0, invert = T)
TCR.removed <- rownames(suv.car.t.integrated.merged.timepoints.premerged.cd4s)[which(rownames(suv.car.t.integrated.merged.timepoints.premerged.cd4s) %!in% TCR.features)]
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- RunPCA(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, npcs = 50, verbose = FALSE, features = TCR.removed)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- RunHarmony(suv.car.t.integrated.merged.timepoints.premerged.cd4s, group.by.vars = c("orig.ident", "timepoint"), project.dim = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- RunUMAP(suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "harmony", dims = 1:50)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- FindNeighbors(suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "harmony", dims = 1:50)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- FindClusters(suv.car.t.integrated.merged.timepoints.premerged.cd4s, graph.name = "RNA_snn", resolution = 0.8)

DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "umap", label = F, pt.size = 1, label.size = 4, cols = mycolors.CD4) + NoLegend() +
  DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "umap", group.by = "genotype", label = F, pt.size = 1, cols = c("black", "red"), label.size = 6) +
  DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "umap", group.by = "timepoint", label = F, pt.size = 1, label.size = 6) + NoLegend()

# Densities with cells as grey points in the background
umap <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "umap"))
ggplot() +
  geom_point(data = umap, aes(x = UMAP_1, y = UMAP_2), color = "grey")

# Densities by genotype and timepoint
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("Mock", "KO"))
ggplot(data = df) +
  geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 40, geom = "polygon") +
  theme_classic() +
  xlim(c(-8, 10)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  facet_grid(. ~ timepoint + genotype) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  )

# Densities by genotype and timepoint
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("Mock", "gSUV"))
ggplot(data = df) +
  geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 40, geom = "polygon") +
  theme_classic() +
  xlim(c(-8, 10)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  facet_grid(. ~ timepoint + genotype) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  )

# Test densities with unique scale for each
df.a <- df %>% filter(timepoint == "D8" & genotype == "Mock")
a <- plot(ggplot(data = df.a) +
  geom_point(data = df.a, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
  theme_classic() +
  xlim(c(-8, 10)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
df.b <- df %>% filter(timepoint == "D8" & genotype == "gSUV")
b <- plot(ggplot(data = df.b) +
  geom_point(data = df.b, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
  theme_classic() +
  xlim(c(-8, 10)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
df.c <- df %>% filter(timepoint == "D28" & genotype == "Mock")
c <- plot(ggplot(data = df.c) +
  geom_point(data = df.c, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
  theme_classic() +
  xlim(c(-8, 10)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
df.d <- df %>% filter(timepoint == "D28" & genotype == "gSUV")
d <- plot(ggplot(data = df.d) +
  geom_point(data = df.d, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
  theme_classic() +
  xlim(c(-8, 10)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
a + b + c + d + plot_layout(ncol = 4) & NoLegend()

# TCR features
"%!in%" <- function(x, y) !("%in%"(x, y))
TCRA <- read.csv("E:/Work/Signatures/TCR_genes/TCRA.txt", sep = "\t")
TCRA.features <- TCRA$Approved.symbol
TCRB <- read.csv("E:/Work/Signatures/TCR_genes/TCRB.txt", sep = "\t")
TCRB.features <- TCRB$Approved.symbol
TCRD <- read.csv("E:/Work/Signatures/TCR_genes/TCRD.txt", sep = "\t")
TCRD.features <- TCRD$Approved.symbol
TCRG <- read.csv("E:/Work/Signatures/TCR_genes/TCRG.txt", sep = "\t")
TCRG.features <- TCRG$Approved.symbol
TCR.features <- c(TCRA.features, TCRB.features, TCRG.features, TCRD.features)
TCR.removed <- rownames(suv.car.t.integrated.merged.timepoints.premerged.cd4s)[which(rownames(suv.car.t.integrated.merged.timepoints.premerged.cd4s) %!in% TCR.features)]

# Heatmap
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- factor(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s), levels = order)
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- factor(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s), levels = c(0:11))
markers <- FindAllMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, assay = "RNA", features = TCR.removed, only.pos = T)
markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 8, wt = avg_log2FC) -> top8
Average <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = unique(top8$gene), return.seurat = T, assays = "RNA")
Average <- ScaleData(object = Average, features = unique(top8$gene))
pdf("Heatmap_average_d8_d28_harmony_sample_CD4s_subcluster.pdf", width = 9, height = 20)
plot(DoHeatmap(Average, features = unique(top8$gene), draw.lines = F, size = 3, group.colors = mycolors.CD4) + scale_fill_viridis_c(option = "inferno"))
dev.off()

Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s$RNA_snn_res.0.2, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint)
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s$RNA_snn_res.0.2, suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype)
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s$RNA_snn_res.0.2

# Saving DEGs
write.csv(markers, "markers_CD4s.txt", quote = F, row.names = F)

# Transfer with NSCLC data
features <- SelectIntegrationFeatures(object.list = c(eleven.tils.cd3.integrated, suv.car.t.integrated.merged.timepoints.premerged.cd4s), nfeatures = 6000)
tumor.anchors <- FindTransferAnchors(
  reference = eleven.tils.cd3.integrated, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s,
  dims = 1:30, query.assay = "RNA", reference.assay = "integrated"
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(eleven.tils.cd3.integrated),
  dims = 1:30
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, group.by = "predicted.id", label = T, repel = T, reduction = "umap", cols = paletteer_d("pals::polychrome", n = 21))
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "prediction.score.max")

# Cell cycle scoring
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- CellCycleScoring(suv.car.t.integrated.merged.timepoints.premerged.cd4s,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes, set.ident = F
)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, group.by = "Phase", label = T, label.size = 6, pt.size = 1.1)

# Nebulosa plotting
p <- plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("TCF7", "IL7R", "SELL", "KLF2"), joint = T, )
p + plot_layout(ncol = 3)

# CytoTRACE
library(CytoTRACE)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, cells = sample(Cells(suv.car.t.integrated.merged.timepoints.premerged.cd4s), 20000))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = Phase == "G1")
suv.car.t.integrated.merged.timepoints.premerged.cd4s[["ident_genotype_timepoint"]] <- paste(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s$RNA_snn_res.0.2, suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype,
  suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint
)

mat <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, assay = "RNA")
mat <- as.data.frame(mat)
results.cyto <- CytoTRACE(mat, ncores = 1)
idents <- as.character(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub$ident_genotype_timepoint)
names(idents) <- names(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub$ident_genotype_timepoint)
emb <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, reduction = "umap"))
plotCytoTRACE(results.cyto, emb = emb, phenotype = idents)

# Monocle3
library(monocle3)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, idents = c("Mono-S100A8", "Mono-FCGR3A", "Mono-TREM1", "Early MAC-ISG15", "Early MAC-CXCR4", "MAC-CXCL2", "MAC-FBP1", "Mo-derived LAM-STAB1", "LAM-APOC1"))
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- RunUMAP(suv.car.t.integrated.merged.timepoints.premerged.cd4s,
  reduction = "harmony",
  dims = 1:50, reduction.name = "UMAP"
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D8" & genotype == "WT")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D28" & genotype == "WT")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D8" & genotype == "KO")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D28" & genotype == "KO")

cds <- as.cell_data_set(suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO.D8)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds@clusters$UMAP$clusters <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO.D8)
cds <- learn_graph(cds, close_loop = F, use_partition = T, learn_graph_control = list(minimal_branch_len = 13))
p <- plot_cells(cds, cell_size = 1, label_cell_groups = F, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = T, trajectory_graph_color = "black", group_label_size = 5, trajectory_graph_segment_size = 1.5, ) + ggtitle("WT D8")
p2 <- plot_cells(cds, cell_size = 1, label_cell_groups = F, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = T, trajectory_graph_color = "black", group_label_size = 5, trajectory_graph_segment_size = 1.5, ) + ggtitle("KO D8")
p3 <- plot_cells(cds, cell_size = 1, label_cell_groups = F, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = T, trajectory_graph_color = "black", group_label_size = 5, trajectory_graph_segment_size = 1.5, ) + ggtitle("WT D28")
p4 <- plot_cells(cds, cell_size = 1, label_cell_groups = F, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = T, trajectory_graph_color = "black", group_label_size = 5, trajectory_graph_segment_size = 1.5, ) + ggtitle("KO D28")
p + p2 + p3 + p4

# Single genes nebulosa
library("Nebulosa")
plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("TOX", "PDCD1", "HAVCR2", "LAG3"))
plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("TCF7", "SELL", "IL7R", "KLF2"))
plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("NFE2L1", "NFE2L3", "JUN"))
plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("PRF1", "IFNG", "TNF", "ITGA1", "ITGA4", "ITGB1", "S1PR1", "CD69", "CTLA4", "NRP1"))

FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("ITGA5", "ITGA6", "ITGAX", "ITGB1", "ITGB2", "ITGAM"), pt.size = 0.8, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm")
plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("S1PR1", "KLF2", "CCR4", "CCR7", "ITGAE", "ITGA4"))

# ReactomeGSA
library(ReactomeGSA)
gsva_result <- analyse_sc_clusters(suv.car.t.integrated.merged.timepoints.premerged.cd4s, verbose = TRUE)
plot_gsva_heatmap(gsva_result, max_pathways = 20, margins = c(6, 25))

# Plot volcano KO vs WT in batch for each cluster
library(EnhancedVolcano)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- NormalizeData(suv.car.t.integrated.merged.timepoints.premerged.cd4s) %>%
  FindVariableFeatures() %>%
  ScaleData()
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s$RNA_snn_res.0.2, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint)
pdf("volcano_CD4_KO_vs_WT_bycluster_timepoint.pdf", width = 9, height = 9)
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s))) {
  degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, ident.1 = "KO", ident.2 = "WT", subset.ident = i, group.by = "genotype", assay = "RNA"))
  try(plot(EnhancedVolcano(degs,
    lab = rownames(degs),
    x = "avg_log2FC",
    y = "p_val_adj",
    drawConnectors = T,
    pCutoff = 0.001,
    FCcutoff = 0.3
  ) + ggtitle(paste("cluster", i, sep = " "))))
}
dev.off()

# Number of DEGs
suv.car.t.integrated.merged.timepoints.premerged.cd4s[["label"]] <- suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype
suv.car.t.integrated.merged.timepoints.premerged.cd4s[["cell_type"]] <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s)
results <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s))) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, ident.1 = "WT", ident.2 = "KO", subset.ident = i, group.by = "genotype", assay = "RNA"))
  results[[i]] <- markers
}
lengths.deg <- c()
for (i in 1:length(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s))) {
  lengths.deg[i] <- print(length(rownames(results[[i]])))
}
names(lengths.deg) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s))
which(names(lengths.deg) %in% augur$AUC$cell_type)
lengths.deg <- as.data.frame(lengths.deg)
lengths.deg[["cell_type"]] <- rownames(lengths.deg)
ggplot(data = lengths.deg, aes(x = cell_type, y = lengths.deg)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_light() +
  RotatedAxis()

# Number of DEGs split for D8
suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D8")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8[["label"]] <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8$genotype
suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8[["cell_type"]] <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8)
results <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8))) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8, ident.1 = "WT", ident.2 = "KO", subset.ident = i, group.by = "genotype", assay = "RNA"))
  results[[i]] <- markers
}
lengths.deg <- c()
for (i in 1:length(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8))) {
  lengths.deg[i] <- print(length(rownames(results[[i]])))
}
names(lengths.deg) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8))
which(names(lengths.deg) %in% augur$AUC$cell_type)
lengths.deg <- as.data.frame(lengths.deg)
lengths.deg[["cell_type"]] <- rownames(lengths.deg)
ggplot(data = lengths.deg, aes(x = cell_type, y = lengths.deg)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_light() +
  RotatedAxis()

# Proportions
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype_timepoint <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint)
suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype_timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype_timepoint, levels = c("Mock D8", "KO D8", "Mock D28", "KO D28"))
suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype, levels = c("WT", "KO"))
suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint, levels = c("D8", "D28"))
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s, plot_type = "prop_fill", group_by = "Phase") + scale_fill_manual(values = mycolors.CD4)

# VP by timepoint / genotype
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4"), pt.size = 0.1, group.by = "genotype_timepoint") + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
  geom_signif(
    comparisons = list(c("WT D8", "KO D8"), c("WT D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 2
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)

VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("TCF7", "SELL", "LEF1", "KLF2", "CCR7", "IL7R"), pt.size = 0.1, group.by = "genotype_timepoint") + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
  geom_signif(
    comparisons = list(c("WT D8", "KO D8"), c("WT D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 2.5
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)


suv.car.t.integrated.merged.timepoints.premerged.cd4s <- DietSeurat(suv.car.t.integrated.merged.timepoints.premerged.cd4s)

# Label transfer of cycling
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- FindVariableFeatures(suv.car.t.integrated.merged.timepoints.premerged.cd4s, nfeatures = 4000)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = Phase == "G1")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref, ident = c("C11-Cycling-CENPF", "C12-Cycling-TOP2A", "C13-Cycling-ATP5F1E"), invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = Phase == "G1", invert = T)
tumor.anchors <- FindTransferAnchors(
  reference = suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query,
  dims = 1:30, reference.assay = "RNA", query.assay = "RNA"
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref),
  dims = 1:30
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, group.by = "predicted.id", label = T, repel = T)
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, "prediction.score.max", pt.size = 0)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$seurat_clusters <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$predicted.id
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$seurat_clusters <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$seurat_clusters, levels = c(0:10))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$genotype <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$genotype, levels = c("WT", "KO"))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$timepoint, levels = c("D8", "D28"))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$genotype_timepoint <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$genotype, suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$timepoint)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$genotype_timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$genotype_timepoint, levels = c("WT D8", "KO D8", "WT D28", "KO D28"))

plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, plot_type = "prop_fill", group_by = "genotype_timepoint") + RotatedAxis() + scale_fill_manual(values = gg_color_hue(7))
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, plot_type = "prop_fill", group_by = "Phase") + RotatedAxis() + scale_fill_manual(values = gg_color_hue(7))

# Confirming that label transfer did work well
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, c("KLF2", "TCF7", "NEAT1", "TMSB10", "IL7R", "SELL"), pt.size = 0, group.by = "predicted.id")

# Boxplots
freq_table <- as.data.frame(prop.table(
  x = table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s), suv.car.t.integrated.merged.timepoints.premerged.cd4s$orig.ident),
  margin = 2
))
colnames(freq_table) <- c("cluster", "replicate", "Freq")
freq_table$condition[which(freq_table$replicate == "Mock.1")] <- "WT"
freq_table$condition[which(freq_table$replicate == "Mock.2")] <- "WT"
freq_table$condition[which(freq_table$replicate == "Mock.3")] <- "WT"
freq_table$condition[which(freq_table$replicate == "Suv.ko.1")] <- "KO"
freq_table$condition[which(freq_table$replicate == "Suv.ko.2")] <- "KO"
freq_table$condition[which(freq_table$replicate == "Suv.ko.3")] <- "KO"
freq_table$condition <- factor(freq_table$condition, levels = c("WT", "KO"))
freq_table$cluster <- factor(freq_table$cluster, levels = order)

pdf("boxplots_suvcart_CD4_WT_vs_KO.pdf", width = 10, height = 5)
ggplot(data = freq_table, aes(x = condition, y = Freq, fill = condition)) +
  geom_boxplot() +
  facet_grid(~cluster) +
  theme_bw() +
  geom_point(pch = 21, position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = c("deepskyblue", "orange1")) +
  theme(strip.text.x = element_text(size = 8, face = "bold", angle = 45)) +
  geom_signif(
    comparisons = list(c("WT", "KO")),
    map_signif_level = F, textsize = 4, test = "wilcox.test", y_position = 0.23
  )
dev.off()

order <- paste(rep(0:11, each = 2), rep(c("WT", "KO"), each = 1), rep(c("D8", "D28"), each = 12))
order <- paste(rep(0:11, each = 2), rep(c("D8", "D28"), each = 12))
order <- paste(rep(0:11, each = 2), rep(c("D8", "D28"), each = 1))
order <- paste(rep(0:11, each = 2), rep(c("WT", "KO"), each = 1))

# Boxplots split by timepoint
suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D8")
freq_table <- as.data.frame(prop.table(
  x = table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8), suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8$orig.ident),
  margin = 2
))
colnames(freq_table) <- c("cluster", "replicate", "Freq")
freq_table$condition[which(freq_table$replicate == "Mock.1")] <- "WT"
freq_table$condition[which(freq_table$replicate == "Mock.2")] <- "WT"
freq_table$condition[which(freq_table$replicate == "Mock.3")] <- "WT"
freq_table$condition[which(freq_table$replicate == "Suv.ko.1")] <- "KO"
freq_table$condition[which(freq_table$replicate == "Suv.ko.2")] <- "KO"
freq_table$condition[which(freq_table$replicate == "Suv.ko.3")] <- "KO"
freq_table$condition <- factor(freq_table$condition, levels = c("WT", "KO"))
order <- paste(c(0:11), "D8")
freq_table$cluster <- factor(freq_table$cluster, levels = order)

pdf("boxplots_suvcart_cd4_WT_vs_KO.pdf", width = 10, height = 5)
ggplot(data = freq_table, aes(x = condition, y = Freq, fill = condition)) +
  geom_boxplot() +
  facet_grid(~cluster) +
  theme_bw() +
  geom_point(pch = 21, position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = c("deepskyblue", "orange1")) +
  theme(strip.text.x = element_text(size = 8, face = "bold", angle = 45)) +
  geom_signif(
    comparisons = list(c("WT", "KO")),
    map_signif_level = F, textsize = 4, test = "wilcox.test", test.args = c(alternative = "two.sided"), y_position = 0.26
  )
dev.off()

# FP for QC
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("CD3D", "cd4A", "cd4B", "CD4"), pt.size = 0.8, order = F, max.cutoff = "q97") & scale_color_paletteer_c("pals::coolwarm") & NoAxes()
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.8, order = F, max.cutoff = "q97") & scale_color_paletteer_c("pals::coolwarm") & NoAxes()

# Splitting cluster 0
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- FindSubCluster(suv.car.t.integrated.merged.timepoints.premerged.cd4s, cluster = "0", resolution = 0.1, graph.name = "RNA_snn")
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, group.by = "sub.cluster", reduction = "umap", label = TRUE, pt.size = 1.3, label.size = 6) + ggtitle("CD4s / Day 8 + 28 / n = 49728")

# Plotting signature scores
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}

# As Featureplots
pdf("Signatures_CD4s_Fraietta.pdf", width = 12, height = 12)
for (i in 1:length(all_sign2)) {
  try(plot(FeaturePlot(
    object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = all_sign2[i], order = T, reduction = "umap",
    pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
  ) + labs(title = all_sign[i], size = 13) + scale_color_paletteer_c("pals::coolwarm")))
}
dev.off()

# Comparison with Chen et al. data
Chen2021_long_response <- readLines("E:/Work/Signatures/Chen_et_al_2021/Long_response_all.txt")
Chen2021_long_response_top50_FC <- readLines("E:/Work/Signatures/Chen_et_al_2021/Long_response_all_top50.txt")
Chen2021_short_response <- readLines("E:/Work/Signatures/Chen_et_al_2021/short_response_all.txt")
fraietta_cr <- readLines("E:/Work/Signatures/Fraietta_et_al_2018/Fraietta_CR.txt")
fraietta_nr <- readLines("E:/Work/Signatures/Fraietta_et_al_2018/Fraietta_NR.txt")

# As violinplots between genotype and time point
all_sign <- c("Chen2021_long_response", "Chen2021_long_response_top50_FC", "Chen2021_short_response")
all_sign <- c("fraietta_cr", "fraietta_nr")
all_sign <- c("sign_ox_ph", "sign_glycolysis", "sign_FA_metabolism_GO", "sign_fatty_acid_metab_KEGG")

pdf("VP_car_t_CD4_Fra_full.pdf", width = 5, height = 6)
for (i in 1:length(all_sign2)) {
  try(plot(VlnPlot(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D8"), "fraietta_nr1", group.by = "genotype", pt.size = 0.1) + labs(title = "Fraietta NR", size = 13) + ylim(-0.5, 0.2) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
    geom_signif(
      comparisons = list(c("WT", "KO")),
      map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.1
    ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)))
}
dev.off()

# FP
FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "fraietta_nr1", split.by = "genotype_timepoint", order = T, reduction = "umap",
  pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")


Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint)
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- factor(x = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s), levels = c("WT D8", "KO D8", "WT D28", "KO D28"))

Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- paste(
  suv.car.t.integrated.merged.timepoints.premerged.cd4s$RNA_snn_res.0.2,
  suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint
)
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- factor(x = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s), levels = order)
order <- paste(rep(0:11, each = 4), rep(c("WT", "KO"), each = 1), rep(c("D8", "D28"), each = 2))

# muscat test
library(muscat)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, ident = c("C11-Cycling-CENPF", "C12-Cycling-TOP2A", "C13-Cycling-ATP5F1E"), invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref$seurat_clusters <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref)
sce <- as.SingleCellExperiment(suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref, assay = "RNA")
sce$id <- paste(sce$timepoint, sce$orig.ident, sep = "_")
muscat.sce <- prepSCE(sce,
  kid = "seurat_clusters", # I have no cluter id
  gid = "genotype_timepoint", # group IDs (ctrl/stim)
  sid = "id", drop = T
)
pb <- aggregateData(muscat.sce,
  assay = "counts", fun = "sum",
  by = c("cluster_id", "sample_id")
)

(pb_mds <- pbMDS(pb))

# sc T cell naive/stemness signature
cd4_stem_sc <- readLines("E:/Work/Signatures/Stemness/cd4_Stem_like_custom_single_cell.txt")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(cd4_stem_sc), name = "cd4_stem_sc", assay = "RNA")

# As Featureplots
FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "CD8_stem_sc1", split.by = "genotype", order = T, reduction = "umap",
  pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")

plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s, "CD8_stem_sc1")

# sc T cell naive/stemness signature
cd4_naive_goldrath <- toupper(readLines("E:/Work/Signatures/Naive/GOLDRATH_NAIVE_VS_EFF_cd4_TCELL_UP.txt"))
cd4_naive_Guo <- readLines("E:/Work/Signatures/Guo_et_al_2018/cd4_LEF1.txt")
CD4_naive_Guo <- readLines("E:/Work/Signatures/Guo_et_al_2018/CD4_CCR7.txt")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(cd4_naive_goldrath), name = "cd4_naive_goldrath", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(cd4_naive_Guo), name = "cd4_naive_Guo", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(CD4_naive_Guo), name = "CD4_naive_Guo", assay = "RNA")

# As Featureplots
p <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "cd4_naive_goldrath1", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "cd4+ naive goldrath") + scale_color_paletteer_c("pals::coolwarm")
p1 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "cd4_naive_Guo1", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "cd4_naive_Guo") + scale_color_paletteer_c("pals::coolwarm")
p2 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "CD4_naive_Guo1", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "CD4_naive_Guo") + scale_color_paletteer_c("pals::coolwarm")
p + p1 + p2 + plot_layout(ncol = 3)

# Fraietta et al. 2018 FGI2C signatures
Fraietta.2C <- readxl::read_xlsx(path = "E:/Work/Signatures/Fraietta_et_al_2018/Fig2C_41591_2018_10_MOESM5_ESM.xlsx")
colnames(Fraietta.2C) <- gsub("[\r\n]", "", colnames(Fraietta.2C))
colnames(Fraietta.2C) <- make.names(colnames(Fraietta.2C))
pdf("Fraietta.2C_SUV_CD4.pdf", height = 10, width = 11)
for (i in colnames(Fraietta.2C)) {
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Fraietta.2C[[i]]), name = i, assay = "RNA")
  # plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, paste0(i, 1), pt.size = 0.8, min.cutoff = "q3", max.cutoff = "q97", order = T) + scale_color_paletteer_c("pals::coolwarm") +ggtitle(i))
}
dev.off()

# Heatmap with signatures
colnames(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data) <- gsub(pattern = 1, replacement = "", x = colnames(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned[["Fraietta"]] <- CreateAssayObject(data = t(x = FetchData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, vars = colnames(Fraietta.2C))))
Average <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = colnames(Fraietta.2C), assays = "Fraietta", slot = "data", return.seurat = T)
Average <- ScaleData(object = Average, features = colnames(Fraietta.2C))
pdf("Heatmap_average_Fraietta_CD4s.pdf", width = 9, height = 12)
plot(DoHeatmap(Average, features = colnames(Fraietta.2C), draw.lines = F, size = 3) + scale_fill_viridis_c(option = "inferno"))
dev.off()

# Reorder idents
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$reordered.idents <- paste(
  Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned),
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint
)
order <- paste(rep(c(
  "C1-Active translation-IL7R", "C2-Stem/memory-KLF2", "C3-Stem/cytotoxic-GNLY", "C4-Resident-ITGA1", "C5-Activated-CD69", "C6-Activated-JUN", "C7-Active translation-RPS27",
  "C8-Effector-KLRB1", "C9-Cytosolic stress-HSPA1B", "C10-Cytosolic stress-LncRNA", "C11-Pro-apopototic-BAX", "C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"
), each = 1), rep(c("Mock", "gSUV"), each = 14), rep(c("D8", "D28"), each = 28))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$reordered.idents <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$reordered.idents, levels = order)

# Redo same with only cluster 2-memory-KLF2
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2$reordered.idents
colnames(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2@meta.data) <- gsub(pattern = 1, replacement = "", x = colnames(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2@meta.data))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2[["Fraietta"]] <- CreateAssayObject(data = t(x = FetchData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, vars = colnames(Fraietta.2C))))
Average <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, features = colnames(Fraietta.2C), assays = "Fraietta", slot = "data", return.seurat = T)
Average <- ScaleData(object = Average, features = colnames(Fraietta.2C))
plot(DoHeatmap(Average, features = colnames(Fraietta.2C), draw.lines = F, size = 4) + scale_fill_viridis_c(option = "inferno") + theme(axis.text.y = element_text(size = 12)))

# SignacX classification
library(SignacX)
labels <- Signac(suv.car.t.integrated.merged.timepoints.premerged.cd4s, num.cores = 14)
celltypes <- GenerateLabels(labels, E = suv.car.t.integrated.merged.timepoints.premerged.cd4s)

suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddMetaData(suv.car.t.integrated.merged.timepoints.premerged.cd4s, metadata = celltypes$CellStates_novel, col.name = "CellStates_novel")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- SetIdent(suv.car.t.integrated.merged.timepoints.premerged.cd4s, value = "CellStates_novel")
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, pt.size = 0.9, label = T)

suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddMetaData(suv.car.t.integrated.merged.timepoints.premerged.cd4s, metadata = celltypes$L8, col.name = "L8")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- SetIdent(suv.car.t.integrated.merged.timepoints.premerged.cd4s, value = "L18")
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, pt.size = 0.9, group.by = "L8")

# Label transfer Guo et al.
tumor.anchors <- FindTransferAnchors(
  reference = guo.seurat, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s,
  dims = 1:30
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = guo.seurat$majorCluster,
  dims = 1:30
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, group.by = "predicted.id", label = T)


# RoE heatmap
library("pheatmap")
suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype_timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype_timepoint, levels = c("Mock D8", "KO D8", "Mock D28", "KO D28"))
clusters.table <- table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s), suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype_timepoint)
chi.prop <- chisq.test(clusters.table)
Roe <- chi.prop$observed / chi.prop$expected
pheatmap(Roe, cluster_rows = F, cluster_cols = F, scale = "column", display_numbers = F, fontsize = 14)

# Roe heatmap split by D8 and D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype, levels = c("Mock", "KO"))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = timepoint == "D8")
clusters.table <- table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8), suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8$genotype)
chi.prop <- chisq.test(clusters.table)
Roe <- chi.prop$observed / chi.prop$expected
pheatmap(Roe, cluster_rows = F, cluster_cols = F, scale = "none", display_numbers = F, fontsize = 14)

# DEenrichR plot
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("KEGG_2021_Human", "WikiPathway_2021_Human", "MSigDB_Hallmark_2020", "GO_Biological_Process_2021", "Elsevier_Pathway_Collection", "HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression")
dbs <- c("MSigDB_Hallmark_2020")
pdf(file = "Pathways_CD4s_DEenrichR_KEGG_2021_Human.pdf", width = 9, height = 8)
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s))) {
  plot(DEenrichRPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, ident.1 = "C5-Memory-KLF2", logfc.threshold = 0.2, balanced = T, assay = "RNA", enrich.database = dbs, max.genes = 300))
}
dev.off()

# DEenrichR plot for WT vs KO
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("KEGG_2021_Human", "WikiPathway_2021_Human", "MSigDB_Hallmark_2020", "GO_Biological_Process_2021", "GO_Molecular_Function_2021", "Elsevier_Pathway_Collection", "HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression")

# ROGUE index by genotype and by cluster
library(ROGUE)
expr <- suv.car.t.integrated.merged.timepoints.premerged.cd4s@assays$RNA@counts
expr <- as.matrix(expr)
rogue.res <- rogue(expr,
  labels = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s),
  samples = suv.car.t.integrated.merged.timepoints.premerged.cd4s$orig.ident,
  platform = "UMI",
  filter = T
)
rogue.plot <- rogue.boxplot(rogue.res)
rogue.plot + RotatedAxis()

# Comparison with Leticia's data
new.symbols <- toupper(leti.integrated.sub3@assays$RNA@counts@Dimnames[[1]])
leti.integrated.sub3@assays$RNA@counts@Dimnames[[1]] <- new.symbols
new.symbols <- toupper(leti.integrated.sub3@assays$RNA@data@Dimnames[[1]])
leti.integrated.sub3@assays$RNA@data@Dimnames[[1]] <- new.symbols
new.symbols <- toupper(leti.integrated.sub3@assays$integrated@data@Dimnames[[1]])
leti.integrated.sub3@assays$integrated@data@Dimnames[[1]] <- new.symbols
Hum_gene_assay <- CreateAssayObject(GetAssayData(leti.integrated.sub3, assay = "integrated"))
leti.integrated.sub3[["hum_gene_assay"]] <- Hum_gene_assay

features <- SelectIntegrationFeatures(object.list = c(leti.integrated.sub3, suv.car.t.integrated.merged.timepoints.premerged.cd4s), features = 5000)
tumor.anchors <- FindTransferAnchors(
  reference = leti.integrated.sub3, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s,
  dims = 1:50, reference.assay = "RNA", query.assay = "RNA", reduction = "cca", features = features
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(leti.integrated.sub3),
  dims = 1:50, weight.reduction = "cca"
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, group.by = "predicted.id", label = T, pt.size = 0.8) + scale_color_brewer(palette = "Set2")

# Cytotoxic signatures
Guo_CD4_GZMA <- readLines("E:/Work/Signatures/Guo_et_al_2018/CD4_GZMA.txt")
Guo_cd4_CX3CR1 <- readLines("E:/Work/Signatures/Guo_et_al_2018/cd4_CX3CR1.txt")
Gueguen_CD4_GZMA <- readLines("E:/Work/Signatures/Gueguen_et_al_2021/CD4-GZMA.txt")
Gueguen_CD8_FCGR3A <- readLines("E:/Work/Signatures/Gueguen_et_al_2021/cd4-FCGR3A.txt")
Biocarta_cytotoxic <- readLines("E:/Work/Signatures/Cytotoxic/BIOCARTA_TCYTOTOXIC_PATHWAY.txt")
GO_BP_cytotoxic <- readLines("E:/Work/Signatures/Cytotoxic/GOBP_T_CELL_MEDIATED_CYTOTOXICITY.txt")

suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(Guo_CD4_GZMA), name = "Guo_CD4_GZMA", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(Guo_cd4_CX3CR1), name = "Guo_cd4_CX3CR1", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(Gueguen_CD4_GZMA), name = "Gueguen_CD4_GZMA", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(Gueguen_CD8_FCGR3A), name = "Gueguen_CD8_FCGR3A", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(Biocarta_cytotoxic), name = "Biocarta_cytotoxic", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(GO_BP_cytotoxic), name = "GO_BP_cytotoxic", assay = "RNA")

# As Featureplots
p <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Guo_CD4_GZMA1", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "Guo_CD4_GZMA") + scale_color_paletteer_c("pals::coolwarm")
p1 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Guo_cd4_CX3CR11", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "Guo_cd4_CX3CR1") + scale_color_paletteer_c("pals::coolwarm")
p2 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Biocarta_cytotoxic1", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "Biocarta_cytotoxic") + scale_color_paletteer_c("pals::coolwarm")
p3 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "GO_BP_cytotoxic1", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "GO_BP_cytotoxic") + scale_color_paletteer_c("pals::coolwarm")
p4 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Gueguen_CD4_GZMA1", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "Gueguen_CD4_GZMA") + scale_color_paletteer_c("pals::coolwarm")
p5 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Gueguen_CD8_FCGR3A1", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "Gueguen_CD8_FCGR3A") + scale_color_paletteer_c("pals::coolwarm")
p + p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 3)

# As Featureplots split by genotype
p <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Guo_CD4_GZMA1", split.by = "genotype", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")
p1 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Guo_cd4_CX3CR11", split.by = "genotype", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")
p2 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Biocarta_cytotoxic1", split.by = "genotype", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")
p3 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "GO_BP_cytotoxic1", split.by = "genotype", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")
p4 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Gueguen_CD4_GZMA1", split.by = "genotype", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")
p5 <- FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = "Gueguen_CD8_FCGR3A1", split.by = "genotype", order = T, reduction = "umap",
  pt.size = 0.7, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm")
p + p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1, nrow = 7)

VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("Guo_CD4_GZMA1", "Guo_cd4_CX3CR11", "Biocarta_cytotoxic1", "GO_BP_cytotoxic1", "Gueguen_CD4_GZMA1", "Gueguen_CD8_FCGR3A1"), pt.size = 0.1) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
  geom_signif(
    comparisons = list(c("WT D8", "KO D8"), c("WT D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 5
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)

# Comparison with Luigia's data
new.symbols <- toupper(suv.nWTp@assays$RNA@counts@Dimnames[[1]])
suv.nWTp@assays$RNA@counts@Dimnames[[1]] <- new.symbols
new.symbols <- toupper(suv.nWTp@assays$RNA@data@Dimnames[[1]])
suv.nWTp@assays$RNA@data@Dimnames[[1]] <- new.symbols

features <- SelectIntegrationFeatures(object.list = c(suv.nWTp, suv.car.t.integrated.merged.timepoints.premerged.cd4s), features = 5000)
tumor.anchors <- FindTransferAnchors(
  reference = suv.nWTp, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s,
  dims = 1:50, reference.assay = "RNA", query.assay = "RNA", reduction = "cca", features = features
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(suv.nWTp),
  dims = 1:50, weight.reduction = "cca"
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, group.by = "predicted.id", label = T, pt.size = 0.8) + scale_color_brewer(palette = "Set2")

# Heatmap of number of DEGs (split by timepoint)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned[["label"]] <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned[["cell_type"]] <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
results.D8 <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))[1:11]) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, logfc.threshold = 0.2, ident.1 = "Mock", ident.2 = "gSUV", subset.ident = i, group.by = "genotype", assay = "RNA"))
  results.D8[[i]] <- markers
  results.D8[[i]] <- results.D8[[i]] %>% filter(p_val < 0.01)
}
lengths.deg.D8 <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))[1:11]) {
  lengths.deg.D8[i] <- print(length(rownames(results.D8[[i]])))
}
names(lengths.deg.D8) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))[1:11]
lengths.deg.D8 <- as.data.frame(lengths.deg.D8)
lengths.deg.D8[["cell_type"]] <- rownames(lengths.deg.D8)

# Same for D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
results.D28 <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))[1:11]) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, logfc.threshold = 0.2, ident.1 = "Mock", ident.2 = "gSUV", subset.ident = i, group.by = "genotype", assay = "RNA"))
  try(results.D28[[i]] <- markers)
  try(results.D28[[i]] <- results.D28[[i]] %>% filter(p_val < 0.01))
}
lengths.deg.D28 <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))[1:11]) {
  lengths.deg.D28[i] <- print(length(rownames(results.D28[[i]])))
}
names(lengths.deg.D28) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))[1:11]
lengths.deg.D28 <- as.data.frame(lengths.deg.D28)
lengths.deg.D28[["cell_type"]] <- rownames(lengths.deg.D28)

lengths.deg.D8 <- cbind(lengths.deg.D8, lengths.deg.D28)
lengths.deg.D8[, c(2, 4)] <- NULL
df <- as.data.frame(lengths.deg.D8)
df <- data.matrix(df)
pheatmap(df, cluster_rows = T, cluster_cols = T, scale = "none", fontsize = 13)

# Total
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[1:11]) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, logfc.threshold = 0.2, ident.1 = "Mock", ident.2 = "gSUV", subset.ident = i, group.by = "genotype", assay = "RNA"))
  try(results.D28[[i]] <- markers)
  try(results.D28[[i]] <- results.D28[[i]] %>% filter(p_val < 0.01))
}
lengths.deg.D28 <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[1:11]) {
  lengths.deg.D28[i] <- print(length(rownames(results.D28[[i]])))
}
names(lengths.deg.D28) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[1:11]
lengths.deg.D28 <- as.data.frame(lengths.deg.D28)
lengths.deg.D28[["cell_type"]] <- rownames(lengths.deg.D28)
df <- as.data.frame(lengths.deg.D28)
df <- data.matrix(df)
df <- df[-1,-1]
pheatmap(df, cluster_rows = T, cluster_cols = T, scale = "none", fontsize = 13)


# Dimplot with cluster names shown
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s$sub.cluster
suv.car.t.integrated.merged.timepoints.premerged.cd4s <- RenameIdents(suv.car.t.integrated.merged.timepoints.premerged.cd4s, 
  "0_0" = "C1-Naivelike-IL7R", "0_1" = "C5-Resident-ITGA1", "1" = "C4-Resident-CD69",
  "2" = "C12-Cycling-TOP2A", "3" = "C3-Stem/memory-KLF2", "4" = "C13-Cycling-ATP5F1E", "5" = "C7-Cytotoxic-GNLY", "6" = "C6-Effector/ Memory-SNHG25",
  "7" = "C11-Cycling-CENPF", "8" = "C8-Activated-HSPA1B", "9" = "C2-Naivelike-RPS27",
  "10" = "C9-Activated-PLCG2", "11" = "C10-Pro-apopototic-BAX"
)
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s) <- factor(x = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s), levels = c(
  "C1-Naivelike-IL7R", "C2-Naivelike-RPS27", "C3-Stem/memory-KLF2", "C4-Resident-CD69", "C5-Resident-ITGA1",
  "C6-Effector/ Memory-SNHG25", "C7-Cytotoxic-GNLY", "C8-Activated-HSPA1B",
  "C9-Activated-PLCG2", "C10-Pro-apopototic-BAX", "C11-Cycling-CENPF", "C12-Cycling-TOP2A", "C13-Cycling-ATP5F1E"
))

# Colors setup
mycolors <- pals::kelly(16)
mycolors.CD4 <- mycolors[4:16]
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, pt.size = 0.9, label = T, repel = T, label.size = 4) + NoLegend() + scale_color_manual(values = mycolors.CD4)

# Selected FP
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("CD8_stem_sc1", "Guo_cd4_CX3CR11", "sign_effector1", "sign_cycling1", "sign_exhausted1", "sign_ICP_complete21"), pt.size = 0.8, min.cutoff = "q1", max.cutoff = "q99", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes()

# VP of cytotoxic and stem in cluster 3 (CD4 and/or cd4s)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, idents = c("C5-Stem/memory-KLF2"))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub$genotype_timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub$genotype_timepoint, levels = c("Mock D8", "KO D8", "Mock D28", "KO D28"))
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, "CD8_stem_sc1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.8
  )
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, "Gueguen_CD8_FCGR3A1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.6
  )
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, "sign_ox_ph1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.7
  )
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, "sign_glycolysis1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.21
  )
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, "sign_FA_metabolism_GO1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.04
  )

# Metabolism single genes
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.sub, "PPA1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 2
  )

# Scatterplot cytotoxic / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "Mock"), feature1 = "Gueguen_CD8_FCGR3A", feature2 = "CD8_stem_sc", group.by = "genotype")) + xlim(c(-0, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "KO"), feature1 = "Gueguen_CD8_FCGR3A", feature2 = "CD8_stem_sc", group.by = "genotype")) + xlim(c(-0, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "Mock"), feature1 = "Gueguen_CD8_FCGR3A", feature2 = "CD8_stem_sc", group.by = "genotype")) + xlim(c(-0, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "KO"), feature1 = "Gueguen_CD8_FCGR3A", feature2 = "CD8_stem_sc", group.by = "genotype")) + xlim(c(-0, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
g <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "Mock"), feature1 = "Gueguen_CD8_FCGR3A", feature2 = "CD8_stem_sc", group.by = "genotype")) + xlim(c(-0, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
h <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "KO"), feature1 = "Gueguen_CD8_FCGR3A", feature2 = "CD8_stem_sc", group.by = "genotype")) + xlim(c(-0, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_sign.pdf", width = 9, height = 20)
plot((c | d) / (e | f) + (g | h))
dev.off()

# Scatterplot cytotoxic / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "Mock"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype"))
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "KO"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype"))
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "Mock"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype"))
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "KO"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype"))
g <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "Mock"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype"))
h <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "KO"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype"))

pdf("FS_CD4_genes.pdf", width = , height = 12)
plot((c | d) / (e | f) + (g | h))
dev.off()

# Computing the % of double positive cells (0.5 threshold)
length(WhichCells(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype == "Mock"), idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], expression = `Gueguen_CD8_FCGR3A_fix` > 0.5 & `CD8_stem_sc_fix` > 0.5)) / length(WhichCells(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype == "Mock", idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3]))) * 100
length(WhichCells(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype == "KO"), idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], expression = `Gueguen_CD8_FCGR3A_fix` > 0.5 & `CD8_stem_sc_fix` > 0.5)) / length(WhichCells(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype == "KO", idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3]))) * 100



suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$CD8_stem_sc_fix <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$CD8_stem_sc
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$Gueguen_CD8_FCGR3A_fix <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$Gueguen_CD8_FCGR3A

# Changing genotype
suv.car.t.integrated.merged.timepoints.premerged.cd4s[["genotype"]] <- suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$orig.ident
suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype == "Mock.1")] <- "Mock"
suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype == "Mock.2")] <- "Mock"
suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype == "Mock.3")] <- "Mock"
suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype == "Suv.ko.1")] <- "KO"
suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype == "Suv.ko.2")] <- "KO"
suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s@meta.data$genotype == "Suv.ko.3")] <- "KO"
suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s$genotype, levels = c("Mock", "KO"))

# tmp
# RoE heatmap
library("pheatmap")
suv.car.t.integrated.merged.timepoints.premerged$genotype_timepoint <- paste(suv.car.t.integrated.merged.timepoints.premerged$genotype, suv.car.t.integrated.merged.timepoints.premerged.cd4s$timepoint)
suv.car.t.integrated.merged.timepoints.premerged$genotype_timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged$genotype_timepoint, levels = c("WT D8", "KO D8", "WT D28", "KO D28"))
clusters.table <- table(Idents(suv.car.t.integrated.merged.timepoints.premerged), suv.car.t.integrated.merged.timepoints.premerged$genotype_timepoint)
chi.prop <- chisq.test(clusters.table)
Roe <- chi.prop$observed / chi.prop$expected
pheatmap(Roe, cluster_rows = F, cluster_cols = F, scale = "none", display_numbers = F, fontsize = 14)

# Label transfer on memory KLF2
suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, ident = c("C3-Stem/memory-KLF2"), invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, ident = c("C3-Stem/memory-KLF2"))
tumor.anchors <- FindTransferAnchors(
  reference = suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query,
  dims = 1:50, reference.assay = "RNA", query.assay = "RNA"
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.ref),
  dims = 1:50
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, group.by = "predicted.id", label = F, repel = T) + scale_color_manual(values = mycolors.CD4[4:16])
suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$predicted.id <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$predicted.id, levels = c(
  "C1-Naivelike-IL7R", "C2-Naivelike-RPS27", "C3-Stem/memory-KLF2", "C4-Resident-CD69", "C5-Resident-ITGA1",
  "C6-Effector/ Memory-SNHG25", "C7-Cytotoxic-GNLY", "C8-Activated-HSPA1B",
  "C9-Activated-PLCG2", "C10-Pro-apopototic-BAX", "C11-Cycling-CENPF", "C12-Cycling-TOP2A", "C13-Cycling-ATP5F1E"
))

suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$seurat_clusters <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.query$predicted.id
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, plot_type = "prop_fill", group_by = "genotype_timepoint") + RotatedAxis() + scale_fill_manual(values = mycolors[4:16])
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.query, "prediction.score.max", pt.size = 0, group.by = "predicted.id")


# Plot all signatures from Caushi et al. 2021
sign_CD4_total_MPR <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/CD4_total_MPR.txt"))
sign_CD4_total_NON_MPR <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/CD4_total_NON_MPR.txt"))

sign_cd4_total_MPR <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/cd4_total_MPR.txt"))
sign_cd4_total_NON_MPR <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/cd4_total_NON_MPR.txt"))
sign_EBV_T_cells <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/EBV_T_cells.txt"))
sign_Influenza_A_T_cells <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/Influenza_A_T_cells.txt"))
sign_MANA_MPR <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/MANA_MPR.txt"))
sign_MANA_NON_MPR <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/MANA_NON_MPR.txt"))
sign_MANA_T_cells <- toupper(readLines("E:/Work/Signatures/Caushi_et_al_2021/MANA_T_cells.txt"))

all_sign <- c(
  "sign_CD4_total_MPR", "sign_CD4_total_NON_MPR", "sign_cd4_total_MPR", "sign_cd4_total_NON_MPR", "sign_EBV_T_cells", "sign_Influenza_A_T_cells", "sign_MANA_MPR",
  "sign_MANA_NON_MPR", "sign_MANA_T_cells"
)
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}

# As Featureplots
pdf("Signatures_Caushi_2021_suv_CD4s.pdf", width = 12, height = 12)
for (i in 1:length(all_sign2)) {
  try(plot(FeaturePlot(
    object = suv.car.t.integrated.merged.timepoints.premerged.cd4s, features = all_sign2[i], order = T, reduction = "umap",
    pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
  ) + labs(title = all_sign[i], size = 13) + scale_color_paletteer_c("pals::coolwarm")))
}
dev.off()

# Nebulosa for key genes
pdf("nebulosa_CD4.pdf", width = 14, height = 13)
plot(plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("ZNF683", "CD27", "IL7R", "TCF7", "KLF2", "LEF1", "SELL", "CCR7", "FGFBP2", "FCGR3A", "GZMK", "GZMB", "GZMA", "GNLY", "KLRG1", "LAG3", "TOX", "ENTPD1", "HAVCR2", "PDCD1"), reduction = "umap") & NoLegend() & NoAxes())
dev.off()

# Proportions with scillus
suv.car.t.integrated.merged.timepoints.premerged.cd4s$seurat_clusters <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s)
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s, plot_type = "prop_fill", group_by = "orig.ident") + RotatedAxis() + scale_fill_manual(values = mycolors.CD4)

# Top genes correlated with a gene/signature
matrix <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s, assay = "RNA")
gene <- as.numeric(matrix["KLF2", ])
correlations <- apply(matrix, 1, function(x) {
  cor(gene, x)
})
correlations <- as.data.frame(correlations)
correlations$genes <- rownames(correlations)
colnames(correlations) <- c("correlation coefficient", "genes")
p <- top_n(correlations, n = 20, wt = `correlation coefficient`) %>%
  arrange(`correlation coefficient`) %>%
  ggplot(., aes(x = reorder(genes, `correlation coefficient`), y = `correlation coefficient`)) +
  geom_col(fill = "deepskyblue") +
  theme_cowplot() +
  ggtitle("Genes correlated with KLF2 - CD4") +
  RotatedAxis() +
  theme(axis.title.x = element_blank())
p1 <- top_n(correlations, n = 20, wt = -`correlation coefficient`) %>%
  arrange(`correlation coefficient`) %>%
  ggplot(., aes(x = reorder(genes, `correlation coefficient`), y = `correlation coefficient`)) +
  geom_col(fill = "salmon") +
  theme_cowplot() +
  ggtitle("Genes anticorrelated with KLF2 - CD4") +
  RotatedAxis() +
  theme(axis.title.x = element_blank())
p / p1
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s, c("KLF2", "ITGA1"), reduction = "umap", ncol = 2, pt.size = 0.8, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm")

# Top genes correlated with a signature
matrix <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s, assay = "RNA")
matrix <- rbind(matrix, suv.car.t.integrated.merged.timepoints.premerged.cd4s$fraietta_cr1)
matrix <- rbind(matrix, suv.car.t.integrated.merged.timepoints.premerged.cd4s$fraietta_nr1)
rownames(matrix)[33539] <- "fraietta_cr1"
rownames(matrix)[33540] <- "fraietta_nr1"
gene <- as.numeric(matrix["fraietta_nr1", ])
correlations <- apply(matrix, 1, function(x) {
  cor(gene, x)
})
correlations <- as.data.frame(correlations)
correlations$genes <- rownames(correlations)
colnames(correlations) <- c("correlation coefficient", "genes")
p <- top_n(correlations, n = 20, wt = `correlation coefficient`) %>%
  arrange(`correlation coefficient`) %>%
  ggplot(., aes(x = reorder(genes, `correlation coefficient`), y = `correlation coefficient`)) +
  geom_col(fill = "deepskyblue") +
  theme_cowplot() +
  ggtitle("Genes correlated with Fraietta NR - CD4") +
  RotatedAxis() +
  theme(axis.title.x = element_blank())
p1 <- top_n(correlations, n = 20, wt = -`correlation coefficient`) %>%
  arrange(`correlation coefficient`) %>%
  ggplot(., aes(x = reorder(genes, `correlation coefficient`), y = `correlation coefficient`)) +
  geom_col(fill = "salmon") +
  theme_cowplot() +
  ggtitle("Genes anticorrelated with Fraietta NR - CD4") +
  RotatedAxis() +
  theme(axis.title.x = element_blank())
p / p1

# Remove CD8s from the CD4 object while keeping them
suv.car.t.integrated.merged.timepoints.premerged.cd4s$old.idents <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = L8 == "T.CD4")
suv.car.t.integrated.merged.timepoints.premerged.extra.cd8s <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s, subset = L8 == "T.CD8")

DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, pt.size = 0.9, label = T, repel = T, label.size = 4) + NoLegend() + scale_color_manual(values = mycolors.CD4)

# Recompute clusters and UMAP
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- RunUMAP(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "harmony", dims = 1:50)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- FindNeighbors(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "harmony", dims = 1:50)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- FindClusters(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, graph.name = "RNA_snn", resolution = 0.25)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, pt.size = 0.9, label = T, repel = T, label.size = 4) + DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, pt.size = 0.9, label = T, repel = T, label.size = 4, group.by = "old.idents")
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("CD4", "CD40LG", "CD8A", "CD8B"), reduction = "umap", ncol = 2, pt.size = 0.8, order = T) & scale_color_paletteer_c("pals::coolwarm") & NoLegend() & NoAxes()



# Deng et al. 2020 signatures
CD4_CRS_high_vs_low <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD4_CRS_high_vs_low.txt"))
CD4_CRS_low_vs_high <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD4_CRS_low_vs_high.txt"))
CD4_PRPD_vs_CR <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD4_PRPD_vs_CR.txt"))
CD4_CR_vs_PRPD <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD4_CR_vs_PRPD.txt"))
CD8_CR_vs_PRPD <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD8_CR_vs_PRPD.txt"))
CD8_CRS_high_vs_low <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD8_CRS_high_vs_low.txt"))
CD8_CRS_low_vs_high <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD8_CRS_low_vs_high.txt"))
CD8_PRPD_vs_CR <- toupper(readLines("E:/Work/Signatures/Deng_et_al_2020/CD8_PRPD_vs_CR.txt"))


all_sign <- c(
  "CD4_CRS_high_vs_low", "CD4_CRS_low_vs_high", "CD4_PRPD_vs_CR", "CD4_CR_vs_PRPD", "CD8_CRS_high_vs_low", "CD8_CRS_low_vs_high",
  "CD8_PRPD_vs_CR", "CD8_CR_vs_PRPD"
)
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}

# As Featureplots
pdf("Signatures_Deng_2020_suv_CD4s.pdf", width = 12, height = 12)
for (i in 1:length(all_sign2)) {
  try(plot(FeaturePlot(
    object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = all_sign2[i], order = T, reduction = "umap",
    pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
  ) + labs(title = all_sign[i], size = 13) + scale_color_paletteer_c("pals::coolwarm")))
}
dev.off()

VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,
  group.by = "genotype_timepoint",
  c("CD4_CRS_high_vs_low1", "CD4_CRS_low_vs_high1", "CD4_PRPD_vs_CR1", "CD4_CR_vs_PRPD1", "CD8_CR_vs_PRPD1"), pt.size = 0.1
) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
  geom_signif(
    comparisons = list(c("WT D8", "KO D8"), c("WT D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 1.7
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)


##
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("WT", "KO"))
df.a <- df %>% filter(timepoint == "D8" & genotype == "WT")
a <- plot(ggplot(data = df.a) +
  geom_point(data = df.a, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
df.b <- df %>% filter(timepoint == "D8" & genotype == "KO")
b <- plot(ggplot(data = df.b) +
  geom_point(data = df.b, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
df.c <- df %>% filter(timepoint == "D28" & genotype == "Mock")
c <- plot(ggplot(data = df.c) +
  geom_point(data = df.c, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
df.d <- df %>% filter(timepoint == "D28" & genotype == "KO")
d <- plot(ggplot(data = df.d) +
  geom_point(data = df.d, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
a + b + c + d + plot_layout(ncol = 4) & NoLegend()

# Rename clusters
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- RenameIdents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,
  "0" = "C1-Active translation-IL7R", "1" = "C5-Activated-CD69",
  "2" = "C4-Resident Memory-ITGA1", "3" = "C2-Stem/memory-KLF2", "4" = "C13-Cycling-TOP2A", "5" = "C14-Cycling/S-PCNA", "6" = "C8-Effector-KLRB1",
  "7" = "C10-Cytosolic stress-LncRNA", "8" = "C12-Cycling/G2M-HSPA1B", "9" = "C9-Cytosolic stress-HSPA1B",
  "10" = "C7-Active translation-RPS27", "11" = "C3-Stem/cytotoxic-GNLY", "12" = "C6-Activated-JUN", "13" = "C11-Pro-apopototic-BAX"
)
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) <- factor(x = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), levels = c(
  "C1-Active translation-IL7R", "C2-Stem/memory-KLF2", "C3-Stem/cytotoxic-GNLY", "C4-Resident Memory-ITGA1", "C5-Activated-CD69", "C6-Activated-JUN", "C7-Active translation-RPS27",
  "C8-Effector-KLRB1", "C9-Cytosolic stress-HSPA1B", "C10-Cytosolic stress-LncRNA", "C11-Pro-apopototic-BAX", "C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"
))
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)

###############################
#### Regenerate all the figures
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, pt.size = 0.8, label = F, cols = mycolors.CD4, repel = T)

# Density
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("Mock", "KO"))
df.a <- df %>% filter(timepoint == "D8" & genotype == "Mock")
a <- plot(ggplot(data = df.a) +
  geom_point(data = df.a, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ) +
  NoAxes() +
  ggtitle("Mock D8") +
  theme(plot.title = element_text(size = 22, face = "bold")))
df.b <- df %>% filter(timepoint == "D8" & genotype == "KO")
b <- plot(ggplot(data = df.b) +
  geom_point(data = df.b, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ) +
  NoAxes() +
  ggtitle("KO D8") +
  theme(plot.title = element_text(size = 22, face = "bold")))
df.c <- df %>% filter(timepoint == "D28" & genotype == "Mock")
c <- plot(ggplot(data = df.c) +
  geom_point(data = df.c, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ) +
  NoAxes() +
  ggtitle("Mock D28") +
  theme(plot.title = element_text(size = 22, face = "bold")))
df.d <- df %>% filter(timepoint == "D28" & genotype == "KO")
d <- plot(ggplot(data = df.d) +
  geom_point(data = df.d, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 15, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-6.5, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ) +
  NoAxes() +
  ggtitle("KO D28") +
  theme(plot.title = element_text(size = 22, face = "bold")))
a + b + c + d + plot_layout(ncol = 2) & NoLegend()

# Heatmap of number of DEGs (split by timepoint)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned[["label"]] <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned[["cell_type"]] <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
results.D8 <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, logfc.threshold = 0.2, ident.1 = "Mock", ident.2 = "KO", subset.ident = i, group.by = "genotype", assay = "RNA"))
  results.D8[[i]] <- markers
  results.D8[[i]] <- results.D8[[i]] %>% filter(p_val < 0.01)
}
lengths.deg.D8 <- c()
for (i in 1:length(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))) {
  lengths.deg.D8[i] <- print(length(rownames(results.D8[[i]])))
}
names(lengths.deg.D8) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))
lengths.deg.D8 <- as.data.frame(lengths.deg.D8)
lengths.deg.D8[["cell_type"]] <- rownames(lengths.deg.D8)

# Same for D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
results.D28 <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, logfc.threshold = 0.2, ident.1 = "Mock", ident.2 = "KO", subset.ident = i, group.by = "genotype", assay = "RNA"))
  try(results.D28[[i]] <- markers)
  try(results.D28[[i]] <- results.D28[[i]] %>% filter(p_val < 0.01))
}
lengths.deg.D28 <- c()
for (i in 1:length(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))) {
  lengths.deg.D28[i] <- print(length(rownames(results.D28[[i]])))
}
names(lengths.deg.D28) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))
lengths.deg.D28 <- as.data.frame(lengths.deg.D28)
lengths.deg.D28[["cell_type"]] <- rownames(lengths.deg.D28)

lengths.deg.D8 <- cbind(lengths.deg.D8, lengths.deg.D28)
lengths.deg.D8[, c(2, 4)] <- NULL
df <- as.data.frame(lengths.deg.D8)
df <- data.matrix(df)
colnames(df) <- c("# of DEGs at D8", "# of DEGs at D28")
pheatmap(df, cluster_rows = T, cluster_cols = T, scale = "none", fontsize = 14)

# Compute the number of shared DEGs with the global signature
markers.WT.KO.CD4 <- FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, logfc.threshold = 0.2, ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA", only.pos = T)

percent.shared <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, logfc.threshold = 0.2, ident.1 = "KO", ident.2 = "Mock", subset.ident = i, group.by = "genotype", assay = "RNA", only.pos = T))
  percent.shared[[i]] <- length(which(rownames(markers.WT.KO.CD4) %in% rownames(markers))) / length(rownames(markers.WT.KO.CD4)) * 100
}
shared.D8 <- c()
for (i in 1:length(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))) {
  shared.D8[i] <- print(percent.shared[[i]])
}
names(shared.D8) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8))
shared.D8 <- as.data.frame(shared.D8)
shared.D8[["cell_type"]] <- rownames(shared.D8)

percent.shared <- c()
for (i in levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))) {
  markers <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, logfc.threshold = 0.2, ident.1 = "KO", ident.2 = "Mock", subset.ident = i, group.by = "genotype", assay = "RNA", only.pos = T))
  percent.shared[[i]] <- length(which(rownames(markers.WT.KO.CD4) %in% rownames(markers))) / length(rownames(markers.WT.KO.CD4)) * 100
  print(percent.shared[[i]])
}
shared.D28 <- c()
for (i in 1:length(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))) {
  shared.D28[i] <- print(percent.shared[[i]])
}
names(shared.D28) <- levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))
shared.D28 <- as.data.frame(shared.D28)
shared.D28[["cell_type"]] <- rownames(shared.D28)

shared.D8 <- cbind(shared.D8, shared.D28)
shared.D8[, c(2, 4)] <- NULL
df <- as.data.frame(shared.D8)
df.CD4 <- data.matrix(df)
colnames(df.CD4) <- c("D8", "D28")
breaksList <- seq(0, 100, by = 1)
pheatmap(df.CD4, cluster_rows = T, cluster_cols = T, scale = "none", fontsize = 14, breaks = breaksList, col = colorRampPalette(c("navy", "white", "firebrick3", "firebrick3"))(100))

# Nebulosa for key genes
plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("CD27", "IL7R", "TCF7", "KLF2", "LEF1", "SELL", "CCR7", "FGFBP2", "FCGR3A", "NKG7", "GZMK", "GZMB", "GZMA", "GNLY", "KLRG1", "LAG3", "TOX", "ENTPD1", "HAVCR2", "PDCD1"), reduction = "umap", ) & NoLegend() & NoAxes()


# DP
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, label = F, pt.size = 0.3, label.size = 6, cols = mycolors.CD4) + NoLegend() +
  DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "genotype", label = F, pt.size = 0.3, label.size = 6, cols = c("black", "red")) +
  DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "timepoint", label = F, pt.size = 0.3, label.size = 6) + NoLegend()

# ROE
# Roe heatmap split by D8 and D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$genotype <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$genotype, levels = c("Mock", "KO"))
clusters.table <- table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8), suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$genotype)
chi.prop <- chisq.test(clusters.table)
Roe <- chi.prop$observed / chi.prop$expected
pheatmap(Roe, cluster_rows = F, cluster_cols = F, scale = "none", display_numbers = F, fontsize = 14)

# Heatmap
markers <- FindAllMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, assay = "RNA", only.pos = T)
markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 8, wt = avg_log2FC) -> top8
Average <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = unique(top8$gene), return.seurat = T, assays = "RNA")
Average <- ScaleData(object = Average, features = unique(top8$gene))
plot(DoHeatmap(Average, features = unique(top8$gene), draw.lines = F, size = 3, group.colors = mycolors.CD4) + scale_fill_viridis_c(option = "inferno"))

# CytoTRACE
library(CytoTRACE)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ident = c("C12-Cycling/G2M-HSPA1B", "C14-Cycling/S-PCNA", "C13-Cycling-TOP2A"), invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, cells = sample(Cells(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub), 15000))
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub$genotype

mat <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, assay = "RNA")
mat <- as.data.frame(mat)
results.cyto <- CytoTRACE(mat, ncores = 1)
idents <- as.character(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub))
names(idents) <- names(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub))
emb <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, reduction = "umap"))
plotCytoTRACE(results.cyto, emb = emb, phenotype = idents, gene = "KLF2")
plotCytoGenes(results.cyto, numOfGenes = 10)

# VP Fraietta
a <- VlnPlot(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8"), "fraietta_cr1", group.by = "genotype", pt.size = 0) + labs(title = "Fraietta CR", size = 13) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
  geom_signif(
    comparisons = list(c("Mock", "KO")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.1
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)
b <- VlnPlot(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8"), "fraietta_nr1", group.by = "genotype", pt.size = 0) + labs(title = "Fraietta NR", size = 13) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
  geom_signif(
    comparisons = list(c("Mock", "KO")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.07
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)
a | b

# Volcano
library(EnhancedVolcano)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- NormalizeData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) %>%
  FindVariableFeatures() %>%
  ScaleData()
degs.klf2 <- FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ident.1 = "KO", ident.2 = "Mock", subset.ident = "C2-Stem/memory-KLF2", group.by = "genotype", assay = "RNA")
plot(EnhancedVolcano(degs.klf2,
  lab = rownames(degs.klf2),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T, xlim = c(-1, 1),
  pCutoff = 0.001, maxoverlapsConnectors = 20,
  FCcutoff = 0.3
)) + ggtitle("CD4-C2-Stem/memory-KLF2 - KO vs Mock")

# Remake Density plot - common scale for Mock & KO within each time point, add cell number per condition
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("Mock", "gSUV"))
df.a <- df %>% filter(timepoint == "D8")
table(df.a$genotype)
labels <- c("Mock - 9322 cells", "gSUV - 9629 cells")
names(labels) <- c("Mock", "gSUV")
a <- plot(ggplot(data = df.a) +
  geom_point(data = df.a, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 20, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
labels <- c("Mock - 14073 cells", "gSUV - 13611 cells")
names(labels) <- c("Mock", "gSUV")
df.b <- df %>% filter(timepoint == "D28")
table(df.b$genotype)
b <- plot(ggplot(data = df.b) +
  geom_point(data = df.b, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-7, 8)) +
  scale_fill_viridis_c(option = "plasma") +
  facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))

a + b + plot_layout(ncol = 1) & NoLegend()

# New michael signatures
Cytotoxic_Effector <- readLines("E:/Work/Signatures/Signatures_Michael_CART/Cytotoxic-Effector_210819.txt")
Dysfunction_Exhaustion <- readLines("E:/Work/Signatures/Signatures_Michael_CART/Dysfunction-Exhaustion_210819.txt")
Stemness_Memory <- readLines("E:/Work/Signatures/Signatures_Michael_CART/Stemness-Memory_210819.txt")
Tissue_homing <- readLines("E:/Work/Signatures/Signatures_Michael_CART/Tissue homing_210927.txt")

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Cytotoxic_Effector), name = "Cytotoxic_Effector", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Dysfunction_Exhaustion), name = "Dysfunction_Exhaustion", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Stemness_Memory), name = "Stemness_Memory", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Tissue_homing), name = "Tissue_homing", assay = "RNA")

# FP
all_sign <- c("Cytotoxic_Effector", "Dysfunction_Exhaustion", "Stemness_Memory", "Tissue_homing")
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}

# As Featureplots
pdf("Signatures_CD4s_Michael.pdf", width = 9, height = 8)
for (i in 1:length(all_sign2)) {
  try(plot(FeaturePlot(
    object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = all_sign2[i], order = T, reduction = "umap",
    pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
  ) + labs(title = all_sign[i], size = 13) + scale_color_paletteer_c("pals::coolwarm")))
}
dev.off()

pdf("Signatures_CD4s_Michael_split_genotype.pdf", width = 16, height = 9)
for (i in 1:length(all_sign2)) {
  try(plot(FeaturePlot(
    object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = all_sign2[i], split.by = "genotype", order = T, reduction = "umap",
    pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
  ) & scale_color_paletteer_c("pals::coolwarm")))
}
dev.off()

# FP single genes
pdf("Genes_CD4s_Michael.pdf", width = 12, height = 10)
plot(FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = c("KLF2", "CCR7", "GNLY", "SELL"), order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm"))
dev.off()

pdf("Genes_CD4s_Michael_genotype.pdf", width = 12, height = 18)
plot(FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = c("KLF2", "CCR7", "GNLY", "SELL"), split.by = "genotype", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm"))
dev.off()

pdf("Genes_CD4s_Michael_genotype_timepoint.pdf", width = 18, height = 18)
plot(FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = c("KLF2", "CCR7", "GNLY", "SELL"), split.by = "genotype_timepoint", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm"))
dev.off()

# Same with sign
pdf("Sign_CD4s_Michael_genotype_timepoint.pdf", width = 18, height = 18)
plot(FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = all_sign2, split.by = "genotype_timepoint", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm"))
dev.off()

# VP with sign
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Cytotoxic_Effector <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Cytotoxic_Effector1
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Dysfunction_Exhaustion <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Dysfunction_Exhaustion1
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Stemness_Memory <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Stemness_Memory1
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Tissue_homing <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Tissue_homing1
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, all_sign, group.by = "genotype_timepoint", ncol = 4, fill.by = "ident", pt.size = 0) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.4
  )

# Volcanos

library(EnhancedVolcano)
pdf("volcano_cd4_KO_vs_WT_total_c2.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset.ident = "C2-Stem/memory-KLF2", ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs), selectLab = c("KLF2","LEF1","IL7R","CCL5","GZMB","S1PR1","CST7","GZMA","KLF3","SELL","TCF7","NKG7","CD27","CCR4","CCR2","ID2"),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T,
  pCutoff = 0.001, title = "c2-KLF2 - SUV39H1-KO vs Mock",
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
pdf("volcano_cd4_KO_vs_WT_D8_c2.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, subset.ident = "C2-Stem/memory-KLF2",ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs), selectLab = c("KLF2","LEF1","IL7R","CCL5","GZMB","S1PR1","CST7","GZMA","KLF3","SELL","TCF7","NKG7","CD27","CCR4","CCR2","ID2"),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T,
  pCutoff = 0.001, title = "c2-KLF2 - D8 - SUV39H1-KO vs Mock",
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
pdf("volcano_cd4_KO_vs_WT_D28_c2.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, subset.ident = "C2-Stem/memory-KLF2", ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs), selectLab = c("KLF2","LEF1","IL7R","CCL5","GZMB","S1PR1","CST7","GZMA","KLF3","SELL","TCF7","NKG7","CD27","CCR4","CCR2","ID2"),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T,
  pCutoff = 0.001, title = "c2-KLF2 - D28 - SUV39H1-KO vs Mock",
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

# Label transfer Guo et al.
tumor.anchors <- FindTransferAnchors(
  reference = guo.seurat, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,
  dims = 1:30
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = guo.seurat$majorCluster,
  dims = 1:30
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "predicted.id", label = F, pt.size = 1, cols = paletteer_d("pals::polychrome", n = 21))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$seurat_clusters <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)
dittoBarPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, var =  "seurat_clusters", group.by = "predicted.id")



# Label Transfer with NSCLC data
tumor.anchors <- FindTransferAnchors(
  reference = eleven.tils.cd3.integrated, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,
  dims = 1:30, query.assay = "RNA", reference.assay = "integrated"
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(eleven.tils.cd3.integrated),
  dims = 1:30
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "predicted.id", label = F, repel = T, cols = paletteer_d("pals::polychrome", n = 21))
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = "prediction.score.max")
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = "prediction.score.max", pt.size = 0)

dittoBarPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, var =  "seurat_clusters", group.by = "predicted.id")
dittoBarPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,group.by =    "seurat_clusters", var ="predicted.id")


# Scatterplot cytotoxic / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "Mock"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "KO"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "Mock"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "KO"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
g <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "Mock"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
h <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "KO"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_cytotox.pdf", width = 12, height = 18)
plot((c | d) / (e | f) + (g | h))
dev.off()

# Scatterplot dysfunction / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "Mock"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "KO"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "Mock"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "KO"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
g <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "Mock"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
h <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "KO"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_dysfunction.pdf", width = 12, height = 18)
plot((c | d) / (e | f) + (g | h))
dev.off()

# Scatterplot tissue homing / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "Mock"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "KO"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "Mock"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "KO"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
g <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "Mock"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
h <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "KO"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_homing.pdf", width = 12, height = 18)
plot((c | d) / (e | f) + (g | h))
dev.off()

# Scatterplot KLF2/GNLY
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "Mock"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype == "KO"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "Mock"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[3], subset = genotype == "KO"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
g <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "Mock"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
h <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[13], subset = genotype == "KO"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_genes.pdf", width = 12, height = 18)
plot((c | d) / (e | f) + (g | h))
dev.off()

##### Same but split by tiempoint

# Scatterplot cytotoxic / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D8"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D8"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D28"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D28"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_cytotox_timepoint.pdf", width = 12, height = 12)
plot((c | d) / (e | f))
dev.off()

# Scatterplot dysfunction / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D8"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D8"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D28"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D28"), feature1 = "Dysfunction_Exhaustion", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_dysf_timepoint.pdf", width = 12, height = 12)
plot((c | d) / (e | f))
dev.off()


# Scatterplot tissue homing / stem-like
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D8"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D8"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D28"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D28"), feature1 = "Tissue_homing", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(-1, 1)) + ylim(c(-1, 1.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_homing_timepoint.pdf", width = 12, height = 12)
plot((c | d) / (e | f))
dev.off()

# Scatterplot KLF2/GNLY
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D8"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D8"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D28"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D28"), feature1 = "GNLY", feature2 = "KLF2", group.by = "genotype")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

pdf("FS_CD4_genes_timepoint.pdf", width = 12, height = 12)
plot((c | d) / (e | f))
dev.off()

# Label transfer on memory KLF2
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ident = c("C2-Stem/memory-KLF2"), invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ident = c("C2-Stem/memory-KLF2"))
tumor.anchors <- FindTransferAnchors(
  reference = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2,
  dims = 1:50, reference.assay = "RNA", query.assay = "RNA"
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref),
  dims = 1:50
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, metadata = predictions)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, group.by = "predicted.id", label = F, repel = T) + scale_color_manual(values = mycolors.CD4[4:16])
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2$predicted.id <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2$predicted.id, levels = c(
  "C1-Active translation-IL7R", "C2-Stem/memory-KLF2", "C3-Stem/cytotoxic-GNLY", "C4-Resident Memory-ITGA1", "C5-Activated-CD69", "C6-Activated-JUN", "C7-Active translation-RPS27",
  "C8-Effector-KLRB1", "C9-Cytosolic stress-HSPA1B", "C10-Cytosolic stress-LncRNA", "C11-Pro-apopototic-BAX", "C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"
))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2$seurat_clusters <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2$predicted.id
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, plot_type = "prop_fill", group_by = "genotype_timepoint") + RotatedAxis() + scale_fill_manual(values = mycolors[c(5, 7:9, 11:14, 16:18)])


# VP?Fraietta
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("fraietta_cr1", "fraietta_nr1"), group.by = "genotype_timepoint", ncol = 2, pt.size = 0) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.1
  )


# VP?Fraietta
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "fraietta_cr1", group.by = "genotype_timepoint", pt.size = 0) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.15
  ) & labs(title = "CD4 - Fraietta CR", size = 13)
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "fraietta_nr1", group.by = "genotype_timepoint", pt.size = 0) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.05
  ) & labs(title = "CD4 - Fraietta NR", size = 13)


# VP integrins
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("ITGA4", "ITGB1"), group.by = "genotype_timepoint", pt.size = 0) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 3.5
  )
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, c("ITGA4", "ITGB1"), group.by = "genotype_timepoint", pt.size = 0) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 3
  )

# TCR features
"%!in%" <- function(x, y) !("%in%"(x, y))
TCRA <- read.csv("E:/Work/Signatures/TCR_genes/TCRA.txt", sep = "\t")
TCRA.features <- TCRA$Approved.symbol
TCRB <- read.csv("E:/Work/Signatures/TCR_genes/TCRB.txt", sep = "\t")
TCRB.features <- TCRB$Approved.symbol
TCRD <- read.csv("E:/Work/Signatures/TCR_genes/TCRD.txt", sep = "\t")
TCRD.features <- TCRD$Approved.symbol
TCRG <- read.csv("E:/Work/Signatures/TCR_genes/TCRG.txt", sep = "\t")
TCRG.features <- TCRG$Approved.symbol
TCR.features <- c(TCRA.features, TCRB.features, TCRG.features, TCRD.features)
TCR.removed <- rownames(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)[which(rownames(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) %!in% TCR.features)]
markers <- FindAllMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, assay = "RNA", features = TCR.removed, only.pos = T)
write.csv(markers, "markers_CD4s.txt", quote = F, row.names = F)
markers <- read.csv(file = "markers_CD4s.txt")
markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

Average <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = unique(top5$gene), return.seurat = T, assays = "RNA")
Average <- ScaleData(object = Average, features = unique(top5$gene))
plot(DoHeatmap(Average, features = unique(top5$gene), draw.lines = F, size = 3, group.colors = mycolors.CD4) + scale_fill_viridis_c(option = "inferno"))

# VP metabolism
a <- VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, "sign_ox_ph1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.7
  ) & labs(title = "Oxphos", size = 13)
b <- VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, "sign_glycolysis1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.21
  ) & labs(title = "Glycolysis", size = 13)
c <- VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, "sign_FA_metabolism_GO1", group.by = "genotype_timepoint", pt.size = 0) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "KO D8"), c("Mock D28", "KO D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.04
  ) & labs(title = "Fatty Acid Metabolism", size = 13)
a | b | c

# Volcano cycling
library(EnhancedVolcano)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ident = c("C14-Cycling/S-PCNA"))
pdf("volcano_CD8_KO_vs_WT_total_c14.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling, ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T,
  pCutoff = 0.001, title = "c14 - SUV39H1-KO vs Mock",
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling, subset = timepoint == "D8")
pdf("volcano_CD8_KO_vs_WT_D8_c14.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling.D8, ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T,
  pCutoff = 0.001, title = "c14 - D8 - SUV39H1-KO vs Mock",
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling, subset = timepoint == "D28")
pdf("volcano_CD8_KO_vs_WT_D28_c14.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cycling.D28, ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T,
  pCutoff = 0.001, title = "c14 - D28 - SUV39H1-KO vs Mock",
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

# Label transfer of cycling
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- FindVariableFeatures(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, nfeatures = 4000)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == "G1")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref, ident = c("C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"), invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == "G1", invert = T)
tumor.anchors <- FindTransferAnchors(
  reference = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query,
  dims = 1:50, reference.assay = "RNA", query.assay = "RNA", reference.reduction = "harmony", reduction = "harmony"
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref),
  dims = 1:50
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, metadata = predictions)
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, "prediction.score.max", pt.size = 0) + scale_color_manual(values = mycolors.CD4)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$seurat_clusters <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id, levels = c(
  "C1-Active translation-IL7R", "C2-Stem/memory-KLF2", "C3-Stem/cytotoxic-GNLY", "C4-Resident Memory-ITGA1", "C5-Activated-CD69", "C6-Activated-JUN", "C7-Active translation-RPS27",
  "C8-Effector-KLRB1", "C9-Cytosolic stress-HSPA1B", "C10-Cytosolic stress-LncRNA", "C11-Pro-apopototic-BAX", "C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"
))
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, group.by = "predicted.id", label = F, repel = T) + scale_color_manual(values = mycolors.CD4)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$genotype_timepoint <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$genotype_timepoint, levels = c("Mock D8", "gSUV D8", "Mock D28", "gSUV D28"))
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, plot_type = "prop_fill", group_by = "genotype_timepoint") + RotatedAxis() + scale_fill_manual(values = mycolors.CD4)
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, plot_type = "prop_fill", group_by = "Phase") + RotatedAxis() + scale_fill_manual(values = mycolors.CD4)


# Remake Density with lines instead of area
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("Mock", "KO"))
df.a <- df %>% filter(timepoint == "D8")
table(df.a$genotype)
labels <- c("Mock - 9322 cells", "KO - 9629 cells")
names(labels) <- c("Mock", "KO")
a <- plot(ggplot(data = df.a) +
  geom_point(data = df.a, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 20, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-7, 8)) +
  scale_fill_gradient(low = "yellow", high = "red") +
  facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))
labels <- c("Mock - 14073 cells", "KO - 13611 cells")
names(labels) <- c("Mock", "KO")
df.b <- df %>% filter(timepoint == "D28")
table(df.b$genotype)
b <- plot(ggplot(data = df.b) +
  geom_point(data = df.b, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
  stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
  theme_classic() +
  xlim(c(-9, 11)) +
  ylim(c(-7, 8)) +
  scale_fill_gradient(low = "yellow", high = "red") +
  facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black",
      face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black",
      face = "bold"
    )
  ))

a + b + plot_layout(ncol = 1) & NoLegend()

# Residency
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("ZNF683", "ITGAE", "ITGA1", "CXCR6"), reduction = "umap", ncol = 2, pt.size = 0.8, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm")
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("ID2", "ID3"), reduction = "umap", ncol = 2, pt.size = 0.8, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm")

# FP with same scale
FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = all_sign, split.by = "genotype_timepoint", order = T, reduction = "umap",
  pt.size = 0.6, min.cutoff = "q1", max.cutoff = "q99"
) & scale_color_paletteer_c("pals::coolwarm") & theme(legend.position = "right")

# Label transfer between CD8 and CD4s
tumor.anchors <- FindTransferAnchors(
  reference = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,
  dims = 1:50
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned),
  dims = 1:50
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, metadata = predictions)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$predicted.id <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$predicted.id, levels = c("C1-Resident Memory-ZNF683", "C2-Stem/memory-KLF2", "C3-Activated-FOS", "C4-Cycling/G2M-TOP2A", "C5-Cycling/S-CENPF"))

DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "predicted.id", label = T, pt.size = 1, cols = mycolors.CD8)
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "predicted.id", features = "prediction.score.max", pt.size = 0, cols = mycolors.CD8)

# Heatmap proportions
library(data.table)
library(gridExtra)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype_timepoint == "Mock D8")
freq_table <- as.data.frame(prop.table(
  x = table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT), suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT@meta.data[, "orig.ident"]),
  margin = 1
))
freq_table <- dcast(freq_table, Var1 ~ Var2)
rownames(freq_table) <- freq_table$Var1
tmp <- rownames(freq_table)
freq_table <- freq_table[, -1]
freq_table <- as.matrix(freq_table)
clusters_correlations <- cor(t(freq_table))
a <- pheatmap(clusters_correlations, scale = "none", fontsize = 12, main = "CD4 correlations - Mock D8")

suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype_timepoint == "KO D8")
freq_table <- as.data.frame(prop.table(
  x = table(suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO@meta.data[, "orig.ident"], Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO)),
  margin = 1
))
freq_table <- dcast(freq_table, Var1 ~ Var2)
rownames(freq_table) <- freq_table$Var1
tmp <- rownames(freq_table)
freq_table <- freq_table[, -1]
freq_table <- as.matrix(freq_table)
clusters_correlations <- cor(freq_table)
b <- pheatmap(clusters_correlations, scale = "none", fontsize = 12, main = "CD4 correlations - KO D8")

suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype_timepoint == "Mock D28")
freq_table <- as.data.frame(prop.table(
  x = table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT), suv.car.t.integrated.merged.timepoints.premerged.cd4s.WT@meta.data[, "orig.ident"]),
  margin = 1
))
freq_table <- dcast(freq_table, Var1 ~ Var2)
rownames(freq_table) <- freq_table$Var1
tmp <- rownames(freq_table)
freq_table <- freq_table[, -1]
freq_table <- as.matrix(freq_table)
clusters_correlations <- cor(t(freq_table))
c <- pheatmap(clusters_correlations, scale = "none", fontsize = 12, main = "CD4 correlations - Mock D28")

suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = genotype_timepoint == "KO D28")
freq_table <- as.data.frame(prop.table(
  x = table(suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO@meta.data[, "orig.ident"], Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.KO)),
  margin = 1
))
freq_table <- dcast(freq_table, Var1 ~ Var2)
rownames(freq_table) <- freq_table$Var1
tmp <- rownames(freq_table)
freq_table <- freq_table[, -1]
freq_table <- as.matrix(freq_table)
clusters_correlations <- cor(freq_table)
d <- pheatmap(clusters_correlations, scale = "none", fontsize = 12, main = "CD4 correlations - KO D28")

plot_list <- list()
plot_list[["p1"]] <- a[[4]]
plot_list[["p2"]] <- b[[4]]
plot_list[["p3"]] <- c[[4]]
plot_list[["p4"]] <- d[[4]]
grid.arrange(arrangeGrob(grobs = plot_list, ncol = 2))

# Exh progenitors signatures
Eff_like_Miller <- toupper(readLines("E:/Work/Signatures/Miller_et_al_2019/eff_like.txt"))
Prog_exh_Miller <- toupper(readLines("E:/Work/Signatures/Miller_et_al_2019/prog_exh.txt"))
Term_exh_Miller <- toupper(readLines("E:/Work/Signatures/Miller_et_al_2019/term_exh.txt"))
Mem_Prec_Yao <- toupper(readLines("E:/Work/Signatures/Exhausted/Yao_et_al_2019/Memory_precursors_Yao.txt"))
Prog_like_Yao <- toupper(readLines("E:/Work/Signatures/Exhausted/Yao_et_al_2019/Progenitor_like_Yao.txt"))
Term_exh_Yao <- toupper(readLines("E:/Work/Signatures/Exhausted/Yao_et_al_2019/terminally_exhausted_Yao.txt"))

all_sign <- c("Eff_like_Miller", "Prog_exh_Miller", "Term_exh_Miller", "Mem_Prec_Yao", "Prog_like_Yao", "Term_exh_Yao")
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}

FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = all_sign2, ncol = 3, pt.size = 0.7, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes() & NoLegend()

# Heatmap proportions
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$seurat_clusters <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, plot_type = "prop_fill", group_by = "timepoint_genotype") + scale_fill_manual(values = mycolors.CD4)

freq_table <- as.data.frame(prop.table(
  x = table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data[, "genotype_timepoint"]),
  margin = 1
))
freq_table <- dcast(freq_table, Var1 ~ Var2)
rownames(freq_table) <- freq_table$Var1
tmp <- rownames(freq_table)
freq_table <- freq_table[, -1]
freq_table <- as.matrix(freq_table)
pheatmap(freq_table, scale = "none", fontsize = 12, main = "Cluster frequencies", cluster_rows = F, cluster_cols = F, color = viridis(50))


# volcano cycling Mock vs cycling KO
# Plot volcano KO vs WT
library(EnhancedVolcano)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.cycling <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == "G1", invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.cycling, subset = timepoint == "D8")
pdf("volcano_CD4_D8_KO_vs_WT_cycling.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.D8, ident.1 = "KO", ident.2 = "Mock", logfc.threshold = 0.1, group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T, title = "D8 - Cycling cells - KO vs Mock",
  pCutoff = 0.001,
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

suv.car.t.integrated.merged.timepoints.premerged.cd4s.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.cycling, subset = timepoint == "D28")
pdf("volcano_CD4_D28_KO_vs_WT_cycling.pdf", width = 9, height = 9)
degs <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.D28, ident.1 = "KO", ident.2 = "Mock", logfc.threshold = 0.1, group.by = "genotype", assay = "RNA"))
try(plot(EnhancedVolcano(degs,
  lab = rownames(degs),
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T, title = "D28 - Cycling cells - KO vs Mock",
  pCutoff = 0.001,
  FCcutoff = 0.3, xlim = c(-1.5, 1.5)
)))
dev.off()

# Scatterplot cytotoxic / stem-like with contours
c <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D8"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()
d <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D8"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()
e <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "Mock D28"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()
f <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[2], subset = genotype_timepoint == "KO D28"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()

g <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[4], subset = genotype_timepoint == "Mock D8"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()
h <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[4], subset = genotype_timepoint == "KO D8"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()
i <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[4], subset = genotype_timepoint == "Mock D28"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()
j <- plot(FeatureScatter(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned))[4], subset = genotype_timepoint == "KO D28"), feature1 = "Cytotoxic_Effector", feature2 = "Stemness_Memory", group.by = "genotype")) + xlim(c(0, 1)) + ylim(c(-0.5, 1)) + geom_hline(yintercept = 0.25) + geom_vline(xintercept = 0.5) + geom_density_2d()


pdf("FS_CD4_cytotox_timepoint.pdf", width = 14, height = 7)
plot((c | d | e | f) / (g | h | i | j))
dev.off()

# Add FP sign TRM

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(sign_memory), name = "Central_memory", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Tissue_Resident <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Central_memory1
pdf(file = "Signatures_CD4.pdf", height = 12, width = 13)
plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ncol = 2, c("Cytotoxic_Effector", "Dysfunction_Exhaustion", "Stemness_Memory", "Tissue_Resident"), pt.size = 0.4, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes())
dev.off()

png(filename = "Signatures_CD4.png", res = 300, height = 3200, width = 3200)
plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ncol = 2, c("Cytotoxic_Effector", "Dysfunction_Exhaustion", "Stemness_Memory", "Tissue_Resident"), pt.size = 0.4, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes())
dev.off()


FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ncol = 1, c("Central_memory1"), pt.size = 0.4, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes() & ggtitle("Central Memory") & theme(title = element_text(hjust = 1, size = 25))


# Plot volcano KO vs WT
library(EnhancedVolcano)
celltype1 <- rownames(degs.D8)
celltype2 <- c("LEF1", "SELL", "TCF7", "S1PR1", "CD27")
colors <- rep("red", 35)

# D8
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
degs.D8 <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.1, assay = "RNA"))
degs.D8.bis <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.3, assay = "RNA"))
colors <- rep("black", 17)
colors[which(rownames(degs.D8.bis) %in%  c("KLF2","LEF1","IL7R","CCL5","GZMB","S1PR1","CST7","GZMA","KLF3","SELL","TCF7","NKG7","CD27","CCR4","CCR2","ID2","MAF","CCR2","CD2"))] <- "red"
try(plot(EnhancedVolcano(degs.D8,
  lab = rownames(degs.D8), labCol = colors,
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T,
  pCutoff = 0.001,
  FCcutoff = 0.3, xlim = c(-1.5, 1.5), title = "CD4 - D8 - KO vs Mock"
)))
write.csv(degs.D8, "degs.D8.csv", quote = F, row.names = T)

# D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
degs.D28 <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.1, assay = "RNA"))
degs.D28.bis <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.25, assay = "RNA"))
colors <- rep("black", 61)
colors[which(rownames(degs.D28.bis) %in% c("KLF2","LEF1","IL7R","CCL5","GZMB","S1PR1","CST7","GZMA","KLF3","SELL","TCF7","NKG7","CD27","CCR4","CCR2","ID2","ITGA1"))] <- "red"
try(plot(EnhancedVolcano(degs.D28,
  lab = rownames(degs.D28), labCol = colors,
  x = "avg_log2FC",
  y = "p_val_adj",
  drawConnectors = T, captionLabSize = 0, labSize = 6, subtitleLabSize = 0, 
  pCutoff = 0.001,
  FCcutoff = 0.25, xlim = c(-1.5, 1.5), title = "CD4 - D28 - gSUV vs Mock"
)))
write.csv(degs.D28, "degs.D28.csv", quote = F, row.names = T)

### CYtotrace Mock vs KO (only for C2/KLF2)
library(CytoTRACE)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ident = c("C2-Stem/memory-KLF2"), invert = F)
# D8
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, subset = timepoint == "D8")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D8) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub$genotype
mat <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D8, assay = "RNA")
mat <- as.data.frame(mat)
results.cyto <- CytoTRACE(mat, ncores = 1)
idents <- as.character(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D8))
names(idents) <- names(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D8))
emb <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D8, reduction = "umap"))
plotCytoTRACE(results.cyto, emb = emb, phenotype = idents, gene = "KLF2")
plotCytoGenes(results.cyto, numOfGenes = 10)

# 28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, subset = timepoint == "D8")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D28) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub$genotype
mat <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D28, assay = "RNA")
mat <- as.data.frame(mat)
results.cyto <- CytoTRACE(mat, ncores = 1)
idents <- as.character(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D28))
names(idents) <- names(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D28))
emb <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.D28, reduction = "umap"))
plotCytoTRACE(results.cyto, emb = emb, phenotype = idents, gene = "KLF2")
plotCytoGenes(results.cyto, numOfGenes = 10)

# Cytotrace
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, cells = sample(Cells(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), 25000))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == "G1")

Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub$reordered.idents
mat <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, assay = "RNA")
mat <- as.data.frame(mat)
results.cyto <- CytoTRACE(mat, ncores = 1)
idents <- as.character(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub))
names(idents) <- names(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub))
"%!in%" <- function(x, y) !("%in%"(x, y))
idents[which(idents %!in% c("C2-Stem/memory-KLF2 Mock D8", "C2-Stem/memory-KLF2 KO D8", "C2-Stem/memory-KLF2 Mock D28", "C2-Stem/memory-KLF2 KO D28"))] <- "other"
idents[which(idents %in% c("C2-Stem/memory-KLF2 Mock D8"))] <- "C2-Mock D8"
idents[which(idents %in% c("C2-Stem/memory-KLF2 Mock D28"))] <- "C2-Mock D28"
idents[which(idents %in% c("C2-Stem/memory-KLF2 KO D8"))] <- "C2-KO D8"
idents[which(idents %in% c("C2-Stem/memory-KLF2 KO D28"))] <- "C2-KO D28"
emb <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, reduction = "umap"))
plotCytoTRACE(results.cyto, emb = emb, phenotype = idents, colors = mycolors.CD4)



# Table of number per mouse
numbers <- table(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$sample)
write.csv(file = "CD4_Cells_per_cluster_per_mouse.csv", x = numbers)

# DEGs C2 / D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
degs.D28.c2.CD4 <- FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, subset.ident ="C2-Stem/memory-KLF2", ident.1 = "KO", ident.2 = "Mock", group.by = "genotype", assay = "RNA")
write.csv(degs.D28.c2.CD4, "degs_D28_c2_CD4.txt", quote = F, row.names = T)

# Heatmap genes per cluster, same genes as nebulosa
genes <- c("IL7R","CD27",  "KLRG1", "TCF7", "KLF2", "LEF1", "SELL", "CCR7", "FGFBP2", "FCGR3A", "NKG7", "GZMK", "GZMB", "GZMA", "GNLY", "LAG3", "ENTPD1", "HAVCR2", "PDCD1", "TOX")
Average <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = genes, return.seurat = T, assays = "RNA")
Average <- ScaleData(object = Average, features = genes)
plot(DoHeatmap(Average, features = genes, draw.lines = F, size = 3, group.colors = mycolors.CD4, ) + scale_fill_viridis_c(option = "inferno")) 

# Send DEGs of all clusters, Mock vs KO TOTAL
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$celltype.condition <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype, sep  = "_")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$celltype.condition
tmp <- c()
for (i in 1:length(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents))){ # or however many clusters you have
  try({
    ident1 <- paste(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents)[i],"Mock", sep = "_")
    ident2 <- paste(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents)[i],"KO", sep = "_")
    tmp.data <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = c(ident1, ident2))
    condition.diffgenes <- FindAllMarkers(tmp.data, min.pct=0.25, logfc.threshold=0.25, only.pos = T)
    tmp <- rbind(tmp, condition.diffgenes)
  })
}
write.csv(tmp, file="DEGs.csv")

# Send DEGs of all clusters, Mock vs KO SPlit by timepoint
## D8
tmp.D8 <- c()
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
for (i in 1:length(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$Idents))){ # or however many clusters you have
  try({
    ident1 <- paste(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$Idents)[i],"Mock", sep = "_")
    ident2 <- paste(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$Idents)[i],"KO", sep = "_")
    tmp.data <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, idents = c(ident1, ident2))
    condition.diffgenes <- FindAllMarkers(tmp.data, min.pct=0.25, only.pos = T)
    tmp.D8 <- rbind(tmp.D8, condition.diffgenes)
  })
}
write.csv(tmp.D8, file="DEGs_D8.csv")

## D28
tmp.D28 <- c()
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
for (i in 1:length(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28$Idents))){ # or however many clusters you have
  try({
    ident1 <- paste(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28$Idents)[i],"Mock", sep = "_")
    ident2 <- paste(levels(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28$Idents)[i],"KO", sep = "_")
    tmp.data <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, idents = c(ident1, ident2))
    condition.diffgenes <- FindAllMarkers(tmp.data, assay = "RNA", only.pos = T)
    tmp.D28 <- rbind(tmp.D28, condition.diffgenes)
  })
}
write.csv(tmp.D28, file="DEGs_D28.csv")

# Rename KO to gSUV
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned[["genotype"]] <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$orig.ident
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype == "Mock.1")] <- "Mock"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype == "Mock.2")] <- "Mock"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype == "Mock.3")] <- "Mock"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype == "Suv.ko.1")] <- "gSUV"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype == "Suv.ko.2")] <- "gSUV"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data$genotype == "Suv.ko.3")] <- "gSUV"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype, levels = c("Mock", "gSUV"))

# New columns for cluster number only
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Active translation-IL7R", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Stem/memory-KLF2", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Stem/cytotoxic-GNLY", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Resident Memory-ITGA1", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Activated-CD69", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Activated-JUN", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Active translation-RPS27", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Effector-KLRB1", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Cytosolic stress-HSPA1B", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Cytosolic stress-LncRNA", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Pro-apopototic-BAX", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Cycling/G2M-HSPA1B", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Cycling-TOP2A", replacement = "")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- gsub(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, pattern = "-Cycling/S-PCNA", replacement = "")

# Rechcking metabolism
all_sign <- c("sign_ox_ph","sign_glycolysis_KEGG","sign_activation", "sign_FA_metabolism_GO", "sign_fatty_acid_metab_KEGG","sign_hypoxia", "sign_core_DNA_repair","sign_er_stress2")
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}

FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ncol = 4, features = c("sign_ox_ph1", "sign_glycolysis_KEGG1","sign_activation1", "sign_FA_metabolism_GO1", "sign_fatty_acid_metab_KEGG1","sign_hypoxia1", "sign_core_DNA_repair1","sign_er_stress21"), 
            pt.size = 0.8, min.cutoff = "q1", max.cutoff = "q99", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes() & NoLegend()


# FP hallmarks
hallmarks_orig <- hallmarks_list
hallmarks <- NULL
hallmarks_name <- NULL
pdf("Featureplot_hallmarks_CD4.pdf", width = 12, height = 12)
for (i in 1:50){
  hallmarks <- hallmarks_orig[,-1:-2][i,]
  names(hallmarks) <- NULL
  hallmarks_name <- hallmarks_orig$V1[i]
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(hallmarks), name = hallmarks_name, ctrl = 20)
  try(plot(VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = paste0(hallmarks_name,"1"), split.by = "genotype_timepoint", pt.size = 0, idents = c("C2-Stem/memory-KLF2"))+ theme(axis.text.x = element_text(hjust = 1, size = 25))))
}
dev.off()

#### GENERATE FIGURES
#####
######
png(filename = "Signatures_CD4.png", res = 300, height = 3200, width = 3200)
plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ncol = 2, c("Cytotoxic_Effector", "Dysfunction_Exhaustion", "Stemness_Memory", "Tissue_Resident"), pt.size = 0.4, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes())
dev.off()

# Nebulosa
genes <- c("BCL2","SATB1","ID3","KLRG1","ZNF683","ID2", "GNLY","GZMK","CX3CR1", "HAVCR2","ENTPD1","S1PR1","RBPJ","FOS","JUN","JUNB","CD69")
genes <- c("TCF7", "KLF2", "LEF1", "CCR7","SELL","CD27", "IL7R", "GZMB", "GZMA", "PDCD1",  "LAG3", "TOX")
genes <- c("BCL2","GZMK")
pdf(file = "Nebulosa_CD4.pdf", height = 3, width = 3)
for (i in 1:length(genes)){
  plot(plot_density(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = genes[i],  reduction = "umap") & NoLegend() & NoAxes())
}
dev.off()

# Heatmap proportions per timepoint
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
freq_table <- as.data.frame(prop.table(
  x = table(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$cluster.number, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8@meta.data[, "genotype"]),
  margin = 1
))
freq_table$Var1 <- factor(x = freq_table$Var1, levels = paste0(rep("C",14),1:14))
freq_table <- dcast(freq_table, Var1 ~ Var2)
rownames(freq_table) <- freq_table$Var1
tmp <- rownames(freq_table)
freq_table <- freq_table[, -1]
freq_table <- as.matrix(freq_table)
pheatmap(freq_table, scale = "none", fontsize = 12, main = "CD4-D8", cluster_rows = F, cluster_cols = F, color = viridis(50))

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
freq_table <- as.data.frame(prop.table(
  x = table(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28$cluster.number, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28@meta.data[, "genotype"]),
  margin = 1
))
freq_table$Var1 <- factor(x = freq_table$Var1, levels = paste0(rep("C",14),1:14))
freq_table <- dcast(freq_table, Var1 ~ Var2)
rownames(freq_table) <- freq_table$Var1
tmp <- rownames(freq_table)
freq_table <- freq_table[, -1]
freq_table <- as.matrix(freq_table)
pheatmap(freq_table, scale = "none", fontsize = 12, main = "CD4-D28", cluster_rows = F, cluster_cols = F, color = viridis(50))

# VP
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")
pdf("VlnPlot_resident.pdf", height = 6, width = 4)
plot(VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "Dysfunction_Exhaustion",pt.size = 0, group.by = "genotype_timepoint", cols = c("darkgrey","firebrick3","black","red"))  + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.4
  ) & labs(title = "Dysfunction Exhaustiont", size = 13))
dev.off()

# Cytotrace
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == "G1")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, cells = sample(Cells(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), 20000))
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub$reordered.idents
mat <- GetAssayData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, assay = "RNA")
mat <- as.data.frame(mat)
results.cyto <- CytoTRACE(mat, ncores = 1)
idents <- as.character(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub))
names(idents) <- names(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub))
"%!in%" <- function(x, y) !("%in%"(x, y))
idents[which(idents %!in% c("C2-Stem/memory-KLF2 Mock D8", "C2-Stem/memory-KLF2 gSUV D8", "C2-Stem/memory-KLF2 Mock D28", "C2-Stem/memory-KLF2 gSUV D28",
                            "C4-Resident Memory-ITGA1 Mock D8", "C4-Resident Memory-ITGA1 gSUV D8", "C4-Resident Memory-ITGA1 Mock D28", "C4-Resident Memory-ITGA1 gSUV D28"))] <- "other"
idents[which(idents %in% c("C2-Stem/memory-KLF2 Mock D8"))] <- "C2-Mock D8"
idents[which(idents %in% c("C2-Stem/memory-KLF2 Mock D28"))] <- "C2-Mock D28"
idents[which(idents %in% c("C2-Stem/memory-KLF2 gSUV D8"))] <- "C2-gSUV D8"
idents[which(idents %in% c("C2-Stem/memory-KLF2 gSUV D28"))] <- "C2-gSUV D28"
idents[which(idents %in% c("C4-Resident Memory-ITGA1 Mock D8"))] <- "C4-Mock D8"
idents[which(idents %in% c("C4-Resident Memory-ITGA1 Mock D28"))] <- "C4-Mock D28"
idents[which(idents %in% c("C4-Resident Memory-ITGA1 gSUV D8"))] <- "C4-gSUV D8"
idents[which(idents %in% c("C4-Resident Memory-ITGA1 gSUV D28"))] <- "C4-gSUV D28"
emb <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, reduction = "umap"))
plotCytoTRACE(results.cyto, emb = emb, phenotype = idents)

# Reimport Cytotrace table
cytotrace.meta <- read.csv("CytoTRACE_plot_table_CD4.txt", sep = "\t", row.names = 1)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- AddMetaData(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, cytotrace.meta)
pdf("CytoTRACE_CD4.pdf", height = 7, width = 7)
plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, features = "CytoTRACE", pt.size = 0.9, order = T) & scale_color_paletteer_c("pals::coolwarm"))
dev.off()

# VlnPlot with CytoTRACE
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, idents = c("C2-Stem/memory-KLF2 Mock D8", "C2-Stem/memory-KLF2 gSUV D8", "C2-Stem/memory-KLF2 Mock D28", "C2-Stem/memory-KLF2 gSUV D28",
                                                                                                                                                                "C4-Resident Memory-ITGA1 Mock D8", "C4-Resident Memory-ITGA1 gSUV D8", "C4-Resident Memory-ITGA1 Mock D28", "C4-Resident Memory-ITGA1 gSUV D28"))
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4) <- factor(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4), levels = c("C2-Stem/memory-KLF2 Mock D8", "C2-Stem/memory-KLF2 gSUV D8", "C2-Stem/memory-KLF2 Mock D28", "C2-Stem/memory-KLF2 gSUV D28",
                                                                                                                                                                                      "C4-Resident Memory-ITGA1 Mock D8", "C4-Resident Memory-ITGA1 gSUV D8", "C4-Resident Memory-ITGA1 Mock D28", "C4-Resident Memory-ITGA1 gSUV D28"))

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4$Idents <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4$Idents <- recode_factor(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4$Idents, 'C2-Stem/memory-KLF2 Mock D8' =  "C2-Mock D8", "C2-Stem/memory-KLF2 gSUV D8"="C2-gSUV D8","C2-Stem/memory-KLF2 Mock D28"="C2-Mock D28",
                                                                                                 "C2-Stem/memory-KLF2 gSUV D28"="C2-gSUV D28", "C4-Resident Memory-ITGA1 Mock D8"="C4-Mock D8",  "C4-Resident Memory-ITGA1 gSUV D8"= "C4-gSUV D8", "C4-Resident Memory-ITGA1 Mock D28"="C4-Mock D28","C4-Resident Memory-ITGA1 gSUV D28"="C4-gSUV D28")
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub.C2.C4, "CytoTRACE", group.by = "Idents", pt.size = 0, cols = rep(c("darkgrey", "firebrick3", "black", "red"), 2))  + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14)) &
  geom_signif(
    comparisons = list(c("C2-Mock D8", "C2-gSUV D8"), c("C2-Mock D28", "C2-gSUV D28"), c("C4-Mock D8", "C4-gSUV D8"), c("C4-Mock D28", "C4-gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.9
  ) & labs(title = "CytoTRACE - CD4 clusters", size = 13)

# Plot on total clusters, comparing conditions
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, "CytoTRACE", group.by = "genotype_timepoint", pt.size = 0, cols = rep(c("darkgrey", "firebrick3", "black", "red"), 2))  + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14)) &
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.9
  ) & labs(title = "CytoTRACE - CD4 - Total clusters", size = 13)


# Redo heatmap with simpler cluster names
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents <- Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents <- recode_factor(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents, "C1-Active translation-IL7R" = "C1", 'C2-Stem/memory-KLF2' ="C2","C3-Stem/cytotoxic-GNLY" = "C3", "C4-Resident Memory-ITGA1"="C4", "C5-Activated-CD69" = "C5",
                                                                                      "C6-Activated-JUN" = "C6", "C7-Active translation-RPS27" = "C7",  "C8-Effector-KLRB1" = "C8", "C9-Cytosolic stress-HSPA1B" ="C9",  "C10-Cytosolic stress-LncRNA"="C10", "C11-Pro-apopototic-BAX"="C11",  "C12-Cycling/G2M-HSPA1B" = "C12","C13-Cycling-TOP2A"="C13", "C14-Cycling/S-PCNA"="C14")  
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents
markers <- read.csv(file = "markers_CD4s.txt")
markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

Average <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = unique(top5$gene), return.seurat = T, assays = "RNA")
Average <- ScaleData(object = Average, features = unique(top5$gene))
plot(DoHeatmap(Average, features = unique(top5$gene), draw.lines = F, size = 3, group.colors = mycolors.CD4) + scale_fill_viridis_c(option = "inferno"))

# Replot FP
pdf("FP_CD4.pdf", height = 18, width = 29)
plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ncol = 4, c("CD8_stem_sc", "Gueguen_CD8_FCGR3A", "sign_effector1", "sign_cycling1", "sign_ICP_negative1", "fraietta_cr1", "fraietta_nr1", "sign_ox_ph1", "sign_glycolysis1", "sign_FA_metabolism_GO1","Tissue_Resident","Tissue_homing1"), pt.size = 1, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes() &  theme(title = element_text(hjust = 1, size = 25)))
dev.off()

# FP Fraietta split by condition
pdf("FP_Fraietta_split_CD4.pdf", height = 6, width = 25)
FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = "fraietta_cr1", split.by = "genotype_timepoint", order = T, reduction = "umap",
  pt.size = 1.2, min.cutoff = "q1", max.cutoff = "q99") & scale_color_paletteer_c("pals::coolwarm") & RestoreLegend()
dev.off()

# Fraietta sign only in c2 - KLF2
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")
pdf("VP_Fraietta_C2_CD4.pdf", height = 7, width = 5)
plot(VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, "fraietta_cr", group.by = "genotype_timepoint", pt.size = 0, cols = c("darkgrey","firebrick3","black","red"))  + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.10) & labs(title = "Fraietta Complete Response", size = 13))& theme(axis.title.x = element_blank()) 
dev.off()

# Is the stem signature enriched in cycling cells

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Phase_condition <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Phase, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype_timepoint, sep = "-")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Phase_condition <- factor(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Phase_condition, levels = c("S-Mock D8", "S-gSUV D8","S-Mock D28", "S-gSUV D28","G2M-Mock D8", "G2M-gSUV D8","G2M-Mock D28", "G2M-gSUV D28"))

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.cycling <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == c("G1"), invert = T)
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.cycling, features = "Stemness_Memory", group.by = "Phase_condition", pt.size = 0, cols = c("darkgrey","firebrick3","black","red","darkgrey","firebrick3","black","red")) + NoLegend() + stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) +
  geom_signif(
    comparisons = list(c("S-Mock D8", "S-gSUV D8"), c("S-Mock D28", "S-gSUV D28"),c("G2M-Mock D8", "G2M-gSUV D8"), c("G2M-Mock D28", "G2M-gSUV D28")),map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.5) & labs(title = "Stemness memory in cycling cells", size = 13)

# DP
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number <- factor(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$cluster.number, levels = paste0(rep("C",14), 1:14))
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "cluster.number", pt.size = 0.9, label = T, repel = T, label.size = 5, label.box = F, cols = mycolors.CD4) + NoLegend() 
p <- DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "cluster.number", pt.size = 1, label = F, repel = T, label.size = 7, label.box = F, cols = mycolors.CD4) & NoLegend() 
LabelClusters(p, id = "cluster.number",  fontface = "bold", color = "black", size = 8, repel = T)

# Markers with D8 and D28 combined
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents_genotype <- paste(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype, sep = "_")
tmp <- c()
for (i in 1:length(levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned)))){ 
  try({
    ident1 <- paste(levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))[13],"Mock", sep = "_")
    ident2 <- paste(levels(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28))[13],"gSUV", sep = "_")
    tmp.data <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, subset = Idents_genotype == c(ident1, ident2))
    Idents(tmp.data) <- tmp.data$Idents_genotype
    condition.diffgenes <- FindMarkers(tmp.data, min.pct=0.25, only.pos = F, ident.1 = ident1, ident.2 = ident2, test.use = "DESeq2")
    tmp <- rbind(tmp, condition.diffgenes)
  })
}
write.csv(tmp, file="DEGs_CD4_merged_timepoints_LR.csv")


# Metafeature test
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- MetaFeature(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = c("KLF2","IL7R","BCL2","SELL","CCR7"), meta.name = "stem", assay = "RNA")
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "stem", min.cutoff = "q1", max.cutoff = "q99", order = T)& scale_color_paletteer_c("pals::coolwarm") & ggtitle("Stem Meta gene")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = c("KLF2","IL7R","BCL2","SELL","CCR7"), name = "stem_ms", assay = "RNA")
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "stem_ms1", min.cutoff = "q1", max.cutoff = "q99", order = T)& scale_color_paletteer_c("pals::coolwarm") & ggtitle("Stem Module Score")

# FindConservedMarkers
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$sample <- paste(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$orig.ident, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, sep = ".")
conserved.markers <- FindConservedMarkers(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ident.1 = "C2-Stem/memory-KLF2", grouping.var = "sample")
DotPlot(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = rownames(conserved.markers)) + RotatedAxis() + scale_color_viridis()

# Rename resident clusters
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- RenameIdents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,
                                                                              "C4-Resident Memory-ITGA1" = "C4-Resident-ITGA1")

# Threshold for cycling prediction score, how does it change the output?
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- FindVariableFeatures(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, nfeatures = 4000)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == "G1")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref, ident = c("C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"), invert = T)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = Phase == "G1", invert = T)
tumor.anchors <- FindTransferAnchors(
  reference = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref, query = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query,
  dims = 1:50, reference.assay = "RNA", query.assay = "RNA"
)
predictions <- TransferData(
  anchorset = tumor.anchors, refdata = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.ref),
  dims = 1:50
)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query <- AddMetaData(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, metadata = predictions)
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, "prediction.score.max", pt.size = 0) + scale_color_manual(values = mycolors.CD4)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id <- factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id, levels = c(
  "C1-Active translation-IL7R", "C2-Stem/memory-KLF2", "C3-Stem/cytotoxic-GNLY", "C4-Resident-ITGA1", "C5-Activated-CD69", "C6-Activated-JUN", "C7-Active translation-RPS27",
  "C8-Effector-KLRB1", "C9-Cytosolic stress-HSPA1B", "C10-Cytosolic stress-LncRNA", "C11-Pro-apopototic-BAX", "C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA"
))
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$seurat_clusters <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, group.by = "predicted.id", label = F, repel = T) + scale_color_manual(values = mycolors.CD4)

# Filter by prediction score
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id[which(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$prediction.score.max < 0.3)] <- "filtered"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$seurat_clusters <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query$predicted.id
plot_stat(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query, plot_type = "prop_fill", group_by = "genotype_timepoint") + RotatedAxis() + scale_fill_manual(values = c(mycolors.CD4[1:6], mycolors.CD4[9:10])) + theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16, color = "black"))

# Do DEGs Venns with Venny 
markers <- FindAllMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, assay = "RNA", features = TCR.removed, only.pos = T)
write.csv(markers, "markers_CD4s.txt", quote = F, row.names = F)

# Volcanos with red gene names D8
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C13-Cycling-TOP2A")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, subset = timepoint == "D28")

degs.D8 <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.1, assay = "RNA"))
degs.D8.bis <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.25, assay = "RNA"))
colors <- rep("black", 128)
colors[which(degs.D8.bis$avg_log2FC > 0 & degs.D8.bis$p_val_adj < 0.001)] <- "salmon"
colors[which(degs.D8.bis$avg_log2FC < 0 & degs.D8.bis$p_val_adj < 0.001)] <- "grey50"
#colors[which(rownames(degs.D8.bis) %in% c("KLF2","LEF1","S1PR1","KLF3","SELL","CD27","TCF7"))] <- "red"
#colors[which(rownames(degs.D8.bis) %in% c("CCL5","ITGA1","CST7","LAG3","ID2","CD2","RBPJ","KLRG1"))] <- "black"
colors[which(rownames(degs.D8.bis) %in% c("KLF2","LEF1","S1PR1","KLF3","SELL","TCF7"))] <- "red"
colors[which(rownames(degs.D8.bis) %in% c("CCL5","ITGA1","CST7","LAG3","CD2","RBPJ"))] <- "black"

try(plot(EnhancedVolcano(degs.D8,
                         lab = rownames(degs.D8), labCol = colors,
                         x = "avg_log2FC",
                         y = "p_val_adj",
                         drawConnectors = T, captionLabSize = 0, labSize = 6, maxoverlapsConnectors = 18, caption = F, legendLabSize = 0, legendPosition = "none", subtitle = " ",
                         pCutoff = 0.001,
                         FCcutoff = 0.25, xlim = c(-3, 3), title = "CD4 - C13/TOP2A - gSUV vs Mock - D28"
)))

# Volcanos with red gene names D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, subset = timepoint == "D28")

degs.D28 <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.1, assay = "RNA"))
degs.D28.bis <- try(FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, ident.1 = "gSUV", ident.2 = "Mock", group.by = "genotype", logfc.threshold = 0.25, assay = "RNA"))
colors <- rep("black", 128)
colors[which(degs.D28.bis$avg_log2FC > 0 & degs.D28.bis$p_val_adj < 0.001)] <- "salmon"
colors[which(degs.D28.bis$avg_log2FC < 0 & degs.D28.bis$p_val_adj < 0.001)] <- "grey50"
colors[which(rownames(degs.D28.bis) %in% c("KLF2","LEF1","S1PR1","KLF3","SELL","CD27","TCF7"))] <- "red"
colors[which(rownames(degs.D28.bis) %in% c("CCL5","ITGA1","CST7","LAG3","ID2","CD2","RBPJ"))] <- "black"

try(plot(EnhancedVolcano(degs.D28,
                         lab = rownames(degs.D28), labCol = colors,
                         x = "avg_log2FC",
                         y = "p_val_adj",
                           drawConnectors = T, captionLabSize = 0, labSize = 6, maxoverlapsConnectors = 18, caption = F, legendLabSize = 0, legendPosition = "none", subtitle = " ",
                         pCutoff = 0.001,
                         FCcutoff = 0.25, xlim = c(-1.1, 1.1), title = "CD4 - C2/KLF2 - gSUV vs Mock - D28"
)))

# Retry pathway enrich with Seurat
library(enrichR)
dbs <- listEnrichrDbs()

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub$genotype
DEenrichRPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.sub, ident.1 = "Mock", ident.2 = "gSUV", assay = "RNA",
              enrich.database = "Azimuth_Cell_Types_2021", max.genes = 500)

# Show total numbers instead of proportions after cycling label transfer
library(reshape2)
tmp <- as.data.frame(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.query@meta.data)
tmp <- tmp %>% group_by(predicted.id) %>% dplyr::count(genotype_timepoint)

ggplot(tmp, aes(x = genotype_timepoint, y = n, fill = factor(predicted.id))) + geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12)) + 
  labs(y = "Number of cells", size = 7) + scale_fill_manual(values = c(mycolors.CD4[1:7], mycolors.CD4[9:10]))+ theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16, color = "black"))

# Reorder Idents
order <- c("C1-Active translation-IL7R", "C2-Stem/memory-KLF2", "C3-Stem/cytotoxic-GNLY", "C4-Resident-ITGA1", "C5-Activated-CD69", "C6-Activated-JUN", "C7-Active translation-RPS27",
  "C8-Effector-KLRB1", "C9-Cytosolic stress-HSPA1B", "C10-Cytosolic stress-LncRNA", "C11-Pro-apopototic-BAX", "C12-Cycling/G2M-HSPA1B", "C13-Cycling-TOP2A", "C14-Cycling/S-PCNA")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) <- factor(x = Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), levels = order)

# Vlnplot per condition for chosen genes
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, c("SELL", "KLF2", "CCL4","CCL5", "CD27", "KLRG1"), pt.size = 0, group.by = "genotype_timepoint") & NoLegend() &
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position =2,
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)

# 4 way representation for DEGs between both genotypes and both timepoints (bulk level and for C2/KLF2)
library(ggrepel)
library(presto)
genes.timepoint <- FindMarkers(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "timepoint", ident.1 = "D28",ident.2 = "D8", assay = "RNA", logfc.threshold = 0, min.pct = 0)
genes.genotype <- FindMarkers(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "genotype", ident.1 = "gSUV",ident.2 = "Mock", assay = "RNA", logfc.threshold = 0, min.pct = 0)

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")
CD4.genes.genotype <- wilcoxauc(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, group_by = "genotype")
CD4.genes.genotype <- CD4.genes.genotype %>% filter(group == "gSUV")
CD4.genes.timepoint <- wilcoxauc(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, group_by = "timepoint")
CD4.genes.timepoint <- CD4.genes.timepoint %>% filter(group == "D28")

tmp <- c()
tmp$"avg_log2FC_timepoint" <- CD4.genes.timepoint[,"logFC"]
tmp$"avg_log2FC_genotype" <- CD4.genes.genotype[,"logFC"]
tmp <- as.data.frame(tmp)
rownames(tmp) <- make.unique(CD4.genes.genotype$feature)
rownames(tmp)
names <- tmp %>% filter(avg_log2FC_timepoint > 0.25 | avg_log2FC_timepoint < -0.25 | avg_log2FC_genotype < -0.25| avg_log2FC_genotype > 0.25)
gene.names <- rownames(names)
ggplot(data = tmp) +
  geom_point(data = tmp, aes(x = avg_log2FC_timepoint, y = avg_log2FC_genotype), size = 1, color = ifelse(tmp$avg_log2FC_timepoint>0.25 | tmp$avg_log2FC_genotype>0.25, "red", "black")) + theme_bw() + 
  geom_text_repel(data = names, aes(x=avg_log2FC_timepoint, y=avg_log2FC_genotype, label = gene.names)) + geom_hline(yintercept=c(-0.25,0.25)) + geom_vline(xintercept=c(-0.25,0.25)) + ggtitle("CD4 DEGs - D8 vs D28 - Mock vs gSUV - C2/KLF2")

# Senescence / apoptosis signature
all_sign <- c("sign_core_DNA_repair", "sign_senescence", "sign_apoptosis", "sign_DNA_repair", "sign_DNA_damage_repair")
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}

# As Featureplots
pdf("Signatures_senescence_apoptosis.pdf", width = 12, height = 12)
for (i in 1:length(all_sign2)) {
  try(plot(FeaturePlot(
    object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = all_sign2[i], order = T, reduction = "umap",
    pt.size = 1.5, min.cutoff = "q1", max.cutoff = "q99"
  ) + labs(title = all_sign[i], size = 13) + scale_color_paletteer_c("pals::coolwarm")))
}
dev.off()
pdf("VP_D8.pdf", width = 12, height = 6)
for (i in all_sign2) {
  try(plot(VlnPlot(subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8"), i, split.by = "genotype", pt.size = 0) + ylim(-0.5, 0.2)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
             geom_signif(
               comparisons = list(c("Mock", "gSUV")),
               map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.1
             )))
}
dev.off()

# DEGs bulk like 
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype
bulk.markers <- FindMarkers(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, logfc.threshold = 0.1, ident.1 = "gSUV", ident.2 = "Mock")
bulk.markers.bis <- FindMarkers(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, logfc.threshold = 0.25, assay = "RNA", ident.1 = "gSUV", ident.2 = "Mock")

colors <- rep("black", 17)
colors[which(bulk.markers.bis$avg_log2FC > 0 & bulk.markers.bis$p_val_adj < 0.01)] <- "salmon"
colors[which(bulk.markers.bis$avg_log2FC < 0 & bulk.markers.bis$p_val_adj < 0.01)] <- "grey50"
colors[which(rownames(bulk.markers.bis) %in% c("KLF2","S1PR1","SELL","KLF3","LEF1"))] <- "red"
colors[which(rownames(bulk.markers.bis) %in% c("RBPJ","GZMK","GZMA","CST7"))] <- "black"

try(plot(EnhancedVolcano(bulk.markers,
                         lab = rownames(bulk.markers), labCol =  colors, 
                         x = "avg_log2FC",
                         y = "p_val_adj",
                         drawConnectors = T, captionLabSize = 0, labSize = 6, maxoverlapsConnectors = 18, caption = F, legendLabSize = 0, legendPosition = "none", subtitle = " ",
                         pCutoff = 0.01,
                         FCcutoff = 0.25, xlim = c(-1.2, 1.2), title = "CD4 - gSUV vs Mock - bulk-like DEGs "
)))

# Sankey plot
plotSankey<-function(seuratObj,idvar=c("varRes.0.3","emt_res.0.3")){
  require(flipPlots)
  message('try install_github("Displayr/flipPlots") if this doesnt work')
  require(dplyr)
  seuratObj@meta.data[,match(idvar,colnames(seuratObj@meta.data))] %>% arrange(.[,1]) %>% group_by_all() %>% summarise(COUNT = n()) ->> my.data
  #my.data<-as.factor(my.data[,1])
  
  SankeyDiagram(my.data[, -grep("COUNT",colnames(my.data))],link.color = "Source",weights = my.data$COUNT,max.categories = 1000)
}
plotSankey(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idvar = c("Phase","genotype_timepoint"))

# Plot signature from Lowery et al. 2022 in batch
Lowery.sign <- readxl::read_xlsx(path = "E:/Work/Signatures/Lowery_et_al_2022/abl5447_Table_S4.xlsx")
colnames(Lowery.sign) <- gsub("[\r\n]", "", colnames(Lowery.sign))
colnames(Lowery.sign) <- make.names(colnames(Lowery.sign))
pdf("Lowery_2022_sign_CD4.pdf", height = 6, width = 7)
for (i in colnames(Lowery.sign)) {
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Lowery.sign[[i]]), name = i, assay = "RNA")
  #plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, paste0(i, 1), pt.size = 0.8, min.cutoff = "q3", max.cutoff = "q97", order = T) + scale_color_paletteer_c("pals::coolwarm") +ggtitle(i))
}
dev.off()

# Same with Ucell
library(UCell)
pdf("Lowery_2022_sign_CD4_UCell.pdf", height = 6, width = 7)
for (i in colnames(Lowery.sign)) {
  suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore_UCell(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Lowery.sign[[i]]), name = i, assay = "RNA")
  plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, paste0( "signature_1",i), pt.size = 0.8, min.cutoff = "q3", max.cutoff = "q97", order = T) + scale_color_paletteer_c("pals::coolwarm") +ggtitle(i))
}
dev.off()

# Scatterplot matrix to show correlations between signatures / genes
library(psych)
mat <- cbind(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Stemness_Memory, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Cytotoxic_Effector, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Tissue_Resident,
              suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Tissue_homing, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$fraietta_cr, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$fraietta_nr)
colnames(mat) <- c("Stemness_Memory","Cytotoxic_Effector","Tissue_Resident", "Tissue_homing", "Fraietta CR", "Fraietta NR")
pdf("scatterplot_matrix_signatures.pdf", width = 9, height = 9)
plot(pairs.panels(mat,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
))
dev.off()

# VP on KLF2 cluster
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Lowery.sign[["Krishna.ACT.Stem.Like"]]), name = "Krishna.ACT.Stem.Like", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Lowery.sign[["Caushi.Stem.like.memory"]]), name = "Caushi.Stem.like.memory", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Lowery.sign[["Caushi.CD8.Stem.like.memory"]]), name = "Caushi.CD8.Stem.like.memory", assay = "RNA")

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2") 
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, c("Krishna.ACT.Stem.Like1"), pt.size = 0, group.by = "genotype_timepoint",  cols = c("darkgrey", "firebrick3", "black", "red")) & theme(axis.title.x = element_blank()) &
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position =0.7,
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) & NoLegend() & ggtitle("Krishna ACT Stem Like")
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, c("Caushi.Stem.like.memory1"), pt.size = 0, group.by = "genotype_timepoint",  cols = c("darkgrey", "firebrick3", "black", "red")) & theme(axis.title.x = element_blank())&
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position =0.6,
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)& NoLegend() & ggtitle("Caushi Stem like memory")
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, c("Caushi.CD8.Stem.like.memory1"), pt.size = 0, group.by = "genotype_timepoint",  cols = c("darkgrey", "firebrick3", "black", "red")) & theme(axis.title.x = element_blank()) &
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position =0.6,
  ) & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95)& NoLegend() & ggtitle("Caushi CD8 Stem like memory")

# VP split per 4 conditions
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.cleaned, 
    c("KLF2"),
    pt.size = 0.01,
    split.by = "genotype_timepoint", y.max = 5,
    cols = c("darkgrey", "firebrick3", "black", "red")
  ) & theme(axis.title.x = element_blank()) &
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1,
    size = 12
  )) & geom_signif(
  comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
  map_signif_level = TRUE,
  textsize = 6,
  test = "wilcox.test",
  y_position = 4.5,
) 

# Selected FP signatures
all_sign <- c("Central_Memory", "CD8_stem_sc", "Gueguen_CD8_FCGR3A", "sign_effector", "sign_cycling","Cytotoxic_Effector", "sign_ICP_negative", "fraietta_cr", "fraietta_nr", 
              "sign_ox_ph","Dysfunction_Exhaustion", "sign_glycolysis", "Caushi.Stem.like.memory1", "Tissue_Resident","Eff_like_Miller1")
all_sign_names <- c("Central Memory", "CD8 stem single-cell", "Gueguen CD8 FCGR3A+", "Effector", "Cycling","Cytotoxic Effector", "Negative ICPs", "Fraietta CR", "Fraietta NR", 
              "Oxydative Phosphorylation","Dysfunction Exhaustion", "Glycolysis", "Caushi Stem-like memory", "Tissue Resident","Eff-like Miller")
all_sign2 <- all_sign
all_sign2 <- paste0(all_sign, "1")
for (i in 1:length(all_sign)) {
  try(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(get(all_sign[i])), name = all_sign[i], assay = "RNA"))
}
gg_Fig <- plot(FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ncol = 5, c("Central_memory1","CD8_stem_sc", "Gueguen_CD8_FCGR3A", "sign_effector1", "sign_cycling1", "Cytotoxic_Effector",
                                                                                                      "sign_ICP_negative1", "fraietta_cr1", "fraietta_nr1", "sign_ox_ph1", "Dysfunction_Exhaustion","sign_glycolysis1","Caushi.Stem.like.memory1", "Tissue_Resident","Eff_like_Miller1"),
                           pt.size = 0.3, min.cutoff = "q1", max.cutoff = "q99", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes()) & theme(title = element_text(size=24))
gg_Fig <- lapply( 1:length(all_sign_names), function(x) { gg_Fig[[x]] + labs(title=all_sign_names[x]) })
gg_Fig <- CombinePlots( gg_Fig, ncol = 5 )
pdf("FP_CD4.pdf", height = 16, width = 28)
plot(gg_Fig)
dev.off()

# Testing using Milo
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)

sce <- as.SingleCellExperiment(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, assay = "RNA")
sce@metadata <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned@meta.data
sce_milo <- Milo(sce)
sce_milo <- buildGraph(sce_milo, k = 20, d = 30)
sce_milo <- buildFromAdjacency(sce_milo, sce@metadata)
sce_milo <- makeNhoods(sce_milo, k = 20, d = 30, refined = TRUE, prop = 0.2)
sce_milo <- countCells(sce_milo, meta.data = data.frame(colData(sce_milo)), sample = "orig.ident")
sce_milo <- calcNhoodDistance(sce_milo, d = 30)

milo.design <- as.data.frame(xtabs(~ genotype + sample, data = milo.meta))
milo.design <- milo.design[milo.design$Freq > 0, ]
milo.res <- testNhoods(sce_milo, design = ~genotype, design.df = milo.design)
head(milo.res)
sce_milo <- buildNhoodGraph(sce_milo)
plotUMAP(sce_milo) + plotNhoodGraphDA(sce_milo, da_results, alpha = 0.05) +
  plot_layout(guides = "collect")

# Remake Density plot CD4
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("Mock", "gSUV"))
df.a <- df %>% filter(timepoint == "D8")
labels <- c("Mock - 2362 cells", "gSUV - 2313 cells")
names(labels) <- c("Mock", "gSUV")
a <- plot(ggplot(data = df.a) +
            geom_point(data = df.a, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 1) +
            stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 20, geom = "polygon") +
            theme_classic()  +
            xlim(c(-9, 10)) +
            ylim(c(-7, 8))+ scale_fill_paletteer_c("pals::coolwarm") +
            facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
            theme(
              strip.text.x = element_text(
                size = 20, color = "black",
                face = "bold"
              ),
              strip.text.y = element_text(
                size = 20, color = "black",
                face = "bold"
              )
            ))
labels <- c("Mock - 493 cells", "gSUV - 391 cells")
names(labels) <- c("Mock", "gSUV")
df.b <- df %>% filter(timepoint == "D28")
b <- plot(ggplot(data = df.b) +
            geom_point(data = df.b, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
            stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
            theme_classic()  +
            xlim(c(-9, 10)) +
            ylim(c(-7, 8))+ scale_fill_paletteer_c("pals::coolwarm")+
            facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
            theme(
              strip.text.x = element_text(
                size = 20, color = "black",
                face = "bold"
              ),
              strip.text.y = element_text(
                size = 20, color = "black",
                face = "bold"
              )
            ))

a + b + plot_layout(ncol = 1) & NoLegend()

# Side by side ratio plots / Fold change between conditions
#suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <-  subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
query.list <- SplitObject(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, split.by = "genotype")
norm.c <- table(Idents(query.list[["Mock"]]))/sum(table(Idents(query.list[["Mock"]])))
norm.q <- table(Idents(query.list[["gSUV"]]))/sum(table(Idents(query.list[["gSUV"]])))
foldchange <- norm.q/norm.c
foldchange <- sort(foldchange,decreasing = T)
tb.m <- melt(foldchange)
colnames(tb.m) <- c("Cluster","Fold_change")
pll <- list()
a <- ggplot(tb.m, aes(x=Cluster, y=Fold_change, fill=Cluster)) + geom_bar(stat="identity") + theme_bw()  +
  scale_fill_manual(values=c(mycolors.CD4[3],mycolors.CD4[6],mycolors.CD4[2],mycolors.CD4[11],mycolors.CD4[1],mycolors.CD4[7],mycolors.CD4[10],mycolors.CD4[9],mycolors.CD4[8],mycolors.CD4[13],mycolors.CD4[12],mycolors.CD4[5],mycolors.CD4[14],mycolors.CD4[4])) + geom_hline(yintercept = 1) + scale_y_continuous(trans='log2', limits = c(0.6,4.3), breaks = seq(-20, 20, by = 1)) + 
  theme(axis.text.x=element_blank(), legend.position="left") + ggtitle("CD4 - gSUV vs Mock - D8") 

#suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <-  subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
query.list <- SplitObject(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, split.by = "genotype")
norm.c <- table(Idents(query.list[["Mock"]]))/sum(table(Idents(query.list[["Mock"]])))
norm.q <- table(Idents(query.list[["gSUV"]]))/sum(table(Idents(query.list[["gSUV"]])))
foldchange <- norm.q/norm.c
foldchange <- sort(foldchange,decreasing = T)
tb.m <- melt(foldchange)
colnames(tb.m) <- c("Cluster","Fold_change")
pll <- list()
b <- ggplot(tb.m, aes(x=Cluster, y=Fold_change, fill=Cluster)) + geom_bar(stat="identity") + theme_bw()  +
  scale_fill_manual(values=c(mycolors.CD4[2],mycolors.CD4[13],mycolors.CD4[3],mycolors.CD4[11],mycolors.CD4[9],mycolors.CD4[5],mycolors.CD4[12],mycolors.CD4[4],mycolors.CD4[8],mycolors.CD4[10],mycolors.CD4[7],mycolors.CD4[1],mycolors.CD4[6],mycolors.CD4[14])) + geom_hline(yintercept = 1) + scale_y_continuous(trans='log2',limits = c(0.6,4.3), breaks = seq(-20, 20, by = 1)) + 
  theme(axis.text.x=element_blank(), legend.position="left") + ggtitle("CD4 - gSUV vs Mock - D28")  + NoLegend()

#suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned.D8 <-  subset(suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned, subset = timepoint == "D8")
query.list <- SplitObject(suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned.D8, split.by = "genotype")
norm.c <- table(Idents(query.list[["Mock"]]))/sum(table(Idents(query.list[["Mock"]])))
norm.q <- table(Idents(query.list[["gSUV"]]))/sum(table(Idents(query.list[["gSUV"]])))
foldchange <- norm.q/norm.c
foldchange <- sort(foldchange,decreasing = T)
tb.m <- melt(foldchange)
colnames(tb.m) <- c("Cluster","Fold_change")
pll <- list()
c <- ggplot(tb.m, aes(x=Cluster, y=Fold_change, fill=Cluster)) + geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c(mycolors.CD8[2],mycolors.CD8[1],mycolors.CD8[4],mycolors.CD8[3],mycolors.CD8[5])) + geom_hline(yintercept = 1) + scale_y_continuous(trans='log2',limits = c(0.8,7.5),breaks = seq(-20, 20, by = 1)) + 
  theme(axis.text.x=element_blank(), legend.position="left") + ggtitle("CD8 - gSUV vs Mock - D8")  

#suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned.D28 <-  subset(suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned, subset = timepoint == "D28")
query.list <- SplitObject(suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned.D28, split.by = "genotype")
norm.c <- table(Idents(query.list[["Mock"]]))/sum(table(Idents(query.list[["Mock"]])))
norm.q <- table(Idents(query.list[["gSUV"]]))/sum(table(Idents(query.list[["gSUV"]])))
foldchange <- norm.q/norm.c
foldchange <- sort(foldchange,decreasing = T)
tb.m <- melt(foldchange)
colnames(tb.m) <- c("Cluster","Fold_change")
pll <- list()
d <- ggplot(tb.m, aes(x=Cluster, y=Fold_change, fill=Cluster)) + geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c(mycolors.CD8[4],mycolors.CD8[5],mycolors.CD8[2],mycolors.CD8[3],mycolors.CD8[1])) + geom_hline(yintercept = 1) + scale_y_continuous(trans='log2',limits = c(0.8,7.5),breaks = seq(-20, 20, by = 1)) + 
  theme(axis.text.x=element_blank(), legend.position="left") + ggtitle("CD8 - gSUV vs Mock - D28")  + NoLegend() 
(a+b) / (c+d)


# Dotplot for shared set of stem genes
DotPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = c("Stemness_Memory","Tissue_Resident","TRM.Milner","Circulating.Milner"), split.by = "genotype", cols = c("red","blue")) + coord_flip() + RotatedAxis() + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) 

# Venn between Tres signature and KLF2 signature
library(ggvenn)
Tres_FC0_5 <- readLines("E:/Work/Signatures/Zheng2022/Tres_FC_0_5.txt")
C2_CD4_D28_gSUV <- readLines("E:/Work/Signatures/Suv39H1_CAR/c2KLF2_CD4_D28_gSUV.txt")

venn.list <- list(Tres_FC0_5, Stemness_Memory, C2_CD4_D28_gSUV)
names(venn.list) <- c("Tres_FC0_5", "Stemness_Memory", "C2_CD4_D28_gSUV")
ggvenn(venn.list)
intersect(venn.list[[1]], intersect(venn.list[[2]], venn.list[[3]]))
intersect(venn.list[[1]],venn.list[[3]])

# Check circulating signature from Milner 2017
sign.circ.milner <- toupper(readLines("E:/Work/Signatures/Milner2017/Core_circulating_Milner.txt"))
sign.trm.milner <- toupper(readLines("E:/Work/Signatures/Milner2017/Core_trm_Milner.txt"))

suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(sign.circ.milner), name = "Circulating.Milner", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Circulating.Milner <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Circulating.Milner1
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "Circulating.Milner", pt.size = 0.4, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes()
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(sign.trm.milner), name = "TRM.Milner", assay = "RNA")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$TRM.Milner <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$TRM.Milner1
FeaturePlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "TRM.Milner", pt.size = 0.4, min.cutoff = "q3", max.cutoff = "q97", order = T) & scale_color_paletteer_c("pals::coolwarm") & NoAxes()

# VPs on Milner signatures
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, "TRM.Milner", group.by = "genotype_timepoint", ncol = 1, fill.by = "ident", pt.size = 0, cols = c("darkgrey","firebrick3","black","red")) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.25
  ) & theme(axis.title.x = element_blank())

# Densities with ratios
library(MASS)
library(reshape2) 
library(scales)
# Calculate the common x and y range for geyser1 and geyser2
xrng = range(c(df.ab$UMAP_1,df.ab$UMAP_2))
yrng = range(c(df.ac$UMAP_1,df.ac$UMAP_2))
# Calculate the 2d density estimate over the common range
df.ab <- df.a %>% filter(genotype == "Mock")
df.ac <- df.a %>% filter(genotype == "gSUV")
d1 = kde2d(df.ab$UMAP_1,df.ab$UMAP_2, lims=c(xrng, yrng), n=200)
d2 = kde2d(df.ac$UMAP_1,df.ac$UMAP_2, lims=c(xrng, yrng), n=200)
# Confirm that the grid points for each density estimate are identical
identical(d1$x, d2$x) # TRUE
identical(d1$y, d2$y) # TRUE
# Calculate the difference between the 2d density estimates
diff12 = d1 
diff12$z = d2$z - d1$z
## Melt data into long format
# First, add row and column names (x and y grid values) to the z-value matrix
rownames(diff12$z) = diff12$x
colnames(diff12$z) = diff12$y
# Now melt it to long format
diff12.m = melt(diff12$z, id.var=rownames(diff12))
names(diff12.m) = c("Duration","Waiting","z")
# Plot difference between geyser2 and geyser1 density
ggplot(diff12.m, aes(Duration, Waiting, z=z, fill=z)) +
  stat_contour(aes(colour=..level..), binwidth=0.001)


# Remake Density plot - common scale for Mock & KO within each time point, add cell number per condition
df <- as.data.frame(Embeddings(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, reduction = "umap"))
df <- cbind(df, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$timepoint, suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
colnames(df) <- c("UMAP_1", "UMAP_2", "timepoint", "genotype")
df$timepoint <- factor(x = df$timepoint, levels = c("D8", "D28"))
df$genotype <- factor(x = df$genotype, levels = c("Mock", "gSUV"))
df.a <- df %>% filter(timepoint == "D8")
table(df.a$genotype)
labels <- c("Mock - 9322 cells", "gSUV - 9629 cells")
names(labels) <- c("Mock", "gSUV")
a <- plot(ggplot(data = df.a) +
            geom_point(data = df.a, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
            stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 20, geom = "polygon") +
            theme_classic() +
            xlim(c(-9, 11)) +
            ylim(c(-7, 8)) +
            scale_fill_viridis_c(option = "plasma") +
            facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
            theme(
              strip.text.x = element_text(
                size = 20, color = "black",
                face = "bold"
              ),
              strip.text.y = element_text(
                size = 20, color = "black",
                face = "bold"
              )
            ))
labels <- c("Mock - 14073 cells", "gSUV - 13611 cells")
names(labels) <- c("Mock", "gSUV")
df.b <- df %>% filter(timepoint == "D28")
table(df.b$genotype)
b <- plot(ggplot(data = df.b) +
            geom_point(data = df.b, aes(x = UMAP_1, y = UMAP_2), color = "lightgrey", size = 0.5) +
            stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level), alpha = (..level..)), bins = 30, geom = "polygon") +
            theme_classic() +
            xlim(c(-9, 11)) +
            ylim(c(-7, 8)) +
            scale_fill_viridis_c(option = "plasma") +
            facet_grid(. ~ genotype, labeller = labeller(genotype = labels)) +
            theme(
              strip.text.x = element_text(
                size = 20, color = "black",
                face = "bold"
              ),
              strip.text.y = element_text(
                size = 20, color = "black",
                face = "bold"
              )
            ))

a + b + plot_layout(ncol = 1) & NoLegend()

# Shared gene program upregulated on all clusters
all.genes.D8 <- read.csv("CD4_DEGs_Mock_vs_KO_D8.csv")
all.genes.D28 <- read.csv("CD4_DEGs_Mock_vs_KO_D28.csv")

# Take genes which are at least in 50% of the clusters
shared.genes.D8 <- all.genes.D8 %>% filter(grepl('KO', cluster))
shared.genes.D8 <- shared.genes.D8 %>% count("gene") %>% filter(freq > 5) %>% arrange(-freq)

shared.genes.D28 <- all.genes.D28 %>% filter(grepl('KO', cluster))
shared.genes.D28 <- shared.genes.D28 %>% count("gene") %>% filter(freq > 5) %>% arrange(-freq)

# Plotting
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <-
  ScaleData(
    suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned,
    features = rownames(
      suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned
    )
  )
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents_genotype <-
  paste(
    suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents,
    suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$genotype)
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents_genotype <-
  factor(x = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned$Idents_genotype,
         levels = paste(rep(c("C1",
                             "C2",
                             "C3",
                             "C4",
                             "C5",
                             "C6",
                             "C7",
                             "C8",
                             "C9",
                             "C10",
                             "C11",
                             "C12",
                             "C13",
                             "C14"), each = 1), rep(c("Mock", "gSUV"), each = 14)))

#Setup KO colors
library(colorspace)
mycolors.CD4.KO <- darken(mycolors.CD4, 0.2)

# D8
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D8")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8$Idents_genotype
Average.D8 <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D8, features = shared.genes.D8$gene, return.seurat = T, assays = "RNA")
Average.D8 <- ScaleData(object = Average.D8, features = shared.genes.D8$gene)
DoHeatmap(Average.D8, features = shared.genes.D8$gene, draw.lines = F, size = 4, group.colors = c(mycolors.CD4,mycolors.CD4.KO)) + scale_fill_viridis_c(option = "inferno") + 
  theme(axis.text.y = element_text(size = 16, colour = "black"))  


library(ggExtra)
ggMarginal(p, type = "histogram", 
           margins = "x")

# D28
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, subset = timepoint == "D28")
Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28) <- suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28$Idents_genotype
Average.D28 <- AverageExpression(object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, features = shared.genes.D28$gene, return.seurat = T, assays = "RNA")
Average.D28 <- ScaleData(object = Average.D28, features = shared.genes.D28$gene)
DoHeatmap(Average.D28, features = shared.genes.D28$gene, draw.lines = F, size = 4, group.colors = c(mycolors.CD4,mycolors.CD4.KO)) + scale_fill_viridis_c(option = "inferno") + 
  theme(axis.text.y = element_text(size = 16, colour = "black"))  

# Test with dotplot
DotPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.D28, features = shared.genes$gene) + scale_color_viridis_c(option = "inferno") + 
  theme(axis.text.y = element_text(size = 16)) + NoLegend() + RotatedAxis()

# UMAP CD3+ with good colors
suv.car.t.integrated.merged.timepoints.premerged <- AddMetaData(suv.car.t.integrated.merged.timepoints.premerged, col.name = "Idents_combined", metadata = c(Idents(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned), Idents(suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned)))
suv.car.t.integrated.merged.timepoints.premerged <- subset(suv.car.t.integrated.merged.timepoints.premerged, subset = Idents_combined %in% NA, invert = T)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged, cols = c(mycolors.CD4, mycolors.CD8), group.by = "Idents_combined")

# Test margins
library(ggExtra)
p <- plot(DoHeatmap(Average.D8, features = shared.genes.D8$gene, draw.lines = F, size = 4, group.colors = c(mycolors.CD4,mycolors.CD4.KO)) + scale_fill_viridis_c(option = "inferno") + 
  theme(axis.text.y = element_text(size = 16, colour = "black")) + NoLegend())

ggMarginal(p, data = as.data.frame(shared.genes.D8),)

# VlnPlot with Tres signatures
Tres.sign <- readLines("E:/Work/Signatures/Zheng2022/Tres_FC_0_5.txt")
Tres.sign <- list(Tres.sign)
names(Tres.sign) <- "Tres.sign"
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore_UCell(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(Tres.sign), assay = "RNA", name ="Tres")
FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = "Tres1_UCell", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "TRES signature (Zhang et al. 2022)") + scale_color_paletteer_c("pals::coolwarm")
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, "Tres1", group.by = "genotype_timepoint", ncol = 1, fill.by = "ident", pt.size = 0, cols = c("darkgrey","firebrick3","black","red")) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 1.5
  ) & theme(axis.title.x = element_blank())
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")

# FC plot with same scale




####
#### REVISIONS


### DNMT3 comparison
library(readxl)
library(UCell)
library(ggsignif)
dnmt3.genes <- read_xlsx(path = "E:/Work/Mnemo_Suv39H1_CART/Signatures/DNMT3 CART/scitranslmed.abh0272_data_file_s1.xlsx", trim_ws = T)
dnmt3.genes <- unlist(dnmt3.genes)

## Add module score for DNMT3 signature
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- AddModuleScore_UCell(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = list(DNMT3 = dnmt3.genes), assay = "RNA", name ="DNMT3", ncores = 6)

#FP
FeaturePlot(
  object = suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, features = "DNMT3DNMT3", order = T, reduction = "umap",
  pt.size = 1, min.cutoff = "q1", max.cutoff = "q99"
) + labs(title = "DNMT3 signature (Prinzing et al. 2021)") + scale_color_paletteer_c("pals::coolwarm")

# VP on KLF2 subset only
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2 <- subset(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, idents = "C2-Stem/memory-KLF2")
VlnPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned.KLF2, "DNMT3DNMT3", group.by = "genotype_timepoint", ncol = 1, fill.by = "ident", pt.size = 0, cols = c("darkgrey","firebrick3","black","red")) & NoLegend() & stat_summary(fun = median, geom = "point", size = 12, colour = "white", shape = 95) &
  geom_signif(
    comparisons = list(c("Mock D8", "gSUV D8"), c("Mock D28", "gSUV D28")),
    map_signif_level = TRUE, textsize = 6, test = "wilcox.test", y_position = 0.44
  ) & theme(axis.title.x = element_blank()) & ggtitle("DNMT3 targets signature")


## Label transfer to detect Tpex
remotes::install_github("carmonalab/ProjecTILs")
library(ProjecTILs)


# Project CD8
ref <- load.reference.map(ref = "E:/Work/Mnemo_Suv39H1_CART/CD8T_human_ref_v1.rds")
suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned <- ProjecTILs.classifier(suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned, ref=ref, filter.cells = F, split.by = "orig.ident")
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned, group.by = "functional.cluster", label = T, cols = ref@misc$atlas.palette, pt.size = 1)
plot.states.radar(ref = ref, query = suv.car.t.integrated.merged.timepoints.premerged.cd8s.cleaned)


# Project CD4
ref <- load.reference.map(ref = "E:/Work/Mnemo_Suv39H1_CART/CD4T_human_ref_v1.rds")
suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned <- ProjecTILs.classifier(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, ref=ref, filter.cells = F)
DimPlot(suv.car.t.integrated.merged.timepoints.premerged.cd4s.cleaned, group.by = "functional.cluster", label = T, cols = ref@misc$atlas.palette, pt.size = 1)



## Redo figure without CR/NR names