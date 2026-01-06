setwd("~/Desktop/multiomics")
library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(dplyr)

data_dir <- "data"

# ============================================================================
# 1) LOAD DATA
# ============================================================================
counts <- Read10X_h5(file.path(data_dir, "filtered_feature_bc_matrix.h5"))
rna_counts <- counts$`Gene Expression`
seurat_obj <- CreateSeuratObject(rna_counts, assay = "RNA")

# Load peaks
peaks <- rtracklayer::import(file.path(data_dir, "atac_peaks.bed"))
seqlevelsStyle(peaks) <- "UCSC"


# Create fragment object
frags <- file.path(data_dir, "atac_fragments.tsv.gz")
frags_obj <- CreateFragmentObject(path = frags)

# Build peak × cell matrix
peak_matrix <- FeatureMatrix(
  fragments = frags_obj,
  features = peaks,
  cells = colnames(seurat_obj),
  sep = c(":", "-")
)

# Create ATAC assay
chrom_assay <- CreateChromatinAssay(
  counts = peak_matrix,
  sep = c(":", "-"),
  fragments = frags_obj
)

# Attach to Seurat object
seurat_obj[["ATAC"]] <- chrom_assay

# Add annotations using direct slot assignment
annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
seurat_obj[["ATAC"]]@annotation <- annotations

# Add seqinfo to fix subsetting issues
genome_info <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
atac_granges <- granges(seurat_obj[["ATAC"]])
seqinfo(atac_granges) <- genome_info[seqlevels(atac_granges)]
seurat_obj[["ATAC"]]@ranges <- atac_granges

print("Initial object:")
print(seurat_obj)

# ============================================================================
# 2) RNA QC AND FILTERING
# ============================================================================
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics before filtering
p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)
print(p1)

# WORKAROUND: Remove ATAC assay temporarily to allow subsetting
atac_assay_backup <- seurat_obj[["ATAC"]]
seurat_obj[["ATAC"]] <- NULL

# Filter RNA cells
cells_to_keep <- colnames(seurat_obj)[
  seurat_obj$nFeature_RNA > 200 & 
    seurat_obj$nFeature_RNA < 5000 & 
    seurat_obj$percent.mt < 10
]

seurat_obj <- seurat_obj[, cells_to_keep]
print(paste("Cells after RNA QC:", ncol(seurat_obj)))

# Manual subsetting of ATAC assay (bypasses subset() function)
atac_counts <- GetAssayData(atac_assay_backup, slot = "counts")
atac_counts_filtered <- atac_counts[, cells_to_keep]

# Get and subset fragments
frags_filtered <- atac_assay_backup@fragments
for (i in seq_along(frags_filtered)) {
  frags_filtered[[i]] <- subset(frags_filtered[[i]], cells = cells_to_keep)
}

# Recreate ChromatinAssay with filtered data
atac_assay_filtered <- CreateChromatinAssay(
  counts = atac_counts_filtered,
  sep = c(":", "-"),
  fragments = frags_filtered
)

# Re-add annotations
atac_assay_filtered@annotation <- annotations

# Re-add seqinfo
atac_granges_filtered <- granges(atac_assay_filtered)
seqinfo(atac_granges_filtered) <- genome_info[seqlevels(atac_granges_filtered)]
atac_assay_filtered@ranges <- atac_granges_filtered

# Add back to Seurat object
seurat_obj[["ATAC"]] <- atac_assay_filtered

print("After RNA filtering:")
print(seurat_obj)


# ============================================================================
# 3) ATAC QC AND FILTERING
# ============================================================================
DefaultAssay(seurat_obj) <- "ATAC"

# Compute QC metrics
seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj, fast = FALSE)

# Visualize ATAC QC
p2 <- VlnPlot(seurat_obj, 
              features = c("nCount_ATAC", "nFeature_ATAC", 
                           "nucleosome_signal", "TSS.enrichment"),
              ncol = 4, pt.size = 0.1)
print(p2)

# WORKAROUND: Remove ATAC assay temporarily again
DefaultAssay(seurat_obj) <- "RNA"
atac_assay_backup2 <- seurat_obj[["ATAC"]]
seurat_obj[["ATAC"]] <- NULL

# Filter ATAC cells
cells_to_keep <- colnames(seurat_obj)[
  seurat_obj$nCount_ATAC > 1000 &
    seurat_obj$nCount_ATAC < 50000 &
    seurat_obj$nFeature_ATAC > 500 &
    seurat_obj$nucleosome_signal < 2 &
    seurat_obj$TSS.enrichment > 1
]

seurat_obj <- seurat_obj[, cells_to_keep]
print(paste("Cells after ATAC QC:", ncol(seurat_obj)))

# Manual subsetting of ATAC assay (bypasses subset() function)
atac_counts2 <- GetAssayData(atac_assay_backup2, slot = "counts")
atac_counts_filtered2 <- atac_counts2[, cells_to_keep]

# Get and subset fragments
frags_filtered2 <- atac_assay_backup2@fragments
for (i in seq_along(frags_filtered2)) {
  frags_filtered2[[i]] <- subset(frags_filtered2[[i]], cells = cells_to_keep)
}

# Recreate ChromatinAssay with filtered data
atac_assay_filtered2 <- CreateChromatinAssay(
  counts = atac_counts_filtered2,
  sep = c(":", "-"),
  fragments = frags_filtered2
)

# Re-add annotations
atac_assay_filtered2@annotation <- annotations

# Re-add seqinfo
atac_granges_filtered2 <- granges(atac_assay_filtered2)
seqinfo(atac_granges_filtered2) <- genome_info[seqlevels(atac_granges_filtered2)]
atac_assay_filtered2@ranges <- atac_granges_filtered2

# Add back to Seurat object
seurat_obj[["ATAC"]] <- atac_assay_filtered2

print("After all QC:")
print(seurat_obj)
print(paste("Cells remaining after all QC:", ncol(seurat_obj)))


# ============================================================================
# 4) RNA PROCESSING
# ============================================================================
DefaultAssay(seurat_obj) <- "RNA"

# Normalize and find variable features
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj)

# PCA and clustering
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
p_elbow <- ElbowPlot(seurat_obj, ndims = 50)
print(p_elbow)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "pca")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, algorithm = 3)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction.name = "umap.rna")

# Visualize RNA clustering
p3 <- DimPlot(seurat_obj, reduction = "umap.rna", label = TRUE) + 
  ggtitle("RNA-based Clustering")
print(p3)

# ============================================================================
# 5) ATAC PROCESSING
# ============================================================================
DefaultAssay(seurat_obj) <- "ATAC"

# Normalization
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 10)
seurat_obj <- RunSVD(seurat_obj)

# Check LSI depth correlation
p_depth <- DepthCor(seurat_obj)
print(p_depth)

# UMAP with ATAC (exclude first component if correlated with depth)
seurat_obj <- RunUMAP(seurat_obj, 
                      reduction = "lsi", 
                      dims = 2:30,
                      reduction.name = "umap.atac")

# Visualize ATAC clustering
p4 <- DimPlot(seurat_obj, reduction = "umap.atac", label = TRUE) + 
  ggtitle("ATAC-based Clustering")
print(p4)


# ============================================================================
# 6) WEIGHTED NEAREST NEIGHBOR (WNN) INTEGRATION
# ============================================================================
seurat_obj <- FindMultiModalNeighbors(
  seurat_obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight"
)

# Build WNN graph
seurat_obj <- RunUMAP(seurat_obj, 
                      nn.name = "weighted.nn", 
                      reduction.name = "wnn.umap",
                      reduction.key = "wnnUMAP_")

seurat_obj <- FindClusters(seurat_obj, 
                           graph.name = "wsnn",
                           algorithm = 3,
                           resolution = 0.5,
                           verbose = FALSE)

# Visualize integrated clustering
p5 <- DimPlot(seurat_obj, reduction = "wnn.umap", label = TRUE, label.size = 6) + 
  ggtitle("WNN-based Clustering")
print(p5)

# Compare modality weights
p6 <- VlnPlot(seurat_obj, features = "RNA.weight", group.by = "seurat_clusters", 
              pt.size = 0.1) + ggtitle("RNA Weight by Cluster")
print(p6)



# ============================================================================
# 7) VISUALIZATION AND EXPLORATION
# ============================================================================

# Plot all three UMAPs side by side
p_rna <- DimPlot(seurat_obj, reduction = "umap.rna", label = TRUE) + 
  ggtitle("RNA") + NoLegend()
p_atac <- DimPlot(seurat_obj, reduction = "umap.atac", label = TRUE) + 
  ggtitle("ATAC") + NoLegend()
p_wnn <- DimPlot(seurat_obj, reduction = "wnn.umap", label = TRUE) + 
  ggtitle("WNN Integration")

p_combined <- p_rna + p_atac + p_wnn
print(p_combined)

# Save this comparison plot
ggsave("umap_comparison.png", p_combined, width = 15, height = 5, dpi = 300)


# ============================================================================
# 8) FIND MARKERS
# ============================================================================

# Find marker genes (RNA)
DefaultAssay(seurat_obj) <- "RNA"
print("Finding RNA markers...")
rna_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
top_markers <- rna_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

print("Top 10 RNA markers per cluster:")
print(top_markers)

# Find differential accessible peaks (ATAC)
DefaultAssay(seurat_obj) <- "ATAC"
print("Finding ATAC markers...")
atac_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC',
  min.pct = 0.05
)

top_peaks <- atac_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

print("Top 10 ATAC peaks per cluster:")
print(head(top_peaks, 50))


# ============================================================================
# 9) VISUALIZE MARKERS
# ============================================================================

# Common immune cell markers to check
marker_genes <- c("CD3D", "CD8A", "CD4", "CD14", "CD79A", "MS4A1", "FCGR3A", "NKG7")
existing_markers <- marker_genes[marker_genes %in% rownames(seurat_obj)]

if (length(existing_markers) > 0) {
  print(paste("Found markers:", paste(existing_markers, collapse = ", ")))
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Feature plots
  p_features <- FeaturePlot(seurat_obj, 
                            features = existing_markers[1:min(6, length(existing_markers))], 
                            reduction = "wnn.umap", 
                            ncol = 3)
  print(p_features)
  tryCatch({
    ggsave("marker_feature_plots.png", p_features, width = 15, height = 10, dpi = 300)
  }, error = function(e) {
    # Use png device directly
    png("marker_feature_plots.png", width = 15, height = 10, units = "in", res = 300)
    print(p_features)
    dev.off()
  })
  
  # Dot plot
  p_dot <- DotPlot(seurat_obj, features = existing_markers) + 
    RotatedAxis()
  print(p_dot)
  ggsave("marker_dotplot.png", p_dot, width = 10, height = 6, dpi = 300)
  
  # Coverage plot for first available marker
  if (length(existing_markers) > 0) {
    DefaultAssay(seurat_obj) <- "ATAC"
    tryCatch({
      p_cov <- CoveragePlot(
        object = seurat_obj,
        region = existing_markers[1],
        extend.upstream = 5000,
        extend.downstream = 5000
      )
      print(p_cov)
      ggsave(paste0("coverage_", existing_markers[1], ".png"), 
             p_cov, width = 12, height = 8, dpi = 300)
    }, error = function(e) {
      message("Coverage plot failed - this is normal if gene has no nearby peaks")
    })
  }
}

# Heatmap of top markers
#Re-scale data to include all genes
DefaultAssay(seurat_obj) <- "RNA"
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

#try the heatmap
top5_markers <- rna_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

p_heatmap <- DoHeatmap(seurat_obj, features = top5_markers$gene) + 
  NoLegend()
print(p_heatmap)
ggsave("marker_heatmap.png", p_heatmap, width = 12, height = 15, dpi = 300)


# ============================================================================
# 10) SAVE RESULTS
# ============================================================================

# Save processed Seurat object
saveRDS(seurat_obj, file = "multiomics_seurat_processed.rds")

# Save marker tables
write.csv(rna_markers, "rna_cluster_markers.csv", row.names = FALSE)
write.csv(atac_markers, "atac_cluster_markers.csv", row.names = FALSE)
write.csv(top_markers, "top10_rna_markers.csv", row.names = FALSE)

# Create summary report
sink("analysis_summary.txt")
cat("==============================================\n")
cat("MULTIOMICS ANALYSIS SUMMARY\n")
cat("==============================================\n\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Final Cell Count:", ncol(seurat_obj), "\n")
cat("Number of Clusters:", length(unique(Idents(seurat_obj))), "\n\n")
cat("RNA Features:", nrow(seurat_obj[["RNA"]]), "\n")
cat("ATAC Features:", nrow(seurat_obj[["ATAC"]]), "\n\n")
cat("Clusters:\n")
print(table(Idents(seurat_obj)))
cat("\n")
cat("==============================================\n")
sink()

# Save session info
writeLines(capture.output(sessionInfo()), "session_info.txt")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n==============================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("==============================================\n\n")
cat("Files saved:\n")
cat("  - multiomics_seurat_processed.rds\n")
cat("  - rna_cluster_markers.csv\n")
cat("  - atac_cluster_markers.csv\n")
cat("  - top10_rna_markers.csv\n")
cat("  - umap_comparison.png\n")
cat("  - marker_feature_plots.png\n")
cat("  - marker_dotplot.png\n")
cat("  - marker_heatmap.png\n")
cat("  - analysis_summary.txt\n")
cat("  - session_info.txt\n\n")

print(seurat_obj)
cat("\nFinal cell count:", ncol(seurat_obj), "\n")
cat("Number of clusters:", length(unique(Idents(seurat_obj))), "\n")
cat("\nTo reload this object later: seurat_obj <- readRDS('multiomics_seurat_processed.rds')\n")

