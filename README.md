# Single-Cell Multiomics Integration: RNA-seq and ATAC-seq Analysis

[![R Version](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.0-green.svg)](https://satijalab.org/seurat/)
[![Signac](https://img.shields.io/badge/Signac-1.16-orange.svg)](https://stuartlab.org/signac/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Integrated analysis pipeline for single-cell multiomics data, combining RNA-seq (gene expression) and ATAC-seq (chromatin accessibility) to identify cell types and their regulatory mechanisms. This project implements Weighted Nearest Neighbor (WNN) integration for superior cell type identification.

![UMAP Comparison](umap_comparison.png)
*Comparison of RNA-only, ATAC-only, and WNN integration clustering*

## Key Features

- **Dual-modality analysis**: Simultaneous RNA and ATAC-seq processing
- **Comprehensive QC**: Multi-metric filtering for both modalities
- **WNN Integration**: State-of-the-art multi-modal clustering algorithm
- **Statistical rigor**: Differential expression/accessibility with FDR correction
- **Visualizations**: UMAPs, heatmaps, and dot plots

## Results Summary

- **Starting cells**: 11,909
- **High-quality cells after QC**: 4,915 (41% pass rate)
- **RNA features analyzed**: 36,601 genes
- **ATAC features analyzed**: 108,377 chromatin accessibility peaks
- **Quality control**: Rigorous filtering based on RNA and ATAC metrics

## Technologies Used

- **R 4.5.2** - Statistical computing
- **Seurat v5** - Single-cell RNA-seq analysis
- **Signac v1.16** - Chromatin accessibility analysis
- **Bioconductor** - Genomic data structures (GenomicRanges, BSgenome, EnsDb)
- **ggplot2** - Data visualization
- **dplyr** - Data manipulation

## Workflow
```
10x Multiome Data (RNA + ATAC)
          ↓
    Quality Control
    • RNA: gene counts, mitochondrial %
    • ATAC: TSS enrichment, nucleosome signal
          ↓
    Normalization & Processing
    • RNA: PCA dimensionality reduction
    • ATAC: LSI dimensionality reduction
          ↓
    WNN Integration
    • Learns optimal weights per cell
    • Combines both modalities
          ↓
    Clustering & Visualization
    • Graph-based clustering
    • UMAP projection
          ↓
    Marker Identification
    • Differential gene expression
    • Differential chromatin accessibility
```

## Installation

### Prerequisites
```r
# Install CRAN packages
install.packages(c("Seurat", "Signac", "ggplot2", "dplyr", "patchwork"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges",
  "GenomeInfoDb",
  "rtracklayer",
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38"
))
```

### Running the Analysis
```r
# Load the processed Seurat object
seurat_obj <- readRDS("multiomics_seurat_processed.rds")

# Explore the data
print(seurat_obj)
DimPlot(seurat_obj, reduction = "wnn.umap", label = TRUE)

# View marker genes
markers <- read.csv("rna_cluster_markers.csv")
head(markers)
```

## Key Results

### 1. Multi-Modal Integration

![UMAP Comparison](umap_comparison.png)

The figure shows three clustering approaches:
- **Left**: RNA-only clustering
- **Center**: ATAC-only clustering  
- **Right**: WNN integrated clustering (optimal)

**Key finding**: WNN integration reveals clearer cell type separation by combining both data modalities.

### 2. Cell Type Markers

![Marker Heatmap](marker_heatmap.png)

Heatmap showing top marker genes for each cluster, revealing distinct expression patterns.

### 3. Marker Gene Tables

**RNA Cluster Markers** (`rna_cluster_markers.csv`):
- Contains all differentially expressed genes per cluster
- Includes statistics: fold-change, p-values, expression percentages

**Top 10 Markers** (`top10_rna_markers.csv`):
- Top 10 marker genes for each cluster
- Useful for rapid cell type identification

**ATAC Cluster Markers** (`atac_cluster_markers.csv`):
- Differentially accessible chromatin peaks per cluster
- Links regulatory elements to cell types

## File Descriptions

| File | Description |
|------|-------------|
| `multiomics_seurat_processed.rds` | Complete analyzed Seurat object with RNA, ATAC, and WNN |
| `rna_cluster_markers.csv` | All differentially expressed genes per cluster |
| `atac_cluster_markers.csv` | All differentially accessible peaks per cluster |
| `top10_rna_markers.csv` | Top 10 marker genes per cluster for quick reference |
| `umap_comparison.png` | Three-way comparison of clustering approaches |
| `marker_heatmap.png` | Heatmap of top marker genes |
| `analysis_summary.txt` | Summary statistics and metadata |
| `session_info.txt` | Complete R package versions used |

## Understanding the Results

### Reading Marker Tables

Each marker gene table contains:

| Column | Meaning |
|--------|---------|
| `cluster` | Cluster ID (0, 1, 2, ...) |
| `gene` | Gene symbol |
| `avg_log2FC` | Log2 fold-change (>1 means 2x higher expression) |
| `pct.1` | % of cells expressing in this cluster |
| `pct.2` | % of cells expressing in other clusters |
| `p_val_adj` | Adjusted p-value (FDR corrected) |

**Example interpretation:**
- Gene with `avg_log2FC = 2.5` is ~5.7x higher in the cluster
- Gene with `pct.1 = 0.90, pct.2 = 0.10` is expressed in 90% of cluster cells vs 10% elsewhere
- `p_val_adj < 0.05` = statistically significant

### Cell Type Identification

Use marker genes to identify cell types:

**Common immune cell markers (if PBMC data):**
- CD3D, CD8A → T cells
- CD4 → Helper T cells
- CD14 → Monocytes
- CD79A, MS4A1 → B cells
- FCGR3A, NKG7 → NK cells

Compare your top markers to:
- [CellMarker Database](http://biocc.hrbmu.edu.cn/CellMarker/)
- [PanglaoDB](https://panglaodb.se/)
- Literature for your tissue type

## Technical Highlights

### Quality Control Metrics

**RNA QC:**
- Genes per cell: >200, <5000
- Mitochondrial percentage: <10%
- **Result**: Filtered 59% of cells

**ATAC QC:**
- Fragments per cell: 1,000-50,000
- TSS enrichment: >1 (signal at promoters)
- Nucleosome signal: <2 (data quality)

### Weighted Nearest Neighbor Integration

**Why WNN is powerful:**
- Learns which data modality (RNA or ATAC) is most informative for each cell
- Some cells are better defined by gene expression
- Others are better defined by chromatin accessibility
- WNN automatically determines optimal weights

**How it works:**
1. Calculate cell similarities in RNA space (PCA)
2. Calculate cell similarities in ATAC space (LSI)
3. Learn per-cell weights for each modality
4. Combine into unified clustering

### Challenges Overcome

1. **Genomic coordinate compatibility**: Resolved Signac ChromatinAssay validation errors
2. **Memory efficiency**: Handled 500M+ data points using sparse matrices
3. **Multi-modal integration**: Successfully combined heterogeneous data types

## Data Source

This analysis uses **10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression** technology, which simultaneously measures:
- Gene expression (RNA-seq) from ~36,000 genes
- Chromatin accessibility (ATAC-seq) from ~100,000 genomic regions
- From the **same individual cells**


**Key references:**
- Hao et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*
- Stuart et al. (2021). Single-cell chromatin state analysis with Signac. *Nature Methods*

## Contact

**venkata pradeep kumar Athota**  
📧 Email: pradeepathota3@gmail.com  

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- 10x Genomics for multiome technology
- Satija Lab for Seurat and WNN algorithm
- Stuart Lab for Signac package
- Bioconductor community

---

⭐ **If you find this useful, please star the repository!**

📊 **See `analysis_summary.txt` for complete analysis details**
