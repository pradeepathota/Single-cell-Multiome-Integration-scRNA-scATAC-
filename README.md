# Single-Cell Multiomics Integration: RNA-seq and ATAC-seq Analysis

[![R Version](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.0-green.svg)](https://satijalab.org/seurat/)
[![Signac](https://img.shields.io/badge/Signac-1.16-orange.svg)](https://stuartlab.org/signac/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Integrated analysis pipeline for single-cell multiomics data, combining RNA-seq (gene expression) and ATAC-seq (chromatin accessibility) to identify cell types and their regulatory mechanisms. This project implements Weighted Nearest Neighbor (WNN) integration for superior cell type identification compared to single-modality approaches.

## Key Features

- **Dual-modality analysis**: Simultaneous RNA and ATAC-seq processing from the same cells
- **Comprehensive QC**: Multi-metric filtering for both modalities (TSS enrichment, nucleosome signal, mitochondrial percentage)
- **WNN Integration**: State-of-the-art algorithm that learns optimal weights for each modality per cell
- **Statistical rigor**: Differential expression and accessibility analysis with FDR correction
- **Publication-quality visualizations**: UMAPs, heatmaps, dot plots, and coverage plots

## Results Summary

- **Starting cells**: 11,909
- **High-quality cells after QC**: 4,915 (41% pass rate)
- **RNA features**: 36,601 genes
- **ATAC features**: 108,377 peaks
- **Cell clusters identified**: Multiple distinct populations
- **Analysis approach**: Integrated multi-modal clustering with marker discovery

## Technologies

| Tool | Version | Purpose |
|------|---------|---------|
| R | 4.5.2 | Statistical computing |
| Seurat | 5.0+ | Single-cell RNA-seq analysis |
| Signac | 1.16+ | Chromatin accessibility analysis |
| Bioconductor | Latest | Genomic data structures |
| ggplot2 | 3.4+ | Data visualization |

## Workflow
```
Raw 10x Multiome Data
        ↓
   Load & Setup
        ↓
RNA QC ← → ATAC QC
   ↓           ↓
  PCA         LSI
   ↓           ↓
   └─────┬─────┘
         ↓
    WNN Integration
         ↓
  UMAP & Clustering
         ↓
  Marker Identification
         ↓
   Visualization
```

## Quick Start

### Prerequisites
```r
# Install required packages
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
# Clone repository
git clone https://github.com/yourusername/single-cell-multiomics-integration.git
cd single-cell-multiomics-integration

# Place your data in data/ directory
# See data/README.md for required files

# Run analysis
setwd("path/to/single-cell-multiomics-integration")
source("scripts/complete_analysis.R")
```

## Data Requirements

This pipeline requires **10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression** data:

**Required files:**
- `filtered_feature_bc_matrix.h5` - Combined RNA and ATAC counts
- `atac_peaks.bed` - Called chromatin accessibility peaks
- `atac_fragments.tsv.gz` - Raw ATAC fragments
- `atac_fragments.tsv.gz.tbi` - Fragment index file

**Download example dataset:**
```bash
# PBMC dataset from 10x Genomics (free, no account needed)
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
# ... (see data/README.md for complete download commands)
```

Place files in `data/` directory before running analysis.

## Analysis Pipeline

### 1. Quality Control

**RNA QC Metrics:**
- Gene detection rate (>200 genes per cell)
- Total RNA counts
- Mitochondrial percentage (<10%)

**ATAC QC Metrics:**
- Fragment counts (1,000-50,000 per cell)
- TSS enrichment (>1, measures signal at promoters)
- Nucleosome signal (<2, measures data quality)

**Result:** 59% of cells filtered, ensuring high-quality dataset.

### 2. Data Processing

**RNA Processing:**
- Normalization (accounts for sequencing depth)
- Feature selection (3,000 most variable genes)
- PCA dimensionality reduction (50 → 30 components)

**ATAC Processing:**
- TF-IDF normalization (standard for accessibility data)
- Feature selection (top accessible peaks)
- LSI dimensionality reduction (excluding component 1)

### 3. Multi-Modal Integration

**Weighted Nearest Neighbor (WNN) Algorithm:**
- Calculates cell similarities in both RNA and ATAC space
- Learns optimal per-cell weights for each modality
- Creates integrated neighbor graph
- Generates unified UMAP embedding

**Why WNN?**
- Some cells are better defined by gene expression
- Others are better defined by chromatin state
- WNN automatically determines optimal balance for each cell

### 4. Cell Type Identification

**Methods:**
- Graph-based Louvain clustering
- UMAP visualization (2D projection)
- Differential expression analysis (Wilcoxon test)
- Differential accessibility analysis (Logistic regression)
- Statistical validation (FDR < 0.05)

### 5. Marker Discovery

**RNA Markers:** Genes specifically expressed in each cluster  
**ATAC Markers:** Regulatory elements specifically accessible in each cluster

## Key Results

### Quality Control

![RNA QC](figures/rna_qc_metrics.png)
*RNA quality metrics showing gene counts, RNA counts, and mitochondrial percentage*

![ATAC QC](figures/atac_qc_metrics.png)  
*ATAC quality metrics showing fragment counts, TSS enrichment, and nucleosome signal*

### Multi-Modal Comparison

![UMAP Comparison](figures/umap_comparison.png)
*Comparison of RNA-only, ATAC-only, and WNN integration clustering*

**Key Finding:** WNN integration reveals clearer cell type separation by optimally combining both data modalities.

### Cell Type Markers

![Marker Heatmap](figures/marker_heatmap.png)
*Top 5 marker genes per cluster showing distinct expression patterns*

![Marker Dotplot](figures/marker_dotplot.png)
*Expression of key cell type markers across all clusters*

## Technical Highlights

### Challenges Overcome

1. **Genomic Coordinate Compatibility**
   - **Issue:** Signac ChromatinAssay validation errors with seqinfo
   - **Solution:** Direct slot assignment and manual seqinfo addition

2. **ChromatinAssay Subsetting**
   - **Issue:** Standard subsetting triggered validation failures
   - **Solution:** Custom workflow reconstructing assay after filtering

3. **Memory Management**
   - **Issue:** Large sparse matrices (>500M data points)
   - **Solution:** Efficient operations on sparse matrix formats

### Performance

- **Runtime:** ~30-45 minutes (standard laptop, 8GB RAM)
- **Peak memory:** ~6GB
- **Scalability:** Tested up to 10,000 cells
- **Efficiency:** Sparse matrix operations preserve memory

## Repository Structure
```
single-cell-multiomics-integration/
├── README.md                          # This file
├── scripts/
│   ├── 01_load_data.R                # Data import
│   ├── 02_qc_filtering.R             # Quality control
│   ├── 03_rna_processing.R           # RNA analysis
│   ├── 04_atac_processing.R          # ATAC analysis
│   ├── 05_wnn_integration.R          # Multi-modal integration
│   ├── 06_marker_identification.R    # Find markers
│   └── 07_visualization.R            # Generate figures
├── results/
│   ├── rna_cluster_markers.csv       # Marker genes
│   ├── atac_cluster_markers.csv      # Differential peaks
│   └── analysis_summary.txt          # Summary statistics
├── figures/                           # All visualizations
└── data/                              # Input data (not tracked)
```

## Output Files

### Results Tables

- `rna_cluster_markers.csv` - Differentially expressed genes per cluster
- `atac_cluster_markers.csv` - Differentially accessible peaks per cluster  
- `top10_rna_markers.csv` - Top 10 marker genes per cluster
- `analysis_summary.txt` - Summary statistics and metadata

### Visualizations

- `umap_comparison.png` - Three-way UMAP comparison (RNA, ATAC, WNN)
- `marker_heatmap.png` - Heatmap of top markers per cluster
- `marker_dotplot.png` - Dot plot of canonical cell type markers
- `marker_feature_plots.png` - Expression patterns on UMAP
- Quality control violin plots

### Processed Objects

- `multiomics_seurat_processed.rds` - Final Seurat object with all analyses

## Interpretation Guide

### Understanding UMAPs

- **Each dot** = one cell
- **Colors** = cluster assignments
- **Proximity** = similarity (close cells are similar)
- **Separation** = distinct cell types/states

### Reading Marker Tables

| Column | Meaning |
|--------|---------|
| `gene` | Gene name |
| `avg_log2FC` | Log2 fold-change (>1 = 2x higher) |
| `pct.1` | % cells expressing in cluster |
| `pct.2` | % cells expressing in other clusters |
| `p_val_adj` | Adjusted p-value (significance) |

### Cluster Annotation

Use marker genes to identify cell types:
- Compare markers to databases (CellMarker, PanglaoDB)
- Look for canonical markers (CD3D = T cells, CD79A = B cells)
- Consider expression patterns and accessibility

## Citation

If you use this pipeline, please cite:
```
[Your Name] (2025). Single-Cell Multiomics Integration Pipeline.
GitHub: https://github.com/yourusername/single-cell-multiomics-integration
```

**Key references:**
- Hao et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*
- Stuart et al. (2021). Single-cell chromatin state analysis with Signac. *Nature Methods*

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add new feature'`)
4. Push to branch (`git push origin feature/improvement`)
5. Create Pull Request

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contact

**[Your Name]**  
📧 Email: your.email@example.com  
💼 LinkedIn: [linkedin.com/in/yourprofile](https://linkedin.com/in/yourprofile)  
🌐 Portfolio: [yourwebsite.com](https://yourwebsite.com)

## Acknowledgments

- 10x Genomics for multiome technology and datasets
- Satija Lab (NYU) for Seurat and WNN algorithm
- Stuart Lab (Stanford) for Signac package
- Bioconductor community for genomic tools

---

⭐ **If you find this project useful, please star the repository!**

📖 **For detailed methods, see [docs/methods.md](docs/methods.md)**

🔬 **For biological interpretation, see [docs/interpretation.md](docs/interpretation.md)**
