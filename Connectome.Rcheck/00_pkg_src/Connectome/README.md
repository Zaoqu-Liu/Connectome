# Connectome

## Single-Cell Connectomics Analysis Toolkit

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/Zaoqu-Liu/Connectome)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Connectome is an R toolkit for analyzing cell-cell connectivity patterns based on ligand-receptor interactions in single-cell RNA sequencing data.

**Maintainer:** [Zaoqu Liu](https://github.com/Zaoqu-Liu)

**Original Author:** Micha Sam Brickman Raredon (Yale University)

---

## Citation

If you use Connectome in your research, please cite:

Raredon, M.S.B., Yang, J., Garritano, J. et al. Computation and visualization of cell–cell signaling topologies in single-cell systems data using Connectome. *Sci Rep* **12**, 4187 (2022). https://doi.org/10.1038/s41598-022-07959-x

---

## Installation

### Prerequisites

```r
install.packages(c("Seurat", "dplyr", "tibble", "data.table", "Matrix", "igraph"))

# For Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("circlize", "ComplexHeatmap"))
```

### Install from GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('Zaoqu-Liu/Connectome')
```

---

## Quick Start

```r
library(Seurat)
library(Connectome)

# Prepare Seurat object
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Create connectome
connectome <- CreateConnectome(
  object = seurat_obj,
  species = 'human',
  LR.database = 'fantom5',
  min.cells.per.ident = 75
)

# Filter and visualize
connectome_filtered <- FilterConnectome(
  connectome,
  min.pct = 0.1,
  min.z = 0.25
)

# Network visualization
NetworkPlot(connectome_filtered)

# Circos plot
CircosPlot(connectome_filtered)
```

---

## Main Functions

### Core Functions

| Function | Description |
|----------|-------------|
| `CreateConnectome()` | Generate connectomic edgelist from Seurat object |
| `FilterConnectome()` | Filter edges by expression, p-value, or other criteria |
| `DifferentialConnectome()` | Compare two connectomes for differential analysis |

### Visualization

| Function | Description |
|----------|-------------|
| `NetworkPlot()` | igraph-based network visualization |
| `CircosPlot()` | Circular chord diagram |
| `CircosDiff()` | Circos plot for differential connectome |
| `EdgeDotPlot()` | Dot plot of edge weights |
| `ModalDotPlot()` | Mode-organized centrality dot plot |
| `Centrality()` | Centrality analysis visualization |

### Analysis

| Function | Description |
|----------|-------------|
| `CompareCentrality()` | Compare centrality across multiple connectomes |
| `SignalScatter()` | Scatter plot of signaling vectors |
| `CellCellScatter()` | Scatter plot of cell-cell interactions |
| `SingleCellConnectome()` | Single-cell resolution connectome |

---

## Parameters Guide

### CreateConnectome Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `species` | - | Species for FANTOM5 database ('human', 'mouse', 'rat', 'pig') |
| `LR.database` | 'fantom5' | Ligand-receptor database to use |
| `min.cells.per.ident` | NULL | Minimum cells per cluster to include |
| `max.cells.per.ident` | NULL | Downsample to this number of cells per cluster |
| `p.values` | TRUE | Calculate statistical significance |
| `include.putative` | TRUE | Include putative L-R pairs |
| `calculate.DOR` | FALSE | Calculate Diagnostic Odds Ratio |

### FilterConnectome Parameters

| Parameter | Description |
|-----------|-------------|
| `min.pct` | Minimum fraction of cells expressing ligand/receptor |
| `min.z` | Minimum z-score threshold |
| `min.exp` | Minimum expression level |
| `max.p` | Maximum p-value threshold |
| `modes.include` | Signaling modes to include |
| `sources.include` | Source cell types to include |
| `targets.include` | Target cell types to include |

---

## Supported Species

The package includes FANTOM5 ligand-receptor databases for:

- Human (*Homo sapiens*)
- Mouse (*Mus musculus*)  
- Rat (*Rattus norvegicus*)
- Pig (*Sus scrofa*)

Custom ligand-receptor lists are also supported via the `custom.list` parameter.

---

## Integration with LIANA

Connectome integrates with [LIANA](https://github.com/saezlab/liana) for comprehensive cell-cell communication analysis:

```r
library(liana)

# Use Connectome method within LIANA framework
liana_res <- liana_wrap(sce, method = "connectome")
```

---

## System Requirements

- **R version:** ≥ 4.0.0
- **Dependencies:** Seurat (≥ 4.0.0), igraph, circlize, ComplexHeatmap
- **Memory:** Scales with dataset size (typically 2-8 GB RAM)
- **OS:** Cross-platform (macOS, Linux, Windows)

---

## License

GPL-3.0

---

## References

1. Raredon, M.S.B., Yang, J., Garritano, J. et al. Computation and visualization of cell–cell signaling topologies in single-cell systems data using Connectome. *Sci Rep* **12**, 4187 (2022). https://doi.org/10.1038/s41598-022-07959-x

2. Ramilowski, J.A. et al. A draft network of ligand–receptor-mediated multicellular signalling in human. *Nat Commun* **6**, 7866 (2015). https://doi.org/10.1038/ncomms8866

3. Raredon, M.S.B. et al. Single-cell connectomic analysis of adult mammalian lungs. *Science Advances* **5**(12), eaaw3851 (2019). https://doi.org/10.1126/sciadv.aaw3851
