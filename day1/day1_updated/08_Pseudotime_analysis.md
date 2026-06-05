---
output:
  html_document:
    keep_md: yes
---


# Introduction

This analysis reconstructs developmental trajectories in the adrenal medulla
using the Jansky et al. (2021) single-cell RNA-seq dataset.

Trajectory inference is performed using Slingshot. Existing UMAP coordinates
provided by the original study are used directly. No cell filtering or
subsetting is performed.

Slingshot is a trajectory inference method that reconstructs developmental paths from single-cell RNA sequencing data. It starts by representing cells in a low-dimensional space, such as UMAP or PCA, where cells with similar gene expression profiles are located close together. Slingshot then connects groups of cells (clusters or annotated cell types) using a minimum spanning tree, which provides a framework for identifying possible developmental relationships. Next, it fits smooth curves through the cells to represent one or more biological trajectories. Finally, each cell is assigned a pseudotime value, which reflects its relative position along a trajectory. Pseudotime does not correspond to real chronological time, but rather to a cell's progression through a biological process such as differentiation, allowing researchers to study how gene expression changes as cells transition between states.

The starting population is defined as Schwann Cell Precursors (SCPs).

# Load libraries


``` r
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(patchwork)

set.seed(1234)
```

# Load data

First, we load the already pre-processed data from Jansky et al 2021.



``` r
seurat_obj <- readRDS("/Users/cbg-mbp-02/Documents/data/janski_2021/janski_2021_seurat_subset_2400_cells.rds")

seurat_obj
```

```
## An object of class Seurat 
## 28422 features across 2400 samples within 1 assay 
## Active assay: RNA (28422 features, 0 variable features)
##  1 layer present: counts
```

# Inspect metadata

We will implement a pseudotime analysis along the already annotated cell types. The calculation
is guided on cluster or cell type information. 



``` r
str(seurat_obj@meta.data)
```

```
## 'data.frame':	2400 obs. of  13 variables:
##  $ cell_id     : chr  "AAACGCTGTCAAAGTA_1" "AAAGGATCAGATCACT_1" "AAAGGATCAGCAGAAC_1" "AACAACCAGTCCTACA_1" ...
##  $ orig.ident  : int  14607 14607 14607 14607 14607 14607 14607 14607 14607 14607 ...
##  $ nCount_RNA  : num  1984 1645 1715 2882 2029 ...
##  $ nFeature_RNA: int  1347 1226 1254 1871 1391 2811 1899 2200 1192 2690 ...
##  $ percent.mt  : num  0.2016 0.0608 0 0.1041 0.0493 ...
##  $ S.Score     : num  0.0879 -0.0594 -0.0209 -0.1035 0.0835 ...
##  $ Phase       : chr  "S" "G1" "G1" "G2M" ...
##  $ cytotrace   : num  0.353 0.203 0.314 0.601 0.317 ...
##  $ pseudotime  : num  17.8 15.6 14 12.9 13.2 ...
##  $ celltype    : chr  "Neuroblasts" "Neuroblasts" "Connecting progenitor cells" "Connecting progenitor cells" ...
##  $ X           : chr  "AAACGCTGTCAAAGTA_1" "AAAGGATCAGATCACT_1" "AAAGGATCAGCAGAAC_1" "AACAACCAGTCCTACA_1" ...
##  $ UMAP_1      : num  5.57 4.8 4.04 4.48 4.52 ...
##  $ UMAP_2      : num  -3.477 -2.193 -0.145 0.913 0.577 ...
```


# Data normalization

Slingshot does not directly require normalized expression values when
operating on existing embeddings, but normalization is useful for
downstream analyses.



``` r
seurat_obj <- NormalizeData(seurat_obj)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

seurat_obj <- ScaleData(
  seurat_obj,
  features = VariableFeatures(seurat_obj)
)

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj), 
  verbose = FALSE
)
```

# Create SingleCellExperiment object

Slingshot uses a different object also used for scRNA-Seq analysis 
from the SingleCellExperiment library.


``` r
sce <- as.SingleCellExperiment(seurat_obj)

sce
```

```
## class: SingleCellExperiment 
## dim: 28422 2400 
## metadata(0):
## assays(2): counts logcounts
## rownames(28422): RP11-34P13.7 AL627309.1 ... RLN3 AC067969.1
## rowData names(0):
## colnames(2400): AAACGCTGTCAAAGTA_1 AAAGGATCAGATCACT_1 ...
##   TTCCGGTTCGGTAGGA_17 TTGTTGTCAGAAATCA_17
## colData names(14): cell_id orig.ident ... UMAP_2 ident
## reducedDimNames(1): PCA
## mainExpName: RNA
## altExpNames(0):
```

As in seurat, the SingleCellExperiment object contains the gene expression
data along with metadata. We will use this object to calculate the trajectories
based in already calculated UMAP coordinates.

# Add existing UMAP coordinates

The dataset already contains UMAP coordinates in the metadata.


``` r
umap_coords <- as.matrix(
  seurat_obj@meta.data[, c("UMAP_1", "UMAP_2")]
)

head(umap_coords)
```

```
##                       UMAP_1     UMAP_2
## AAACGCTGTCAAAGTA_1  5.574948 -3.4768778
## AAAGGATCAGATCACT_1  4.800121 -2.1926302
## AAAGGATCAGCAGAAC_1  4.037976 -0.1445995
## AACAACCAGTCCTACA_1  4.479241  0.9133858
## AACAAGAAGGACTTCT_1  4.515640  0.5770715
## AACAAGATCTGCGTCT_1 -8.895806  1.5624598
```

``` r
reducedDims(sce)$UMAP <- umap_coords
```

# Visualize annotated cell types



``` r
ggplot(
  seurat_obj@meta.data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = celltype
  )
) +
  geom_point(size = 0.5) +
  theme_classic()
```

![](08_Pseudotime_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

# Run Slingshot

The slingshot() function performs trajectory inference using the cell annotations and low-dimensional representation of the data. In this example, Slingshot uses the celltype labels to define groups of biologically similar cells and the existing UMAP coordinates (reducedDim = "UMAP") to determine how these groups are arranged in the developmental landscape. The argument start.clus = "SCPs" specifies Schwann Cell Precursors (SCPs) as the starting population, allowing Slingshot to infer trajectories that originate from this early cell state and progress toward more differentiated cell types. The resulting object contains the inferred lineages as well as pseudotime values for each cell.


``` r
sce <- slingshot(
  sce,
  clusterLabels = "celltype",  # use annotated cell types as groups
  reducedDim = "UMAP",         # use existing UMAP coordinates
  start.clus = "SCPs"          # define SCPs as the starting population
)
```

# Visualize inferred trajectories


``` r
celltype_factor <- as.factor(colData(sce)$celltype)

plot(
  reducedDims(sce)$UMAP,
  col = brewer.pal(
    min(length(levels(celltype_factor)), 12),
    "Set3"
  )[as.numeric(celltype_factor)],
  pch = 16,
  asp = 1,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  main = "Slingshot trajectory",
  xlim = c(-15, 15),
  ylim = c(-15, 15)
)

lines(
  SlingshotDataSet(sce),
  lwd = 3,
  col = "black",
  xlim = c(-15, 15),
)
```

![](08_Pseudotime_analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

## Inspecting and subsetting trajectories


``` r
slingLineages(sce)
```

```
## $Lineage1
## [1] "SCPs"                        "Bridge"                     
## [3] "Connecting progenitor cells" "Chromaffin cells"           
## [5] "Late chromaffin cells"      
## 
## $Lineage2
## [1] "SCPs"                        "Bridge"                     
## [3] "Connecting progenitor cells" "Neuroblasts"                
## [5] "Late neuroblasts"           
## 
## $Lineage3
## [1] "SCPs"      "Late SCPs"
```

## Extracting pseudotime 

seudotime <- slingPseudotime(sce) extracts the inferred developmental ordering of each cell from a Slingshot trajectory, assigning a numeric value that reflects where each cell sits along a lineage path from early to late states.


``` r
pseudotime_df <- slingPseudotime(sce) %>% as.data.frame()

seurat_obj$pt_lineage1 <- pseudotime_df$Lineage1
```

The NA points in pt_lineage1 are cells for which Slingshot (or whatever generated pseudotime) did not assign a pseudotime value for that specific lineage.

In practice, that usually means one of two things:

the cell is not part of the inferred trajectory/lineage 2 (it belongs to another branch), or
the algorithm couldn’t place it along that lineage because it lies outside the fitted path or isn’t connected to it.



``` r
seurat_obj@meta.data %>%
    ggplot(aes(UMAP_1, UMAP_2, colour = pt_lineage1)) +
      geom_point() +
      viridis::scale_color_viridis() +
      theme_classic()
```

![](08_Pseudotime_analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

We can visualize the smooth fitted trajectory as follows.


``` r
plot(
  reducedDims(sce)$UMAP,
  col = "grey80",
  pch = 16
)

lines(
    slingCurves(sce)[["Lineage1"]],
    lwd = 4,
    col = "red"
)
```

![](08_Pseudotime_analysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


# Saving results



``` r
pseudo_time_analysis <- seurat_obj@meta.data 

write.csv(
  pseudo_time_analysis,
  file = "slingshot_pseudotime.csv",
  row.names = FALSE
)
```


[Previous Chapter (Profiling cells)](./07-Profiling_cells.md)|
