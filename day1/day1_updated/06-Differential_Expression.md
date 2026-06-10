---
output:
  html_document:
    keep_md: yes
---






# 6. Differential Expression Analysis

The main advantage of using scRNA-seq technologies is the possibility of 
assessing cell type specificity and heterogeneity, which is not possible using bulk assays. 

We expect that some of the identified clusters in the UMAP might correspond
to distinct cell types or cell states (for example different metabolic states, etc...). The assignation of cell types identities is not always
straightforward and is one of the biggest challenge of single-cell analysis. 

Cell type profiling is generally done by assessing the expression of *marker genes*. 
This task can be done manually by inspecting known markers in dimensional reduced data
projected in UMAP or tSNE (for example classifical cell type markers such as CD3 for T-cells). It can also be done in a automatic manner scoring 
cells using gene signatures, which are lists of marker genes. Scores are generally
based in the median expression of all the markers. Using scores has the advantage
or reducing bias due to the arbitrary selection of markers.

Another approach is to determine genes which are differentially expressed between clusters, and then trying to understand what these genes are. Finding differential expressed markers (or differentially expressed genes, DEGs) is important for cluster profiling and
identification. We will use the `FindAllMarkers()` function, which performs
a statistical test comparing the distribution of gene expression values for 
each gene separately comparing one assigned cell type cluster (in this case 
the `seurat_clusters` column) *vs* the rest of cells. 

First, we will set up the column used to define the clusters using the 
function `Idents()`. 


``` r
Idents(ad_medula_filtered) <- ad_medula_filtered$seurat_clusters
```

Now we can calculate the DEGs. 
There are several parameters for `FindAllMarkers()`, we will discuss
`logfc.threshold`, `min.pct` and `min.cells.feature` that corresponds to the threshold of gene expression fold change, the minimum percentage of cells expressing the marker and the minimum of cells expressing (counts > 0) the feature. These parameters 
are used to filter out genes prior calculating DEGs. Lowering the values of these parameters will increase the sensibility of the method at the expense of increasing computation time.



``` r
ad_medula.degs <- FindAllMarkers(ad_medula_filtered, 
                            logfc.threshold = 1, 
                            min.pct = 0.05, 
                            min.cells.feature = 10, 
                            verbose = FALSE)
```


The output `ad_medula.degs` consist of a data frame contanning the DEGs with
p-vales, p-adjusted values and log fold change values for each gene as 
we can see next:



``` r
head(ad_medula.degs)
```

```
##                  p_val avg_log2FC pct.1 pct.2     p_val_adj cluster     gene
## SLC24A2  1.103268e-232   4.983958 0.631 0.031 2.765231e-228       0  SLC24A2
## CCSER1   3.469903e-223   1.702038 0.990 0.488 8.696964e-219       0   CCSER1
## AGBL4    1.053599e-189   1.610500 0.963 0.476 2.640739e-185       0    AGBL4
## TH       1.021532e-184   2.278201 0.854 0.281 2.560368e-180       0       TH
## ST18     7.427605e-146   2.297766 0.739 0.212 1.861655e-141       0     ST18
## KIAA1244 8.635505e-144   2.750580 0.646 0.143 2.164403e-139       0 KIAA1244
```


Differential gene expression are usually represented graphically as **volcano plots**, in which each gene is represented using the log-fold change between one cluster and the rest (x-axis) and the statistical significance of this logFC (-log10(Pvalue) in y-axis).

We can make a volcano plot using `ggplot`:


``` r
library(ggplot2)
library(dplyr)         ## for handling data frames
library(ggrepel)

ad_medula.degs %>%
  filter(cluster=='0') %>%
  filter(p_val_adj > 0) %>%
  arrange(desc(abs(avg_log2FC))) %>%       ## Arranging genes by FC
  mutate(highlight=ifelse(-log10(p_val_adj)>280, TRUE, FALSE)) %>% ## highlighting top FC markers
  mutate(gene_label=ifelse(highlight==TRUE, gene, '')) %>% ## Adding labels for top markers
  ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj),
             colour=highlight,
             label=gene_label)) +         ## adding labels for top markers
      geom_point() +
      geom_text_repel() +
      theme_bw()
```

![](06-Differential_Expression_files/figure-html/vulcano_plot-1.png)<!-- -->


## Exercises

> What can you learn from this plot, if you look at how the genes are distributed? Can you make the equivalent volcano plot for the other clusters?


[Previous Chapter (Cluster Visualization)](./05-Cluster_visualization.md)|
[Next Chapter (Profiling cells)](./07-Profiling_cells.md)
