---
output:
  html_document:
    keep_md: yes
---






# 5. Cell clustering


Detection of groups or cluster of cells is an important task in scRNA-seq 
analysis. These groups could represent different cell states or cell types. 

## Defining clusters

Seurat implements a clustering method based in KNN graphs and 
community detection using the Louvain algorithm. An important parameter
for clustering is the *resolution* which can be set to increase/reduce the granularity of the clusters.

This method can be implemented by using the functions `FindNeighbors()` and
`FindClusters()` as follows:



``` r
ad_medula_filtered <- FindNeighbors(ad_medula_filtered)
ad_medula_filtered <- FindClusters(ad_medula_filtered, 
                            resolution = 0.1, 
                            verbose = FALSE)
```


The results of the clustering are stored in the Seurat metadata slot, which
can be accessed as a simple data.frame using the `$` operator. A column vector
containing each cluster for each cell is name `seurat_cluster` as shown next:




``` r
head(ad_medula_filtered$seurat_clusters)
```

```
## AAACGCTGTCAAAGTA_1 AAAGGATCAGATCACT_1 AAAGGATCAGCAGAAC_1 AACAACCAGTCCTACA_1 
##                  1                  1                  2                  2 
## AACAAGAAGGACTTCT_1 AACAAGATCTGCGTCT_1 
##                  2                  3 
## Levels: 0 1 2 3
```




There are 3 different clusters, labeled from 0 to 2 and stored like a factor.
We can plot a frequency table of the number of cells assigned to each cluster
by the algorithm.




``` r
table(ad_medula_filtered$seurat_clusters) 
```

```
## 
##   0   1   2   3 
## 624 570 557 536
```


So, 624 cells were assigned to the cluster 0.

## Exercise

> Try different different parameters for the clustering. For example, `k.param` in the *FindNeighbors()* function and higher levels of resolution. How do these 2 parameters influence the number of clusters?


## Cluster visualization


Transformations like PCA, tSNE or UMAP are used to project multidimensional
data into 2D or 3D representations that can be visualized at the expense
of the lose of information. tSNE and UMAP transformations aims to preserve
global relations between sample points. We will use UMAPs to visualize the
scRNA-Seq data.


### UMAP

We can use the `RunUMAP` function to calculate the UMAP transformation. The 
calculation of a UMAP projection can intensive computationally and is 
usually carried out on already dimensional reduced data using, for example,
PCA. The RunUMAP from seurat by default will use the PCA reduced data, the
parameter `dims` sets the number of dimensions, PCs that should be used, as
we saw we can use 7 PCs which are the ones in which there is more variability.




``` r
ad_medula_filtered <- RunUMAP(ad_medula_filtered, 
                       dims = 1:7, 
                       verbose = FALSE)
```

After the calculation of the UMAP we can visualize it using the function
`DimPlot()`.




``` r
DimPlot(ad_medula_filtered)
```

![](05-Cluster_visualization_files/figure-html/unnamed-chunk-5-1.png)<!-- -->



## Exercise


> Perform a UMAP visualization using different sets
of parameters.


[Previous Chapter (Normalization & Dim. reduction)](./04-Dimensional_reduction.md)|
[Next Chapter (Differential expression)](./06-Differential_Expression.md)
