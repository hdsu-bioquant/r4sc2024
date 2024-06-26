---
output:
  html_document:
    keep_md: yes
---





# Dimensional Reduction



The size of scRNA-Seq matrices can be huge and for this reason techniques to reduce the dimensionality
of this data are used. Here, we will use PCA, a very common techniques for dimension
reduction and visualization.

We will run a PCA using the already calculated top 1000 HVGs using the function `RunPCA()`.


```r
pbmc.filtered <- RunPCA(pbmc.filtered, 
                        features = VariableFeatures(pbmc.filtered))
```

We can assess the dimensionality, a measure of the complexity, by using an
elbow plot of the standard deviation for each principal component (PC)
from the PCA.

We will use the function `ElbowPlot()`.



```r
ElbowPlot(pbmc.filtered)
```

<img src="04-Normalization_and_Dimensional_Reduction_files/figure-html/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

The PC components in a PCA reflects corresponds to the directions in which
more variability is observed. These PCs are ranked by using the eigenvalues
of the covariance matrix. We can the plot a Elbow or joystick plot of the 
standard deviation and the rank of each PC. Top ranked PCs are expected to 
have higher values of variability and then to gradually decrease. So, we
can use the elbow plot representation to keep PCs from the top to the bottom
until we do not see further variability changes, in these case we can use 
the number of PC equal to 7.


## Quizzes

### Quiz 1

<details>
<summary> <b>Which command(s) can be used to extract the PCA matrix from the seurat object?</b>
<br>

1. <code>pca <- Embeddings(pbmc.filtered, reduction = 'pca')</code>
2. <code>pca <- pbmc.filtered@reductions$pca@cell.embeddings</code>
3. <code>pca <- pbmc.filtered@pca$reductions@cell.embeddings</code>
</summary>
TIP: Use str(pbmc.filtered) to explore the slots present in the seurat object. Two options
are correct.
</details>

### Quiz 2

<details>
<summary>
Look at the feature loadings of the principal components. Which are the genes with the highest/lowest loadings for PC1?

1. highest: MTMR11,NTMT1,DGAT1,CCL4L1,LDHB lowest: OSBP,DFFB,EFCAB2,RNF170,POMZP3 
2. highest: MALAT1,IL32,LTB,CD3D,LDHB lowest: S100A9,CST3,S100A8,LYZ,FTL 
3. highest: C2CD2L,C16orf58,CENPL,ZNF181,RABL5 lowest: INO80D,CHST2,OAF,LAGAP2YZ,ADRB2 

</summary>
<code>
> sort(pbmc.filtered@reductions$pca@feature.loadings[,1])[1:5]
    S100A9       CST3     S100A8        LYZ        FTL 
-0.1380440 -0.1375347 -0.1338024 -0.1336727 -0.1324367 
> sort(pbmc.filtered@reductions$pca@feature.loadings[,1],decreasing=TRUE)[1:5]
    MALAT1       IL32        LTB       CD3D       LDHB 
0.11474983 0.07876368 0.07786466 0.07604930 0.06500361 
</code>
</details>

## Exercise

> Perform the PCA analysis on the `pbmc200.seurat` object!


[Previous Chapter (Feature selection)](./03-Feature_selection.md)|
[Next Chapter (Clustering)](./05-Cluster_visualization.md)


