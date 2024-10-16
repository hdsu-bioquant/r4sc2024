---
output:
  html_document:
    keep_md: yes
---





# 4. Dimensional Reduction



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


## Exercises

> Manipulation of PCA embedded projection

 * Extract PCA embedding matrix
 * Make a PCA plot showing the first 2 principal components 
 




