# Feature selection


Because of the sparsity in the sequencing data many genes or features are almost no expressed.
Additionally, some genes are constantly expressed across cells. These features are then probably
not playing any function in cells and on the other hand can just add noise and unnecessary complexity
to further analysis. Then, it's usual to remove genes with very low variability and to select
only top highly variable genes (HVG).

We will use the function `FindVariableFeatures()` to calculate the top most variable genes.
The parameter nfeatures is used to set the number of top selected genes. We set to the top
1000 features.


```r
pbmc.filtered <- FindVariableFeatures(pbmc.filtered, nfeatures = 1000)
```


We can access to the top 1000 variable features using the VariableFeatures function. In the next
chunk we display the top first 6 (head) of this set. 



```r
head(VariableFeatures(pbmc.filtered))
```

```
## [1] "S100A9" "LYZ"    "S100A8" "FTL"    "GNLY"   "FTH1"
```


In the next scatter plot we can see the average expression *vs* the standardized variance for each feature.
Genes in red are the selected HVG.



```r
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc.filtered)
plot1 <- LabelPoints(plot = plot1, 
                     points = head(VariableFeatures(pbmc.filtered),
                                   10), 
                     repel = TRUE)
plot1 
```

<img src="03-Feature_selection_files/figure-html/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

For further analysis we will use only HVGs. 

## Quizzes

Repeat this analysis using the `pbmc200.seurat` object.

<details> 
Normalize the gene expression values from the Seurat object.
Calculate variances manually from the matrix. Sort genes based on variances in 
decreasing order and show top 6 genes.
<summary> Find the top 6 genes with the highest variance in descending order
<br>
a) HLA-DRA, CST3, S100A8, NKG7, S100A9, LYZ
<br>
b) NKG7, S100A8, S100A9, LYZ, CST3, HLA-DRA
<br>
c) HLA-DRA, LYZ, NKG7, S100A9, TYROBP, CST3 
</summary>
<br>
<b>Answer:</b>
<br>

<code>
pbmc200.seurat<- NormalizeData(pbmc200.seurat) %>%
                      ScaleData()
norm.exp <- GetAssayData(pbmc200.seurat, slot = 'data')
</code>

<b>Calculation of the variance in genes</b>
<code>std.devs <- apply(norm.exp, 1, var)</code>

<b>Showing the top 6 genes with highest variance</b>
<code>head(sort(std.devs, decreasing = T))</code>

</details>

[Previous Chapter (Quality control)](./02-Quality_control.md)|
[Next Chapter (Normalization & Dim. reduction)](./04-Normalization_and_Dimensional_Reduction.md)