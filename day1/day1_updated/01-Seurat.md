---
output:
  html_document:
    keep_md: yes
---




# 1. Standard Preprocessing using Seurat


Some standard steps are usually carried out in scRNA-Seq prior to further analysis as QC, dimensional
reduction and marker visualization. Here, we will use the Seurat R package to perform these steps which
is increasingly becoming the most popular tool, however, there are some other options as SingleCellExperiment
in R and scanpy available for python. First, we need to define a Seurat object.








## Creating Seurat object

We create a Seurat object using the `CreateSeuratObject` function as follows. Use the `?function` helper
in R to get information about the parameters that are need to be provided to the function.

 * counts - a count matrix. It can be a matrix, sparse matrix or dataframe.
 * project - A single character string. Correspond to a arbitrary name to label the object.
 * assay - A single character string. An arbitrary name that usually is assigned to label the
            type of sequencing information, for example, RNA, spliced RNA, ATAC, etc...
 * min.cells - An integer. Indicates a threshold of the number of cells for which a feature was
              recorded. Cells below that threshold will be filtered out.
 * min.features - An integer. Similar to min cells but for the number of features of a cell.
 



``` r
ad_medula_seurat <- CreateSeuratObject(
  counts = counts, 
  project = 'human_ad_medulla', 
  assay = 'RNA', 
  min.cells = 1, 
  min.features = 1
)
Idents(ad_medula_seurat) <- 'human_ad_medulla'
```

The variable `ad_medula_seurat` now contains the Seurat object that we can feed into the package.
If we print the variable we get information about the number of genes and cells.



``` r
ad_medula_seurat
```

```
## An object of class Seurat 
## 25064 features across 2400 samples within 1 assay 
## Active assay: RNA (25064 features, 0 variable features)
##  1 layer present: counts
```



## Exploring the Seurat object


Seurat objects can be seen as a container of different features. At this step it contains
our gene expression matrix, but in addition it can store metadata, processed data,
information from different assays, for example, scATACSeq, scCITESeq or unspliced transcripts.

We can explore the seurat object using the `$` to explore its *metadata* in combination with the tab 
key. For example, during the creation of the seurat object the number of counts quality metric
is calculated and added to the metadata. We can explore this metric by accessing the metadata
as follows.



``` r
ad_medula_seurat$nCount_RNA %>% head
```

```
## AAACGCTGTCAAAGTA_1 AAAGGATCAGATCACT_1 AAAGGATCAGCAGAAC_1 AACAACCAGTCCTACA_1 
##           2673.072           2563.070           2593.677           3183.051 
## AACAAGAAGGACTTCT_1 AACAAGATCTGCGTCT_1 
##           2718.732           3705.433
```


 We can do the same with the `@` operator to explore the different *slots*. For example,
 we can extract the original count matrix that we used to create the seurat object as follows:
 
 


``` r
ad_medula_seurat@assays$RNA$counts[1:5, 1:5]
```

```
## 5 x 5 sparse Matrix of class "dgCMatrix"
##               AAACGCTGTCAAAGTA_1 AAAGGATCAGATCACT_1 AAAGGATCAGCAGAAC_1
## RP11-34P13.7                   .                  .                  .
## AL627309.1                     .                  .                  .
## AP006222.2                     .                  .                  .
## RP4-669L17.10                  .                  .                  .
## RP5-857K21.2                   .                  .                  .
##               AACAACCAGTCCTACA_1 AACAAGAAGGACTTCT_1
## RP11-34P13.7                   .                  .
## AL627309.1                     .                  .
## AP006222.2                     .                  .
## RP4-669L17.10                  .                  .
## RP5-857K21.2                   .                  .
```


## Extracting expression values 


Next, we visualize gene counts to see its behavior. We take a look at the expression of the 
house keeping gene ACTIN Beta and plot an histogram of count values. We will use the
the function `FetchData` which is used to extract values from selected features in the Seurat
object and then plot it using an histogram.




``` r
actin <- FetchData(ad_medula_seurat, vars = 'ACTB')
hist(actin$ACTB)
```

![](01-Seurat_files/figure-html/unnamed-chunk-7-1.png)<!-- -->





## Quizzes

<br>
<details>
<summary> Find and display the metadata in the seurat object
<br>
TIP: You can have 
a look at the [documentation](https://github.com/satijalab/seurat/wiki/Seurat#object-information) 
of the seurat objects from the GitHub Wiki.
</summary>
<br>
<b>Answer:</b>
<br>

```r
ad_medula_seu_filtered@meta.data %>% head
```

</details> 



## Exercises

<blockquote>
Create a Seurat object 

The file in the follwing URL:

`https://raw.githubusercontent.com/caramirezal/caramirezal.github.io/master/bookdown-minimal/data/pbmc_10X_250_cells.tsv` 

contains 250 cells downsampled from the 10x PBMC data and stored in tsv format

 * Load the count matrix in tsv format
 
 * Create a Seurat object using the count matrix
 
 * How many features and cells are present in the count matrix?

</blockquote>



[Previous Chapter (Seurat)](./01-Seurat.md)|
[Next Chapter (Feature selection)](./03-Feature_selection.md)





