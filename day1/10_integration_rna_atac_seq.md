---
title: "Single Cell ATAC-Seq Integration"
author: "Ramirez C and Herrmann C"
date: "2024-06-10"
output:
  html_document: 
    keep_md: yes
---


In this section we start using the PBMC data already preprocessed from the last 
section. Please, make sure not to skip any step.









## Gene Activity Analysis


In scATAC-seq analysis, gene activity analysis aims to infer gene expression potential from chromatin accessibility data, providing a proxy for transcriptional activity in the absence of direct RNA measurements. Since scATAC-seq profiles accessible regulatory regions, we can estimate the activity of genes by measuring the accessibility of their promoters and nearby regulatory elements, such as enhancers. This is often achieved by aggregating accessibility counts from peaks located within or near a gene body, typically within a defined distance from the transcription start site (TSS). These aggregated signals can be used to build a "gene activity matrix," which resembles an expression matrix in scRNA-seq data and enables comparative analysis. Gene activity analysis is especially valuable in multi-modal or integration workflows, as it allows researchers to connect chromatin accessibility profiles to transcriptional states, facilitating cell type annotation and regulatory network inference in single-cell studies.

The following function calculates gene activities.



``` r
## This step takes some time to run
gene.activities <- GeneActivity(pbmc)
```

The following step performs a normalization of the calculated gene activities.


``` r
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
rm(gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
```

Then, we can visualize specific gene activity in clusters by inspecting
the UMAP projection.



``` r
DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
```

<img src="10_integration_rna_atac_seq_files/figure-html/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

## Label transfer from RNA to ATAC


Since scRNA-seq captures gene expression and scATAC-seq captures regulatory region 
accessibility, label transfer leverages the shared biological information between 
the two modalities. This is typically achieved by first clustering the scRNA-seq data 
(the reference) and assigning cell-type labels based on expression profiles. Then, a 
correspondence is established between the RNA and ATAC data by integrating both datasets 
in a shared low-dimensional space, often using anchor-based methods or common latent spaces.
The scRNA-seq labels can then be transferred to the scATAC-seq cells based on their 
proximity to scRNA-seq cells in this space. This approach enables researchers to 
annotate cell types in scATAC-seq data even in the absence of direct transcriptional 
information, improving the interpretation of chromatin accessibility landscapes 
and supporting integrative analyses across modalities.

First, we load the reference scRNA-Seq dataset, that you can download 
[here](https://figshare.com/s/5e0afc577026e2e7413a).






``` r
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("/home/carlos/Documents/r4sc2024/day2/data/pbmc_1k_v3.rds")
pbmc_rna
```

```
## An object of class Seurat 
## 19089 features across 1000 samples within 1 assay 
## Active assay: RNA (19089 features, 3000 variable features)
##  3 layers present: counts, data, scale.data
##  3 dimensional reductions calculated: pca, tsne, umap
```

Using the shared low-dimensional space, algorithms (such as those in Seurat) identify anchors by matching cells in the scRNA-seq dataset to similar cells in the scATAC-seq dataset. Anchors are defined based on the similarity of cell profiles and are typically selected using nearest-neighbor or mutual nearest neighbor (MNN) methods, ensuring each anchor reflects genuine similarity rather than noise.



``` r
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)
```

Once anchors are established, cell-type labels from the scRNA-seq dataset can be transferred to the scATAC-seq cells based on anchor assignments. Cells that share an anchor are likely to represent the same cell type or state, enabling accurate annotation of cell types in scATAC-seq data.



``` r
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
```



``` r
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2
```

<img src="10_integration_rna_atac_seq_files/figure-html/umap_plots_label_transfer-1.png" style="display: block; margin: auto;" />

## Identification of cell type chromosomic features

First, we filter out under-represented cell types since is difficult 
to state valid conclusions based in few cells.



``` r
predicted_id_counts <- table(pbmc$predicted.id)

# Identify the predicted.id values that have more than 20 cells
major_predicted_ids <- names(predicted_id_counts[predicted_id_counts > 20])
pbmc <- pbmc[, pbmc$predicted.id %in% major_predicted_ids]
```


In the following step the function FindMarkers() is employed to identify 
differentially accessible peaks between two specified cell types: CD4 Naive 
and CD14+ Monocytes. The Wilcoxon test is set as the default statistical method 
for comparing the accessibility profiles of these two cell types, and a minimum 
percentage threshold of 0.1 is applied to ensure that only peaks present in at least 
10% of the cells in either group are considered. Finally, the results of the 
differential accessibility analysis are displayed using the head function, 
which shows the top entries of the resulting data frame da_peaks.



``` r
# change cell identities to the per-cell predicted labels
Idents(pbmc) <- pbmc$predicted.id

# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

## Setting cell types to compare
cell_type_1 <- "B cell progenitor"
#cell_type_1 <- "pDC"
cell_type_2 <- "pre-B cell"

# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = cell_type_1,
  ident.2 = cell_type_2,
  test.use = 'wilcox',
  min.pct = 0.1
)

tail(da_peaks)
```

```
##                          p_val avg_log2FC pct.1 pct.2 p_val_adj
## chr8-133209008-133209914     1 -0.2021179 0.111  0.10         1
## chr9-107340919-107341821     1 -0.3892266 0.111  0.10         1
## chr9-131178035-131178937     1 -0.4869067 0.111  0.10         1
## chr11-64449391-64450120      1 -0.1083792 0.167  0.16         1
## chr3-13548722-13549622       1 -0.1467612 0.167  0.16         1
## chr10-30455881-30456731      1 -0.1524558 0.139  0.14         1
```

Now we can visualize cell type specific peaks as follows.


``` r
plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c(cell_type_1, cell_type_2)
)
plot2 <- FeaturePlot(
  object = subset(pbmc,
                  predicted.id %in% c(cell_type_1, cell_type_2)
                  ),
  features = rownames(da_peaks)[1],
  pt.size = 2
)

plot1 | plot2
```

<img src="10_integration_rna_atac_seq_files/figure-html/diff_peaks_umap-1.png" style="display: block; margin: auto;" />



## Interpreting identified peaks

In order to drive conclusions about the identified cell type specific open peaks
it's very useful to identify which genomic locations are closed to that
positions. The `ClosestFeature` is used to annotate peaks to closest annotated
features in chromosome references.



``` r
open_cell_type_1 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cell_type_2 <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cell_type_1 <- ClosestFeature(pbmc, regions = open_cell_type_1)
closest_genes_cell_type_2 <- ClosestFeature(pbmc, regions = open_cell_type_2)
```




``` r
head(closest_genes_cell_type_1)
```

```
##                           tx_id gene_name         gene_id   gene_biotype type
## ENSE00001239301 ENST00000576742     TRPV3 ENSG00000167723 protein_coding exon
## ENST00000524347 ENST00000524347      SGCD ENSG00000170624 protein_coding  gap
## ENST00000672146 ENST00000672146     OPRL1 ENSG00000125510 protein_coding  utr
## ENST00000492665 ENST00000492665    ZBTB20 ENSG00000181722 protein_coding  gap
## ENST00000270162 ENST00000270162      SIK1 ENSG00000142178 protein_coding  utr
## ENST00000182527 ENST00000182527     TRAM2 ENSG00000065308 protein_coding  utr
##                           closest_region             query_region distance
## ENSE00001239301    chr17-3557676-3557995    chr17-3558417-3559246      421
## ENST00000524347 chr5-156393866-156508600 chr5-156402205-156403126        0
## ENST00000672146  chr20-64098887-64100633  chr20-64099915-64100843        0
## ENST00000492665 chr3-114389106-114500351 chr3-114427942-114428782        0
## ENST00000270162  chr21-43414483-43416741  chr21-43352428-43353343    61139
## ENST00000182527   chr6-52497408-52503196   chr6-52497518-52498336        0
```




``` r
head(closest_genes_cell_type_2)
```

```
##                           tx_id gene_name         gene_id   gene_biotype type
## ENST00000340722 ENST00000340722     TCL1B ENSG00000213231 protein_coding  utr
## ENSE00003713479 ENST00000517566      OXR1 ENSG00000164830 protein_coding exon
## ENST00000395153 ENST00000395153     DACT1 ENSG00000165617 protein_coding  cds
## ENST00000652626 ENST00000652626   GUCY1B1 ENSG00000061918 protein_coding  cds
## ENSE00003630226 ENST00000603223     GFOD1 ENSG00000145990 protein_coding exon
## ENST00000498049 ENST00000498049   C2orf76 ENSG00000186132 protein_coding  gap
##                           closest_region             query_region distance
## ENST00000340722  chr14-95691931-95692628  chr14-95695641-95696543     3012
## ENSE00003713479 chr8-106359476-106359636 chr8-106269756-106270665    88810
## ENST00000395153  chr14-58638203-58638547  chr14-58637398-58638262        0
## ENST00000652626 chr4-155759141-155759143 chr4-155758521-155759415        0
## ENSE00003630226   chr6-13486638-13487662   chr6-13526432-13527329    38769
## ENST00000498049 chr2-119339972-119366789 chr2-119355536-119356453        0
```


## Visualization of peaks

For peaks visualization we can use built-in function like `CoveragePlot`.
This function generates coverage plots that display the number of fragments 
(or reads) mapped to specified peaks or genomic features, such as promoters 
or enhancers, in individual cells or aggregated across a group of cells. 
By providing a graphical representation of the accessibility landscape, 
CoveragePlot helps identify patterns of chromatin accessibility that may 
correlate with gene activity or regulatory mechanisms. 
The function can also highlight differences in accessibility between different 
cell types or conditions, making it a valuable tool for interpreting the 
functional significance of chromatin dynamics in single-cell ATAC-seq analyses. 
Additionally, users can customize the plot to focus on specific genes or genomic 
regions of interest, facilitating a targeted exploration of the data.

Here, we visualize a peak identified to be open in CD4 Naive T cells. The peak
is located by looking at a gene and annotated peaks to it. Here, we use the
CD4 gene. 


``` r
## Sorting samples by similarity
pbmc <- SortIdents(pbmc)

# find DA peaks overlapping gene of interest
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cell_type_1), 
                                      LookupGeneCoords(pbmc, "TRAM2"))

CoveragePlot(
  object = pbmc,
  region = "TRAM2",
  region.highlight = regions_highlight,
  extend.upstream = 100,
  extend.downstream = 100, 
  idents = c(cell_type_1, cell_type_2),
)
```

<img src="10_integration_rna_atac_seq_files/figure-html/vis_peaks-1.png" style="display: block; margin: auto;" />

It's interesting to play with the parameters for plotting. We can zoom in/out
chromosomic regions by extending bases up/downstream to the selected.





## Exercises 

**Exercise 1**

<br>

<details>
<summary> How many peaks are differentially open the previous
comparison.
TIP: Check at the dimensions of the table of differential peaks
after taking the statistically significant peaks.
</summary>

<b>Answer:</b>
<br>
<tt> There are 14 differentially open and statistically significant peaks.
You can check this by running.
``` 
da_peaks %>% dplyr::filter(p_val_adj<0.05)
```
</tt>
<br>
</details> 

<br>


**Exercise 2**

<br>

<details>
<summary> Check the same comparing now pre-B and pDC cell types.
TIP: Re-run the Differential peak analysis.
</summary>

<b>Answer:</b>
<br>
<tt> There are 428 differentially open and statistically significant peaks.
You can check this by running.

``` 
da_peaks %>% dplyr::filter(p_val_adj<0.05)
```

<br>

Is this expected? Try visualizing some differential peak in a coverage
plot as before for those cell types.

</tt>
<br>
</details> 

<br>

**Exercise 3**

<details>
<summary> Use the code day 2 to create a vulcano plot of the
differentially open peaks.
</summary>

</details> 

<br>



``` r
sessionInfo()
```

```
## R version 4.4.1 (2024-06-14)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 22.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Berlin
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] ensembldb_2.28.1        AnnotationFilter_1.28.0 GenomicFeatures_1.56.0 
##  [4] AnnotationDbi_1.66.0    Biobase_2.64.0          AnnotationHub_3.12.0   
##  [7] BiocFileCache_2.12.0    dbplyr_2.5.0            hdf5r_1.3.11           
## [10] patchwork_1.3.0         ggplot2_3.5.1           GenomicRanges_1.56.2   
## [13] GenomeInfoDb_1.40.1     IRanges_2.38.0          S4Vectors_0.42.0       
## [16] BiocGenerics_0.50.0     Seurat_5.1.0            SeuratObject_5.0.2     
## [19] sp_2.1-4                Signac_1.14.0          
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22            splines_4.4.1              
##   [3] later_1.3.2                 BiocIO_1.14.0              
##   [5] bitops_1.0-9                filelock_1.0.3             
##   [7] tibble_3.2.1                polyclip_1.10-7            
##   [9] rpart_4.1.19                XML_3.99-0.17              
##  [11] fastDummies_1.7.4           lifecycle_1.0.4            
##  [13] globals_0.16.3              lattice_0.21-8             
##  [15] MASS_7.3-60                 backports_1.5.0            
##  [17] magrittr_2.0.3              Hmisc_5.2-0                
##  [19] plotly_4.10.4               sass_0.4.9                 
##  [21] rmarkdown_2.29              jquerylib_0.1.4            
##  [23] yaml_2.3.8                  httpuv_1.6.15              
##  [25] sctransform_0.4.1           spam_2.11-0                
##  [27] spatstat.sparse_3.1-0       reticulate_1.39.0          
##  [29] cowplot_1.1.3               pbapply_1.7-2              
##  [31] DBI_1.2.3                   RColorBrewer_1.1-3         
##  [33] abind_1.4-5                 zlibbioc_1.50.0            
##  [35] Rtsne_0.17                  purrr_1.0.2                
##  [37] biovizBase_1.52.0           RCurl_1.98-1.16            
##  [39] nnet_7.3-19                 VariantAnnotation_1.50.0   
##  [41] rappdirs_0.3.3              GenomeInfoDbData_1.2.12    
##  [43] ggrepel_0.9.6               irlba_2.3.5.1              
##  [45] listenv_0.9.1               spatstat.utils_3.1-1       
##  [47] goftest_1.2-3               RSpectra_0.16-1            
##  [49] spatstat.random_3.3-2       fitdistrplus_1.2-1         
##  [51] parallelly_1.39.0           DelayedArray_0.30.1        
##  [53] leiden_0.4.3.1              codetools_0.2-19           
##  [55] RcppRoll_0.3.1              tidyselect_1.2.1           
##  [57] UCSC.utils_1.0.0            farver_2.1.2               
##  [59] base64enc_0.1-3             matrixStats_1.3.0          
##  [61] spatstat.explore_3.3-3      GenomicAlignments_1.40.0   
##  [63] jsonlite_1.8.8              Formula_1.2-5              
##  [65] progressr_0.15.0            ggridges_0.5.6             
##  [67] survival_3.5-5              tools_4.4.1                
##  [69] ica_1.0-3                   Rcpp_1.0.12                
##  [71] glue_1.7.0                  SparseArray_1.4.8          
##  [73] gridExtra_2.3               xfun_0.45                  
##  [75] MatrixGenerics_1.16.0       dplyr_1.1.4                
##  [77] withr_3.0.0                 BiocManager_1.30.23        
##  [79] fastmap_1.2.0               fansi_1.0.6                
##  [81] digest_0.6.34               R6_2.5.1                   
##  [83] mime_0.12                   colorspace_2.1-0           
##  [85] scattermore_1.2             tensor_1.5                 
##  [87] dichromat_2.0-0.1           spatstat.data_3.1-2        
##  [89] RSQLite_2.3.7               utf8_1.2.4                 
##  [91] tidyr_1.3.1                 generics_0.1.3             
##  [93] data.table_1.15.4           rtracklayer_1.64.0         
##  [95] S4Arrays_1.4.1              httr_1.4.7                 
##  [97] htmlwidgets_1.6.4           uwot_0.2.2                 
##  [99] pkgconfig_2.0.3             gtable_0.3.5               
## [101] blob_1.2.4                  lmtest_0.9-40              
## [103] XVector_0.44.0              htmltools_0.5.8.1          
## [105] dotCall64_1.2               ProtGenerics_1.36.0        
## [107] scales_1.3.0                png_0.1-8                  
## [109] spatstat.univar_3.1-1       knitr_1.47                 
## [111] rstudioapi_0.17.1           rjson_0.2.23               
## [113] reshape2_1.4.4              checkmate_2.3.2            
## [115] nlme_3.1-162                curl_5.2.1                 
## [117] cachem_1.1.0                zoo_1.8-12                 
## [119] stringr_1.5.1               BiocVersion_3.19.1         
## [121] KernSmooth_2.23-22          parallel_4.4.1             
## [123] miniUI_0.1.1.1              foreign_0.8-82             
## [125] restfulr_0.0.15             pillar_1.9.0               
## [127] grid_4.4.1                  vctrs_0.6.5                
## [129] RANN_2.6.2                  promises_1.3.0             
## [131] xtable_1.8-4                cluster_2.1.4              
## [133] htmlTable_2.4.3             evaluate_0.23              
## [135] cli_3.6.3                   compiler_4.4.1             
## [137] Rsamtools_2.20.0            rlang_1.1.4                
## [139] crayon_1.5.3                future.apply_1.11.3        
## [141] labeling_0.4.3              plyr_1.8.9                 
## [143] stringi_1.8.4               viridisLite_0.4.2          
## [145] deldir_2.0-4                BiocParallel_1.38.0        
## [147] munsell_0.5.1               Biostrings_2.72.1          
## [149] lazyeval_0.2.2              spatstat.geom_3.3-3        
## [151] Matrix_1.7-1                BSgenome_1.72.0            
## [153] RcppHNSW_0.6.0              bit64_4.5.2                
## [155] future_1.34.0               KEGGREST_1.44.1            
## [157] shiny_1.9.1                 highr_0.11                 
## [159] SummarizedExperiment_1.34.0 ROCR_1.0-11                
## [161] igraph_2.1.1                memoise_2.0.1              
## [163] bslib_0.7.0                 fastmatch_1.1-4            
## [165] bit_4.5.0
```
