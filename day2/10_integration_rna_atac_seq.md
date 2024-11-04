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

First, we load the reference scRNA-Seq dataset, that we previously worked with.


``` r
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("/Users/cbg-mbp-02/Documents/data/r2sc2024/pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)
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

# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks)
```

```
##                                  p_val avg_log2FC pct.1 pct.2    p_val_adj
## chr1-59814285-59815204    3.410189e-50   7.140323 0.617 0.006 2.726276e-45
## chr12-119988529-119989427 1.577804e-43   4.602373 0.617 0.032 1.261376e-38
## chr7-142808582-142809442  1.662968e-42   7.953730 0.525 0.003 1.329460e-37
## chr17-82126437-82127355   3.934464e-41  12.420762 0.500 0.000 3.145407e-36
## chr8-11016117-11017015    8.371607e-40   4.974866 0.567 0.029 6.692681e-35
## chr6-96924166-96925086    1.137904e-37   3.576527 0.633 0.067 9.096973e-33
```

Now we can visualize cell type specific peaks as follows.


``` r
plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
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
open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)
```




``` r
head(closest_genes_cd4naive)
```

```
##                           tx_id gene_name         gene_id   gene_biotype type
## ENST00000455990 ENST00000455990     HOOK1 ENSG00000134709 protein_coding  cds
## ENSE00002206071 ENST00000397558    BICDL1 ENSG00000135127 protein_coding exon
## ENST00000632998 ENST00000632998     PRSS2 ENSG00000275896 protein_coding  utr
## ENST00000665763 ENST00000665763    CCDC57 ENSG00000176155 protein_coding  gap
## ENST00000416569 ENST00000416569      XKR6 ENSG00000171044 protein_coding  gap
## ENST00000544166 ENST00000544166    KLHL32 ENSG00000186231 protein_coding  utr
##                            closest_region              query_region distance
## ENST00000455990    chr1-59815118-59815180    chr1-59814285-59815204        0
## ENSE00002206071 chr12-119989869-119990297 chr12-119988529-119989427      441
## ENST00000632998  chr7-142774509-142774564  chr7-142808582-142809442    34017
## ENST00000665763   chr17-82101867-82127691   chr17-82126437-82127355        0
## ENST00000416569    chr8-10924831-11200575    chr8-11016117-11017015        0
## ENST00000544166    chr6-96924620-96925026    chr6-96924166-96925086        0
```




``` r
head(closest_genes_cd14mono)
```

```
##                           tx_id gene_name         gene_id   gene_biotype type
## ENSE00001389095 ENST00000340607     PTGES ENSG00000148344 protein_coding exon
## ENST00000336600 ENST00000336600  C6orf223 ENSG00000181577 protein_coding  utr
## ENST00000606214 ENST00000606214    TBC1D7 ENSG00000145979 protein_coding  gap
## ENST00000554237 ENST00000554237     VASH1 ENSG00000071246 protein_coding  gap
## ENST00000373313 ENST00000373313      MAFB ENSG00000204103 protein_coding  utr
## ENST00000610832 ENST00000610832      KLF4 ENSG00000136826 protein_coding  utr
##                           closest_region             query_region distance
## ENSE00001389095 chr9-129752887-129753042 chr9-129776921-129777829    23878
## ENST00000336600   chr6-44003127-44007612   chr6-44058502-44059239    50889
## ENST00000606214   chr6-13267836-13305061   chr6-13302519-13303438        0
## ENST00000554237  chr14-76763131-76769962  chr14-76768047-76768963        0
## ENST00000373313  chr20-40688851-40689236  chr20-40688870-40689743        0
## ENST00000610832 chr9-107489168-107489766 chr9-107489471-107490352        0
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
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), 
                                      LookupGeneCoords(pbmc, "CD4"))

CoveragePlot(
  object = pbmc,
  region = "CD4",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)
```

<img src="10_integration_rna_atac_seq_files/figure-html/vis_peaks-1.png" style="display: block; margin: auto;" />

It's interesting to play with the parameters for plotting. We can zoom in/out
chromosomic regions by extending bases up/downstream to the selected



``` r
# find DA peaks overlapping gene of interest
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), 
                                      LookupGeneCoords(pbmc, "CD4"))

CoveragePlot(
  object = pbmc,
  region = "CD4",
  region.highlight = regions_highlight,
  extend.upstream = 50000,
  extend.downstream = 50000
)
```

<img src="10_integration_rna_atac_seq_files/figure-html/vis_peaks_alt-1.png" style="display: block; margin: auto;" />





``` r
sessionInfo()
```

```
## R version 4.4.0 (2024-04-24)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sonoma 14.4.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: Europe/Berlin
## tzcode source: internal
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
## [13] GenomeInfoDb_1.40.1     IRanges_2.38.1          S4Vectors_0.42.1       
## [16] BiocGenerics_0.50.0     Seurat_5.1.0            SeuratObject_5.0.2     
## [19] sp_2.1-4                Signac_1.14.0          
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22            splines_4.4.0              
##   [3] later_1.3.2                 BiocIO_1.14.0              
##   [5] bitops_1.0-8                filelock_1.0.3             
##   [7] tibble_3.2.1                polyclip_1.10-7            
##   [9] rpart_4.1.23                XML_3.99-0.17              
##  [11] fastDummies_1.7.3           lifecycle_1.0.4            
##  [13] globals_0.16.3              lattice_0.22-6             
##  [15] MASS_7.3-61                 backports_1.5.0            
##  [17] magrittr_2.0.3              Hmisc_5.2-0                
##  [19] plotly_4.10.4               sass_0.4.9                 
##  [21] rmarkdown_2.27              jquerylib_0.1.4            
##  [23] yaml_2.3.10                 httpuv_1.6.15              
##  [25] sctransform_0.4.1           spam_2.10-0                
##  [27] spatstat.sparse_3.1-0       reticulate_1.38.0          
##  [29] cowplot_1.1.3               pbapply_1.7-2              
##  [31] DBI_1.2.3                   RColorBrewer_1.1-3         
##  [33] abind_1.4-5                 zlibbioc_1.50.0            
##  [35] Rtsne_0.17                  presto_1.0.0               
##  [37] purrr_1.0.2                 biovizBase_1.52.0          
##  [39] RCurl_1.98-1.16             nnet_7.3-19                
##  [41] VariantAnnotation_1.50.0    rappdirs_0.3.3             
##  [43] GenomeInfoDbData_1.2.12     ggrepel_0.9.5              
##  [45] irlba_2.3.5.1               listenv_0.9.1              
##  [47] spatstat.utils_3.0-5        goftest_1.2-3              
##  [49] RSpectra_0.16-2             spatstat.random_3.3-1      
##  [51] fitdistrplus_1.2-1          parallelly_1.38.0          
##  [53] DelayedArray_0.30.1         leiden_0.4.3.1             
##  [55] codetools_0.2-20            RcppRoll_0.3.1             
##  [57] tidyselect_1.2.1            UCSC.utils_1.0.0           
##  [59] farver_2.1.2                base64enc_0.1-3            
##  [61] matrixStats_1.3.0           spatstat.explore_3.3-1     
##  [63] GenomicAlignments_1.40.0    jsonlite_1.8.8             
##  [65] Formula_1.2-5               progressr_0.14.0           
##  [67] ggridges_0.5.6              survival_3.7-0             
##  [69] tools_4.4.0                 ica_1.0-3                  
##  [71] Rcpp_1.0.13                 glue_1.7.0                 
##  [73] SparseArray_1.4.8           gridExtra_2.3              
##  [75] xfun_0.46                   MatrixGenerics_1.16.0      
##  [77] dplyr_1.1.4                 withr_3.0.1                
##  [79] BiocManager_1.30.23         fastmap_1.2.0              
##  [81] fansi_1.0.6                 digest_0.6.36              
##  [83] R6_2.5.1                    mime_0.12                  
##  [85] colorspace_2.1-1            scattermore_1.2            
##  [87] tensor_1.5                  dichromat_2.0-0.1          
##  [89] spatstat.data_3.1-2         RSQLite_2.3.7              
##  [91] utf8_1.2.4                  tidyr_1.3.1                
##  [93] generics_0.1.3              data.table_1.15.4          
##  [95] rtracklayer_1.64.0          S4Arrays_1.4.1             
##  [97] httr_1.4.7                  htmlwidgets_1.6.4          
##  [99] uwot_0.2.2                  pkgconfig_2.0.3            
## [101] gtable_0.3.5                blob_1.2.4                 
## [103] lmtest_0.9-40               XVector_0.44.0             
## [105] htmltools_0.5.8.1           dotCall64_1.1-1            
## [107] ProtGenerics_1.36.0         scales_1.3.0               
## [109] png_0.1-8                   spatstat.univar_3.0-0      
## [111] knitr_1.48                  rstudioapi_0.16.0          
## [113] rjson_0.2.21                reshape2_1.4.4             
## [115] checkmate_2.3.2             nlme_3.1-165               
## [117] curl_5.2.1                  cachem_1.1.0               
## [119] zoo_1.8-12                  stringr_1.5.1              
## [121] BiocVersion_3.19.1          KernSmooth_2.23-24         
## [123] parallel_4.4.0              miniUI_0.1.1.1             
## [125] foreign_0.8-87              restfulr_0.0.15            
## [127] pillar_1.9.0                grid_4.4.0                 
## [129] vctrs_0.6.5                 RANN_2.6.1                 
## [131] promises_1.3.0              xtable_1.8-4               
## [133] cluster_2.1.6               htmlTable_2.4.3            
## [135] evaluate_0.24.0             cli_3.6.3                  
## [137] compiler_4.4.0              Rsamtools_2.20.0           
## [139] rlang_1.1.4                 crayon_1.5.3               
## [141] future.apply_1.11.2         labeling_0.4.3             
## [143] plyr_1.8.9                  stringi_1.8.4              
## [145] viridisLite_0.4.2           deldir_2.0-4               
## [147] BiocParallel_1.38.0         munsell_0.5.1              
## [149] Biostrings_2.72.1           lazyeval_0.2.2             
## [151] spatstat.geom_3.3-2         Matrix_1.7-0               
## [153] BSgenome_1.72.0             RcppHNSW_0.6.0             
## [155] bit64_4.0.5                 future_1.34.0              
## [157] KEGGREST_1.44.1             shiny_1.9.1                
## [159] highr_0.11                  SummarizedExperiment_1.34.0
## [161] ROCR_1.0-11                 igraph_2.0.3               
## [163] memoise_2.0.1               bslib_0.8.0                
## [165] fastmatch_1.1-4             bit_4.0.5
```
