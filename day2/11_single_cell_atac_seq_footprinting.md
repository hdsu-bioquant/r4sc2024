---
title: "Single Cell ATAC-Seq Motif Footprinting"
author: "Ramirez C and Herrmann C"
date: "2024-06-10"
output:
  html_document: 
    keep_md: yes
---










 In single-cell ATAC-seq (scATAC-seq) footprinting analysis, the goal is to achieve fine-grained identification of transcription factor (TF) binding sites at the single-cell level. In ATAC-seq, footprinting analysis leverages the fact that TFs bound to DNA protect the underlying sequence from transposase insertion, creating “footprints” or regions with fewer fragments. These footprints help pinpoint active regulatory elements and offer insights into TF binding activity in specific cell types or states.

Using Signac, an R package for single-cell chromatin data, scATAC-seq footprinting can be combined with clustering and motif enrichment analyses to map out how chromatin accessibility changes in relation to TF binding across individual cells. This approach is highly valuable for identifying TF regulatory networks, understanding cell identity, and exploring the mechanisms behind cellular differentiation and response. Signac’s streamlined functions make it possible to visualize these footprints and link them to TF motifs, ultimately revealing how specific TFs contribute to gene regulation in various cellular contexts.




``` r
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(AnnotationHub)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
```






## Loading motifs information

First we use the function `getMatrixSet()` to retireve motifs from the JASPAR database, which are collections of transcription factor (TF) binding motifs. In the context of single-cell ATAC-seq, these motifs represent sequence patterns that TFs typically recognize and bind to in the genome.
getMatrixSet() is often used in combination with JASPAR2020 or other motif libraries, allowing you to pull a set of motifs for specific TFs or a group of TFs, based on various parameters such as species or collection type.

Then we add motif information to the Seurat object in Signac by associating the retrieved motif matrices (like those from getMatrixSet()) with accessible chromatin regions in the single-cell ATAC-seq dataset.
Using motifmatchr under the hood, it scans these regions for sequence matches to the motifs, effectively identifying where specific TFs may bind in the genome. This information is essential for downstream motif enrichment and footprinting analyses, as it lets you explore how TF binding sites are distributed across different clusters or conditions.



``` r
# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
pbmc <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
```


## Performing Footprint


A footprint Analysis Workflow in ATAC-Seq performs the following steps: (i) Motif Matching: 
Known TF binding motifs are mapped onto the peaks to predict where each TF could theoretically bind based on sequence patterns. (ii) Footprint Detection: For each motif-matched region, the read coverage is examined closely to detect a footprint pattern (low coverage within the motif, with higher coverage flanking it). This pattern is taken as a sign of TF binding. And (iii) Aggregation and Scoring: To
enhance signal detection, footprints are often aggregated across many cells and regions, allowing for more robust detection. Various scoring metrics are then applied to rank the likelihood of TF binding.


In a TN5 insertion enrichment vs. distance from motif plot, the goal is to visualize how transposase activity is influenced by transcription factor (TF) binding near a motif site, which helps reveal footprint patterns in ATAC-seq data. ATAC-seq uses the TN5 transposase enzyme to cut open regions of chromatin and insert sequencing adapters. This insertion activity provides a proxy for chromatin accessibility.
The following footprint plot can be used to evaluate
when a TF binds to a motif in DNA, it blocks TN5 insertion at the motif’s exact site while leaving adjacent regions accessible. This creates a pattern in the sequencing data: lower insertion frequencies (coverage dips) at the motif site with higher insertion frequencies just outside this region.



``` r
# gather the footprinting information for sets of motifs
pbmc <- Footprint(
  object = pbmc,
  motif.name = c("CEBPA", "EBF1", "SPI1"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)


# plot the footprint data for each group of cells
p2 <- PlotFootprint(pbmc, 
                    features = c("CEBPA", "EBF1", "SPI1"), 
                    group.by = "predicted.id", 
                    show.expected = FALSE, 
                    repel = TRUE, 
                    idents = c("pre-B cell",
                               "CD14+ Monocytes",
                               "CD16+ Monocytes",
                               "pDC",
                               "CD4 Naive"), 
                    label.top = 5)
p2 + patchwork::plot_layout(ncol = 1)
```

<img src="11_single_cell_atac_seq_footprinting_files/figure-html/footprinting-1.png" style="display: block; margin: auto;" />


Do you notice any difference in binding for the TFs shown above. 

* Is there a differential binding for CEBPA in mononuclear compared to lymphocytes?

* Can we say something about EBF1?

* Is SPI1 involved in the maturation of PBMCs?




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
##  [1] ensembldb_2.28.1                  AnnotationFilter_1.28.0          
##  [3] GenomicFeatures_1.56.0            AnnotationDbi_1.66.0             
##  [5] Biobase_2.64.0                    BSgenome.Hsapiens.UCSC.hg38_1.4.5
##  [7] BSgenome_1.72.0                   rtracklayer_1.64.0               
##  [9] BiocIO_1.14.0                     Biostrings_2.72.1                
## [11] XVector_0.44.0                    TFBSTools_1.42.0                 
## [13] JASPAR2020_0.99.10                motifmatchr_1.26.0               
## [15] AnnotationHub_3.12.0              BiocFileCache_2.12.0             
## [17] dbplyr_2.5.0                      hdf5r_1.3.11                     
## [19] patchwork_1.3.0                   ggplot2_3.5.1                    
## [21] GenomicRanges_1.56.2              GenomeInfoDb_1.40.1              
## [23] IRanges_2.38.1                    S4Vectors_0.42.1                 
## [25] BiocGenerics_0.50.0               Seurat_5.1.0                     
## [27] SeuratObject_5.0.2                sp_2.1-4                         
## [29] Signac_1.14.0                    
## 
## loaded via a namespace (and not attached):
##   [1] ProtGenerics_1.36.0         matrixStats_1.3.0          
##   [3] spatstat.sparse_3.1-0       bitops_1.0-8               
##   [5] DirichletMultinomial_1.46.0 httr_1.4.7                 
##   [7] RColorBrewer_1.1-3          backports_1.5.0            
##   [9] tools_4.4.0                 sctransform_0.4.1          
##  [11] utf8_1.2.4                  R6_2.5.1                   
##  [13] lazyeval_0.2.2              uwot_0.2.2                 
##  [15] withr_3.0.1                 gridExtra_2.3              
##  [17] progressr_0.14.0            cli_3.6.3                  
##  [19] spatstat.explore_3.3-1      fastDummies_1.7.3          
##  [21] labeling_0.4.3              sass_0.4.9                 
##  [23] spatstat.data_3.1-2         readr_2.1.5                
##  [25] ggridges_0.5.6              pbapply_1.7-2              
##  [27] Rsamtools_2.20.0            foreign_0.8-87             
##  [29] R.utils_2.12.3              dichromat_2.0-0.1          
##  [31] parallelly_1.38.0           rstudioapi_0.16.0          
##  [33] RSQLite_2.3.7               generics_0.1.3             
##  [35] gtools_3.9.5                ica_1.0-3                  
##  [37] spatstat.random_3.3-1       dplyr_1.1.4                
##  [39] GO.db_3.19.1                Matrix_1.7-0               
##  [41] fansi_1.0.6                 abind_1.4-5                
##  [43] R.methodsS3_1.8.2           lifecycle_1.0.4            
##  [45] yaml_2.3.10                 SummarizedExperiment_1.34.0
##  [47] SparseArray_1.4.8           Rtsne_0.17                 
##  [49] grid_4.4.0                  blob_1.2.4                 
##  [51] promises_1.3.0              crayon_1.5.3               
##  [53] pwalign_1.0.0               miniUI_0.1.1.1             
##  [55] lattice_0.22-6              cowplot_1.1.3              
##  [57] annotate_1.82.0             KEGGREST_1.44.1            
##  [59] pillar_1.9.0                knitr_1.48                 
##  [61] rjson_0.2.21                future.apply_1.11.2        
##  [63] codetools_0.2-20            fastmatch_1.1-4            
##  [65] leiden_0.4.3.1              glue_1.7.0                 
##  [67] spatstat.univar_3.0-0       data.table_1.15.4          
##  [69] vctrs_0.6.5                 png_0.1-8                  
##  [71] spam_2.10-0                 gtable_0.3.5               
##  [73] poweRlaw_0.80.0             cachem_1.1.0               
##  [75] xfun_0.46                   S4Arrays_1.4.1             
##  [77] mime_0.12                   pracma_2.4.4               
##  [79] survival_3.7-0              RcppRoll_0.3.1             
##  [81] fitdistrplus_1.2-1          ROCR_1.0-11                
##  [83] nlme_3.1-165                bit64_4.0.5                
##  [85] filelock_1.0.3              RcppAnnoy_0.0.22           
##  [87] bslib_0.8.0                 irlba_2.3.5.1              
##  [89] rpart_4.1.23                KernSmooth_2.23-24         
##  [91] Hmisc_5.2-0                 colorspace_2.1-1           
##  [93] seqLogo_1.70.0              DBI_1.2.3                  
##  [95] nnet_7.3-19                 tidyselect_1.2.1           
##  [97] bit_4.0.5                   compiler_4.4.0             
##  [99] curl_5.2.1                  htmlTable_2.4.3            
## [101] DelayedArray_0.30.1         plotly_4.10.4              
## [103] checkmate_2.3.2             scales_1.3.0               
## [105] caTools_1.18.2              lmtest_0.9-40              
## [107] rappdirs_0.3.3              stringr_1.5.1              
## [109] digest_0.6.36               goftest_1.2-3              
## [111] spatstat.utils_3.0-5        rmarkdown_2.27             
## [113] base64enc_0.1-3             htmltools_0.5.8.1          
## [115] pkgconfig_2.0.3             MatrixGenerics_1.16.0      
## [117] highr_0.11                  fastmap_1.2.0              
## [119] rlang_1.1.4                 htmlwidgets_1.6.4          
## [121] UCSC.utils_1.0.0            shiny_1.9.1                
## [123] farver_2.1.2                jquerylib_0.1.4            
## [125] zoo_1.8-12                  jsonlite_1.8.8             
## [127] BiocParallel_1.38.0         R.oo_1.26.0                
## [129] VariantAnnotation_1.50.0    RCurl_1.98-1.16            
## [131] magrittr_2.0.3              Formula_1.2-5              
## [133] GenomeInfoDbData_1.2.12     dotCall64_1.1-1            
## [135] munsell_0.5.1               Rcpp_1.0.13                
## [137] reticulate_1.38.0           stringi_1.8.4              
## [139] zlibbioc_1.50.0             MASS_7.3-61                
## [141] plyr_1.8.9                  parallel_4.4.0             
## [143] listenv_0.9.1               ggrepel_0.9.5              
## [145] deldir_2.0-4                CNEr_1.40.0                
## [147] splines_4.4.0               tensor_1.5                 
## [149] hms_1.1.3                   igraph_2.0.3               
## [151] spatstat.geom_3.3-2         RcppHNSW_0.6.0             
## [153] reshape2_1.4.4              TFMPvalue_0.0.9            
## [155] BiocVersion_3.19.1          XML_3.99-0.17              
## [157] evaluate_0.24.0             biovizBase_1.52.0          
## [159] BiocManager_1.30.23         tzdb_0.4.0                 
## [161] httpuv_1.6.15               RANN_2.6.1                 
## [163] tidyr_1.3.1                 purrr_1.0.2                
## [165] polyclip_1.10-7             future_1.34.0              
## [167] scattermore_1.2             xtable_1.8-4               
## [169] restfulr_0.0.15             RSpectra_0.16-2            
## [171] later_1.3.2                 viridisLite_0.4.2          
## [173] tibble_3.2.1                memoise_2.0.1              
## [175] GenomicAlignments_1.40.0    cluster_2.1.6              
## [177] globals_0.16.3
```




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
##  [1] ensembldb_2.28.1                  AnnotationFilter_1.28.0          
##  [3] GenomicFeatures_1.56.0            AnnotationDbi_1.66.0             
##  [5] Biobase_2.64.0                    BSgenome.Hsapiens.UCSC.hg38_1.4.5
##  [7] BSgenome_1.72.0                   rtracklayer_1.64.0               
##  [9] BiocIO_1.14.0                     Biostrings_2.72.1                
## [11] XVector_0.44.0                    TFBSTools_1.42.0                 
## [13] JASPAR2020_0.99.10                motifmatchr_1.26.0               
## [15] AnnotationHub_3.12.0              BiocFileCache_2.12.0             
## [17] dbplyr_2.5.0                      hdf5r_1.3.11                     
## [19] patchwork_1.3.0                   ggplot2_3.5.1                    
## [21] GenomicRanges_1.56.2              GenomeInfoDb_1.40.1              
## [23] IRanges_2.38.1                    S4Vectors_0.42.1                 
## [25] BiocGenerics_0.50.0               Seurat_5.1.0                     
## [27] SeuratObject_5.0.2                sp_2.1-4                         
## [29] Signac_1.14.0                    
## 
## loaded via a namespace (and not attached):
##   [1] ProtGenerics_1.36.0         matrixStats_1.3.0          
##   [3] spatstat.sparse_3.1-0       bitops_1.0-8               
##   [5] DirichletMultinomial_1.46.0 httr_1.4.7                 
##   [7] RColorBrewer_1.1-3          backports_1.5.0            
##   [9] tools_4.4.0                 sctransform_0.4.1          
##  [11] utf8_1.2.4                  R6_2.5.1                   
##  [13] lazyeval_0.2.2              uwot_0.2.2                 
##  [15] withr_3.0.1                 gridExtra_2.3              
##  [17] progressr_0.14.0            cli_3.6.3                  
##  [19] spatstat.explore_3.3-1      fastDummies_1.7.3          
##  [21] labeling_0.4.3              sass_0.4.9                 
##  [23] spatstat.data_3.1-2         readr_2.1.5                
##  [25] ggridges_0.5.6              pbapply_1.7-2              
##  [27] Rsamtools_2.20.0            foreign_0.8-87             
##  [29] R.utils_2.12.3              dichromat_2.0-0.1          
##  [31] parallelly_1.38.0           rstudioapi_0.16.0          
##  [33] RSQLite_2.3.7               generics_0.1.3             
##  [35] gtools_3.9.5                ica_1.0-3                  
##  [37] spatstat.random_3.3-1       dplyr_1.1.4                
##  [39] GO.db_3.19.1                Matrix_1.7-0               
##  [41] fansi_1.0.6                 abind_1.4-5                
##  [43] R.methodsS3_1.8.2           lifecycle_1.0.4            
##  [45] yaml_2.3.10                 SummarizedExperiment_1.34.0
##  [47] SparseArray_1.4.8           Rtsne_0.17                 
##  [49] grid_4.4.0                  blob_1.2.4                 
##  [51] promises_1.3.0              crayon_1.5.3               
##  [53] pwalign_1.0.0               miniUI_0.1.1.1             
##  [55] lattice_0.22-6              cowplot_1.1.3              
##  [57] annotate_1.82.0             KEGGREST_1.44.1            
##  [59] pillar_1.9.0                knitr_1.48                 
##  [61] rjson_0.2.21                future.apply_1.11.2        
##  [63] codetools_0.2-20            fastmatch_1.1-4            
##  [65] leiden_0.4.3.1              glue_1.7.0                 
##  [67] spatstat.univar_3.0-0       data.table_1.15.4          
##  [69] vctrs_0.6.5                 png_0.1-8                  
##  [71] spam_2.10-0                 gtable_0.3.5               
##  [73] poweRlaw_0.80.0             cachem_1.1.0               
##  [75] xfun_0.46                   S4Arrays_1.4.1             
##  [77] mime_0.12                   pracma_2.4.4               
##  [79] survival_3.7-0              RcppRoll_0.3.1             
##  [81] fitdistrplus_1.2-1          ROCR_1.0-11                
##  [83] nlme_3.1-165                bit64_4.0.5                
##  [85] filelock_1.0.3              RcppAnnoy_0.0.22           
##  [87] bslib_0.8.0                 irlba_2.3.5.1              
##  [89] rpart_4.1.23                KernSmooth_2.23-24         
##  [91] Hmisc_5.2-0                 colorspace_2.1-1           
##  [93] seqLogo_1.70.0              DBI_1.2.3                  
##  [95] nnet_7.3-19                 tidyselect_1.2.1           
##  [97] bit_4.0.5                   compiler_4.4.0             
##  [99] curl_5.2.1                  htmlTable_2.4.3            
## [101] DelayedArray_0.30.1         plotly_4.10.4              
## [103] checkmate_2.3.2             scales_1.3.0               
## [105] caTools_1.18.2              lmtest_0.9-40              
## [107] rappdirs_0.3.3              stringr_1.5.1              
## [109] digest_0.6.36               goftest_1.2-3              
## [111] spatstat.utils_3.0-5        rmarkdown_2.27             
## [113] base64enc_0.1-3             htmltools_0.5.8.1          
## [115] pkgconfig_2.0.3             MatrixGenerics_1.16.0      
## [117] highr_0.11                  fastmap_1.2.0              
## [119] rlang_1.1.4                 htmlwidgets_1.6.4          
## [121] UCSC.utils_1.0.0            shiny_1.9.1                
## [123] farver_2.1.2                jquerylib_0.1.4            
## [125] zoo_1.8-12                  jsonlite_1.8.8             
## [127] BiocParallel_1.38.0         R.oo_1.26.0                
## [129] VariantAnnotation_1.50.0    RCurl_1.98-1.16            
## [131] magrittr_2.0.3              Formula_1.2-5              
## [133] GenomeInfoDbData_1.2.12     dotCall64_1.1-1            
## [135] munsell_0.5.1               Rcpp_1.0.13                
## [137] reticulate_1.38.0           stringi_1.8.4              
## [139] zlibbioc_1.50.0             MASS_7.3-61                
## [141] plyr_1.8.9                  parallel_4.4.0             
## [143] listenv_0.9.1               ggrepel_0.9.5              
## [145] deldir_2.0-4                CNEr_1.40.0                
## [147] splines_4.4.0               tensor_1.5                 
## [149] hms_1.1.3                   igraph_2.0.3               
## [151] spatstat.geom_3.3-2         RcppHNSW_0.6.0             
## [153] reshape2_1.4.4              TFMPvalue_0.0.9            
## [155] BiocVersion_3.19.1          XML_3.99-0.17              
## [157] evaluate_0.24.0             biovizBase_1.52.0          
## [159] BiocManager_1.30.23         tzdb_0.4.0                 
## [161] httpuv_1.6.15               RANN_2.6.1                 
## [163] tidyr_1.3.1                 purrr_1.0.2                
## [165] polyclip_1.10-7             future_1.34.0              
## [167] scattermore_1.2             xtable_1.8-4               
## [169] restfulr_0.0.15             RSpectra_0.16-2            
## [171] later_1.3.2                 viridisLite_0.4.2          
## [173] tibble_3.2.1                memoise_2.0.1              
## [175] GenomicAlignments_1.40.0    cluster_2.1.6              
## [177] globals_0.16.3
```
