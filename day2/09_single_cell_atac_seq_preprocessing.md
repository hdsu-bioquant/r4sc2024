---
title: "Single Cell ATAC-Seq Pre-processing"
author: "Ramirez C and Herrmann C"
date: "2024-06-10"
output:
  html_document: 
    keep_md: yes
---





In this tutorial we introduce an standard pipeline of Single Cell ATAC-Seq
analysis. In order to run the steps described below please download the
files from the following 
[link](https://figshare.com/articles/dataset/Single_Cell_RNA_ATAC_Seq_integration/27331188).

In genomics, data integration generally refers to combining different types of biological data to create a more comprehensive understanding of complex biological systems. Integrating multiple data types—such as genomic, transcriptomic, and epigenomic data—allows researchers to capture diverse molecular layers that influence cellular function, uncovering relationships and regulatory mechanisms that might not be evident from a single dataset.

Specifically, in the context of ATAC-seq (chromatin accessibility data) and RNA-seq (transcriptomic data), integration seeks to connect chromatin accessibility (which regions of DNA are open and potentially active) with gene expression (which genes are being transcribed and expressed as mRNA). This integration helps to build a bridge between gene regulation and gene activity.

ATAC- and RNA integration is a useful approach than can be used for addressing: 
(i) the link between Open Chromatin Regions to gene expression, 
(ii) inferring Gene Regulatory Networks interactions 
(iii) performing cell type annotation 
(iv) analyzing dynamical trajectories.


## Pre-processing workflow


### Installations

For the pipeline, we will use the Signac R library. Please, make sure that
you have installed the following libraries by running the following commands.



``` r
install.packages('hdf5r')      ## More instructions about installing 
                                ## https://github.com/hhoeflin/hdf5r
install.packages('ggplot2')
install.packages('patchwork')
install.packages('Seurat')
install.packages('Signac') 


if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("AnnotationHub")
```


Then, we load the libraries. 



``` r
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(AnnotationHub)
```


### Loading data

The input formats for the Signac pipeline might vary. However, there are at 
least to important files which are necessary. 

i) Peaks per cell matrix. Containing the read counts per peak and per cell
as a square matrix. 

ii) Fragment file. A list of all recorded fragment reads including peaks which
are not mapped to peaks. These reads might be necessary for computing some 
QC metrics and further analysis.

Optionally, cell annotations might be also added. 

In the following lines we load the peaks counts matrix, fragment files and 
metadata. In the same way as we did for scRNA-Seq we create a Seurat
object to store the data.



``` r
## Loading the counts matrix
counts <- Read10X_h5(filename = "/Users/cbg-mbp-02/Documents/data/r2sc2024/atac_pbmc_1k_nextgem_filtered_peak_bc_matrix.h5")

## Loading metadata
metadata <- read.csv(
  file = "/Users/cbg-mbp-02/Documents/data/r2sc2024/atac_pbmc_1k_nextgem_singlecell.csv",
  header = TRUE,
  row.names = 1
)

## Loading reads-fragments
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "/Users/cbg-mbp-02/Documents/data/r2sc2024/atac_pbmc_1k_nextgem_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

## Adding the data to a Seurat object
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
```

The counts matrix can be accessed in the following path:


``` r
pbmc[['peaks']]$counts[1:5, 1:5]
```

```
## 5 x 5 sparse Matrix of class "dgCMatrix"
##                    AAACGAATCGCATAAC-1 AAACGAATCTGTGTGA-1 AAACTCGAGAGGAACA-1
## chr1-9767-10657                     .                  .                  .
## chr1-191273-192086                  .                  .                  .
## chr1-267563-268458                  .                  .                  2
## chr1-629540-630395                  .                  .                  .
## chr1-633600-634489                  .                  .                  .
##                    AAACTCGAGCCTGTAT-1 AAACTCGAGCTGAGGT-1
## chr1-9767-10657                     .                  .
## chr1-191273-192086                  .                  .
## chr1-267563-268458                  .                  .
## chr1-629540-630395                  .                  .
## chr1-633600-634489                  .                  .
```


``` r
granges(pbmc)
```

```
## GRanges object with 79973 ranges and 0 metadata columns:
##             seqnames        ranges strand
##                <Rle>     <IRanges>  <Rle>
##       [1]       chr1    9767-10657      *
##       [2]       chr1 191273-192086      *
##       [3]       chr1 267563-268458      *
##       [4]       chr1 629540-630395      *
##       [5]       chr1 633600-634489      *
##       ...        ...           ...    ...
##   [79969] KI270726.1   27104-27985      *
##   [79970] KI270713.1     3938-4823      *
##   [79971] KI270713.1   21364-22265      *
##   [79972] KI270713.1   29610-30522      *
##   [79973] KI270713.1   36938-37826      *
##   -------
##   seqinfo: 34 sequences from an unspecified genome; no seqlengths
```

In the next step we keep only the peaks which are located in conventional 
chromosomes. 



``` r
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]
```

Next, we annotate the peaks to genomic regions. The AnnotationHub R library 
provides access to genome databales like ENSEMBL and UCSC resources. 


``` r
ah <- AnnotationHub()

# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")
```

```
## AnnotationHub with 1 record
## # snapshotDate(): 2024-04-30
## # names(): AH75011
## # $dataprovider: Ensembl
## # $species: Homo sapiens
## # $rdataclass: EnsDb
## # $rdatadateadded: 2019-05-02
## # $title: Ensembl 98 EnsDb for Homo sapiens
## # $description: Gene and protein annotations for Homo sapiens based on Ensem...
## # $taxonomyid: 9606
## # $genome: GRCh38
## # $sourcetype: ensembl
## # $sourceurl: http://www.ensembl.org
## # $sourcesize: NA
## # $tags: c("98", "AHEnsDbs", "Annotation", "EnsDb", "Ensembl", "Gene",
## #   "Protein", "Transcript") 
## # retrieve record with 'object[["AH75011"]]'
```

``` r
ensdb_v98 <- ah[["AH75011"]]


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(pbmc) <- annotations
```


## Quality control

In single-cell ATAC-seq (scATAC-seq) data, quality control is essential to ensure meaningful interpretation of chromatin accessibility at the single-cell level. Several metrics are commonly used to assess data quality. First, the nucleosome banding pattern provides insight into chromatin structure by identifying whether fragments align in characteristic periodic patterns, reflecting nucleosome-bound and nucleosome-free regions. This pattern helps distinguish high-quality cells with well-defined chromatin states. TSS enrichment score measures the accumulation of fragments near transcription start sites (TSS), with high scores indicating open, accessible regions around gene promoters, characteristic of active chromatin. Additionally, the total number of fragments in peaks is evaluated to confirm that sufficient reads fall within identified accessible regions (peaks), a sign of effective enrichment in regulatory regions. Lastly, the ratio of reads in genomic blacklist regions is checked; these blacklist regions are typically artifactual and non-informative, so a low proportion of reads mapping to them suggests cleaner, high-quality data with minimal background noise. 


In the following chunk these QC metrics are calculated.



``` r
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)

# add fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# add blacklist ratio
pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)
```

During quality control it is advisable to use jointly use different
QC metrics in order to decide which cells to retain and what to filter.
In the next plot the total number of counts in cells and signal around
TSS are shown as an scatter density plot.



``` r
DensityScatter(pbmc, 
               x = 'nCount_peaks', 
               y = 'TSS.enrichment', 
               log_x = TRUE, 
               quantiles = TRUE)
```

<img src="09_single_cell_atac_seq_preprocessing_files/figure-html/scatter_density_counts_per_peak-1.png" style="display: block; margin: auto;" />


The red lines represent 5, 10, 90 and 95 percentiles. So, for example most of the cells
are in a range of ~ 15k and 45k total reads in peaks (nCount_peaks).


Now, we show a nucleosome bad patterning. A healthy good quality ATAC-Seq pattern 
should show a periodicity like the following one.



``` r
FragmentHistogram(object = pbmc, 
                  group.by = "orig.ident")
```

<img src="09_single_cell_atac_seq_preprocessing_files/figure-html/nucleosome_pattern-1.png" style="display: block; margin: auto;" />


The following violin plots show the separate distributions of the QC metrics. 



``` r
VlnPlot(
  object = pbmc,
  group.by = "orig.ident",
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
```

<img src="09_single_cell_atac_seq_preprocessing_files/figure-html/vilin_plots-1.png" style="display: block; margin: auto;" />

We can use this information to filter in/out cells. Filtering usually requires 
visual inspection in order to set up thresholds. Outliers in the distributions
might represent doublets, empty droplets, dying cells or artefactual barcodes. 


``` r
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)
pbmc
```

```
## An object of class Seurat 
## 79945 features across 878 samples within 1 assay 
## Active assay: peaks (79945 features, 0 variable features)
##  2 layers present: counts, data
```

In the Signac pipeline for scATAC-seq data analysis, RunTFIDF and RunSVD are critical functions used for data normalization and dimensionality reduction. RunTFIDF applies a Term Frequency-Inverse Document Frequency (TF-IDF) transformation to the cell-by-peak matrix, where each entry reflects how often a particular peak is accessible in a cell, adjusted for how frequently it appears across all cells. This normalization enhances the signal of cell-specific accessibility patterns and reduces the influence of universally accessible regions. Following this, RunSVD (Singular Value Decomposition) reduces the dimensionality of the TF-IDF matrix, extracting key components that capture the variation in chromatin accessibility across cells. This enables clustering and visualization of cells based on chromatin accessibility profiles, facilitating identification of cell types and states within scATAC-seq data.



``` r
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
```

In the process of dimensional reduction, peaks data is transformed into a 
different features/components space (like in Principal Component Analysis. PCA). 
In the next plot the correlation between each components of the
dimensional reduction and the sequencing depth. Normally we avoid to drive
any conclusion based in features correlated to the sequencing process since
this might be just an artefactual and for that reason we remove this effect
in the visualization step by removing this dimensions.



``` r
DepthCor(pbmc)
```

<img src="09_single_cell_atac_seq_preprocessing_files/figure-html/cor_lsi-1.png" style="display: block; margin: auto;" />


As in the Gene Expression data we can project our data now in two dimensions
in order to visualize the cells. The next



``` r
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()
```

<img src="09_single_cell_atac_seq_preprocessing_files/figure-html/umap-1.png" style="display: block; margin: auto;" />



**QUIZ 1**

<br>

<details>
<summary> Is the clustering biased by the sequencing depth? 
TIP: Evaluate this possibility by plotting a UMAP or violin plot
by cluster.
</summary>

<b>Answer:</b>
<br>
<tt> 
``` 
VlnPlot(object = pbmc, group.by = "seurat.clusters", features = 'nCount_peaks', pt.size = 0.1)
```
</tt>
<br>
</details> 

<br>

**QUIZ 2**

<br>

<details>
<summary> What happen if we don't remove the dimension correlated to 
sequencing depth? 
TIP: Add the first dimension to the `RunUMAP()` function and run 
the clustering again.
</summary>

<b>Answer:</b>
<br>
<tt> 
``` 
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 1:30)
```
</tt>
<br>
</details> 

<br>


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
##  [35] Rtsne_0.17                  purrr_1.0.2                
##  [37] biovizBase_1.52.0           RCurl_1.98-1.16            
##  [39] nnet_7.3-19                 VariantAnnotation_1.50.0   
##  [41] rappdirs_0.3.3              GenomeInfoDbData_1.2.12    
##  [43] ggrepel_0.9.5               irlba_2.3.5.1              
##  [45] listenv_0.9.1               spatstat.utils_3.0-5       
##  [47] goftest_1.2-3               RSpectra_0.16-2            
##  [49] spatstat.random_3.3-1       fitdistrplus_1.2-1         
##  [51] parallelly_1.38.0           DelayedArray_0.30.1        
##  [53] leiden_0.4.3.1              codetools_0.2-20           
##  [55] RcppRoll_0.3.1              tidyselect_1.2.1           
##  [57] UCSC.utils_1.0.0            farver_2.1.2               
##  [59] base64enc_0.1-3             matrixStats_1.3.0          
##  [61] spatstat.explore_3.3-1      GenomicAlignments_1.40.0   
##  [63] jsonlite_1.8.8              Formula_1.2-5              
##  [65] progressr_0.14.0            ggridges_0.5.6             
##  [67] survival_3.7-0              tools_4.4.0                
##  [69] ica_1.0-3                   Rcpp_1.0.13                
##  [71] glue_1.7.0                  SparseArray_1.4.8          
##  [73] gridExtra_2.3               xfun_0.46                  
##  [75] MatrixGenerics_1.16.0       dplyr_1.1.4                
##  [77] withr_3.0.1                 BiocManager_1.30.23        
##  [79] fastmap_1.2.0               fansi_1.0.6                
##  [81] digest_0.6.36               R6_2.5.1                   
##  [83] mime_0.12                   colorspace_2.1-1           
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
## [105] dotCall64_1.1-1             ProtGenerics_1.36.0        
## [107] scales_1.3.0                png_0.1-8                  
## [109] spatstat.univar_3.0-0       knitr_1.48                 
## [111] rstudioapi_0.16.0           rjson_0.2.21               
## [113] reshape2_1.4.4              checkmate_2.3.2            
## [115] nlme_3.1-165                curl_5.2.1                 
## [117] cachem_1.1.0                zoo_1.8-12                 
## [119] stringr_1.5.1               BiocVersion_3.19.1         
## [121] KernSmooth_2.23-24          parallel_4.4.0             
## [123] miniUI_0.1.1.1              foreign_0.8-87             
## [125] restfulr_0.0.15             pillar_1.9.0               
## [127] grid_4.4.0                  vctrs_0.6.5                
## [129] RANN_2.6.1                  promises_1.3.0             
## [131] xtable_1.8-4                cluster_2.1.6              
## [133] htmlTable_2.4.3             evaluate_0.24.0            
## [135] cli_3.6.3                   compiler_4.4.0             
## [137] Rsamtools_2.20.0            rlang_1.1.4                
## [139] crayon_1.5.3                future.apply_1.11.2        
## [141] labeling_0.4.3              plyr_1.8.9                 
## [143] stringi_1.8.4               viridisLite_0.4.2          
## [145] deldir_2.0-4                BiocParallel_1.38.0        
## [147] munsell_0.5.1               Biostrings_2.72.1          
## [149] lazyeval_0.2.2              spatstat.geom_3.3-2        
## [151] Matrix_1.7-0                BSgenome_1.72.0            
## [153] RcppHNSW_0.6.0              bit64_4.0.5                
## [155] future_1.34.0               KEGGREST_1.44.1            
## [157] shiny_1.9.1                 highr_0.11                 
## [159] SummarizedExperiment_1.34.0 ROCR_1.0-11                
## [161] igraph_2.0.3                memoise_2.0.1              
## [163] bslib_0.8.0                 fastmatch_1.1-4            
## [165] bit_4.0.5
```
