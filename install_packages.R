## Course R for single-cell (R4SC)

## you can run this script to install the required R packages
## simply hit the Run button in the RStudio interface!

## For day 1, we need the following packages

## installing RColorBrewer
install.packages('RColorBrewer')

## install tidyverse
install.packages("tidyverse")

## install pheatmap to make nice heatmaps
install.packages("pheatmap")

## install vioplot to make violin plots
install.packages("vioplot")

## For day 2, we need some additional packages

install.packages('rmarkdown')
install.packages('Seurat')
install.packages('ggplot2')
install.packages('ggrepel')
install.packages('gridExtra')

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("destiny")
BiocManager::install("Biobase")
