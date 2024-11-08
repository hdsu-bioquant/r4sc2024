# Course R for single-cell analysis (R4SC 2024)


Welcome to the **R for single-cell analysis** This workshop is meant for individuals with little previous knowledge R. 

The course will run over 2 days **(Monday, 11.11 and Tuesday, 12.11)** from 9am - 12am and 2pm - 5pm.


******
## Tutors

* Carl Herrmann, [Bioinformatics group](https://www.hdsu.org/) IPMB and BioQuant, Heidelberg University (carl.herrmann@uni-heidelberg.de)
* Carlos Ramirez, [Bioinformatics group](https://www.hdsu.org/) IPMB and BioQuant, Heidelberg University (carlos.ramirez@bioquant.uni-heidelberg.de)

********

## Schedule and practical infos

The course will take place at the [Institute for Pharmacy and Molecular Biotechnology (IPMB), Im Neuenheimer Feld 364, 69120 Heidelberg](https://maps.app.goo.gl/56pMBgnT7HnGQesT9). The computer room in in the 5th floor.

* Monday: 10h - 13h and 14h - 17h
* Tuesday: 9h - 13h and 14h - 16h


********

## Is this course for me?

In this two-day course, we want to give you an **introduction to the analysis of single-cell RNA and ATAC-seq** which will guide you through the first steps of this kind of analysis!

### IMPORTANT NOTE! 

 While we will start at a very basic level, we would **strongly encourage absolute beginners**, who have never ever worked with R, to complete a very simple online R intro course on DataCamp (["Introduction to R"](https://learn.datacamp.com/courses/free-introduction-to-r)), which will give you the very basic first concepts on what R is, and how to do some very simple operations with it.
 We will send you a link so that you can freely register to DataCamp and follow this course. 

********

## Technical pre-requisites

Every participant will work on her/his own laptop. The easiest way to work with R is using the **RStudio** interface.
Please install RStudio Desktop prior to the start of the course:

1. First install R for your operating system; you will find the correct version [on this website](https://cran.rstudio.com/) 
2. Once R is installed, you can install the RStudio Desktop version, which you find [here](https://www.rstudio.com/products/rstudio/download/#download)
3. You need to install a couple of R packages; you can download the 
[following script](./dependencies.R). Load it into RStudio, and then hit the *Run* button at the top to execute it. It should run smoothly!

Please check that you can open RStudio without error message!


## Resources

Here is a list of usefull resources if you want to perform basic analysis with R

* [Base R cheatsheet](https://github.com/rstudio/cheatsheets/blob/main/base-r.pdf)
* [Basic plots in R](http://www.sthda.com/english/wiki/r-base-graphs)
* [Seurat tutorials](https://satijalab.org/seurat/articles/get_started.html)

## Preliminary: Intro to basic R

Here are the links to the slides for the introduction to R

* Day 0 : [introduction](./r4sc_intro.pdf)
* Day 0 : [R markdown](./r4sc_markdown.pdf)
* Day 0 : [Data types](./r4sc_datatypes.pdf)
* Day 0 : [Cleanup](./r4sc_cleanup.pdf)
* Day 0 : [Plots](./r4sc_plots.pdf)
* Day 0 : [Statistical tests](./r4sc_test.pdf)

These Markdown files contain some exercises in R

* Part 0 : [Objectives](./day0/00_Objectives.md)
* Part 1 : [Rstudio](./day0/01_rstudio.md)
* Part 2 : [Dataframes](./day0/02_dataframe.md)
* Part 3 : [Data cleanup](./day0/03_cleanup.md)
* Part 4 : [Plotting](./day0/04_plotting.md)
* Part 5 : [Hypothesis tests](./day0/05_test.md)



## Practical parts

### Day 1: a simple single-cell RNA-seq analysis workflow

On the first day, we will go through a step by step simple analysis of a small scRNA-seq dataset using the Seurat toolkit. **Don't expect to be able to carry a full scRNA-seq analysis after this!** This is meant to give you an idea of a typical workflow rather.

Slides:

* Day 1 : [Introduction to single-cell analysis](https://docs.google.com/presentation/d/15N_4US7Z-1RgQmsHEXHQkbQerku00HFuafntAicYLfY/edit?usp=sharing)


Tutorials:

* Part 0 : [Initial steps](./day1/index.md)
* Part 1 : [Introduction to Seurat](./day1/01-Seurat.md)
* Part 2 : [Quality control](./day1/02-Quality_control.md)
* Part 3 : [Feature selection](./day1/03-Feature_selection.md)
* Part 4 : [Dimensional reduction](./day1/04-Normalization_and_Dimensional_Reduction.md)
* Part 5 : [Cluster_visualization - UMAP](./day1/05-Cluster_visualization.md)
* Part 6 : [Differential expression](./day1/06-Differential_Expression.md)
* Part 7 : [Profiling cells](./day1/07-Profiling_cells.md)

### Day 2: analysis of scATAC-seq and integration

* Part 1 : [Preprocessing and QC](./day2/09_single_cell_atac_seq_preprocessing.md)
* Part 2 : [Integration scRNA/scATAC](./day2/10_integration_rna_atac_seq.md)
* Part 2 : [Motif analysis](./day2/11_single_cell_atac_seq_footprinting.md)

