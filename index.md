# Course R for single-cell analysis (R4SC 2026)


Welcome to the **R for single-cell analysis** This workshop is meant for individuals with little previous knowledge R. 

The course will run over 2 days **(Thursday, 11.06 and Friday, 12.06)** 


******
## Tutors

* Carl Herrmann, [Bioinformatics group](https://www.hdsu.org/) IPMB and BioQuant, Heidelberg University (carl.herrmann@uni-heidelberg.de)
* Carlos Ramirez, [Bioinformatics group](https://www.hdsu.org/) IPMB and BioQuant, Heidelberg University (carlos.ramirez@bioquant.uni-heidelberg.de)

********

## Schedule and practical infos

The course will take place at the Medizinisches Lehrgebäude. 

* Thursday: 10h - 13h and 14h - 17h
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


**IF YOU HAVE FINISHED DAY1 and DAY2 early: You can redo the analysis steps of Day 1 using this 10x dataset from glioblastoma (brain tumor), which contains 2000 cells. Download the files [here](https://drive.google.com/drive/folders/1PwWkTIBjewgZpo5nHJK4W2dQsvYuylGg?usp=sharing), and rerun all the steps of day1!**


On the first day, we will go through a step by step simple analysis of a small scRNA-seq dataset using the Seurat toolkit. **Don't expect to be able to carry a full scRNA-seq analysis after this!** This is meant to give you an idea of a typical workflow rather.

Slides:

* Day 1 : [Introduction to single-cell analysis](https://docs.google.com/presentation/d/1NI2XDHyUHCa7SV7ktgpwW41V9INUWEr3gBfsIQIaKI4/edit?usp=sharing)


Tutorials:

* Part 0 : [Initial steps](./day1/index.md)
* Part 1 : [Introduction to Seurat](./day1/01-Seurat.md)
* Part 2 : [Quality control](./day1/02-Quality_control.md)
* Part 3 : [Feature selection](./day1/03-Feature_selection.md)
* Part 4 : [Dimensional reduction](./day1/04-Dimensional_reduction.md)
* Part 5 : [Cluster_visualization - UMAP](./day1/05-Cluster_visualization.md)
* Part 6 : [Differential expression](./day1/06-Differential_Expression.md)
* Part 7 : [Profiling cells](./day1/07-Profiling_cells.md)
* Part 8 : [Pseudotime analysis](./day1/08_Pseudotime_analysis.md)


### Day 2: analysis of scATAC-seq and integration

* Part 1 : [Preprocessing and QC](./day2/09_single_cell_atac_seq_preprocessing.md)
* Part 2 : [Integration scRNA/scATAC and label transfer](./day2/10_integration.md)
* Part 3 : [Motif analysis](./day2/11_single_cell_atac_seq_footprinting.md)

