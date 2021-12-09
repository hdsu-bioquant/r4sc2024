# IRTG Course - Introduction to R for genomics

![circos](./circos.png)

Welcome to the **IRTG Course - Introduction to R** This workshop is meant for individuals with little previous knowledge R. 

The course will run over 2 days **(Wednesday, 8.12 and Thursday, 9.12)** from 10 am - 12.30 pm and 1.30 pm - 5.30 pm.


******
## Tutors

* Carl Herrmann, [Health Data Science Unit](https://www.hdsu.org/) Medical Faculty Heidelberg and BioQuant (carl.herrmann@bioquant.uni-heidelberg.de)
* Carlos Ramirez, [Health Data Science Unit](https://www.hdsu.org/) Medical Faculty Heidelberg and BioQuant (carlos.ramirez@bioquant.uni-heidelberg.de)


********

## Is this course for me?

In this two-day course, we want to give you an **introduction to working with R in simple data analysis tasks**; you will learn the basic principles of reading in a data table, doing some descriptive statistics, making nice plots.
On the second day, we will focus on a simple single-cell analysis workflow, which will guide you through the first steps of this kind of analysis!

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
[following script](./install_packages.R). Load it into RStudio, and then hit the *Run* button at the top to execute it. It should run smoothly!

Please check that you can open RStudio without error message!


## Slides

Here are the links to the slides

* Day 1 : [introduction](./irtg2021_intro.pdf)
* Day 1 : [R markdown](./irtg2021_rmarkdown.pdf)
* Day 1 : [Data types](./irtg2021_datatypes.pdf)
* Day 1 : [Cleanup](./irtg2021_cleanup.pdf)
* Day 1 : [Plots](./irtg2021_plots.pdf)
* Day 1 : [Statistical tests](./irtg2021_tests.pdf)

* Day 2 : [Introduction to single-cell analysis](https://docs.google.com/presentation/d/1DSC6gUIbO6PzrqLCt1jp-sIx1U31TvMdDGgKdhohCIY/edit?ts=60c8bafb#slide=id.gdf238a40cf_0_5)
## Practical parts

**Please document your progress in this [Google Sheet](https://docs.google.com/spreadsheets/d/1rFcWJJD-qOqeRWZvhqPEqMCt_ddtinvdTlLPl2Syomw/edit?usp=sharing)**

### Day 1: General introduction - (almost) first steps in R!                                        

On the first day, we will guide you through the first steps of working with R, from reading data to exploratory analysis and basic statistics.

* Part 0 : [Objectives](./day1/00_Objectives.md)
* Part 1 : [Rstudio](./day1/01_rstudio.md)
* Part 2 : [Dataframes](./day1/02_dataframe.md)
* Part 3 : [Data cleanup](./day1/03_cleanup.md)
* Part 4 : [Plotting](./day1/04_plotting.md)
* Part 5 : [Hypothesis tests](./day1/05_test.md)

### Day 2: a simple single-cell RNA-seq analysis workflow

On the second day, we will go through a step by step simple analysis of a small scRNA-seq dataset using the Seurat toolkit. **Don't expect to be able to carry a full scRNA-seq analysis after this!** This is meant to give you an idea of a typical workflow rather.

* Part 0 : [Initial steps](./day1/index.md)
* Part 1 : [Introduction to Seurat](./day2/01-Seurat.md)
* Part 2 : [Quality control](./day2/02-Quality_control.md)
* Part 3 : [Feature selection](./day2/03-Feature_selection.md)
* Part 4 : [Normalization](./day2/04-Normalization_and_Dimensional_Reduction.md)
* Part 5 : [Cluster_visualization](./day2/05-Cluster_visualization.md)
* Part 6 : [Differential expressoin](./day2/06-Differential_Expression.md)
* Part 7 : [Profiling cells](./day2/07-Profiling_cells.md)

You can check your results using [this form](https://forms.gle/QhSRRSLd9PpjP8NL8)

*********
## Organisation

The course will be online only! 
* Lectures will be over **Zoom** (we will send the link via email prior to the course)
* We will use **Discord channels** for the practical sessions (register [using this link](https://discord.gg/RDaYPKPQ))

