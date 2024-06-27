# Course R for single-cell analysis (R4SC 2024)


Welcome to the **R for single-cell analysis** This workshop is meant for individuals with little previous knowledge R. 

The course will run over 2 days **(Thursday, 27.06 and Friday, 28.06)** from 9am - 12am and 2pm - 5pm.


******
## Tutors

* Carl Herrmann, [Bioinformatics group](https://www.hdsu.org/) IPMB and BioQuant, Heidelberg University (carl.herrmann@uni-heidelberg.de)
* Carlos Ramirez, [Bioinformatics group](https://www.hdsu.org/) IPMB and BioQuant, Heidelberg University (carlos.ramirez@bioquant.uni-heidelberg.de)


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


## Resources

Here is a list of usefull resources if you want to perform basic analysis with R

* [Base R cheatsheet](https://github.com/rstudio/cheatsheets/blob/main/base-r.pdf)
* [Basic plots in R](http://www.sthda.com/english/wiki/r-base-graphs)
* [Seurat tutorials](https://satijalab.org/seurat/articles/get_started.html)

## Slides

Here are the links to the slides

* Day 1 : [introduction](./r4sc_intro.pdf)
* Day 1 : [R markdown](./r4sc_markdown.pdf)
* Day 1 : [Data types](./r4sc_datatypes.pdf)
* Day 1 : [Cleanup](./r4sc_cleanup.pdf)
* Day 1 : [Plots](./r4sc_plots.pdf)
* Day 1 : [Statistical tests](./r4sc_test.pdf)

* Day 2 : [Introduction to single-cell analysis](https://docs.google.com/presentation/d/15N_4US7Z-1RgQmsHEXHQkbQerku00HFuafntAicYLfY/edit?usp=sharing)
## Practical parts

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

* Part 0 : [Initial steps](./day2/index.md)
* Part 1 : [Introduction to Seurat](./day2/01-Seurat.md)
* Part 2 : [Quality control](./day2/02-Quality_control.md)
* Part 3 : [Feature selection](./day2/03-Feature_selection.md)
* Part 4 : [Dimensional reduction](./day2/04-Normalization_and_Dimensional_Reduction.md)
* Part 5 : [Cluster_visualization - UMAP](./day2/05-Cluster_visualization.md)
* Part 6 : [Differential expression](./day2/06-Differential_Expression.md)
* Part 7 : [Profiling cells](./day2/07-Profiling_cells.md)
* Part 8 : [Diffusion analysis](./day2/08-Intro_to_pseudotime_analysis.md)

