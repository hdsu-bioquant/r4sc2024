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
## Slides

Here are the links to the slides

* Day 1 : [introduction](./irtg2021_intro.pdf)
* Day 1 : [R markdown](./irtg2021_rmarkdown.pdf)
* Day 1 : [Data types](./irtg2021_datatypes.pdf)
* Day 1 : [Statistical tests](./irtg2021_tests.pdf)

* Day 2 : [Introduction to single-cell analysis](https://docs.google.com/presentation/d/1DSC6gUIbO6PzrqLCt1jp-sIx1U31TvMdDGgKdhohCIY/edit?ts=60c8bafb#slide=id.gdf238a40cf_0_5)
## Practical parts


**You need to install a couple of R packages; you can download the 
[following script](./install_packages.R). Load it into RStudio, and then hit the *Run* button at the top to execute it. It should run smoothly!**

**Please document your progress in this [Google Sheet](https://docs.google.com/spreadsheets/d/1rFcWJJD-qOqeRWZvhqPEqMCt_ddtinvdTlLPl2Syomw/edit?usp=sharing)**

### Day 1: General introduction - (almost) first steps in R!                                        


On the first day, we will guide you through the first steps of working with R, from reading data to exploratory analysis and basic statistics.

0. [Goals of Day 1](./day1/00_Objectives.md)
1. [Getting started with RStudio](./day1/01_rstudio.md)
2. [Reading in a data table](./day1/02_dataframe.md)
3. [Cleaning the dataset](./day1/03_cleanup.md)
4. [Making plots](./day1/04_plotting.md)
5. [Statistical tests](./day1/05_test.md)
                                                   

### Day 2: a simple single-cell RNA-seq analysis workflow

On the second day, we will go through a step by step simple analysis of a small scRNA-seq dataset using the Seurat toolkit. **Don't expect to be able to carry a full scRNA-seq analysis after this!** This is meant to give you an idea of a typical workflow rather.

0. [Intro on scRNA-seq analysis](./day2/00.md)
1. [Creating a SEURAT object](./day2/01.md)
2. [Basic QC and normalization ](./day2/02.md)
3. [Feature selection](./day2/03.md)
4. [Normalization and Dimensionality reduction](./day2/04.md)
5. [Cluster visualization](./day2/05.md)
6. [Differential gene expression analysis](./day2/06.md)
7. [Profiling cells](./day2/07.md)




*********
## Organisation

The course will be online only! 
* Lectures will be over **Zoom** (we will send the link via email prior to the course)
* We will use **Discord channels** for the practical sessions (register [using this link](https://discord.gg/gPXJDukGfQ))

**********
## Technical pre-requisites

Every participant will work on her/his own laptop. The easiest way to work with R is using the **RStudio** interface.
Please install RStudio Desktop prior to the start of the course:

1. first install R for your operating system; you will find the correct version [on this website](https://cran.rstudio.com/) 
2. once R is installed, you can isntall the RStudio Desktop version, which you find [here](https://www.rstudio.com/products/rstudio/download/#download)

Please check that you can open RStudio without error message!