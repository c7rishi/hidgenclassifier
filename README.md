# `hidgenclassifier`: An R package implementing methodologies described in  "Mining Mutation Contexts across the Genome to Map Tumor Site of Origin" by Chakraborty et al.

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)


# Overview
`hidgenclassifier` is an R package implementing Bayesian hierarchical hidden genome classifier for cancer sites developed in "Mining Mutation Contexts across the Genome to Map Tumor Site of Origin" by Chakraborty, Martin, Guan, Begg and Shen (2021). It provides various pre-processing, fitting, and post-processing functions that collectively simplify handling of genomic datasets for use in the classifier, facilitate training of the hidden genome model,  compute predicted cancer type probabilities of new tumors based on trained models, and aid rigorous quantification of predictor effects (via odds ratios) in fitted models. The repository also includes an interactive html version of one of the figures (namely, Figure 1) displayed in the main manuscript.


# Repo Contents

- [R](./R): `R` package code.
- [data](./data): filtered subsets of TCGA whole-exome and MSK-IMPACT targeted cancer gene panel sequencing datasets used in the analysis presented in the manuscript.
- [man](./man): package manual for help in R session.
- [src](./src): C++ source codes implementing various computation-heavy back-end functions. 
- [vignettes](./vignettes): `R` vignettes for R session html help pages.
- [figures](./figures): Interactive `.html` version of Figure 1 in the main manuscript.


# System Requirements

## Hardware Requirements

The package `hidgenclassifier` can be run on a standard computer with 2 GB of RAM. For optimal performance we recommend a computer with specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core

The installation-times noted in the following are  from a computer with the recommended specs (16 GB RAM, 4 cores@3.3 GHz) and internet of speed 100 Mbps.

## Software Requirements

### OS Requirements

The GitHub development version of `hidgenclassifier` has been tested on *Linux* and *Windows* operating systems as follows:
Linux: CentOS Linux release 7.8.2003 (Core) 
Windows: Windows 10


### R version

The package `hidgenclassifier` depends on R v3.5.0 or newer. See the installation notes on the [R project homepage](https://www.r-project.org/) for details on how to install the latest version of R.

### R build tools

`hidgenclassifier` contains source C++ codes, and thus requires the necessary C++ compilers to be pre-installed. This, for example, can be ensured in Windows computers by installing [Rtools](https://cran.r-project.org/bin/windows/Rtools/). See the [CRAN manual on installing R packages](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages) for more details on installing source R packages on various platforms.



# Installation Guide

## Installing Bioconductor dependencies

`hidgenclassifier` depends on a number of Bioconductor packages. To install these dependencies run the following commands in R:
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(
  c("SomaticSignatures",
    "VariantAnnotation",
    "IRanges",
    "BSgenome.Hsapiens.UCSC.hg19")
)
```



## Installing `hidgenclassifier`


Once the Bioconductor dependencies are all installed, the easiest way to install `hidgenclassifier` from GitHub is via `R` package [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html). Run the following commands in `R` to install `devtools`, if it is not already installed:
```{r}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
```

Then install `hidgenclassifier` as follows:
```{r}
devtools::install_github("c7rishi/hidgenclassifier", build_vignettes = TRUE)
```


## Typical install time

`hidgenclassifier` depends on a number of R packages (both on CRAN and on Bioconductor), and installing them all from scratch on a Windows computer using binary packages take about 10 minutes. Install time on Linux where binary sources are not available can be substantially longer (~30 minutes). If the dependencies are all installed, `hidgenclassifier` takes about 1 minute to install. 


# Demo

After installation, a vignette illustrating an analysis of the publicly available MSK-IMPACT dataset (contained in the package) can be accessed by entering the following in the R console:
```{r}
vignette("impact_anlaysis", package = "hidgenclassifier")
```
A rendered copy of the vignette from this repository can be found [here](https://htmlpreview.github.io/?https://github.com/c7rishi/hidgenclassifier/blob/master/vignettes/impact_anlaysis.html).


# Interactive html Version of Figure 1 in the Article

An interactive html version of Figure 1 in the article is included in this repository (inside [figures](./figures)). A rendered copy of the html figure is available [here](https://htmlpreview.github.io/?https://github.com/c7rishi/hidgenclassifier/blob/master/figures/interactive-figure-1.html).   
