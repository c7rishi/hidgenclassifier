# hidgenclassifier
An R package for Bayesian hierarchical hidden genome classification. Provides functionalilites for pre-processing genomic datasets to be used in the classifier, training the hidden genome classifier, and predicting cancer classes of new tumors based on a trained model. 


# How to Install


The easiest way to install `hidgenclassifier` is via `R` package [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html). Run the following commands in `R` to install `devtools`, if it is not already installed:
```{r}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
```
Then install `hidgenclassifier` as follows:
```{r}
devtools::install_github("c7rishi/hidgenclassifier")
```

Note that `hidgenclassifier` contains source C++ codes which require compilation, and thus requires the necessary C++ compilers to be pre-installed. This, for example, can be ensured in Windows computers by installing [Rtools](https://cran.r-project.org/bin/windows/Rtools/). See the [CRAN manual on installing R packages](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages) for more details on installing source R packages on various platforms.
