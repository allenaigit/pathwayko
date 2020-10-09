# PathwayKO

**PathwayKO**: An **R** package for **knock-out** pathway enrichment analysis

PathwayKO is an R package for knock-out pathway enrichment analysis, excavating the true positive KO KEGG pathways that contain and are impacted by a KO gene at the system-level. It enables ROC curve-based statistics analysis for assessing the performance of methods in terms of AUC (area under ROC curve), partial AUCs, Youden's best p-value threshold, specificity, sensitivity, accuracy, precision and recall. It is flexible to incorporate custom methods and currently integrates the state-of-the-art SPIA, ROntoTools (PE and pDIS), PADOG, GSA, GESA and SAFE. A benchmark dataset of mouse 10 KO GEO datasets is embedded for demo. The PathwayKO package provides a novel solution to the knock-out pathway enrichment analysis.

## Installation

This is an R package. So you should first have R installed on your OS. [(here's how)](https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-R-be-installed_003f).

To install this package from source code (and in this case, source code from Github), [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html) should be installed once R is in place. Also, you need to have [BiocManager](http://bioconductor.org/install/) installed so that all R package dependencies are kept consistent under BioConductor. PathwayKO was built under **R version 4.0.2** and **BioConductor version 3.11**, older versions of R and BioConductor may still work but were not tested.

Afterwards having installed R and BioConductor, launch your R session, and use the following commands to install PathwayKO ("**pathwayko**" is the package name identifier used in R):

```R
library(devtools)
# If installing from Github
devtools::install_github("allenaigit/pathwayko") 
```

If there is any missing dependencies, you will be prompted during installation and your installation may fail. Just update/install any missing packages following the prompt then try again. You may need more than one try to take care of all dependencies.

**Note**: under \*nix systems, default path for R package installation is "/usr/local/lib/R" which usually requires root access to write into it. Alternatively you can install packages under your home directory (e.g. "\~/R/"), which makes such packages only available to you and not other users.

## Usage

```R
library(pathwayko)
preprocess()
pathwayko()
```

**Note**: more details for usages are available in pathwayko/inst/docs/Users_manual.pdf

## Issue

All issues and/or suggestions should be directed to [Issues](https://github.com/allenaigit/pathwayko/issues) page, and please be sure to attach sufficient descriptions and materials to reproduce any issue you may encounter.

## License
[CC BY-NC-ND 4.0](http://creativecommons.org/licenses/by-nc-nd/4.0/)
