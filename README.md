# RDA-forest
## application of gradient forest for genotype-environment association analysis of genetic distance matrices

Citation: K. Black, J. P. Rippe, and M. Matz 2022. “Environmental Drivers of Genetic Adaptation in Florida Corals.” europepmc.org/article/ppr/ppr579436. (preprint; in revision for Molecular Ecology)

RDA-forest is designed to identify associations between principal coordinates of a distance matrx *Y* and a matrix of potential explanatory variables *X*. While the metod is applicable to any problem that can be formulated in that way, in the example here *Y* is the identity-by-state (IBS) matrix derived from reduced-representation sequencing (i.e., RAD and its flavors), and *X* is a matrix of environmental variables measured for samples in *Y*. Unlike most other GEA approaches, we are not trying to find specific loci associated with environment. Instead, we leverage the polygenic nature of adaptation and look for clusters, extensions, and bumps in the multivariate cloud of points (samples) defined by the genetic distance matrix that can be explained by certain environmental variables. The most fundamental underlying assumption is that adaptation to specific conditions make organisms slightly more similar to each other, genome-wide, than to their peers living in a different environment. The goal is to identify the environmental gradients that drive this subtle genetic structure.

> NOTE on sampling design for GEA: It is necessary to maximize the number of sampled locations, with just a few (or even just one) sample per location. I mean, instead of sampling 50 individuals from 2 locations (a typical design for Fst-based genome scanning), sample 2 individuals from 50 locations.

## Overview ##

The method relies on `gradientForest` package in R, which is extension of the random forest approach to multiple response variables. Its main advantags over "previous generation" regression-based GEA methods are:
- identifies linear and all sorts of non-linear relationships
- automatically accounts for all possible interactions between predictors
- uses proper cross-validation to compute importance of predictors
The idea here is to apply `gradientForest` to principal coordinates of a distance matrix.
There are two novel ideas implemented within the main `RDAforest` function here:
- It uses "spatial bootstrap" to ensure that detected associations are not due to the chance configuration of the principal coordinates. The analysis is run N times (N = 25 is reasonable) each time rotating the cloud of points (defined by the distance matrix) in a random direction and re-forming the principal coordinates to be orthogonal to that direction. The rotation procedure is essentially redundancy analysis (RDA) using random predictor, hence the name "RDA forest".
- To remove predictors that are not by theselves important but are correlated with important ones, we run not one but two gradient forest analyses on each spatial bootstrap replicate, with different values of `mtry`. `mtry` is the number of randomly chosen predictors to find the next split in the tree. Simulations show (Strobl et al 2008) that with higher `mtry` setting, predictors that are "standing-in" for the truly important ones decrease in importance because there is a higher chance that the actual important predictor will also be picked. We compare importances of each variable at two `mtry` values, 0.25 **p* and 0.375**p*, where *p* is the total number of predictors (`mtry` is *p*/3 by default). We then discard variables that show diminishing importance at higher `mtry` in more than half of all spatial bootstrap replicates.
There are only two R functions in this repository: 


## Suggested readings:
- [short and sweet intro into decision trees and random forest](https://towardsdatascience.com/understanding-random-forest-58381e0602d2)
- `death00_classification_regression_trees.pdf` : great intro into regression trees, for ecological data
- `ellis_gradient_forest.pdf` : dense but very comprehensive outline of the gradient forest approach
-  vignette for the `gradientForest` package
- `strobl_conditional_permutation.pdf` : deals with the problem of correlated predictors, contains simulations which are the basis for predictor selection method used here.

## Installing *gradientForest* package in R  

Installing `gradientForest` is more involved than a typical R package since it must be compiled from source. 

First, install **devtools**. 
```R
install.packages("devtools")
```
This may require additional installations outside R. Hopefully they will happen automatically, if not, see [here](https://www.r-project.org/nosvn/pandoc/devtools.html).
Then:
```R
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
```
If this one fails, chances are, you need to install **gfortran** first, FOR YOUR SYSTEM from here:
https://gcc.gnu.org/wiki/GFortranBinaries or, for Mac, https://github.com/fxcoudert/gfortran-for-macOS/releases. If there is no precompiled **gfortran** for your combination of processor and OS, choose the one for the closest OS for your processor.
On a Mac you might also need to point your Rstudio compiler to the location **gfortran** is installed at, by creating/modifying the file *~/.R/Makevars*. The following spell in **Terminal** should work:
```sh
cd
mkdir .R
echo "FC = /usr/local/bin/gfortran
F77 = /usr/local/gfortran
FLIBS = -L/usr/local/gfortran/lib" >> ~/.R/Makevars
```
To check if everything was intalled correctly, do this in Rstudio and see if the package is loaded without errors.
```R
library(gradientForest)
```

