# RDA-forest
## association analysis of distance matrices using gradient forest approach

Citation: K. Black, J. P. Rippe, and M. Matz 2022. “Environmental Drivers of Genetic Adaptation in Florida Corals.” europepmc.org/article/ppr/ppr579436. (preprint; in revision for Molecular Ecology)

RDA-forest looks for associations between principal coordinates of a distance matrx *Y* and a matrix of potential explanatory variables *X*. **Our example here is genotype-environment association (GEA) analysis** where *Y* is the identity-by-state (IBS) matrix derived from reduced-representation sequencing (i.e., RAD and its flavors), and *X* is a matrix of environmental variables associated with samples in *Y*. Unlike most other GEA approaches, we are not trying to find specific loci associated with environment. Instead, we leverage the polygenic nature of adaptation and look for clusters, extensions, and bumps in the multivariate cloud of points (samples) defined by the genetic distance matrix that can be explained by certain environmental variables. Our underlying assumption is that adaptation to specific conditions make organisms slightly more similar to each other, genome-wide, than to their peers living in a different environment. The goal is to identify the environmental gradients that are driving this subtle genetic structure.

> NOTE on sampling design for GEA: It is necessary to maximize the number of sampled locations, with just a few (maybe just one) sample per location. I mean, instead of sampling 50 individuals from 2 locations (a typical design for Fst-based genome scanning), sample 2 individuals from 50 locations.

The method relies on 'gradientForest' package for R, which is extension of random forest to multiple response variables. It has been developed to analyze ecological community data, but is perfectly fine for genetic or gene expression data.  

## Suggested readings (if you never worked gradient forest or random forest):
- [short and sweet intro into decision trees and random forest](https://towardsdatascience.com/understanding-random-forest-58381e0602d2)
- `death00_classification_regression_trees.pdf` : great intro into regression trees, for ecological data
- `ellis_gradient_forest.pdf` : dense but very comprehensive outline of the gradient forest approach
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

