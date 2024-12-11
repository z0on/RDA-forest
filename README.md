# RDA-forest
## Using random forest to analyze shapes of multivariate datasets

> Citation: Matz, Mikhail V., and Kristina L. Black. 2024. “RDAforest: Identifying Environmental Drivers of Polygenic Adaptation.” bioRxiv. https://doi.org/10.1101/2024.10.21.619525.

RDA forest is a way to detect associations between principal components of a response matrx *Y* and a matrix of potential explanatory variables *X*. Essentially, the method looks for clusters, extensions, and bumps in the multivariate cloud of data points that can be explained by any combination of variables in *X* (including all sorts of non-linear dependencies and multi-way interactions). We call this approach RDAforest, to reflect the fact that it has the same purpose as redundancy analysis (RDA) - to find associations between highly dimensional data and multiple predictor variables - except RDA-forest relies on more versatile RF instead of linear regressions. 

### Application to Genotype-Environment Association (GEA) analysis

In the example here *Y* is the matrix of genetic distances between individuals and *X* is a matrix of environmental variables measured for all individuals in *Y*. Unlike most other GEA methods, we are not trying to find specific loci associated with environment, but want to identify the environmental variables that drive polygenic local adaptation. Polygenic adaptation [alters the genome-wide covariance structure](https://doi.org/10.1111/2041-210X.13722), which means that similarly adapted organisms become slighly more similar to each other genetically than to their peers adapting to a different environment. Our method aims to identify environmental parameters driving this subtle pattern of genetic similarity, captured by the leading principal components (PCs) of the genetic distance matrix. We also provide solutions to visualize the risk of future genetic maladaptation ("genetic offset"), plan assisted gene flow interventions, and identify most suitable locations for a given individual (based on its genotype).

### Overview

The method relies on `gradientForest` package in R, which is extension of the random forest approach to multiple response variables. Its main advantags over regression- and model-based methods are:
- it identifies all sorts of non-linear and non-monotonous relationships as well as linear ones;
- it automatically accounts for all possible interactions between predictors;
- it handles correlated predistors properly, using ![conditional permutation](strobl18_conditional_permutation_mtry.pdf) to determine their importance;
- it uses cross-validation to compute importance of predictors, so what it reports is the actual predictive power of the model for a completely new set of data.

In addition, there are two novel ideas in our RDA-forest method:
- **Jackknifing**: We use "ordination jackknife" procedure, which rebuilds the ordination multiple times based on a subset (default 0.8 of total) of all samples and reruns the analysis. This models the uncertainty of PCs determination.
- **Mtry-based variable selection**: To detect and discard predictors that are not by themselves important but are correlated with an important one, we use the `mtry`criterion. `mtry` is the number of randomly chosen predictors to find the next split in the tree. With higher `mtry` there is a higher chance that the actual driver is chosen together with the non-influential correlated variable and is then used for the split. As a result, the correlated variable is used less often, which drives its importance down (Strobl et al 2008). So, we fit two models with different `mtry` settings to each ordination jackknife replicate. Predictors consistently showing diminished raw importance at the higher `mtry` setting are then discarded.

### Installation 

The RDA-forest functions come in the form of an R package, `RDAforest_2.6.4.tar.gz`. To install it, run this in Rstudio
```R
install.packages("/path/to/downloaded/file/RDAforest_2.6.4.tar.gz")
library(RDAforest)
```

The package depends on `vegan`, `dplyr`, and `gradientForest`. Installing `gradientForest` is more involved than a typical R package since it must be compiled from source. 

First, install `devtools`. 
```R
install.packages("devtools")
```
This may require additional installations outside R. Hopefully they will happen automatically, if not, see [here](https://www.r-project.org/nosvn/pandoc/devtools.html).
Then:
```R
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
```
If this one fails, chances are, you need to install `gfortran` first, FOR YOUR SYSTEM from here:
https://gcc.gnu.org/wiki/GFortranBinaries or, for Mac, https://github.com/fxcoudert/gfortran-for-macOS/releases. If there is no precompiled `gfortran` for your combination of processor and OS, choose the one for the closest OS for your processor.
On a Mac you might also need to point your Rstudio compiler to the location `gfortran` is installed at, by creating/modifying the file *~/.R/Makevars*. The following spell in `Terminal` should work:
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
### RDAforest functions
All functions have documentation accessible as usual by asking `?functionName` in R, for example to see what are the necessary arguments and what does the function return.

#### Main functions:
- **`ordinationJackknife`** : Runs `gradientForest` on *nreps* jackknife replicates, each time rebuilding the ordination withholding a fraction of datapoints (default 0.2). Forms two kinds of predictions for the supplied *newX* data, averaging over replicates - turnover curves from gradient forest and straight-up random forest predictions. These can be used for plotting adaptive neighborhoods (`plot_adaptation` function). Makes sure no predictions are made beyond the numeric range of predictors used to fit the model, but can allow a bit of extension specified as fraction of the original range (option *extra*).
- **`importance_RDAforest`** : Recalculates R2-based importances stored in the `gradientForest` model into proportion of variation attributable to each predictor (takes into account eigenvalues of the ordination that was analyzed).
- **`mtrySelJack`** : Performs variable selection based on *mtry* criterion: variables that are not important by themselves, only correlated with the actually important ones, begin to lose importance at higher *mtry* setting (Strobl et al 2018). The function runs *nreps*  `ordinationJackknife` replicates, fitting two `gradientForest` models with different *mtry* settings. It then selects variables that do not decrease in importance at higher *mtry*. Can be made more allowing (i.e. retain more predictors) by using lower *mtry* values (add options `mintry=3, maxtry=6`) or lower *prop.positive.cutoff*. Also discards variables whose importance is less than *importance.cutoff* fraction of the best predictor (default 0.1). Default importance scaling is "GF", like in `gradientForest` - average R2 across all predicted variables.
- **`plot_adaptation`** : Plots first two or three PCs of the supplied matrix (*[result of `ordinationJackknife`]$predictions.direct*), then plots geographical map of adaptation colored according to these PCs. Can color points by their continuous values along PCs, or by cluster they fall into (with *nclust* argument). Clustering is done using function `cluster::clara`. Can use turnover curves (*[result of `ordinationJackknife`]$predictions.turnover*) for clustering (option *clustering.guide*) and merge clusters if they are too similar accoring to *[result of `ordinationJackknife`]$predictions.direct*. In simulations this generates less noisy clustering than clustering random forest predictions straight up. Cluster merging is controlled by *cluster.merge.threshold* parameter, which is specified as the fraction of the maximal observed between-cluster distance (default 0.333).
- **`plot_nice_map`** : plots a map of adaptive neighborhoods (same as `plot_adaptation`) in Universal Transverse Merkator (UTM) coordinates. Can plot two series of points: the colored raster of adaptive neighborhoods, and (optionally) second series of points assumed to be locations of actual samples used for modeling. 
- **`gen_offset_oj`** : computes genetic offset (measure of anticipated risk of future maladaptation) based on two series of predictions generated by `ordinationJackknife`, one for present-day and another for the future.
- **`env_mismatch`** : calculates maladaptation across the landscape for a given set of genetic PCs. Can be used for assisted gene flow planning (find locations that best match required future genetic PCs now), or for finding most suitable environment for an individual based on its genotype.

#### Minor/accessory functions:
- **`dummify`** : Turns a dataframe containing numerical and categorical predictors into fully numerical.
- **`sum_up_importances`** : Sums up importances of original factors that were dummified using `dummify`.
- **`plot_gf_turnovers`** : Plots turnover curves for a `gradientForest` model. Wrapper for `plot.gradientForest` *plot.type="Cumulative.Importance"*
- **`plot_turnover`** : Plots turnover curve for the specific *X* predictor, based on *[result of `ordinationJackknife`]$predictions.turnover*. 
- **`makeGF_simple`** : A simple wrapper for the `gradientForest` function, uses straight-up response matrix *Y*.
- **`makeGF`** : Runs gradient forest analysis on an ordination object made by `vegan::capscale` or `vegan::rda`.
- **`predict_rf`** : (To be used directly on *Y* data, without ordination or jackknifing) Predicts *Y* values based on `extendedForest::randomForest`. These values can be used for plotting adaptation with `plot_adaptation` instead of *[result of `ordinationJackknife`]$predictions.direct*. Handles modeled and predicted ranges like `ordinationJackknife`. 
- **`predict_gf`** : (To be used directly on *Y* data, without ordination or jackknifing) Generates turnover curves with `gradientForest`. These should **not** be used plotting adaptation but are good for clustering points into adaptive neighborhoods in *plot_adaptation* with option *clustering.guide* (instead of *[result of `ordinationJackknife`]$predictions.turnover*). Handles modeled and predicted ranges like `ordinationJackknife`.
- **`gcd.dist`** : calculates great circle distances based on longitude and latitude, returns adjusted coordinates and the distance matrix.
- **`latlon2UTM`**, **`epsg.maker`**, **`bw_choose`** ,**`gen_offset`**,**`adapt_scale`** - various accessory functions.

### Example analysis: [North American Wolves](https://rpubs.com/cmonstr/1242635)
Download Rmarkdown script `RDAforest_wolves.Rmd` and the dataset `wolf_v4.RData` to replicate this.

### Suggested readings
- [short and sweet intro into decision trees and random forest](https://towardsdatascience.com/understanding-random-forest-58381e0602d2)
- `death00_classification_regression_trees.pdf` : great intro into regression trees, for ecological data
- `ellis_gradient_forest.pdf` : dense but very comprehensive outline of the gradient forest approach
-  [vignette for the `gradientForest` package](https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf)
- `strobl_conditional_permutation_mtry.pdf` : deals with the problem of correlated predictors, contains simulations which are the basis for predictor selection method used here.
- `bay18_songbird_science.pdf` : example application of gradient forest to predict bird adaptation and vulnerability to climate change

