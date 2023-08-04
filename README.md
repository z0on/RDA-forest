# RDA-forest
association analysis of distance matrices using gradient forest approach

Citation: K. Black, J. P. Rippe, and M. Matz 2022. “Environmental Drivers of Genetic Adaptation in Florida Corals.” europepmc.org/article/ppr/ppr579436. (preprint; in revision for Molecular Ecology)

RDA-forest looks for associations between principal coordinates of a distance matrx **Y** and a matrix of potential explanatory variables **X**. Our example here is genotype-environment association (GEA) analysis where **Y** is the identity-by-state (IBS) matrix derived from reduced-representation sequencing (i.e., RAD and its flavors), and **X** is a matrix of environmental variables associated with samples in **Y**. Unlike most other GEA approaches, we are not trying to find specific loci associated with environment. Instead, we leverage the polygenic nature of adaptation and look for clusters, extensions, and bumps in the multivariate cloud of points (samples) defined by the genetic distance matrix that can be explained by certain environmental variables. In other words, our underlying assumption is that adaptation to specific conditions make organisms slightly more similar to each other, genome-wide, than to their peers living in a different environment. The goal is to identify the environmental gradients that are driving this subtle genetic structure.

> NOTE on sampling design for GEA: It is necessary to maximize the number of sampled locations, with just a few (maybe just one) sample per location. I mean, instead of sampling 50 individuals from 2 locations (a typical design for Fst-based genome scanning), sample 2 individuals from 50 locations.


