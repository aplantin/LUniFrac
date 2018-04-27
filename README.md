## LUniFrac
Calculates a Longitudinal UniFrac dissimilarity matrix. 

LUniFrac is a two-stage UniFrac-type metric that calculates (normalized) changes in the microbiota between two time points, then compares these changes between individuals. For unweighted LUniFrac, the changes across time are defined specifically as changes in taxon presence or absence, i.e., acquisition or loss between time points. For generalized LUniFrac, the changes are defined via normalized changes in abundance (similar to fold-changes in abundance, but the denominator is average proportion, not proportion at time 1).  

These dissimilarities may be used in any downstream analysis requiring a beta-diversity summary. 

## Installation instructions: 

    library(devtools)
    install_github("aplantin/LUniFrac")
    library(LUniFrac) 
