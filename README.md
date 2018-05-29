## LUniFrac
Calculates a Longitudinal UniFrac dissimilarity matrix. 

LUniFrac is a two-stage UniFrac-type metric that calculates (normalized) changes in the microbiota between two time points, then compares these changes between individuals. For unweighted LUniFrac, the changes across time are defined specifically as changes in taxon presence or absence, i.e., acquisition or loss between time points. For generalized LUniFrac, the changes are defined via normalized changes in abundance (similar to fold-changes in abundance, but the denominator is average proportion, not proportion at time 1).  

These dissimilarities may be used in any downstream analysis requiring a beta-diversity summary. 

## Installation instructions: 

    library(devtools)
    install_github("aplantin/LUniFrac")
    library(LUniFrac) 

## Sample usage with simulated data  

    library(ape)
    ntaxa = 5; nsubj = 10
    sim.tree = rtree(n=ntaxa)  # simulated rooted phylogenetic tree
    sim.otu <- matrix(runif(2*ntaxa*nsubj, 0, 100), nrow = 2*nsubj)
    for (i in 1:nrow(sim.otu)) { sim.otu[i, ] <- sim.otu[i, ] / sum(sim.otu[i, ]) }
    
    # Row names must match samples in metadata
    rownames(sim.otu) <- paste("Subj", rep(1:nsubj, each = 2), "_Time", rep(1:2, nsubj), sep = "")  
    
    # Column names must match tree tip labels
    colnames(sim.otu) <- paste("t", 1:ntaxa, sep = "")   
    
    # Metadata includes sample ID (unique for each sample), 
    #    subject ID (unique for each subject, so observed twice in a dataset), 
    #    and time point (must take exactly two unique values). 
    sim.meta <- data.frame(sampleID = rownames(sim.otu), 
                           subjID = rep(paste("Subj", 1:nsubj, sep = ""), each = 2), 
                           time = rep(c(1,2), nsubj))
                           
    # Calculate LUniFrac dissimilarity 
    Ds <- LUniFrac(otu.tab = sim.otu, tree = sim.tree, alpha = c(0, 0.5, 1), metadata = sim.meta)

    # Each distance/dissimilarity matrix is a slice of the returned array
    D_unweighted <- Ds[,,"d_UW"]  
    D_gen_a0.5 <- Ds[,,"d_0.5"]
    
