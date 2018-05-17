#' LUniFrac
#'
#' Longitudinal UniFrac distances for comparing changes in
#' microbial communities across 2 time points.
#'
#' Based in part on Jun Chen & Hongzhe Li (2012), GUniFrac.
#'
#' Computes difference between time points and then calculates
#' difference of these differences, resulting in a dissimilarity
#' matrix that can be used as-is or converted to a kernel
#' (similarity) matrix for hypothesis testing with MIRKAT
#' and related methods.
#'
#' @param otu.tab OTU count table, containing 2*n rows (samples) and q columns (OTUs)
#' @param tree Rooted phylogenetic tree of R class "phylo"
#' @param alpha Parameter controlling weight on abundant lineages. The same weight is used within a subjects as between subjects.
#' @param metadata Data frame with three columns: sample identifiers (must match row names of otu.tab), subject identifiers (n unique values), and time point (variable with two unique levels).
#' @return Returns a (K+1) dimensional array containing the longitudinal UniFrac dissimilarities with the K specified alpha values plus the unweighted distance. The unweighted dissimilarity matrix may be accessed by result[,,"d_UW"], and the generalized dissimilarities by result[,,"d_A"] where A is the particular choice of alpha.
#' @importFrom ape is.rooted drop.tip rtree
#' @examples
#' ## Example: simulated tree and OTU table
#' ntaxa = 5; nsubj = 10
#' library(ape)
#' sim.tree = rtree(n=ntaxa)  # simulated rooted phylogenetic tree
#' sim.otu <- matrix(runif(2*ntaxa*nsubj, 0, 100), nrow = 2*nsubj)
#' for (i in 1:nrow(sim.otu)) { sim.otu[i, ] <- sim.otu[i, ] / sum(sim.otu[i, ]) }
#' rownames(sim.otu) <- paste("Subj", rep(1:nsubj, each = 2), "_Time", rep(1:2, nsubj), sep = "")  # Row names must match samples in metadata
#' colnames(sim.otu) <- paste("t", 1:ntaxa, sep = "")   # Column names must match tree tip labels
#' sim.meta <- data.frame(sampleID = rownames(sim.otu), subjID = rep(paste("Subj", 1:nsubj, sep = ""), each = 2), time = rep(c(1,2), nsubj))
#' Ds <- LUniFrac(otu.tab = sim.otu, tree = sim.tree, alpha = c(0, 0.5, 1), metadata = sim.meta)
#' D_unweighted <- Ds[,,"d_UW"]   ## access individual array elements
#' D_gen_a0.5 <- Ds[,,"d_0.5"]
#'
#' @export
#'
LUniFrac <- function (otu.tab, tree, alpha = c(0, 0.5, 1), metadata) {
  if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")

  # Convert into proportions (per sample)
  otu.tab <- as.matrix(otu.tab)
  row.sum <- rowSums(otu.tab)
  otu.tab <- otu.tab / row.sum
  n <- nrow(otu.tab)

  # Add row names if needed
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste("comm", 1:n, sep="_")
  }

  # Check OTU name consistency
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names
         in the OTU table and the tree should match." )
  }

  # Get the subtree if tree contains more OTUs
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }

  # Reorder the otu.tab matrix if the OTU orders are different
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]

  ntip <- length(tip.label)
  nbr <- nrow(tree$edge)        # number of branches = 2*(ntip - 1)
  edge <- tree$edge             # edges entering a node (1 through (ntip - 1))
  edge2 <- edge[, 2]            # edges leaving a node (1 through 2*(ntip -1))
  br.len <- tree$edge.length    # branch lengths, corresponds to edge2

  #  Accumulate OTU proportions up the tree
  cum <- matrix(0, nbr, n)							# Branch abundance matrix
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
    node <- edge[tip.loc, 1]						# Assume the direction of edge
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }

  ### *** Up to this point, same as GUniFrac -- now we change *** ###
  ### Step 1: calculate within-subject distance data
  metadata[,3] <- as.numeric(as.factor(metadata[,3]))
  cum.t1 <- cum[, which(metadata[,3] == 1)]
  colnames(cum.t1) <- metadata[,2][which(metadata[,3] == 1)]

  cum.t2 <- cum[, which(metadata[,3] == 2)]
  colnames(cum.t2) <- metadata[,2][which(metadata[,3] == 2)]
  cum.t2 <- cum.t2[, colnames(cum.t1)]

  if (any(colnames(cum.t1) != colnames(cum.t2))) {
    stop("Same set of samples is not present at both time points!")
  }

  cum.avg <- Reduce("+", list(cum.t1, cum.t2))/2

  # Construct the returning array
  # d_UW: unweighted
  dimname3 <- c(paste("d", alpha, sep="_"), "d_UW")
  lunifracs <- array(NA, c(ncol(cum.t1), ncol(cum.t1), length(alpha) + 1),
                    dimnames=list(colnames(cum.t1), colnames(cum.t1), dimname3))
  for (i in 1:(length(alpha)+1)){
    for (j in 1:ncol(cum.t2)){
      lunifracs[j, j, i] <- 0
    }
  }

  # Calculate within-subject distances (for generalized/weighted)
  cum.diff <- matrix(0, nrow = nrow(cum), ncol = ncol(cum.t2))
  for (i in 1:ntip) {
    cum.diff[i,] <- (cum.t2[i,] - cum.t1[i,]) / (cum.t2[i,] + cum.t1[i,])
  }
  cum.diff[is.na(cum.diff)] <- 0

  # Calculate within-subject distances (for unweighted)
  cum.diff.uw <- matrix(0, nrow = nrow(cum), ncol = ncol(cum.t2))
  for (i in 1:ntip) {
    ## this is the change in V2
    cum.diff.uw[i,] <- as.numeric(cum.t2[i,] > 0) - as.numeric(cum.t1[i,] > 0)
  }

  ### Step 2: calculate distances based on within-subject summaries
  for (i in 2:ncol(cum.diff)) {
    for (j in 1:(i-1)) {
      d1 <- cum.diff[, i]
      d2 <- cum.diff[, j]
      avg1 <- cum.avg[, i]
      avg2 <- cum.avg[, j]

      ind <- which((abs(d1) + abs(d2)) != 0)
      d1 <- d1[ind]
      d2 <- d2[ind]
      avg1 <- avg1[ind]
      avg2 <- avg2[ind]
      br.len2 <- br.len[ind]

      diff <- abs(d2 - d1) / 2 

      # Generalized UniFrac distance
      for(k in 1:length(alpha)){
        w <- br.len2 * (avg1 + avg2)^alpha[k]
        lunifracs[i, j, k] = lunifracs[j, i, k] = sum(diff * w) / sum(w)
      }

      #	Unweighted UniFrac Distance
      d1 <- cum.diff.uw[, i]
      d2 <- cum.diff.uw[, j]
      diff <- abs(d2 - d1)/2

      # only branches with some change contribute
      ind <- which((abs(cum.diff.uw[, i]) + abs(cum.diff.uw[,j])) != 0)
      if (length(ind) > 0) {
        diff <- diff[ind]
        br.len2 <- br.len[ind]
        lunifracs[i, j, (k + 1)] = lunifracs[j, i, (k + 1)] = sum(br.len2*diff) / sum(br.len2)
      } else {
        lunifracs[i, j, (k + 1)] = lunifracs[j, i, (k + 1)] = 0
      }
    }
  }
  return(lunifracs)
}

