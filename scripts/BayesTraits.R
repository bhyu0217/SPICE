library(posterior)
library(tidybayes)
library(coda)
library(btw)
library(tidytree)
library(ape)
library(tidyverse)
library(patchwork)
library(parallel)


AddNodeTagsForBayesTraits <- function(tree, outFile = NULL) {
  # This function creates BayesTraits "AddTag" and "AddNode" lines for each internal node in a tree,
  # along with the tip labels of that node's descendant clade.
  #
  # INPUT:
  #   tree    : A phylo object (ape format)
  #   outFile : (optional) file path to save the commands. If NULL, no file is written.
  #
  # OUTPUT:
  #   A character vector of lines suitable for a BayesTraits command file.
  #   Each pair of lines corresponds to a single node:
  #     1) "AddTag T<NODE_ID> [list of tip labels]"
  #     2) "AddNode <NODE_ID> T<NODE_ID>"
  #
  # EXPLANATION:
  #   - subtrees(tree) creates a list of subtrees, one subtree per internal node.
  #   - We loop through each node, generate the "AddTag" line that includes the node's descendant tips,
  #     and then create the "AddNode" line referencing that tag.
  #   - Finally, we format and optionally write the lines to outFile.

  # -----------------------------------------------------
  # 1) Generate all subtrees, one subtree per internal node
  # -----------------------------------------------------
  subtr <- subtrees(tree)
  # subtr is a list where each element is a subtree (phylo) corresponding to a node in 'tree'.

  # -----------------------------------------------------
  # 2) Build the matrix (MAT) that will store commands
  #    We have 2 rows per node (AddTag, AddNode),
  #    and columns for "command" + tip labels
  # -----------------------------------------------------
  Nnode <- tree$Nnode
  MAT <- matrix(nrow = (Nnode * 2), ncol = (2 + length(tree$tip.label)))

  # -----------------------------------------------------
  # 3) Write the AddTag / AddNode lines to the matrix
  #    Each node has 2 lines:
  #       Row1: "AddTag T<NODE_ID>"
  #       Row2: "AddNode <NODE_ID> T<NODE_ID>"
  # -----------------------------------------------------
  for (i in seq_len(Nnode)) {
    MAT[(i*2 - 1), 1] <- paste0("AddTag ", "T", i)
    MAT[(i*2), 1]     <- paste0("AddNode ", i, " T", i)
  }

  # -----------------------------------------------------
  # 4) Loop over each node's subtree tips, placing them on the "AddTag" line
  #    The AddTag line is (i*2 - 1) in the matrix. We shift by +1 in columns
  #    so the first column is the command, subsequent columns hold tip labels.
  # -----------------------------------------------------
  for (j in seq_len(Nnode)) {
    tips_j <- subtr[[j]]$tip.label
    for (i in seq_along(tips_j)) {
      MAT[(j*2 - 1), (i + 1)] <- tips_j[i]
    }
  }

  # -----------------------------------------------------
  # 5) Convert each row of MAT into a single tab-delimited string
  #    We also remove NAs for columns not used.
  # -----------------------------------------------------
  formatted_MAT <- apply(MAT, 1, function(row) {
    # Remove NA entries
    row <- row[!is.na(row)]
    # Join by tab
    paste(row, collapse = "\t")
  })

  # -----------------------------------------------------
  # 6) If outFile is provided, write the lines to that file
  # -----------------------------------------------------
  if (!is.null(outFile)) {
    writeLines(formatted_MAT, outFile)
  }

  # -----------------------------------------------------
  # 7) Return the character vector of commands
  # -----------------------------------------------------
  return(formatted_MAT)
}

# This function writes a BayesTraits command file for MCMC with MultiState model.
# Then it calls BayesTraits from R via system().
bayestraits_mcmc_cmd <- function(treefile, traitfile, outputdir, outprefix, logfile, iters=1000000, burnin=200000, sample=1000) {
  cmd_lines <- c(
    "1",	# Use discrete MultiState model
    "2",	# Bayesian MCMC analysis
    paste("Iterations", format(iters, scientific = FALSE)),
    paste("Burnin", format(burnin, scientific = FALSE)),
    "HyperPriorAll exp 0 10",	# Adjust to your data
    "Stones 10 1000",	# if you want stepping stones for marginal likelihood
    paste("Sample", format(sample, scientific = FALSE)),
    AddNodeTagsForBayesTraits(tree, paste0(outputdir, outprefix, "_AddNodes.txt")), 
    paste("LogFile", logfile),
    "Run"
  )
  return(cmd_lines)
}
