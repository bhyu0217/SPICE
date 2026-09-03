library(posterior)
library(tidybayes)
library(coda)
library(btw)
library(tidytree)
library(ape)
library(tidyverse)
library(patchwork)
library(parallel)
source("./scripts/BayesTraits.R")


run_bayestraits_mcmc <- function(treefile, traitfile, outputdir, outprefix, iters=1000000, burnin=200000, sample=1000) {
  if (!dir.exists(outputdir)) {
	  dir.create(outputdir, recursive=TRUE)
	  message("Folder created: ", outputdir)
  }

  # 1) Create command file
  cmdfile <- paste0(outputdir, "/", outprefix, "_cmd.txt")
  logfile <- paste0(outputdir, "/", outprefix)
  cmd_lines <- bayestraits_mcmc_cmd(treefile, traitfile, outputdir, outprefix, logfile, iters, burnin, sample)
  writeLines(cmd_lines, cmdfile)

  # 2) Run BayesTraits
  BT <- "/c4/home/bhyu0217/build/BayesTraitsV4.1.3-Linux/BayesTraitsV4"
  cmd <- paste(BT, treefile, traitfile, "<", cmdfile)
  cat(cmd)
  system(cmd)
}

args <- commandArgs(trailingOnly = TRUE)
clone_dir <- args[1]
clone_id <- args[2]
trait_file <- args[3]
trait_order <- args[4]
mcmc_chains <- as.integer(args[5])
iterations <- as.integer(args[6])
burnin <- as.integer(args[7])
log_sample_period <- as.integer(args[8])

# Load in subclone tree
#clone_path <- paste0(clone_dir, "Clone/", clone_id, "/")
clone_path <- paste0(clone_dir, "/")
treefile <- paste0(clone_path, clone_id, ".reformatted.nex")
tree <- read.nexus(treefile)

# Tests to see whether the tree is binary, should be true
is.binary(tree)

# Read in trait data
trait <- read.table(trait_file, header=FALSE)
colnames(trait) <- c("CellBarcode","Label")
trait <- trait[trait$CellBarcode %in% tree$tip.label, ]
trait_ordered <- trait[match(tree$tip.label, trait$CellBarcode), ]
rownames(trait_ordered) <- NULL

# Is the tree tip order the same as the trait order
sum(tree$tip.label == trait[,1]) == nrow(trait)

desired_order <- strsplit(trait_order, ",")[[1]]
trait_ordered$Label <- factor(trait_ordered$Label, levels = desired_order)
trait_ordered$Label <- as.numeric(trait_ordered$Label) - 1

traitfile <- paste0(clone_path, clone_id, ".trait.tsv")
write.table(trait_ordered, traitfile,
	    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Run BayesTraits
mcmc_ids <- paste0("MCMC", seq_len(mcmc_chains))
results <- mclapply(mcmc_ids, function(id) {
  run_bayestraits_mcmc(treefile, traitfile, paste0(clone_path, id, "/"), id, iterations, burnin, log_sample_period)},
  mc.cores = 3
  )
