library(posterior)
library(tidybayes)
library(coda)
library(btw)
library(tidytree)
library(ape)
library(tidyverse)
library(patchwork)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
clone_dir <- args[1]
clone_id <- args[2]
mcmc_chains <- as.integer(args[3])
eff_thr <- as.double(args[4])
psrf_thr <- as.double(args[5])

clone_path <- paste0(clone_dir, "/")
trait_file <- paste0(clone_path, clone_id, ".trait.tsv")
df <- read.table(trait_file, header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("label", "state")
num_categories <- length(unique(df$state))

# read in output
mcmc1_log <- parse_log(paste0(clone_path, "MCMC1/MCMC1.Log.txt"))
mcmc2_log <- parse_log(paste0(clone_path, "MCMC2/MCMC2.Log.txt"))
mcmc3_log <- parse_log(paste0(clone_path, "MCMC3/MCMC3.Log.txt"))

# grab sampling period
mcmc1_sampling <- parse_number(mcmc1_log$options[which(str_detect(mcmc1_log$options, 'Sample Period'), TRUE)])
mcmc2_sampling <- parse_number(mcmc2_log$options[which(str_detect(mcmc2_log$options, 'Sample Period'), TRUE)])
mcmc3_sampling <- parse_number(mcmc3_log$options[which(str_detect(mcmc3_log$options, 'Sample Period'), TRUE)])

# create mcmc objects for compatibility with coda
mcmc1 <- as_draws_matrix(mcmc1_log$results)

mcmc1 <- mcmc1_log$results %>%
  janitor::clean_names() %>%
  select(., -tree_no, -iteration) %>%
  mcmc(start = min(mcmc1_log$results$Iteration), end = max(mcmc1_log$results$Iteration), thin = mcmc1_sampling)
mcmc2 <- mcmc2_log$results %>%
  janitor::clean_names() %>%
  select(., -tree_no, -iteration) %>%
  mcmc(start = min(mcmc2_log$results$Iteration), end = max(mcmc2_log$results$Iteration), thin = mcmc2_sampling)
mcmc3 <- mcmc3_log$results %>%
  janitor::clean_names() %>%
  select(., -tree_no, -iteration) %>%
  mcmc(start = min(mcmc3_log$results$Iteration), end = max(mcmc3_log$results$Iteration), thin = mcmc3_sampling)

# combine chains
res_mcmc <- mcmc.list(mcmc1, mcmc2, mcmc3)

pdf(paste0(clone_path, clone_id, ".traceplots.pdf"), width=10, height=6)
par(mfrow = c(3,5))
if (num_categories == 4) {
	traceplot(
	  res_mcmc[,c('lh', 'q01', 'q02', 'q03', 'q10', 'q12', 'q13', 'q20', 'q21', 'q23', 'q30', 'q31', 'q32')]
	)
} else if (num_categories == 3) {
	traceplot(
	  res_mcmc[,c('lh', 'q01', 'q02', 'q10', 'q12', 'q20', 'q21')]
	)
} else if (num_categories == 2) {
	traceplot(
          res_mcmc[,c('lh', 'q01', 'q10')]
	)
}
dev.off()

# calculate effective sample sizes
ess <- effectiveSize(res_mcmc)
#summary(ess)

# extract the median from summary(ess)
ess_sum <- summary(ess)
# 'summary(ess)' should be a named numeric vector with element "Median"
ess_median <- unname(ess_sum["Median"])

cat("[ESS] median =", ess_median, " | threshold =", eff_thr, "\n")

if (is.na(ess_median) || ess_median < eff_thr) {
  message(
    sprintf(
      "[ESS] The median ESS (%.0f) is below the threshold (%.0f). ",
      ess_median, eff_thr
    ),
    "It appears the 'effective_size_threshold' was not met. ",
    "Consider increasing the number of MCMC iterations (and/or thinning/chains)."
  )
} else {
  message("[ESS] Effective size threshold satisfied.")
}

# calculate gelman diagnostic
gelman <- gelman.diag(res_mcmc, autoburnin = FALSE, multivariate = FALSE)
#summary(gelman$psrf)

lower_bound <- 0.99 * psrf_thr
upper_bound <- 1.01 * psrf_thr

# Safely compute medians for Point estimate and Upper CI
psrf_mat <- gelman$psrf
if (!all(c("Point est.", "Upper C.I.") %in% colnames(psrf_mat))) {
  stop("gelman$psrf must have columns 'Point est.' and 'Upper C.I.'.")
}

psrf_point_median <- median(psrf_mat[, "Point est."], na.rm = TRUE)
psrf_upper_median <- median(psrf_mat[, "Upper C.I."], na.rm = TRUE)

cat("[PSRF] median(Point est.) =", psrf_point_median, "\n")
cat("[PSRF] median(Upper C.I.) =", psrf_upper_median, "\n")
cat("[PSRF] allowed range =", lower_bound, "to", upper_bound, "\n")

# Check both medians
point_in_range <- !is.na(psrf_point_median) &&
  psrf_point_median >= lower_bound &&
  psrf_point_median <= upper_bound

upper_in_range <- !is.na(psrf_upper_median) &&
  psrf_upper_median >= lower_bound &&
  psrf_upper_median <= upper_bound

if (point_in_range && upper_in_range) {
  message("[PSRF] PSRF threshold satisfied (both Point est. and Upper C.I. within acceptable range).")
} else {
  message(
    sprintf(
      "[PSRF] PSRF median values out of range. Point est. = %.4f, Upper C.I. = %.4f (allowed [%.2f, %.2f]). ",
      psrf_point_median, psrf_upper_median, lower_bound, upper_bound
    ),
    "This suggests chains may not have mixed sufficiently. ",
    "Consider increasing MCMC iterations and/or adjusting burn-in or number of chains."
  )
}

# grab draws
if (num_categories == 4) {
  d_mcmc <- tidy_draws(res_mcmc) %>%
    select(., -lh, -q01, -q02, -q03, -q10, -q12, -q13, -q20, -q21, -q23, -q30, -q31, -q32, -root_p_0, -root_p_1, -root_p_2, -root_p_3) %>%
    pivot_longer(starts_with('x'), names_to = 'var', values_to = 'prob') %>%
    filter(!is.na(prob)) %>%
    separate(col = var, c('node', 'p', 'state'), sep = '_') %>%
    select(-p) %>%
    mutate(node = gsub('x', '', node))
} else if (num_categories == 3) {
  d_mcmc <- tidy_draws(res_mcmc) %>%
    select(., -lh, -q01, -q02, -q10, -q12, -q20, -q21, -root_p_0, -root_p_1, -root_p_2) %>%
    pivot_longer(starts_with('x'), names_to = 'var', values_to = 'prob') %>%
    filter(!is.na(prob)) %>%
    separate(col = var, c('node', 'p', 'state'), sep = '_') %>%
    select(-p) %>%
    mutate(node = gsub('x', '', node))
} else if (num_categories == 2) {
  d_mcmc <- tidy_draws(res_mcmc) %>%
    select(., -lh, -q01, -q10, -root_p_0, -root_p_1) %>%
    pivot_longer(starts_with('x'), names_to = 'var', values_to = 'prob') %>%
    filter(!is.na(prob)) %>%
    separate(col = var, c('node', 'p', 'state'), sep = '_') %>%
    select(-p) %>%
    mutate(node = gsub('x', '', node))
}

# plot the distributions for each chain
#p <- ggplot(d_mcmc, aes(prob)) +
  #geom_line(aes(col = state, group = interaction(.chain, state)), stat = 'density') +
  #facet_wrap(~node, scales = 'free') +
  #theme_bw()# +
  #scale_color_manual(values = c('black', 'blue'))

#ggsave(paste0(clone_path, "ASR.pdf"), units = 'cm', width = 20, height = 16, plot = p)

# Create a data frame from d_mcmc
df <- data.frame(
  node = paste0("T", d_mcmc$node),
  state = paste0("state", d_mcmc$state),
  prob = d_mcmc$prob
  )

# Create an empty dataframe to store results
unique_nodes <- unique(df$node)  # Get unique node values
df_wide <- data.frame(node = unique_nodes)

# Extract probability values for each state and add them as separate columns
for (state in unique(df$state)) {
  df_wide[[state]] <- df$prob[df$state == state][match(df_wide$node, df$node)]
}

# Fix column types
df_wide <- df_wide %>%
  mutate(across(starts_with("state"), as.numeric))

# Find the most probable state and its probability correctly
df_wide <- df_wide %>%
  rowwise() %>%
  mutate(
    anc_state = names(pick(starts_with("state")))[which.max(c_across(starts_with("state")))],  # Get correct column name
    anc_prob = max(c_across(starts_with("state")), na.rm = TRUE)  # Extract max probability
  ) %>%
  ungroup()

# Print the result
print(df_wide)

# Save the result as a tab-separated file
write.table(df_wide, paste0(clone_path, clone_id, "_ASE.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
