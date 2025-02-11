# Single-cell Plasticity Inference and Clonal Evolution (SPICE)

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Quick Start](#quickstart)
- [Subclone-Specific Mutational Signature Analysis](#mutsignatures)
- [Citation](#citation)
- [Contact Information](#contact)

## <a name="introduction"></a> Introduction
**SPICE** is a package that leverages single-cell data to classify subclonal populations and assess cellular plasticity via phylogenetic analysis. It delineates the intricate subclonal architecture of tumors and reveals evolutionary pathways that drive aggressive phenotypes and treatment resistance. Moreover, SPICE investigates the associated epigenetic regulations and mutational signatures, providing valuable insights into cancer progression and recurrence while paving the way for effective precision oncology strategies.

SPICE is composed of three main modules.

  **1) Somatic SNV Matrix Construction:** SPICE uses Monopogen (https://github.com/KChen-lab/Monopogen.git) output files as input to extract somatic SNVs from single-cell data. A hard filtering process selects high-confidence SNVs while removing variants with very low allele frequencies that could compromise phylogenetic tree construction. Additionally, cells with insufficient sequencing depth are excluded, resulting in a robust cell-by-variant matrix.

  **2) Phylogenetic Inference and Subclone Classification:** SPICE employs IQ-TREE2 (https://github.com/iqtree/iqtree2.git), a maximum likelihood-based tool, to infer phylogenetic trees. It determines the optimal number of subclones by cutting branches with branch support values, as calculated using UFBoot and SH-aLRT tests. This approach reliably classifies trustworthy subclones.

  **3) Ancestral State Estimation and Cellular Plasticity Evaluation:** For each subclone’s phylogenetic tree, SPICE maps the given cell states and uses Bayesian MCMC-based ancestral state estimation from BayesTraits (https://github.com/AndrewPMeade/BayesTraits-Release.git) to determine the parent node’s cell state. Based on these results, cellular plasticity is computed for each subclone, and its significance is assessed via permutation testing.

## <a name="installation"></a> Installation

### Dependencies
```python
python
R
iqtree2
bayestraits
```

python packages
```python
pandas
```

R packages
```python
ape
btw
coda
phangorn
phytools
posterior
tidybayes
tidytree
tidyverse
patchwork
readr
dplyr
tidyr
ggtree
ggplot2
```

### Installation
Install the current stable PyPi version of SPICE:
```
$ pip install SPICE
```

## <a name="quickstart"></a> Quick Start

### Somatic SNV Matrix Construction

Using monopogen somatic SNV calling files as input

```python
usage: python SPICE.py filter [-h] [--depth_total DEPTH_TOTAL] [--depth_ref DEPTH_REF]
                              [--depth_alt DEPTH_ALT] [--svm_pos_score SVM_POS_SCORE]
                              [--ldrefine_merged_score LDREFINE_MERGED_SCORE]
                              [--baf_alt BAF_ALT] [--min_alt_cells MIN_ALT_CELLS]
                              [--min_snvs_per_cell MIN_SNVS_PER_CELL]
                              [--threads NTHREADS]
                              input_directory output_directory sample_id cell_barcode

optional arguments:
  --depth_total
  --depth_ref
  --depth_alt
  --svm_pos_score
  --ldrefine_merged_score
  --baf_alt
  --min_alt_cells
  --min_snvs_per_cell
  --threads
```

### Phylogenetic Inference and Subclone Classification

Detects and uses the output files from the **filter** module as input.

```python
usage: python SPICE.py phylogeny [-h] [--include_failed_chisq {true,false}] [--model MODEL]
                                 [--uf_bootstrap_replicates UF_BOOTSTRAP_REPS] [--sh_alrt_replicates SH_ALRT_REPS]
                                 [--uf_support_threshold UF_SUPPORT_THRESHOLD] [--sh_support_threshold SH_SUPPORT_THRESHOLD]
                                 [--branch_cut_min BRANCH_CUT_MIN] [--branch_cut_max BRANCH_CUT_MAX]
                                 [--branch_cut_step BRANCH_CUT_STEP] [--min_tips MIN_TIPS]
                                 [--threads NTHREADS]
                                 output_directory sample_id

optional arguments:
  --include_failed_chisq         Determines whether to include cells that do not pass the IQTREE2 composition chi-square test.
  --model                        Specifies the model selection option for IQTREE2.
  --uf_bootstrap_replicates      Number of replicates (≥1000) for ultrafast bootstrap analysis.
  --sh_alrt_replicates           Number of replicates (≥1000) to perform the SH-like approximate likelihood ratio test (SH-aLRT).
  --uf_support_threshold         Branch support threshold value to be applied if ultrafast bootstrap is performed.
  --sh_support_threshold         Branch support threshold value to be applied if the SH-aLRT is performed.
  --branch_cut_min               Minimum value for the branch-length cutting range.
  --branch_cut_max               Maximum value for the branch-length cutting range.
  --branch_cut_step              Step size for the branch-length cutting range
  --min_tips                     Threshold for the minimum number of tips in the subclonal phylogenetic tree.
  --threads
```

### Ancestral State Estimation

Detects and uses the output files from the **phylogeny** module as input.

```python
usage: python SPICE.py ancestry [-h] [--mcmc_chains MCMC_CHAINS] [--discrete_states DISCRETE_STATES]
                                [--iterations ITERATIONS] [--burnin BURNIN]
                                [--rate_prior RATE_PRIOR] [--stepping_stones STEPPING_STONES]
                                [--log_sample_period LOG_SAMPLE_PERIOD]
                                [--effective_size_threshold EFFECTIVE_SIZE_THRESHOLD]
                                [--psrf_threshold PSRF_THRESHOLD]
                                output_directory sample_id cell_state

optional arguments:
  --mcmc_chains                 Number of MCMC chains to run.
  --discrete_states             Number of discrete multistates to be used for ASE analysis.
  --iterations                  Total number of iterations for the MCMC.
  --burnin                      Number of initial iterations to discard as burn-in.
  --rate_prior                  Prior value for the substitution rates.
  --stepping_stones             Number of stepping stones used for marginal likelihood estimation.
  --log_sample_period           Sample period (in iterations) for log output.
  --effective_size_threshold    The effective size threshold used to assess MCMC convergence.
  --psrf_threshold              The Gelman diagnostic PSRF threshold for evaluating MCMC convergence.
```

### Cellular Plasticity Evaluation

Detects and uses the output files from the **ancestry** module as input.

```python
usage: python SPICE.py plasticity [-h] [--perm_replicates PERM_REPLICATES]
                                  [--significance_direction {greater,less,two-sided}]
                                  output_directory sample_id

optional arguments:
  --perm_replicates               Number of permutation replicates to perform.
  --significance_direction        Specifies the test direction for calculating statistical significance.
```

## <a name="mutsignatures"></a> Subclone-Specific Mutational Signature Analysis
SPICE enables the calculation of mutational signatures (COSMIC) for individual subclones using a predefined set of signatures, thereby defining subclones linked to specific mutational patterns.

## <a name="citation"></a> Citation

Yu, B., Rosatti, E., Jain, A., Diaz, A. (2025)

## <a name="contact"></a> Contact Information
For any questions or to report issues, please contact Bohyeon Yu at bohyeon.yu@ucsf.edu.
