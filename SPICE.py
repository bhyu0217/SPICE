#!/usr/bin/env python3
"""
SPICE: Single-cell Plasticity Inference and Clonal Evolution

This module provides a command-line interface with four subcommands:
- filter
- phylogeny
- ancestry
- plasticity
"""

import os
import argparse
import sys
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any

from scripts.monopogen_filter import *
from scripts.convert_matrix_to_fasta import *

# ------------------------- Subcommand Implementations -------------------------
def run_filter(args: argparse.Namespace) -> None:
    """
    Filter step entrypoint.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments containing input/output paths and thresholds.
    """
    input_directory = args.input_directory
    output_directory = args.output_directory
    sample_id = getattr(args, "sample_id", None) or getattr(args, "prefix", None)
    cell_barcode = args.cell_barcode
    threads = args.threads

    print("[SPICE:filter] Input directory :", input_directory)
    print("[SPICE:filter] Output directory:", output_directory)
    print("[SPICE:filter] Sample ID       :", sample_id)

    # Validate input directory
    if not os.path.isdir(input_directory):
        print(f"Error: The input directory '{input_directory}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory, exist_ok=True)
            print(f"Created output directory: '{output_directory}'")
        except Exception as e:
            print(f"Error: Failed to create output directory '{output_directory}': {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Output directory '{output_directory}' already exists.")

    # Define thresholds with their comparison direction for clearer reporting
    thresholds = {
        "Depth_ref": (args.depth_ref, ">="),
        "Depth_alt": (args.depth_alt, ">="),
        "SVM_pos_score": (args.svm_pos_score, ">="),
        "LDrefine_merged_score": (args.ldrefine_merged_score, ">="),
        "BAF_alt": (args.baf_alt, "<="),
        "min_alt_cells_per_snv": (args.min_alt_cells_per_snv, ">="),
        "min_snvs_per_cell": (args.min_snvs_per_cell, ">="),
    }

    print("\nThresholds for filtering putative SNVs:")
    for key, (value, op) in thresholds.items():
        print(f"  {key}: {op} {value}")

    # ------------------------ Step 1: Check required files ---------------------
    print("\nChecking for required files...")
    try:
        missing = check_monopogen_variants_files(input_directory)
    except NameError:
        print("Error: 'check_monopogen_variants_files' is not implemented. Please provide this helper.", file=sys.stderr)
        sys.exit(2)

    if not missing:
        print("All required files are present for chromosomes 1 through 22.")
    else:
        print("Missing files detected:")
        for chr_key, files in missing.items():
            print(f"{chr_key}:")
            for file in files:
                print(f"  - {file}")
        print("\nPlease ensure all required files are present before proceeding.", file=sys.stderr)
        sys.exit(1)  # Exit if any files are missing

    # -------- Step 2: Load initial cell barcode and index information ----------
    print("\nLoading initial cell barcode and index information from chr1.cell_snv.cellID.csv...")
    try:
        cell_info_dict = load_initial_cell_info(input_directory)
    except NameError:
        print("Error: 'load_initial_cell_info' is not implemented. Please provide this helper.", file=sys.stderr)
        sys.exit(2)

    # --- Step 3: Load filtered cell barcodes from all chromosomes --------------
    print("\nLoading filtered cell barcodes from all chromosomes and identifying common cells...")
    try:
        common_cells = load_filtered_cells(input_directory)
    except NameError:
        print("Error: 'load_filtered_cells' is not implemented. Please provide this helper.", file=sys.stderr)
        sys.exit(2)

    # -------- Step 4: Filter initial info to retain only common cells ----------
    print("\nFiltering initial cell info to retain only common cells...")
    try:
        filtered_cell_info = filter_common_cells(cell_info_dict, common_cells)
    except NameError:
        print("Error: 'filter_common_cells' is not implemented. Please provide this helper.", file=sys.stderr)
        sys.exit(2)

    # --- Step 5: Save filtered common cell barcodes/index information ----------
    print("\nSaving filtered common cell barcodes and index information...")
    try:
        filtered_cells_df = save_filtered_cell_info(filtered_cell_info, sample_id, output_directory)
    except NameError:
        print("Error: 'save_filtered_cell_info' is not implemented. Please provide this helper.", file=sys.stderr)
        sys.exit(2)

    # ---------------- Step 6: Process putativeSNVs.csv files -------------------
    print("\nProcessing putativeSNVs.csv files for all chromosomes...")
    try:
        process_putative_snvs(input_directory, sample_id, output_directory, thresholds)
    except NameError:
        print("Error: 'process_putative_snvs' is not implemented. Please provide this helper.", file=sys.stderr)
        sys.exit(2)

    # ---------------- Step 7: Run R merging script --------------
    print("\nMerging results via R script (monopogen_merge.R)...")
    cmd = ["Rscript", "scripts/monopogen_merge.R", input_directory, output_directory, str(sample_id)]
    try:
        subprocess.run(cmd, check=True)
        print("R merge script completed successfully.")
    except FileNotFoundError:
        print("Warning: 'Rscript' not found or 'scripts/monopogen_merge.R' missing. Skipping this step.", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error executing R script: {e}", file=sys.stderr)
        sys.exit(1)

    # ---------------- Step 8: Mutation filter --------------
    print("\nRunning mutation filter R script (mutation_filter.R)...")
    cmd = ["Rscript", "scripts/mutation_filter.R", output_directory, str(sample_id), cell_barcode, str(thresholds["min_alt_cells_per_snv"][0]), str(thresholds["min_snvs_per_cell"][0]), str(threads)]
    try:
        subprocess.run(cmd, check=True)
        print("Mutation filter R script completed successfully.")
    except FileNotFoundError:
        print("Warning: 'Rscript' not found or 'scripts/mutation_filter.R' missing. Skipping this step.", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error executing mutation_filter.R: {e}", file=sys.stderr)
        sys.exit(1)

    # ------------------------- Step 9: Convert to FASTA ------------------------
    print("\nConverting somatic SNV matrix to FASTA...")
    try:
        convert_snv_matrix_to_fasta(output_directory, sample_id)
    except NameError:
        print("Error: 'convert_snv_matrix_to_fasta' is not implemented. Please provide this helper.", file=sys.stderr)
        sys.exit(2)

    print("\nAll processes completed successfully.")


def run_phylogeny(args: argparse.Namespace) -> None:
    """
    Phylogeny step entrypoint.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments including IQ-TREE2 configuration and branch cutting grid.
    """
    print("[SPICE:phylogeny] Starting phylogeny with arguments:")
    for k, v in vars(args).items():
        print(f"  - {k}: {v}")


def run_ancestry(args: argparse.Namespace) -> None:
    """
    Ancestry step entrypoint.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments for MCMC-based discrete multistate ancestry analysis.
    """
    print("[SPICE:ancestry] Starting ancestry with arguments:")
    for k, v in vars(args).items():
        print(f"  - {k}: {v}")


def run_plasticity(args: argparse.Namespace) -> None:
    """
    Plasticity step entrypoint.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments including permutation replicates and test direction.
    """
    print("[SPICE:plasticity] Starting plasticity with arguments:")
    for k, v in vars(args).items():
        print(f"  - {k}: {v}")


# ------------------------------ Helper Functions ------------------------------
def str_to_bool_strict(val: str) -> bool:
    """
    Convert a string to boolean with strict 'true'/'false' semantics.

    Parameters
    ----------
    val : str
        Expected to be 'true' or 'false' (case-insensitive).

    Returns
    -------
    bool
        Parsed boolean value.

    Raises
    ------
    argparse.ArgumentTypeError
        If value is not one of {'true', 'false'}.
    """
    lv = val.lower()
    if lv in {"true", "t", "1", "yes", "y"}:
        return True
    if lv in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("Expected a boolean string: 'true' or 'false'.")


def ensure_dir(path: Path) -> None:
    """Create the directory if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)


# --------------------------------- CLI Parser ---------------------------------
def build_parser() -> argparse.ArgumentParser:
    """
    Build the top-level parser and subparsers to match the provided usage text.
    """
    parser = argparse.ArgumentParser(
        prog="python SPICE.py",
        description="SPICE command-line interface",
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", metavar="")

    # ------------------------------- filter -----------------------------------
    p_filter = subparsers.add_parser(
        "filter",
        help="Filter somatic SNVs and cells using quality and model thresholds",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Mandatory positional arguments
    p_filter.add_argument("input_directory", help="Path to the monopogen somatic variants calling output folder")
    p_filter.add_argument("output_directory", help="Path to the directory where outputs will be saved")
    p_filter.add_argument("prefix", help="Identifier to prefix output filenames")
    p_filter.add_argument("cell_barcode", help="File containing cell barcodes to be used in the analysis")

    # Optional arguments
    p_filter.add_argument("--depth_ref", type=int, default=5,
                          help="Minimum threshold for the number of cells supporting the reference allele")
    p_filter.add_argument("--depth_alt", type=int, default=5,
                          help="Minimum threshold for the number of cells supporting the alternative allele")
    p_filter.add_argument("--svm_pos_score", type=float, default=0.1,
                          help="Minimum threshold from the Monopogen SVM module")
    p_filter.add_argument("--ldrefine_merged_score", type=float, default=0.25,
                          help="Minimum threshold from the Monopogen LD refinement module")
    p_filter.add_argument("--baf_alt", type=float, default=0.5,
                          help="Maximum threshold for the alternative allele frequency (BAF)")
    p_filter.add_argument("--min_alt_cells_per_snv", type=int, default=5,
                          help="Minimum number of cells that must support a mutated allele")
    p_filter.add_argument("--min_snvs_per_cell", type=int, default=5,
                          help="Minimum number of somatic SNVs that must be supported")
    p_filter.add_argument("--threads", type=int, default=1,
                          help="Number of threads to use")
    p_filter.set_defaults(func=run_filter)

    # ------------------------------ phylogeny ---------------------------------
    p_phy = subparsers.add_parser(
        "phylogeny",
        help="Infer phylogeny with IQ-TREE2 and apply support/branch filters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_phy.add_argument("--include_failed_chisq", type=str_to_bool_strict, default=False,
                       choices=[True, False],
                       help="Determines whether to include cells that do not pass the IQTREE2 composition chi-square test")
    p_phy.add_argument("--model", type=str, default="AUTO",
                       help="Specifies the model selection option for IQTREE2")
    p_phy.add_argument("--uf_bootstrap_replicates", type=int, default=1000,
                       help="Number of replicates (≥1000) for ultrafast bootstrap analysis")
    p_phy.add_argument("--sh_alrt_replicates", type=int, default=1000,
                       help="Number of replicates (≥1000) to perform the SH-like approximate likelihood ratio test (SH-aLRT)")
    p_phy.add_argument("--uf_support_threshold", type=int, default=90,
                       help="Branch support threshold value to be applied if ultrafast bootstrap is performed")
    p_phy.add_argument("--sh_support_threshold", type=int, default=75,
                       help="Branch support threshold value to be applied if the SH-aLRT is performed")
    p_phy.add_argument("--branch_cut_min", type=float,
                       help="Minimum value for the branch-length cutting range")
    p_phy.add_argument("--branch_cut_max", type=float,
                       help="Maximum value for the branch-length cutting range")
    p_phy.add_argument("--branch_cut_step", type=float,
                       help="Step size for the branch-length cutting range")
    p_phy.add_argument("--min_tips", type=int, default=50,
                       help="Threshold for the minimum number of tips in the subclonal phylogenetic tree")
    p_phy.add_argument("--threads", type=int, default=1,
                       help="Number of threads to use")

    # Required positional arguments
    p_phy.add_argument("output_directory", help="")
    p_phy.add_argument("sample_id", help="")
    p_phy.set_defaults(func=run_phylogeny)

    # ------------------------------ ancestry ----------------------------------
    p_anc = subparsers.add_parser(
        "ancestry",
        help="Estimate ancestral states / ASE with MCMC and convergence diagnostics",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_anc.add_argument("--mcmc_chains", type=int, help="Number of MCMC chains to run.")
    p_anc.add_argument("--discrete_states", type=int, help="Number of discrete multistates to be used for ASE analysis.")
    p_anc.add_argument("--iterations", type=int, help="Total number of iterations for the MCMC.")
    p_anc.add_argument("--burnin", type=int, help="Number of initial iterations to discard as burn-in.")
    p_anc.add_argument("--rate_prior", type=float, help="Prior value for the substitution rates.")
    p_anc.add_argument("--stepping_stones", type=int, help="Number of stepping stones used for marginal likelihood estimation.")
    p_anc.add_argument("--log_sample_period", type=int, help="Sample period (in iterations) for log output.")
    p_anc.add_argument("--effective_size_threshold", type=float, help="The effective size threshold used to assess MCMC convergence.")
    p_anc.add_argument("--psrf_threshold", type=float, help="The Gelman diagnostic PSRF threshold for evaluating MCMC convergence.")
    # Required positional arguments
    p_anc.add_argument("output_directory", help="")
    p_anc.add_argument("sample_id", help="")
    p_anc.add_argument("cell_state", help="")
    p_anc.set_defaults(func=run_ancestry)

    # ------------------------------ plasticity --------------------------------
    p_pl = subparsers.add_parser(
        "plasticity",
        help="Quantify plasticity with permutation testing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_pl.add_argument("--perm_replicates", type=int, help="Number of permutation replicates to perform.")
    p_pl.add_argument("--sig_direction", choices={"greater", "less", "two-sided"},
                      help="Specifies the test direction for calculating statistical significance.")
    # Required positional arguments
    p_pl.add_argument("output_directory", help="")
    p_pl.add_argument("sample_id", help="")
    p_pl.set_defaults(func=run_plasticity)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    """
    Program entrypoint.

    Parameters
    ----------
    argv : list[str], optional
        Argument vector; defaults to sys.argv[1:] when None.

    Returns
    -------
    int
        Exit status code.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if not hasattr(args, "func"):
        # No subcommand provided
        parser.print_help()
        return 2

    # Ensure output directory exists for subcommands that define it
    if hasattr(args, "output_directory") and args.output_directory:
        ensure_dir(Path(args.output_directory))

    # Dispatch to the selected subcommand
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
