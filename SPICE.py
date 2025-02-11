import os
import argparse
import sys
import subprocess
from scripts.monopogen_filter import *
from scripts.convert_matrix_to_fasta import *

def parse_arguments():
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Check for required monopogen somatic variants calling output files "
                    "for chromosomes 1 through 22 in a specified input directory, "
                    "process cell barcode information, filter putative SNVs, and save results to an output directory."
    )
    
    parser.add_argument(
        "input_directory",
        type=str,
        help="Path to the monopogen somatic variants calling output folder."
    )
    
    parser.add_argument(
        "output_directory",
        type=str,
        help="Path to the directory where filtered results will be saved."
    )
    
    parser.add_argument(
        "sample_id",
        type=str,
        help="Sample identifier to prefix output filenames."
    )
    
    # Optional arguments for thresholds
    parser.add_argument(
        "--depth_total",
        type=int,
        default=10,
        help="Minimum Depth_total threshold (default: 10)."
    )
    
    parser.add_argument(
        "--depth_ref",
        type=int,
        default=5,
        help="Minimum Depth_ref threshold (default: 5)."
    )
    
    parser.add_argument(
        "--depth_alt",
        type=int,
        default=5,
        help="Minimum Depth_alt threshold (default: 5)."
    )
    
    parser.add_argument(
        "--svm_pos_score",
        type=float,
        default=0.1,
        help="Minimum SVM_pos_score threshold (default: 0.1)."
    )
    
    parser.add_argument(
        "--ldrefine_merged_score",
        type=float,
        default=0.25,
        help="Minimum LDrefine_merged_score threshold (default: 0.25)."
    )
    
    parser.add_argument(
        "--baf_alt",
        type=float,
        default=0.5,
        help="Maximum BAF_alt threshold (default: 0.5)."
    )
    
    return parser.parse_args()

def main():
    """
    Main function to execute the file checking and processing.
    """
    # Parse command-line arguments
    args = parse_arguments()
    input_directory = args.input_directory
    output_directory = args.output_directory
    sample_id = args.sample_id
    
    # Define thresholds
    thresholds = {
        'Depth_total': args.depth_total,
        'Depth_ref': args.depth_ref,
        'Depth_alt': args.depth_alt,
        'SVM_pos_score': args.svm_pos_score,
        'LDrefine_merged_score': args.ldrefine_merged_score,
        'BAF_alt': args.baf_alt
    }
    
    print("Thresholds for filtering putative SNVs:")
    for key, value in thresholds.items():
        if key == "BAF_alt":
            print(f"  {key}: <= {value}")
        else:
            print(f"  {key}: >= {value}")
    
    # Check if the input directory exists
    if not os.path.isdir(input_directory):
        print(f"Error: The input directory '{input_directory}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
            print(f"Created output directory '{output_directory}'.")
        except Exception as e:
            print(f"Error: Failed to create output directory '{output_directory}'. Exception: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Output directory '{output_directory}' already exists.")
    
    # Step 1: Check for missing files
    print("\nChecking for required files...")
    missing = check_monopogen_variants_files(input_directory)
    
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
    
    # Step 2: Load initial cell barcode and index information
    print("\nLoading initial cell barcode and index information from chr1.cell_snv.cellID.csv...")
    cell_info_dict = load_initial_cell_info(input_directory)
    
    # Step 3: Load filtered cell barcodes from all chromosomes and find common cells
    print("\nLoading filtered cell barcodes from all chromosomes and identifying common cells...")
    common_cells = load_filtered_cells(input_directory)
    
    # Step 4: Filter the initial cell info to retain only common cells
    print("\nFiltering initial cell info to retain only common cells...")
    filtered_cell_info = filter_common_cells(cell_info_dict, common_cells)
    
    # Step 5: Save the filtered common cell barcodes and index information to the output directory
    print("\nSaving filtered common cell barcodes and index information...")
    filtered_cells_df = save_filtered_cell_info(filtered_cell_info, sample_id, output_directory)
    
    # Step 6: Process putativeSNVs.csv files
    print("\nProcessing putativeSNVs.csv files for all chromosomes...")
    process_putative_snvs(input_directory, sample_id, output_directory, thresholds)
    
    # Step 7
    cmd = ["Rscript", "scripts/monopogen_merge.R", input_directory, output_directory, sample_id]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing R script: {e}")

    # Step 8

    # Step 9
    print("\nConverting somatic SNV matrix to fasta...")
    convert_snv_matrix_to_fasta(output_directory, sample_id)

    print("\nAll processes completed successfully.")

if __name__ == "__main__":
    main()
