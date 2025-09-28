import os
import glob
import sys
import pandas as pd

class SNV:
    """
    A class to represent a Single Nucleotide Variant (SNV).
    """
    def __init__(self, chr, pos, Ref_allele, Alt_allele, Depth_total, Depth_ref,
                 Depth_alt, SVM_pos_score, LDrefine_twoLoci_score,
                 LDrefine_trioLoci_score, LDrefine_merged_score, BAF_alt):
        self.chr = chr
        self.pos = pos
        self.Ref_allele = Ref_allele
        self.Alt_allele = Alt_allele
        self.Depth_total = Depth_total
        self.Depth_ref = Depth_ref
        self.Depth_alt = Depth_alt
        self.SVM_pos_score = SVM_pos_score
        self.LDrefine_twoLoci_score = LDrefine_twoLoci_score
        self.LDrefine_trioLoci_score = LDrefine_trioLoci_score
        self.LDrefine_merged_score = LDrefine_merged_score
        self.BAF_alt = BAF_alt

    def passes_filter(self, thresholds):
        """
        Check if the SNV passes the specified thresholds.
        
        Parameters:
            thresholds (dict): A dictionary containing threshold values.
        
        Returns:
            bool: True if SNV passes all thresholds, False otherwise.
        """
        #if self.Depth_total < thresholds['Depth_total']:
        if self.Depth_total < sum([thresholds['Depth_ref'][0], thresholds['Depth_alt'][0]]):
            return False
        if self.Depth_ref < thresholds['Depth_ref'][0]:
            return False
        if self.Depth_alt < thresholds['Depth_alt'][0]:
            return False
        if self.SVM_pos_score < thresholds['SVM_pos_score'][0]:
            return False
        if self.LDrefine_merged_score < thresholds['LDrefine_merged_score'][0]:
            return False
        if self.BAF_alt > thresholds['BAF_alt'][0]:
            return False
        return True

def check_monopogen_variants_files(directory):
    """
    Check for the presence of required monopogen somatic variants calling output files 
    for chromosomes 1 through 22 in the specified directory.
    
    Parameters:
        directory (str): Path to the monopogen somatic variants calling output folder.
    
    Returns:
        dict: A dictionary with chromosome numbers as keys and lists of missing file patterns as values.
              If all files are present for a chromosome, the list will be empty.
    """
    # Define the required file patterns for each chromosome
    required_files_patterns = [
        "chr{chr}.cell_snv.*.csv",          # e.g., chr1.cell_snv.cellA.csv
        "chr{chr}.cell_snv.*.filter.csv",   # e.g., chr1.cell_snv.cellA.filter.csv
        "chr{chr}.cell_snv.mat.gz",         # e.g., chr1.cell_snv.mat.gz
        "chr{chr}.cell_snv.snvID.csv",      # e.g., chr1.cell_snv.snvID.csv
        "chr{chr}.putativeSNVs.csv"         # e.g., chr1.putativeSNVs.csv
    ]
    
    # Initialize a dictionary to store missing files per chromosome
    missing_files = {}
    
    # Iterate over chromosomes 1 to 22
    for chr_num in range(1, 23):
        # Initialize a list to keep track of missing files for the current chromosome
        chr_missing = []
        
        # Iterate over each required file pattern
        for pattern in required_files_patterns:
            # Format the pattern with the current chromosome number
            formatted_pattern = pattern.format(chr=chr_num)
            
            # Create the full path with wildcard for cellID where applicable
            search_pattern = os.path.join(directory, formatted_pattern)
            
            # Use glob to find files matching the pattern
            matched_files = glob.glob(search_pattern)
            
            # If no files match the pattern, consider it missing
            if not matched_files:
                # Replace '*' with '<cellID>' for clarity in the report
                missing_pattern = formatted_pattern.replace('*', '<cellID>')
                chr_missing.append(missing_pattern)
        
        # If there are any missing files for the current chromosome, add them to the dictionary
        if chr_missing:
            missing_files[f"chr{chr_num}"] = chr_missing
    
    return missing_files

def load_initial_cell_info(directory):
    """
    Load initial cell barcode and index information from chr1.cell_snv.cellID.csv.
    
    Parameters:
        directory (str): Path to the monopogen somatic variants calling output folder.
    
    Returns:
        dict: Dictionary with cell barcodes as keys and their indices as values.
    """
    # Path to the chr1.cell_snv.cellID.csv file
    initial_file = os.path.join(directory, "chr1.cell_snv.cellID.csv")
    
    if not os.path.isfile(initial_file):
        print(f"Error: Initial cell info file '{initial_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Read the CSV file
    initial_df = pd.read_csv(initial_file)
    
    # Check if required columns are present
    if not {'cell', 'id', 'index'}.issubset(initial_df.columns):
        print(f"Error: The file '{initial_file}' does not contain the required columns.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded initial cell info from '{initial_file}'. Total cells: {initial_df.shape[0]}")
    
    # Convert to dictionary for easy lookup
    cell_info_dict = dict(zip(initial_df['cell'], initial_df['index']))
    
    return cell_info_dict

def load_filtered_cells(directory):
    """
    Load filtered cell barcodes from all chromosomes and identify common cell barcodes.
    
    Parameters:
        directory (str): Path to the monopogen somatic variants calling output folder.
    
    Returns:
        set: Set of common cell barcodes present in all chromosomes.
    """
    # List to hold sets of cell barcodes from each chromosome
    cell_sets = []
    
    for chr_num in range(1, 23):
        # Construct the filename pattern for the filtered cell barcodes
        filter_file_pattern = os.path.join(directory, f"chr{chr_num}.cell_snv.*.filter.csv")
        matched_files = glob.glob(filter_file_pattern)
        
        if not matched_files:
            print(f"Warning: No filter files found for chr{chr_num}. Skipping this chromosome.", file=sys.stderr)
            continue
        
        # Assuming only one filter file per chromosome
        filter_file = matched_files[0]
        
        # Read the filter CSV file
        try:
            filter_df = pd.read_csv(filter_file)
        except Exception as e:
            print(f"Error: Failed to read '{filter_file}'. Exception: {e}", file=sys.stderr)
            sys.exit(1)
        
        # Check if 'cell' column exists
        if 'cell' not in filter_df.columns:
            print(f"Error: The file '{filter_file}' does not contain the 'cell' column.", file=sys.stderr)
            sys.exit(1)
        
        # Extract the set of cell barcodes
        cell_set = set(filter_df['cell'])
        cell_sets.append(cell_set)
        
        #print(f"Loaded filtered cells from '{filter_file}'. Number of cells: {len(cell_set)}")
    
    if not cell_sets:
        print("Error: No filtered cell barcode files were loaded.", file=sys.stderr)
        sys.exit(1)
    
    # Find the intersection of all cell sets
    common_cells = set.intersection(*cell_sets)
    
    print(f"Number of common cells across all chromosomes: {len(common_cells)}")
    
    return common_cells

def filter_common_cells(cell_info_dict, common_cells):
    """
    Retain only the common cell barcodes and their index information.
    
    Parameters:
        cell_info_dict (dict): Dictionary with cell barcodes as keys and indices as values.
        common_cells (set): Set of common cell barcodes.
    
    Returns:
        dict: Filtered dictionary with common cell barcodes and their indices.
    """
    # Filter the dictionary to retain only common cells
    filtered_dict = {cell: idx for cell, idx in cell_info_dict.items() if cell in common_cells}
    
    print(f"Number of cells after filtering to common cells: {len(filtered_dict)}")
    
    return filtered_dict

def save_filtered_cell_info(filtered_dict, sample_id, output_directory):
    """
    Save the filtered cell barcode and index information to a CSV file.
    
    Parameters:
        filtered_dict (dict): Dictionary containing filtered cell info.
        sample_id (str): Sample identifier to prefix the output filename.
        output_directory (str): Path to the directory where output files will be saved.
    """
    output_file = os.path.join(output_directory, f"{sample_id}.cellID.filter.csv")
    
    # Convert dictionary to DataFrame
    filtered_df = pd.DataFrame(list(filtered_dict.items()), columns=['cell', 'index'])
    
    # If 'id' is required, it can be included from the initial data if needed
    # For now, only 'cell' and 'index' are saved as per the dictionary
    
    try:
        filtered_df.to_csv(output_file, index=False)
        print(f"Filtered cell info saved to '{output_file}'.")
    except Exception as e:
        print(f"Error: Failed to save filtered cell info to '{output_file}'. Exception: {e}", file=sys.stderr)
        sys.exit(1)
    
    return filtered_df

def process_putative_snvs(directory, sample_id, output_directory, thresholds):
    """
    Load, filter, and save putative SNVs from all chromosomes.
    
    Parameters:
        directory (str): Path to the monopogen somatic variants calling output folder.
        sample_id (str): Sample identifier to prefix the output filename.
        output_directory (str): Path to the directory where output files will be saved.
        thresholds (dict): Dictionary containing threshold values for filtering.
    
    Returns:
        None
    """
    all_filtered_snvs = []
    total_snvs = 0
    total_filtered_snvs = 0
    
    for chr_num in range(1, 23):
        # Construct the filename pattern for putativeSNVs.csv
        snv_file_pattern = os.path.join(directory, f"chr{chr_num}.putativeSNVs.csv")
        matched_files = glob.glob(snv_file_pattern)
        
        if not matched_files:
            print(f"Warning: No putativeSNVs.csv file found for chr{chr_num}. Skipping this chromosome.", file=sys.stderr)
            continue
        
        snv_file = matched_files[0]
        
        # Read the putativeSNVs.csv file
        try:
            snv_df = pd.read_csv(snv_file)
        except Exception as e:
            print(f"Error: Failed to read '{snv_file}'. Exception: {e}", file=sys.stderr)
            sys.exit(1)
        
        # Check if required columns exist
        required_columns = ['chr', 'pos', 'Ref_allele', 'Alt_allele',
                            'Depth_total', 'Depth_ref', 'Depth_alt',
                            'SVM_pos_score', 'LDrefine_twoLoci_score',
                            'LDrefine_trioLoci_score', 'LDrefine_merged_score',
                            'BAF_alt']
        
        if not set(required_columns).issubset(snv_df.columns):
            print(f"Error: The file '{snv_file}' does not contain all required columns.", file=sys.stderr)
            sys.exit(1)
        
        #print(f"\nProcessing SNVs from '{snv_file}'...")
        
        # Replace 'NA' with appropriate values, e.g., NaN
        snv_df.replace('NA', pd.NA, inplace=True)
        
        # Convert relevant columns to numeric, coerce errors to NaN
        numeric_columns = ['Depth_total', 'Depth_ref', 'Depth_alt',
                           'SVM_pos_score', 'LDrefine_twoLoci_score',
                           'LDrefine_trioLoci_score', 'LDrefine_merged_score',
                           'BAF_alt']
        snv_df[numeric_columns] = snv_df[numeric_columns].apply(pd.to_numeric, errors='coerce')
        
        # Instantiate SNV objects and apply filters
        for _, row in snv_df.iterrows():
            snv = SNV(
                chr=row['chr'],
                pos=row['pos'],
                Ref_allele=row['Ref_allele'],
                Alt_allele=row['Alt_allele'],
                Depth_total=row['Depth_total'] if pd.notna(row['Depth_total']) else 0,
                Depth_ref=row['Depth_ref'] if pd.notna(row['Depth_ref']) else 0,
                Depth_alt=row['Depth_alt'] if pd.notna(row['Depth_alt']) else 0,
                SVM_pos_score=row['SVM_pos_score'] if pd.notna(row['SVM_pos_score']) else 0,
                LDrefine_twoLoci_score=row['LDrefine_twoLoci_score'] if pd.notna(row['LDrefine_twoLoci_score']) else 0,
                LDrefine_trioLoci_score=row['LDrefine_trioLoci_score'] if pd.notna(row['LDrefine_trioLoci_score']) else 0,
                LDrefine_merged_score=row['LDrefine_merged_score'] if pd.notna(row['LDrefine_merged_score']) else 0,
                BAF_alt=row['BAF_alt'] if pd.notna(row['BAF_alt']) else 0
            )
            total_snvs += 1
            if snv.passes_filter(thresholds):
                all_filtered_snvs.append(snv)
                total_filtered_snvs += 1
        
        #print(f"Total SNVs in chr{chr_num}: {len(snv_df)}")
        #print(f"Filtered SNVs in chr{chr_num}: {total_filtered_snvs}")
    
    # After processing all chromosomes, save the filtered SNVs
    print("\nSaving filtered SNVs...")
    filtered_snvs_data = [{
        'chr': snv.chr,
        'pos': snv.pos,
        'Ref_allele': snv.Ref_allele,
        'Alt_allele': snv.Alt_allele,
        'Depth_total': snv.Depth_total,
        'Depth_ref': snv.Depth_ref,
        'Depth_alt': snv.Depth_alt,
        'SVM_pos_score': snv.SVM_pos_score,
        'LDrefine_twoLoci_score': snv.LDrefine_twoLoci_score,
        'LDrefine_trioLoci_score': snv.LDrefine_trioLoci_score,
        'LDrefine_merged_score': snv.LDrefine_merged_score,
        'BAF_alt': snv.BAF_alt
    } for snv in all_filtered_snvs]
    
    filtered_snvs_df = pd.DataFrame(filtered_snvs_data)
    
    output_file = os.path.join(output_directory, f"{sample_id}.SNVs.filter.csv")
    
    try:
        filtered_snvs_df.to_csv(output_file, index=False)
        print(f"Filtered SNVs saved to '{output_file}'.")
    except Exception as e:
        print(f"Error: Failed to save filtered SNVs to '{output_file}'. Exception: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nTotal SNVs processed: {total_snvs}")
    print(f"Total SNVs after filtering: {total_filtered_snvs}")
