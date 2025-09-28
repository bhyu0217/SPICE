import pandas as pd

# Define IUPAC code mappings
IUPAC_CODES = {
    'AA': 'A',
    'AC': 'M',
    'AG': 'R',
    'AT': 'W',
    'CC': 'C',
    'CG': 'S',
    'CT': 'Y',
    'GG': 'G',
    'GT': 'K',
    'TT': 'T',
}

def get_iupac_code(ref, alt):
    """
    Returns the IUPAC code for a given pair of reference and alternate alleles.
    
    Parameters:
    - ref (str): Reference allele (e.g., 'A')
    - alt (str): Alternate allele (e.g., 'G')
    
    Returns:
    - str: IUPAC code representing the allele combination
    """
    # Convert alleles to uppercase
    ref = ref.upper()
    alt = alt.upper()
    
    # Create a sorted pair to ensure consistency
    alleles = ''.join(sorted([ref, alt]))
    
    return IUPAC_CODES.get(alleles, 'N')  # Return 'N' if combination is not defined

def map_genotype(genotype, ref, alt):
    """
    Converts a genotype string into a nucleotide based on reference and alternate alleles.
    
    Parameters:
    - genotype (str): Genotype string (e.g., '0/0', '1/0', '0/1', '1/1', '2/0', '0/2', '3/1', etc.)
    - ref (str): Reference allele (e.g., 'A')
    - alt (str): Alternate allele (e.g., 'G')
    
    Returns:
    - str: Mapped nucleotide ('A', 'G', 'R', 'M', 'W', 'S', 'Y', 'K', 'C', 'T', 'N')
    """
    if genotype == '0/0':
        return 'N'  # No allele present
    
    parts = genotype.split('/')
    if len(parts) != 2:
        return 'N'  # Unexpected genotype format
    
    try:
        ref_count = int(parts[0])
        alt_count = int(parts[1])
    except ValueError:
        return 'N'  # Non-integer genotype values
    
    if ref_count > 0 and alt_count == 0:
        return ref  # Support for reference allele
    elif alt_count > 0 and ref_count == 0:
        return alt  # Support for alternate allele
    else:
        # Both support; use IUPAC code
        #print(genotype, ref, alt)
        return get_iupac_code(ref, alt)

def convert_snv_matrix_to_fasta(output_directory, sample_id):
    """
    Converts an SNV matrix CSV file into FASTA format and saves it.
    
    Parameters:
    - merged_snv_mat_path (str): Path to the merged SNV matrix CSV file
    - output_fasta_path (str): Path where the FASTA file will be saved
    """
    # 1. Load data
    print("Loading SNV matrix...")
    merged_snv_mat_path = output_directory + sample_id + '.SNV_mat.filter.csv'
    df = pd.read_csv(merged_snv_mat_path, index_col=0)
    print(f"Matrix size: {df.shape}")
    
    # 2. Extract Ref and Alt alleles from Variant_ID
    print("Extracting Ref and Alt alleles from Variant_ID...")
    variant_info = df.index.to_series().str.split(':').apply(lambda x: pd.Series({
        'chr': x[0],
        'pos': x[1],
        'Ref': x[2],
        'Alt': x[3]
    }))
    df = df.join(variant_info)
    print("Ref and Alt allele extraction complete.")
    
    # 3. Merge variant information with the matrix
    print("Merging variant information with the matrix...")
    # Uncomment the following line if you need to drop 'chr' and 'pos' columns
    # df = df.drop(columns=['chr', 'pos'])
    
    # 4. Verify and subset cell barcodes in the matrix
    print("Verifying and subsetting cell barcodes in the matrix...")
    cell_barcodes_in_matrix = df.columns.difference(['chr', 'pos', 'Ref', 'Alt'])
    print(f"Number of cell barcodes present in the matrix: {len(cell_barcodes_in_matrix)}")
    
    # 5. Generate FASTA sequences
    print("Generating FASTA sequences...")
    fasta_entries = []
    total_cells = len(cell_barcodes_in_matrix)
    #for idx, cell in enumerate(cell_barcodes_in_matrix, 1):
    for cell in cell_barcodes_in_matrix:
        # Extract genotype, Ref, and Alt for each variant in the cell
        genotypes = df[cell]
        refs = df['Ref']
        alts = df['Alt']

        nucleotides = []
        for variant_id, genotype in genotypes.items():
            parts = variant_id.split(':')
            if len(parts) == 4:
                chrom, pos, ref, alt = parts
                nucleotides.append(map_genotype(genotype, ref, alt))
        
        # Generate nucleotide sequence based on genotype
        #nucleotides = genotypes.apply(lambda g, r, a: map_genotype(g, r, a), args=(refs, alts))
        sequence = ''.join(nucleotides)
        
        # Add to FASTA format
        fasta_entries.append(f">{cell}\n{sequence}\n")
        
        # Print progress every 100 cells or on the last cell
        #if idx % 100 == 0 or idx == total_cells:
            #print(f"Completed: {idx}/{total_cells} cells")
    
    # 6. Save to FASTA file
    output_fasta_path = output_directory + sample_id + '.SNV_mat.filter.fasta'
    print(f"Saving to FASTA file: {output_fasta_path}")
    with open(output_fasta_path, 'w') as fasta_file:
        for entry in fasta_entries:
            fasta_file.write(entry + '\n')
    print("FASTA file saved successfully.")
