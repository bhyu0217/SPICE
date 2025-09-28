import os
import glob

def extract_failed_sequences(directory):
    """
    Searches for all '*.fasta.log' files in the given directory, reads each file,
    and extracts sequence names that are marked as "failed" in the composition analysis.
    
    Assumptions:
      - Each log file is a text file.
      - The composition section of the log includes lines where one of the columns
        indicates "passed" or "failed".
      - The sequence name is assumed to be the second column in such lines.
      
    Parameters:
      directory (str): The path to the directory containing the '*.fasta.log' files.
    
    Returns:
      A sorted list of unique sequence names (str) that have "failed" composition.
    """
    failed_seqs = set()  # Use a set to store unique sequence names
    
    # Create a glob pattern for all '*.fasta.log' files in the directory
    pattern = os.path.join(directory, "*.fasta.log")
    log_files = glob.glob(pattern)
    
    # Loop over each log file found
    for log_file in log_files:
        with open(log_file, 'r') as f:
            for line in f:
                # Remove any leading/trailing whitespace
                line = line.strip()
                # Skip empty lines
                if not line:
                    continue
                # Check if the line contains the word "failed"
                if "failed" in line:
                    # Split the line into columns (assuming whitespace separation)
                    parts = line.split()
                    # Check that there are enough parts (we expect at least 5 columns)
                    if len(parts) == 5:
                        # Assume the sequence name is in the second column (index 1)
                        seq_name = parts[1]
                        failed_seqs.add(seq_name)
    
    # Convert the set to a sorted list and return
    return sorted(list(failed_seqs))

def filter_fasta_by_barcode(output_directory, sample_id, failed_seqs):
    """
    Filter entries in the failed_sequences to keep only those whose header
    matches a barcode in the barcode_file.
    """
    # Open the FASTA file and the output FASTA file
    fasta_file = output_directory + sample_id + '.SNV_mat.filter.fasta'
    output_fasta = output_directory + sample_id + '.SNV_mat.filter.comp.fasta'
    with open(fasta_file, 'r') as fin, open(output_fasta, 'w') as fout:
        keep_entry = False  # flag to indicate if current FASTA entry should be kept

        for line in fin:
            line = line.rstrip('\n')
            # If this line starts a new FASTA entry (header line):
            if line.startswith('>'):
                # The header is everything after '>'
                header = line[1:].strip()
                # Check if the header (barcode) is in the set
                if not header in failed_seqs:
                    keep_entry = True
                    fout.write(line + '\n')  # write the header line
                else:
                    keep_entry = False
            else:
                # If this line is part of a sequence, write only if keep_entry is True
                if keep_entry:
                    fout.write(line + '\n')
