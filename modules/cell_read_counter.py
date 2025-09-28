import sys
import pysam
from collections import Counter
from tqdm import tqdm
from pathlib import Path


def load_barcodes(barcode_file):
    """
    Load barcodes from a text file into a set.

    Parameters
    ----------
    barcode_file : str
        Path to the file containing valid barcodes (one per line).

    Returns
    -------
    set
        A set of barcodes.
    """
    with open(barcode_file, "r") as f:
        return {line.strip() for line in f}

def count_reads_per_barcode(bam_file, barcodes):
    """
    Count reads per cell barcode from a BAM file.

    Parameters
    ----------
    bam_file : str
        Path to the aligned BAM file.
    barcodes : set
        Set of valid barcodes to count.

    Returns
    -------
    Counter
        A Counter dictionary mapping barcode -> read count.
    """
    read_counts = Counter()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in tqdm(bam.fetch(until_eof=True), desc="Processing BAM"):
            # Only process reads with a cell barcode ("CB" tag)
            if read.has_tag("CB"):
                barcode = read.get_tag("CB")
                if barcode in barcodes:
                    read_counts[barcode] += 1
    return read_counts

def save_counts(read_counts, output_file):
    """
    Save read counts per barcode to a TSV file.

    Parameters
    ----------
    read_counts : Counter
        Dictionary of barcode -> read count.
    output_file : str
        Path to the output TSV file.
    """
    sorted_read_counts = sorted(read_counts.items(), key=lambda x: x[1], reverse=True)

    with open(output_file, "w") as out:
        out.write("cell\tcount\n")
        for barcode, count in sorted_read_counts:
            out.write(f"{barcode}\t{count}\n")

    print(f"Output saved to {output_file}")

def main(bam_file, barcode_file, output_path, prefix):
    """
    Main function to process BAM file and count reads per cell barcode.
    """
    Path(output_path).mkdir(parents=True, exist_ok=True)
    output_file = Path(output_path) / f"{prefix}.cell_read_counts.tsv"

    barcodes = load_barcodes(barcode_file)
    read_counts = count_reads_per_barcode(bam_file, barcodes)
    save_counts(read_counts, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python cell_read_counter.py <bam_file> <barcode_file> <output_path> <prefix>")
        sys.exit(1)

    bam_file = sys.argv[1]
    barcode_file = sys.argv[2]
    output_path = sys.argv[3]
    prefix = sys.argv[4]

    main(bam_file, barcode_file, output_path, prefix)
