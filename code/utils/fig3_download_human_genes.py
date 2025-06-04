#!/usr/bin/env python3
"""
This script automatically downloads a bulk FASTA file of human cDNA sequences (DNA for transcripts)
from Ensembl, filters it to retain one DNA sequence per protein-coding gene (yielding ~20,000 sequences),
and splits the filtered records into individual gene FASTA files.

Usage:
    python download_human_genes.py [options]

Optional arguments:
    --bulk_url      URL to download the bulk human cDNA FASTA file.
                    (Default: Ensembl GRCh38 cDNA file)
    --bulk_file     Local filename for the downloaded bulk FASTA (default: bulk_human_genes_dna.fa.gz)
    --gene_dir      Directory to store individual gene FASTA files (default: human_genes_dna)
    --verbose       Enable verbose output.

Dependencies:
    - Biopython
"""

import argparse
import glob
import gzip
import os
import sys
import urllib.request

from Bio import SeqIO

# Default URL: Ensembl release 104 cDNA file for Homo sapiens (GRCh38 DNA sequences with annotations)
DEFAULT_BULK_URL = "ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
DEFAULT_BULK_FILE = "bulk_human_genes_dna.fa.gz"
DEFAULT_GENE_DIR = "human_genes_dna"

def download_bulk_file(url, local_filename, verbose=False):
    """Download the bulk FASTA file if it does not already exist."""
    if os.path.exists(local_filename):
        if verbose:
            print(f"[INFO] Bulk file '{local_filename}' already exists. Skipping download.")
        return
    if verbose:
        print(f"[INFO] Downloading bulk FASTA file from {url} ...")
    try:
        urllib.request.urlretrieve(url, local_filename)
        if verbose:
            print(f"[INFO] Download completed. Saved as '{local_filename}'.")
    except Exception as e:
        print(f"[ERROR] Error downloading file: {e}")
        sys.exit(1)

def split_bulk_fasta_unique_genes(bulk_filename, output_dir, verbose=False):
    """
    Parse the bulk FASTA file (which contains many transcripts) and select one transcript
    per gene (only for protein-coding genes) by extracting the gene ID from the header.

    Writes each unique gene record to an individual FASTA file in the output directory.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        if verbose:
            print(f"[INFO] Created directory: {output_dir}")

    # Support for gzipped files
    if bulk_filename.endswith(".gz"):
        open_func = lambda f: gzip.open(f, "rt")
    else:
        open_func = open

    gene_records = {}
    count_total = 0
    if verbose:
        print(f"[INFO] Parsing bulk FASTA file '{bulk_filename}' for unique protein-coding gene records...")
    with open_func(bulk_filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            count_total += 1
            # Look for the 'gene:' field in the description.
            gene_id = None
            for token in record.description.split():
                if token.startswith("gene:"):
                    gene_id = token.split("gene:")[1]
                    break
            if gene_id and "gene_biotype:protein_coding" in record.description:
                if gene_id not in gene_records:
                    gene_records[gene_id] = record
                    if verbose:
                        print(f"[DEBUG] Selected gene {gene_id} from transcript: {record.id}")
    if verbose:
        print(f"[INFO] Parsed {count_total} transcript records; selected {len(gene_records)} unique protein-coding gene records.")

    for gene_id, record in gene_records.items():
        safe_gene_id = "".join(c if c.isalnum() or c in ('-', '_') else "_" for c in gene_id)
        output_path = os.path.join(output_dir, f"{safe_gene_id}.fa")
        with open(output_path, "w") as out_handle:
            SeqIO.write(record, out_handle, "fasta")
        if verbose:
            print(f"[DEBUG] Wrote gene {gene_id} to file '{output_path}'.")

    print(f"[INFO] Wrote {len(gene_records)} individual gene FASTA files in '{output_dir}'.")

def main():
    parser = argparse.ArgumentParser(
        description="Download a bulk cDNA FASTA file, filter to one DNA sequence per protein-coding gene (≈20,000 sequences), and split into individual gene files."
    )
    parser.add_argument("--bulk_url", default=DEFAULT_BULK_URL,
                        help=f"URL to download bulk human cDNA FASTA file (default: {DEFAULT_BULK_URL}).")
    parser.add_argument("--bulk_file", default=DEFAULT_BULK_FILE,
                        help=f"Local filename for downloaded bulk file (default: {DEFAULT_BULK_FILE}).")
    parser.add_argument("--gene_dir", default=DEFAULT_GENE_DIR,
                        help=f"Directory to store individual gene FASTA files (default: {DEFAULT_GENE_DIR}).")
    parser.add_argument("--cleanup", action="store_true",
                        help="Remove both bulk file and gene directory after writing to output_file (implies --keep_bulk=False).")
    parser.add_argument("--keep_bulk", action="store_true",
                        help="Keep the downloaded bulk file after processing (default: delete).")
    parser.add_argument("--output_file",
                        help="Concatenate all individual FASTA files into a single output file before cleanup.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()

    verbose = args.verbose

    print("[INFO] Starting the download_human_genes pipeline.")
    print("[INFO] Step 1: Downloading bulk FASTA file...")
    download_bulk_file(args.bulk_url, args.bulk_file, verbose)

    print("[INFO] Step 2: Filtering and splitting bulk FASTA to obtain one transcript per protein-coding gene...")
    split_bulk_fasta_unique_genes(args.bulk_file, args.gene_dir, verbose)

    # Collect individual gene files to report the count
    gene_files = []
    for pattern in ["*.fasta", "*.fa"]:
        gene_files.extend(glob.glob(os.path.join(args.gene_dir, pattern)))

    gene_count = len(gene_files)
    print(f"[INFO] Created {gene_count} individual gene files in '{args.gene_dir}'.")

    # If an output file is specified, concatenate all gene files
    if args.output_file and gene_count > 0:
        try:
            with open(args.output_file, 'w') as outfile:
                for gene_file in gene_files:
                    with open(gene_file) as infile:
                        outfile.write(infile.read())
                        # Ensure there's a blank line between records
                        outfile.write('\n')
            if verbose:
                print(f"[INFO] Concatenated all gene files into '{args.output_file}'")
        except Exception as e:
            print(f"[WARNING] Failed to create concatenated output file: {e}")

    # Delete the bulk file unless --keep_bulk is specified
    if (not args.keep_bulk or args.cleanup) and os.path.exists(args.bulk_file):
        try:
            os.remove(args.bulk_file)
            if verbose:
                print(f"[INFO] Deleted bulk file '{args.bulk_file}'")
        except Exception as e:
            print(f"[WARNING] Failed to delete bulk file '{args.bulk_file}': {e}")

    # Delete the gene directory if cleanup is specified
    if args.cleanup and os.path.exists(args.gene_dir):
        try:
            import shutil
            shutil.rmtree(args.gene_dir)
            if verbose:
                print(f"[INFO] Deleted gene directory '{args.gene_dir}'")
        except Exception as e:
            print(f"[WARNING] Failed to delete gene directory '{args.gene_dir}': {e}")

    print("[INFO] Pipeline completed.")

if __name__ == "__main__":
    main()
