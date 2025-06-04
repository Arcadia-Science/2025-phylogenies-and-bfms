#!/bin/bash
set -e

# Script to process a directory of FASTA files with MMSeqs2 clustering
# Based on fig4_cluster_sample_families.sh clustering logic

# Default parameters
COV=0.8
FASTA_DIR=""
OUTPUT_DIR="clustering_results"

# Usage function
usage() {
    echo "Usage: $0 -d <fasta_directory> [-c coverage] [-o output_directory]"
    echo "  -d: Directory containing FASTA files (required)"
    echo "  -c: Coverage threshold (default: 0.8)"
    echo "  -o: Output directory (default: clustering_results)"
    echo ""
    echo "Example: $0 -d /path/to/fasta/files -c 0.8 -o my_results"
    exit 1
}

# Parse command-line arguments
while getopts "d:c:o:h" opt; do
  case $opt in
    d) FASTA_DIR="$OPTARG" ;;
    c) COV="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    h) usage ;;
    *) echo "Invalid option. Use -h for help." >&2
       exit 1 ;;
  esac
done

# Check required arguments
if [[ -z "$FASTA_DIR" ]]; then
    echo "Error: FASTA directory (-d) is required."
    usage
fi

# Validate FASTA directory exists
if [[ ! -d "$FASTA_DIR" ]]; then
    echo "Error: FASTA directory '$FASTA_DIR' does not exist."
    exit 1
fi

# Check if MMSeqs2 is installed
if ! command -v mmseqs &> /dev/null; then
    echo "Error: MMSeqs2 is not installed. Please install it first."
    echo "Installation instructions: https://github.com/soedinglab/MMseqs2"
    exit 1
fi

# Find FASTA files in the directory
fasta_files=($(find "$FASTA_DIR" -name "*.fasta" -o -name "*.fa" -o -name "*.fas"))

if [[ ${#fasta_files[@]} -eq 0 ]]; then
    echo "Error: No FASTA files found in directory '$FASTA_DIR'"
    echo "Looking for files with extensions: .fasta, .fa, .fas"
    exit 1
fi

echo "Found ${#fasta_files[@]} FASTA files in '$FASTA_DIR'"
echo "Using coverage threshold: $COV"
echo "Output directory: $OUTPUT_DIR"

# Define sequence identity thresholds
SEQ_IDS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

# Create output directories
mkdir -p "$OUTPUT_DIR"

# Create summary header
summary_file="$OUTPUT_DIR/clustering_summary.csv"
echo "family,seq_identity,num_clusters,total_sequences" > "$summary_file"

echo "Starting clustering analysis..."
echo "Sequence identity thresholds: ${SEQ_IDS[*]}"

# Process each family at each sequence identity threshold
for fasta in "${fasta_files[@]}"; do
    # Extract base name without path and extension
    base=$(basename "$fasta")
    base="${base%.*}"  # Remove extension
    
    echo ""
    echo "Processing family: $base"
    echo "File: $fasta"
    
    # Create directory for this family
    family_dir="$OUTPUT_DIR/$base"
    mkdir -p "$family_dir"
    
    # Calculate total sequences
    total_seqs=$(grep -c ">" "$fasta" 2>/dev/null || echo "0")
    
    if [[ $total_seqs -eq 0 ]]; then
        echo "  Warning: No sequences found in $fasta, skipping..."
        continue
    fi
    
    echo "  Total sequences: $total_seqs"
    
    # Run MMSeqs2 clustering at each threshold
    for seq_id in "${SEQ_IDS[@]}"; do
        echo "  Clustering at ${seq_id} identity threshold..."
        
        # Create temporary directory for this run
        tmp_dir="$family_dir/tmp_${seq_id}"
        mkdir -p "$tmp_dir"
        
        # Run MMSeqs2 clustering with error handling
        if mmseqs createdb "$fasta" "$tmp_dir/seq_db" &>/dev/null; then
            if mmseqs cluster "$tmp_dir/seq_db" "$tmp_dir/cluster_db" "$tmp_dir/tmp" \
                --min-seq-id "${seq_id}" -c "${COV}" --cov-mode 0 &>/dev/null; then
                
                # Extract cluster results
                mmseqs createtsv "$tmp_dir/seq_db" "$tmp_dir/seq_db" "$tmp_dir/cluster_db" \
                    "$family_dir/${base}_${seq_id}.tsv" &>/dev/null
                
                # Create FASTA files for representative sequences
                mmseqs createsubdb "$tmp_dir/cluster_db" "$tmp_dir/seq_db" "$tmp_dir/cluster_rep" &>/dev/null
                mmseqs convert2fasta "$tmp_dir/cluster_rep" "$family_dir/${base}_${seq_id}.rep.fasta" &>/dev/null
                
                # Calculate statistics
                if [[ -f "$family_dir/${base}_${seq_id}.tsv" ]]; then
                    num_clusters=$(cut -f1 "$family_dir/${base}_${seq_id}.tsv" | sort -u | wc -l | tr -d ' ')
                else
                    num_clusters="0"
                fi
                
                # Add to summary file
                echo "$base,$seq_id,$num_clusters,$total_seqs" >> "$summary_file"
                
                echo "    Clusters: $num_clusters"
            else
                echo "    Error: MMSeqs2 clustering failed"
                echo "$base,$seq_id,ERROR,${total_seqs}" >> "$summary_file"
            fi
        else
            echo "    Error: MMSeqs2 database creation failed"
            echo "$base,$seq_id,ERROR,${total_seqs}" >> "$summary_file"
        fi
        
        # Clean up temporary files
        rm -rf "$tmp_dir"
    done
    
    echo "  Completed clustering for $base"
done

echo ""
echo "Clustering analysis complete!"
echo "Results saved in: '$OUTPUT_DIR'"
echo "Summary statistics: '$summary_file'"
echo ""

# Display summary statistics
if [[ -f "$summary_file" ]]; then
    echo "Summary of results:"
    echo "Total families processed: $(tail -n +2 "$summary_file" | cut -d',' -f1 | sort -u | wc -l)"
    echo "Total clustering runs: $(tail -n +2 "$summary_file" | wc -l)"
    echo "Failed runs: $(grep -c "ERROR" "$summary_file" || echo "0")"
    
    echo ""
    echo "First few lines of summary:"
    head -6 "$summary_file"
    
    if [[ $(wc -l < "$summary_file") -gt 6 ]]; then
        echo "..."
        echo "(See full results in $summary_file)"
    fi
fi