#!/usr/bin/env python3
"""
Extract all Newick trees from the Compara EMF file and save them to separate files.
"""

import os
import gzip
import argparse
from tqdm import tqdm

def sanitize_filename(name):
    """
    Sanitize a filename to ensure it's valid across all operating systems.
    
    Args:
        name (str): The original filename
        
    Returns:
        str: A sanitized filename
    """
    # Replace problematic characters with underscores
    import re
    # Remove characters that are invalid in filenames
    sanitized = re.sub(r'[\\/*?:"<>|]', '_', name)
    # Replace spaces with underscores
    sanitized = sanitized.replace(' ', '_')
    # Limit length of filename
    if len(sanitized) > 100:
        sanitized = sanitized[:100]
    return sanitized

def extract_protein_family_name(seq_ids, species_seqs):
    """
    Extract a meaningful name for the protein family based on sequence IDs or annotations.
    
    Args:
        seq_ids (list): List of Ensembl protein IDs
        species_seqs (list): List of (species, seq_id) tuples
    
    Returns:
        str: A meaningful name for the protein family
    """
    # Try to find gene names in the species_seqs data
    gene_names = []
    
    for parts in species_seqs:
        # The gene name is the last column (if it's not NULL)
        if len(parts) >= 8 and parts[-1] != "NULL":
            gene_names.append(parts[-1])
    
    # If we found gene names, use the most common one
    if gene_names:
        from collections import Counter
        name_counts = Counter(gene_names)
        most_common_name = name_counts.most_common(1)[0][0]
        return sanitize_filename(f"{most_common_name}")
    
    # If no gene names, try to extract protein family info from the sequence IDs
    
    # Try to find human (ENSP) IDs first
    human_ids = [sid for sid in seq_ids if sid.startswith('ENSP')]
    if human_ids:
        return sanitize_filename(f"human_{human_ids[0]}")
    
    # If no human IDs, look for mouse (ENSMUSP)
    mouse_ids = [sid for sid in seq_ids if sid.startswith('ENSMUSP')]
    if mouse_ids:
        return sanitize_filename(f"mouse_{mouse_ids[0]}")
    
    # If no human or mouse, look for any mammalian protein
    mammal_ids = [sid for sid in seq_ids if any(sid.startswith(prefix) for prefix in 
                 ['ENSP', 'ENSMUSP', 'ENSRNOP', 'ENSCAFP', 'ENSBTAP'])]
    if mammal_ids:
        return sanitize_filename(f"mammal_{mammal_ids[0]}")
    
    # If no mammalian IDs, use the first ID
    if seq_ids:
        return sanitize_filename(f"protein_{seq_ids[0]}")
    
    # If all else fails, use a generic name with the family number
    return "protein_family_unknown"

def extract_trees(emf_file, output_dir, max_families=0, save_seq_ids=False):
    """
    Extract Newick trees from the EMF file and save to separate files.
    
    Args:
        emf_file (str): Path to the EMF file
        output_dir (str): Output directory to save trees
        max_families (int): Maximum number of families to process (0 means all)
        save_seq_ids (bool): Whether to save sequence IDs for each tree
    """
    # Check if file is gzipped
    is_gzipped = emf_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    read_mode = 'rt' if is_gzipped else 'r'
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    family_count = 0
    current_seq_ids = []
    current_species_seqs = []
    processing_tree = False
    
    print(f"Extracting trees from {emf_file} to {output_dir}")
    
    with open_func(emf_file, read_mode) as f:
        # Count total lines to set up progress bar
        if max_families == 0:
            print("Counting families in file...")
            total_families = 0
            for line in f:
                if line.strip() == "//":
                    total_families += 1
            f.seek(0)  # Reset file pointer
        else:
            total_families = max_families
        
        with tqdm(total=total_families, desc="Extracting trees") as pbar:
            for line in f:
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Sequence lines start with SEQ
                if line.startswith("SEQ"):
                    parts = line.split()
                    if len(parts) >= 3:  # Make sure we have enough parts
                        species = parts[1]    # Species name
                        seq_id = parts[2]     # The Ensembl protein ID
                        
                        # Store sequence ID for tree naming
                        current_seq_ids.append(seq_id)
                        
                        # Store all sequence data for protein family name extraction
                        current_species_seqs.append(parts[1:])
                
                # DATA line indicates the tree will follow on the next line
                elif line == "DATA":
                    processing_tree = True
                    continue
                
                # Process the Newick tree line
                elif processing_tree:
                    processing_tree = False
                    tree_line = line
                    
                    # Generate a protein family name based on the sequence IDs
                    protein_family = extract_protein_family_name(current_seq_ids, current_species_seqs)
                    
                    # Use the protein family name for the file name
                    # Also include a numeric ID to ensure uniqueness
                    tree_file = os.path.join(output_dir, f"{protein_family}_{family_count:05d}.nh")
                    
                    with open(tree_file, 'w') as tf:
                        tf.write(tree_line + "\n")
                    
                    # Save sequence IDs if requested
                    if save_seq_ids:
                        seq_file = os.path.join(output_dir, f"{protein_family}_{family_count:05d}_seqs.txt")
                        with open(seq_file, 'w') as sf:
                            for seq_id in current_seq_ids:
                                sf.write(seq_id + "\n")
                
                # Family separator
                elif line == "//":
                    family_count += 1
                    pbar.update(1)
                    
                    # Clear sequence data for the next family
                    current_seq_ids = []
                    current_species_seqs = []
                    
                    # Check if we've reached the maximum number of families
                    if max_families > 0 and family_count >= max_families:
                        break
    
    print(f"Extracted {family_count} trees to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Extract Newick trees from Compara EMF file")
    parser.add_argument("--emf_file", default="data/Compara.113.protein_default.nh.emf.gz",
                        help="Path to the Compara EMF file")
    parser.add_argument("--output_dir", default="trees",
                        help="Directory to save extracted trees")
    parser.add_argument("--max_families", type=int, default=0,
                        help="Maximum number of families to process (0 means all)")
    parser.add_argument("--save_seq_ids", action="store_true",
                        help="Save sequence IDs for each tree")
    
    args = parser.parse_args()
    
    extract_trees(args.emf_file, args.output_dir, args.max_families, args.save_seq_ids)

if __name__ == "__main__":
    main()