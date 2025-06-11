# Utility script for biological data analysis
# Shared packages and functions across figure scripts

# Load all packages used across figure scripts
library(ape)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(arcadiathemeR)
library(pracma)
library(httr)
library(jsonlite)
library(parallel)
library(pbapply)
library(dplyr)
library(tidyr)
library(scales)

# Function to calculate hill diversity (from fig2)
normalize_hill_q1 <- function(tree, epsilon = 1e-10) {
  branches <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),
                                              function(x, y) which(y == x),
                                              y = tree$edge[, 2])],
                      tree$tip.label) + epsilon
  p_branch <- branches / sum(branches)
  entropy <- -sum(p_branch * log(p_branch))
  hill_div <- exp(entropy)
  
  n_tips <- length(tree$tip.label)
  faith_pd <- sum(branches)
  
  # Normalizations
  norm_by_tips <- hill_div / n_tips
  norm_by_pd <- hill_div / faith_pd
  
  return(list(
    Hill_diversity = hill_div,
    Norm_by_tips = norm_by_tips,
    Norm_by_PD = norm_by_pd,
    n_tips = n_tips,
    entropy = entropy,
    tree = tree
  ))
}

# Function to retrieve evo2 likelihoods for a fasta file
get_evo2 <- function(fasta_file) {
  # Extract gene ID from filename
  gene_id <- basename(tools::file_path_sans_ext(fasta_file))
  temp_result <- file.path(TEMP_DIR, paste0(gene_id, ".result"))
  
  tryCatch({
    # Read the FASTA file
    content <- readLines(fasta_file)
    
    # Extract sequence (only non-header lines)
    sequence_lines <- content[!grepl("^>", content)]
    sequence <- paste(sequence_lines, collapse = "")
    sequence <- substr(sequence, 1, MAX_SEQ_LENGTH)  # Limit sequence length
    
    # Skip if empty sequence
    if (nchar(sequence) == 0) {
      writeLines(paste0(gene_id, ",ERROR_EMPTY_SEQUENCE"), temp_result)
      return(data.frame(gene_id = gene_id, mean_likelihood = "ERROR_EMPTY_SEQUENCE"))
    }
    
    # Prepare API request
    request_data <- list(
      sequence = sequence,
      num_tokens = 8,
      top_k = 1,
      enable_sampled_probs = TRUE
    )
    
    # Make API request
    response <- httr::POST(
      url = URL,
      body = jsonlite::toJSON(request_data, auto_unbox = TRUE),
      httr::add_headers(
        'Content-Type' = 'application/json',
        'Authorization' = paste('Bearer', API_KEY),
        'nvcf-poll-seconds' = '30'
      ),
      httr::timeout(TIMEOUT)
    )
    
    # Parse response
    response_json <- httr::content(response, "parsed")
    
    # Extract probabilities
    if (!is.null(response_json$sampled_probs)) {
      probs <- unlist(response_json$sampled_probs)
      mean_val <- mean(probs, na.rm = TRUE)
      
      # Write result to file
      writeLines(paste0(gene_id, ",", mean_val), temp_result)
      return(data.frame(gene_id = gene_id, mean_likelihood = mean_val))
    } else {
      writeLines(paste0(gene_id, ",ERROR_NO_PROBS"), temp_result)
      return(data.frame(gene_id = gene_id, mean_likelihood = "ERROR_NO_PROBS"))
    }
  }, error = function(e) {
    # Handle any errors
    error_msg <- substr(as.character(e), 1, 50)
    writeLines(paste0(gene_id, ",ERROR_", error_msg), temp_result)
    return(data.frame(gene_id = gene_id, mean_likelihood = paste0("ERROR_", error_msg)))
  })
}

# Function to calculate Hill's diversity for a single column (from fig4)
hill_diversity <- function(column, q) {
  freqs <- table(column) / length(column)
  
  if (q == 1) {
    # Exponential of Shannon entropy
    diversity <- exp(-sum(freqs * log(freqs), na.rm = TRUE))
  } else if (q == 0) {
    # Species richness
    diversity <- length(freqs)
  } else {
    # General Hill number
    diversity <- (sum(freqs^q))^(1 / (1 - q))
  }
  return(diversity)
}

# Normalization function (from fig4)
normalize_hill <- function(D_q, S) {
  return((D_q - 1) / (S - 1))
}