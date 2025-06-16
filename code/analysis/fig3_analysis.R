library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils/utils.R")))

# Download
system('wget -O data/homo_sapiens_gene_ages.csv https://chenzxlab.hzau.edu.cn/static/GenOrigin/download/Age_Inference/Homo_sapiens.csv')

# Load gene age data
all_ages <- read.csv(here("data/homo_sapiens_gene_ages.csv"))

print(paste("Gene age data loaded:", nrow(all_ages), "entries"))
print("First few rows of gene age data:")
print(head(all_ages))
print("Column names:")
print(colnames(all_ages))

# Clean gene age formatting (remove '>' symbols)
all_ages$gene_age <- gsub(">", "", all_ages$gene_age)

print("Gene age range after cleaning:")
age_numeric <- as.numeric(all_ages$gene_age)
print(paste("Min age:", min(age_numeric, na.rm=TRUE), "MYA"))
print(paste("Max age:", max(age_numeric, na.rm=TRUE), "MYA"))
print(paste("Median age:", median(age_numeric, na.rm=TRUE), "MYA"))

system('python3 code/utils/fig3_download_human_genes.py --gene_dir data/human_genes_dna')

# Configuration 
# Insert API key below
API_KEY <- NULL
URL <- "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate"
OUTPUT_DIR <- here("data/evo2_output")
TEMP_DIR <- file.path(OUTPUT_DIR, "temp")
MAX_PARALLEL <- min(parallel::detectCores() - 1, 15)  # Use up to 15 cores or max available - 1
MAX_SEQ_LENGTH <- 300  # Shortened sequence length for faster API responses
TIMEOUT <- 60  # Timeout for API requests in seconds

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)

# Set working directory to project root
ROOT_DIR <- getwd()

# Get maximum number of parallel connections supported by R session
print(paste("Using", MAX_PARALLEL, "cores for parallel processing"))

# Define the process_file function for parallel processing
process_file <- function(fasta_file) {
  return(get_evo2(fasta_file))
}

# Find all FASTA files
fasta_files <- list.files(path = here("data/human_genes_dna"), 
                          pattern = "\\.fa$", 
                          full.names = TRUE)

print(paste("Found", length(fasta_files), "FASTA files to process"))

# Option to limit the number of files for testing
# Comment out the next line for full processing
#fasta_files <- fasta_files[1:50]  # Process only the first 50 files for testing

print(paste("Will process", length(fasta_files), "files"))

# Process files in parallel
start_time <- Sys.time()

# Set up parallel cluster
cl <- makeCluster(MAX_PARALLEL)

# Export variables and functions to the cluster
clusterExport(cl, c("API_KEY", "URL", "TEMP_DIR", "MAX_SEQ_LENGTH", "TIMEOUT", "get_evo2", "process_file"))

# Load libraries on all workers
clusterEvalQ(cl, {
    library(httr)
    library(jsonlite)
    library(tools)
})

# Process files with progress bar
results <- pblapply(fasta_files, process_file, cl = cl)

# Stop the cluster
stopCluster(cl)

# Combine results into a data frame
results_df <- do.call(rbind, results)

# Calculate success and error counts
success_count <- sum(!grepl("ERROR", results_df$mean_likelihood))
error_count <- sum(grepl("ERROR", results_df$mean_likelihood))

elapsed_time <- Sys.time() - start_time
print(paste("Processing complete in", round(as.numeric(elapsed_time), 2), "seconds!"))
print(paste("Successfully processed:", success_count, "files"))
print(paste("Errors:", error_count, "files"))

# Save combined results to CSV
print("Combining results...")

write.csv(results_df, here("mean_likelihoods.csv"), row.names = FALSE)
print(paste("Results saved to", here("mean_likelihoods.csv")))

likelihoods <- read.csv("~/Documents/Research/github/2025-evo2-gene-age/output/mean_likelihoods.csv")

# Clean up gene identifiers (extract base gene ID before underscore)
original_ids <- likelihoods$gene_id[1:5]
likelihoods$gene_id <- unlist(lapply(strsplit(likelihoods$gene_id, "_"), function(x) x[1]))
cleaned_ids <- likelihoods$gene_id[1:5]

print("Gene ID cleaning example:")
print("Original IDs:")
print(original_ids)
print("Cleaned IDs:")
print(cleaned_ids)

print(paste("Total likelihood entries after cleaning:", length(unique(likelihoods$gene_id))))

# Identify genes present in both datasets
genes <- intersect(all_ages$ensembl_gene_id, likelihoods$gene_id)

print(paste("Genes in age dataset:", length(unique(all_ages$ensembl_gene_id))))
print(paste("Genes in likelihood dataset:", length(unique(likelihoods$gene_id))))
print(paste("Genes in both datasets:", length(genes)))

# Calculate overlap percentages
age_overlap <- round(100 * length(genes) / length(unique(all_ages$ensembl_gene_id)), 1)
lik_overlap <- round(100 * length(genes) / length(unique(likelihoods$gene_id)), 1)
print(paste("Overlap with age data:", age_overlap, "%"))
print(paste("Overlap with likelihood data:", lik_overlap, "%"))

# Match genes to their evolutionary ages
ages <- as.numeric(all_ages$gene_age[match(genes, all_ages$ensembl_gene_id)])

# Match genes to their likelihood values
liks <- as.numeric(likelihoods$mean_likelihood[match(genes, likelihoods$gene_id)])

print("Matched data summary:")
print(paste("Valid age values:", sum(!is.na(ages))))
print(paste("Valid likelihood values:", sum(!is.na(liks))))
print(paste("Complete pairs:", sum(!is.na(ages) & !is.na(liks))))

# Calculate rolling averages across evolutionary time
bins <- list()
window_size <- 250
step_size <- 1

print(paste("Calculating rolling windows with", window_size, "MYA windows, step size:", step_size))
print("Processing time windows...")

for (i in seq(1, 1500, step_size)) {
  if (i %% 100 == 1) print(paste("Processing window starting at", i, "MYA"))
  
  # Extract log-likelihood values for 250 MYA windows
  x <- log(liks[ages > i & ages < (i + window_size)]) * -1
  bins[[as.character(i)]] <- x
}

print(paste("Created", length(bins), "time windows"))

# Show window sizes
window_sizes <- unlist(lapply(bins, length))
print(paste("Window size range:", min(window_sizes), "to", max(window_sizes), "genes"))
print(paste("Median window size:", median(window_sizes), "genes"))

# Calculate median and standard error for each time window
roll_m <- unlist(lapply(bins, function(x) median(x, na.rm = TRUE)))
roll_error_m <- unlist(lapply(bins, function(x) pracma::std_err(na.omit(x))))

print("Rolling statistics summary:")
print(paste("Valid median values:", sum(!is.na(roll_m))))
print(paste("Valid error values:", sum(!is.na(roll_error_m))))

print("Median negative log-likelihood range:")
print(paste("Min:", round(min(roll_m, na.rm=TRUE), 3)))
print(paste("Max:", round(max(roll_m, na.rm=TRUE), 3)))
print(paste("Mean:", round(mean(roll_m, na.rm=TRUE), 3)))

print("Standard error range:")
print(paste("Min:", round(min(roll_error_m, na.rm=TRUE), 4)))
print(paste("Max:", round(max(roll_error_m, na.rm=TRUE), 4)))
print(paste("Mean:", round(mean(roll_error_m, na.rm=TRUE), 4)))

# Prepare data frame for plotting
val <- data.frame(x = seq(1, 1500, 1),
                 y = roll_m,
                 error = roll_error_m)

# Remove any rows with missing values
val_clean <- val[complete.cases(val), ]

# Create time series plot with confidence bands
p <- ggplot(data = val,
       aes(x = x,
           y = y,
           ymin = y - roll_error_m,
           ymax = y + roll_error_m)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  theme_arcadia() +
  xlab("Years before present (MYA)") +
  ylab("Negative log likelihood")

# Display the plot
print(p)

print("Time series plot generated")

# Calculate correlation between time and likelihood
correlation <- cor(val$x, val$y, use = "complete.obs")

print(paste("Correlation between time (MYA) and negative log-likelihood:", round(correlation, 4)))

# Additional correlation statistics
cor_test <- cor.test(val$x, val$y, use = "complete.obs")
print(paste("P-value:", format.pval(cor_test$p.value, digits = 3)))
print(paste("95% Confidence interval:", 
           round(cor_test$conf.int[1], 4), "to", round(cor_test$conf.int[2], 4)))

