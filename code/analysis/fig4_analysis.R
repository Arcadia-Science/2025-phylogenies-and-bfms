library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils/utils.R")))

# Get list of FASTA files containing multiple sequence alignments
files <- list.files(here('data/example_pfam_families/'))
files <- files[grep("fasta", files)]

print(paste("Total FASTA files found:", length(files)))
print("First 10 MSA files:")
print(head(files, 10))

if (length(files) == 0) {
  stop("No FASTA files found in the current directory")
}

# Initialize storage for MSA Hill's diversity results
msa_hill <- list()

print("Processing MSA files for diversity analysis...")
print("This may take several minutes depending on file sizes")

# Process each MSA file
for (i in 1:length(files)) {
  if (i %% 10 == 0 || i == 1) print(paste("Processing file", i, "of", length(files), ":", files[i]))
  
  # Read multiple sequence alignment from FASTA format
  msa <- readBStringSet(paste(here('data/example_pfam_families/'), files[i], sep = ''))
  
  if (length(msa) > 0) {
    # Convert alignment to matrix format for column-wise analysis
    msa_mat <- as.matrix(msa)
    
    # Set Hill's diversity parameter (q=1 for exponential Shannon entropy)
    q <- 1
    
    # Calculate Hill's diversity for each alignment column
    column_diversities <- apply(msa_mat, 2, hill_diversity, q = q)
    
    # Store both column-wise and mean diversity values
    msa_hill[[files[i]]] <- list(all_hill = column_diversities,
                                mean_hill = mean(column_diversities))
  } else {
    print(paste("Warning: Empty MSA file:", files[i]))
  }
}

print("MSA diversity analysis complete!")
print(paste("Processed", length(msa_hill), "MSA files"))

# Remove file extensions from names for consistency
original_names <- names(msa_hill)[1:min(5, length(msa_hill))]
names(msa_hill) <- gsub(".fasta", "", names(msa_hill))
cleaned_names <- names(msa_hill)[1:min(5, length(msa_hill))]

print("Name cleaning example:")
print("Original names:")
print(original_names)
print("Cleaned names:")
print(cleaned_names)

# Show diversity statistics before normalization
raw_diversities <- unlist(lapply(msa_hill, function(x) x$mean_hill))
print("Raw Hill's diversity statistics:")
print(summary(raw_diversities))

# Normalize diversity values (assuming 20 amino acids max)
norm_hill <- lapply(msa_hill, function(x) normalize_hill(x$mean_hill, 20))

print("Normalized Hill's diversity statistics:")
print(summary(unlist(norm_hill)))

# Remove last element (likely incomplete/problematic)
if (length(norm_hill) > 1) {
  print(paste("Removing last element:", names(norm_hill)[length(norm_hill)]))
  norm_hill <- norm_hill[-length(norm_hill)]
}

print(paste("Final MSA dataset size:", length(norm_hill), "families"))

system('./code/utils/fig4_cluster_fasta_directory.sh -d data/example_pfam_families/ -o data/example_pfam_families/clustering_results/')

# Load MMseqs clustering results
res <- read.csv(here("data/example_pfam_families/clustering_results/clustering_summary.csv"),
               header = TRUE)

print(paste("Clustering data loaded:", nrow(res), "entries"))
print("Column names:")
print(colnames(res))
print("First few rows:")
print(head(res))

# Show unique families
print(paste("Number of protein families:", length(unique(res$family))))
print("Families represented:")
print(unique(res$family)[1:min(10, length(unique(res$family)))])

# Split results by protein family
res <- split(res, res$family)

print(paste("Data split into", length(res), "families"))

# Calculate clustering proportions for each family
for (i in 1:length(res)) {
  res[[i]]$prop <- res[[i]]$num_clusters / res[[i]]$total_sequences
}

print("Clustering proportions calculated")

# Show example proportions
if (length(res) > 0) {
  print("Example clustering proportions for first family:")
  print(res[[1]][c("seq_identity", "num_clusters", "total_sequences", "prop")])
}

# Remove last element for consistency
if (length(res) > 1) {
  print(paste("Removing last element:", names(res)[length(res)]))
  res <- res[-length(res)]
}

print(paste("Final clustering dataset size:", length(res), "families"))

# Show proportion statistics across all families
all_props <- unlist(lapply(res, function(x) x$prop))
print("Clustering proportion statistics across all families:")
print(summary(all_props))

# Calculate cumulative clustering proportions
props <- list()
for (i in 1:length(res)) {
  props[[i]] <- sort(res[[i]]$num_clusters)
}
props <- do.call(rbind, props)
props <- colSums(props)
props <- props / props[length(props)]

print("Cumulative clustering proportions:")
print(props)
print(paste("Range:", round(min(props), 3), "to", round(max(props), 3)))

# Calculate slopes of clustering trends vs sequence identity
slopes <- lapply(res, function(x) {
  if (nrow(x) > 1 && var(x$seq_identity) > 0) {
    coef(lm(x$prop ~ x$seq_identity))[2]
  } else {
    NA
  }
})

# Remove NA values
slopes <- slopes[!is.na(slopes)]

print("Slope analysis results:")
print(paste("Number of families with valid slopes:", length(slopes)))
print("Slope statistics:")
print(summary(unlist(slopes)))

# Display sorted slopes
sorted_slopes <- sort(unlist(slopes))
print("Sorted slopes (first 10):")
print(head(sorted_slopes, 10))

# Match MSA diversity data with clustering results
print("Matching datasets...")
print(paste("MSA families:", length(names(norm_hill))))
print(paste("Clustering families:", length(names(res))))

# Find common families
common_families <- intersect(names(norm_hill), names(res))
print(paste("Common families:", length(common_families)))

if (length(common_families) > 0) {
  print("First few common families:")
  print(head(common_families))
} else {
  print("Warning: No common families found between datasets")
}

# Extract data for matching families
x <- res[match(names(norm_hill), names(res))]
x <- unlist(lapply(x, function(y) {
  if (!is.null(y) && nrow(y) >= 9) {
    y$prop[9]  # Extract proportion at 9th identity level
  } else {
    NA
  }
}))

z <- unlist(slopes[match(names(norm_hill), names(slopes))])
y <- unlist(norm_hill)

print("Extracted data summary:")
print(paste("Clustering proportions (x):", sum(!is.na(x)), "valid values"))
print(paste("Slopes (z):", sum(!is.na(z)), "valid values"))
print(paste("Hill diversities (y):", sum(!is.na(y)), "valid values"))

# Show complete cases
complete_cases <- sum(!is.na(x) & !is.na(z) & !is.na(y))
print(paste("Complete cases for analysis:", complete_cases))

# Prepare data for cumulative clustering plot
df <- data.frame(sim = seq(0.1, 0.9, 0.1), props = props)

p1 <- ggplot(data = df,
            aes(x = sim, y = props)) +
  xlab("Sequence similarity") +
  ylab("Proportion of clusters") +
  geom_line() +
  theme(legend.position = "none") +
  ylim(c(0, 1)) +
  theme_arcadia()

print(p1)
print("Plot 1 created: Overall clustering trends")

# Prepare data for individual family clustering curves
df <- data.frame(props = unlist(lapply(res, function(x) c(sort(x$prop), 1))),
                variable = rep(1:length(res), each = length(seq(0.1, 1, 0.1))))
df$sim <- rep(seq(0.1, 1, 0.1), length(res))
df$color <- rep(unlist(lapply(res, function(x) min(x$prop))),
               each = length(seq(0.1, 1, 0.1)))

print(paste("Family clustering data prepared:", nrow(df), "points across", length(res), "families"))

p2 <- ggplot(data = df,
            aes(x = sim,
                y = props,
                group = variable,
                color = color)) +
  scale_color_gradientn(colours = arcadia_gradient_palette("verde")$colors) +
  xlab("Sequence similarity") +
  ylab("Proportion of clusters") +
  geom_line(show.legend = FALSE) +
  theme(legend.position = "none") +
  ylim(c(0, 1)) +
  theme_arcadia()

print(p2)
print("Plot 2 created: Family-specific clustering patterns")

# Prepare data for slope vs Hill's diversity correlation
# Use matched data to ensure color values are available
color_values <- unlist(lapply(res[match(names(norm_hill), names(res))], function(x) {
  if (!is.null(x)) {
    min(x$prop)
  } else {
    NA
  }
}))

df <- data.frame(slopes = z,
                hill = y,
                color = color_values)

# Remove rows with missing data
df <- df[complete.cases(df), ]

print(paste("Correlation data prepared:", nrow(df), "complete cases"))

if (nrow(df) > 1) {
  p3 <- ggplot(data = df,
              aes(x = slopes,
                  y = hill,
                  color = color)) +
    scale_color_gradientn(colours = arcadia_gradient_palette("verde")$colors) +
    xlab("Slope") +
    ylab("Hill's diversity (normalized)") +
    geom_point(size = 2, show.legend = FALSE) +
    theme(legend.position = "none") +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    theme_arcadia()
  
  print(p3)
  print("Plot 3 created: Diversity vs clustering relationship")
  
  # Calculate correlation
  correlation <- cor(df$slopes, df$hill, use = "complete.obs")
  print(paste("Correlation between slopes and Hill's diversity:", round(correlation, 4)))
} else {
  print("Insufficient data for correlation plot")
}

# Combine and display plots
if (exists("p2") && exists("p3")) {
  grid.arrange(p2, p3, nrow = 1)
  print("Combined plots displayed")
} else {
  print("Not all plots available for combination")
}

