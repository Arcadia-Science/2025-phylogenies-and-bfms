library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils/utils.R")))

system('mkdir data/trees/')
system('wget -P data/ https://ftp.ensembl.org/pub/release-114/emf/ensembl-compara/homologies/Compara.114.protein_default.nh.emf.gz')

system('python3 code/utils/fig2_extract_trees.py --emf_file data/Compara.114.protein_default.nh.emf.gz --output_dir data/trees/')

# Set working directory
setwd(here('data/trees/'))

# List all tree files (exclude taxonomy files)
files <- list.files()

print(paste("Total tree files found:", length(files)))
print("First 10 files:")
print(head(files, 10))

# Initialize storage for Hill diversity values
hill_values <- list()

# Progress bar for tree processing
pb <- txtProgressBar(
  min = 1,
  max = length(files),
  style = 3,
  width = 100,
  char = "."
)

print("Processing trees...")

# Calculate Hill diversity for each tree
for (i in 1:length(files)) {
  # Update progress indicator
  setTxtProgressBar(pb, i)
  
  # Load phylogenetic tree
  tree <- ape::read.tree(files[i])
  
  # Calculate normalized Hill diversity metrics
  hill_values[[as.character(files[i])]] <- normalize_hill_q1(tree)
}

close(pb)
print("\nProcessing complete!")
print(paste("Processed", length(hill_values), "trees"))

#Calculate sum of effective sequence count as a proportion of total
sum_eff = sum(unlist(lapply(hill_values, function(x) x$Hill_diversity)))
sum_total = sum(unlist(lapply(hill_values, function(x) length(x$tree$tip.label))))

print(sum_eff)
print(sum_total)
print(sum_eff/sum_total)

# Filter trees with at least 100 tips for statistical robustness
hill_values_filter <- hill_values[lapply(hill_values, function(x) length(x$tree$tip.label)) >= 100]

print(paste("Trees before filtering:", length(hill_values)))
print(paste("Trees after filtering (≥100 tips):", length(hill_values_filter)))

# Show distribution of tree sizes
tree_sizes <- unlist(lapply(hill_values_filter, function(x) length(x$tree$tip.label)))
print(paste("Tree size range:", min(tree_sizes), "to", max(tree_sizes), "tips"))
print(paste("Median tree size:", median(tree_sizes), "tips"))

# Prepare data for diversity plot
val <- data.frame(x = 1:length(hill_values_filter),
                 y = sort(unlist(lapply(hill_values_filter, function(x) x$Norm_by_tips))))

print("Hill's diversity statistics:")
print(summary(val$y))

# Add row names for tracking
rownames(val) <- names(sort(unlist(lapply(hill_values_filter, function(x) x$Norm_by_tips))))

# Calculate coefficient of variation for branch lengths
br_variance <- list()

print("Calculating branch length variance...")
for (i in 1:length(hill_values_filter)) {
  if (is.null(hill_values_filter[[i]]$tree$edge.length) == FALSE) {
    
    # Extract tree for analysis
    tree <- hill_values_filter[[i]]$tree
    
    # Get all branch lengths (simplified from terminal-only)
    br <- tree$edge.length
    
    # Calculate coefficient of variation (CV = sd/mean)
    br_variance[[names(hill_values_filter)[i]]] <- sd(br) / mean(br)
  } else {
    br_variance[[names(hill_values_filter)[i]]] <- NA
  }
}

print("Branch length variance calculation complete!")

# Match branch variance order to Hill's diversity values
br_variance <- unlist(br_variance[match(rownames(val), names(br_variance))])
y <- val$y[match(names(br_variance), rownames(val))]

# Prepare data for variance vs diversity plot
val2 <- data.frame(x = y, y = br_variance)

# Remove NA values
val2 <- val2[complete.cases(val2), ]

print("Final dataset statistics:")
print(paste("Families with diversity data:", nrow(val)))
print(paste("Families with both diversity and variance data:", nrow(val2)))
print("Variance statistics:")
print(summary(val2$y))

# Load required plotting libraries
library(ggplot2)
library(gridExtra)

# Plot 1: Hill's diversity distribution across protein families
p1 <- ggplot(data = val,
            aes(x = x,
                y = y,
                color = y)) +
  scale_color_gradientn(colours = arcadia_gradient_palette("magma")$colors) +
  xlab("Protein family #") +
  ylab("Hills diversity (normalized)") +
  geom_point(size = 2, show.legend = FALSE) +
  theme(legend.position = "none") +
  theme_arcadia()

print("Plot 1 created: Hill's diversity distribution")

# Plot 2: Relationship between diversity and branch length variance
p2 <- ggplot(data = val2,
            aes(x = x,
                y = y,
                color = x)) +
  scale_color_gradientn(colours = arcadia_gradient_palette("magma")$colors) +
  xlab("Hills diversity (normalized)") +
  ylab("Variance (mean normalized)") +
  geom_point(size = 2, show.legend = FALSE) +
  theme(legend.position = "none") +
  theme_arcadia()

print("Plot 2 created: Diversity vs variance relationship")

# Display combined plots
grid.arrange(p1, p2, nrow = 1)

print("Combined plots displayed")

# Plot representative trees at different diversity levels
eff_counts <- sort(unlist(lapply(hill_values_filter, function(x) x$Norm_by_tips)))

print("Diversity levels for example trees:")
print(paste("Minimum:", round(min(eff_counts), 3)))
print(paste("0.2 level:", round(eff_counts[which.min(abs(eff_counts - 0.2))], 3)))
print(paste("0.4 level:", round(eff_counts[which.min(abs(eff_counts - 0.4))], 3)))
print(paste("0.6 level:", round(eff_counts[which.min(abs(eff_counts - 0.6))], 3)))
print(paste("Maximum:", round(max(eff_counts), 3)))

# Set up panel layout for example trees
par(mfrow = c(1, 5))

# Plot trees at different normalized diversity levels
plot(hill_values_filter[names(eff_counts)[1]][[1]]$tree,
     show.tip.label = FALSE, main = "0")
plot(hill_values_filter[names(eff_counts)[which.min(abs(eff_counts - 0.2))]][[1]]$tree,
     show.tip.label = FALSE, main = "0.2")
plot(hill_values_filter[names(eff_counts)[which.min(abs(eff_counts - 0.4))]][[1]]$tree,
     show.tip.label = FALSE, main = "0.4")
plot(hill_values_filter[names(eff_counts)[which.min(abs(eff_counts - 0.6))]][[1]]$tree,
     show.tip.label = FALSE, main = "0.6")
plot(hill_values_filter[names(eff_counts)[which.min(abs(eff_counts - 1))]][[1]]$tree,
     show.tip.label = FALSE, main = "0.824")

print("Example phylogenies displayed across diversity gradient")

# Generate summary statistics
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total protein families processed:", length(hill_values), "\n")
cat("Families with ≥100 tips:", length(hill_values_filter), "\n")
cat("Families with complete data:", nrow(val2), "\n")
cat("\nHill's Diversity Range:", round(min(val$y), 3), "to", round(max(val$y), 3), "\n")
cat("Branch Length Variance Range:", round(min(val2$y, na.rm=TRUE), 3), "to", round(max(val2$y, na.rm=TRUE), 3), "\n")

# Calculate correlation between diversity and variance
correlation <- cor(val2$x, val2$y, use="complete.obs")
cat("\nCorrelation (diversity vs variance):", round(correlation, 3), "\n")

