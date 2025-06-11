library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils/utils.R")))

# Load cox1 tree from the specified data source
# Tree from: https://www.sciencedirect.com/org/science/article/pii/S2534970821000028#aep-e-component-id12
cox1 <- read.tree(here("data/cox1.tre"))

# Display basic tree information
print(paste("Number of tips:", Ntip(cox1)))
print(paste("Number of nodes:", Nnode(cox1)))

# Extract major plant and animal phyla
plant_cox1 <- keep.tip(cox1, cox1$tip.label[grep("Streptophyta", cox1$tip.label)])
nonplant_cox1 <- keep.tip(cox1, cox1$tip.label[grep("Chordata", cox1$tip.label)])

print(paste("Plant (Streptophyta) tips:", Ntip(plant_cox1)))
print(paste("Animal (Chordata) tips:", Ntip(nonplant_cox1)))

# Remove outlier branches based on 95th percentile threshold
L <- quantile(nonplant_cox1$edge.length, 0.95)
tip_edges <- which(nonplant_cox1$edge[, 2] <= Ntip(nonplant_cox1))
too_long_edges <- tip_edges[nonplant_cox1$edge.length[tip_edges] > L]
tips_2_drop <- nonplant_cox1$tip.label[nonplant_cox1$edge[too_long_edges, 2]]
nonplant_cox1 <- drop.tip(nonplant_cox1, tips_2_drop)

print(paste("Removed", length(tips_2_drop), "outlier tips from animal tree"))
print(paste("Final animal tree tips:", Ntip(nonplant_cox1)))

# Remove outlier branches based on 95th percentile threshold
L <- quantile(plant_cox1$edge.length, 0.95)
tip_edges <- which(plant_cox1$edge[, 2] <= Ntip(plant_cox1))
too_long_edges <- tip_edges[plant_cox1$edge.length[tip_edges] > L]
tips_2_drop <- plant_cox1$tip.label[plant_cox1$edge[too_long_edges, 2]]
plant_cox1 <- drop.tip(plant_cox1, tips_2_drop)

print(paste("Removed", length(tips_2_drop), "outlier tips from plant tree"))
print(paste("Final plant tree tips:", Ntip(plant_cox1)))

# Set up plotting parameters for side-by-side plots
par(mfrow = c(1, 2),
    mar = c(1, 1, 1, 1))

# Plot animal tree (Chordata) in pink
plot(nonplant_cox1,
     show.tip.label = FALSE,
     edge.color = "#F898AE",
     edge.width = 2,
     type = "fan")
title("Animal (Chordata) Cox1 Tree", line = -1)

# Plot plant tree (Streptophyta) in green
plot(plant_cox1,
     show.tip.label = FALSE,
     edge.color = "#3B9886",
     edge.width = 2,
     type = "fan")
title("Plant (Streptophyta) Cox1 Tree", line = -1)

