# NMDS Analysis Script

# 1. Set working directory ----
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# 2. Load Packages ----
packs <- c("ggplot2", "dplyr", "vegan", "tidyverse", "broom")
for (q in 1:length(packs)) {
  if (!require(packs[q], character.only = TRUE)) {
    install.packages(packs[q])
    require(packs[q])
  }
}

# 3. Read in the multi-visit dataset (superdata) ----
superdata <- read.csv("Output/07_cleaned_multivisit_data.csv")

# 4. Wrangle Data ----
# Identify species columns (4-letter bird codes) using regex to exclude non-species columns
species_columns <- colnames(superdata)[grep("^[A-Z]{4}$", colnames(superdata))]
non_species_columns <- setdiff(colnames(superdata), species_columns)

# Sort species columns alphabetically
species_columns_sorted <- sort(species_columns)

# Arrange the DataFrame so species columns are sorted alphabetically
wide_data <- superdata[, c(non_species_columns, species_columns_sorted)]

# Calculate total richness across all species
wide_data$richness <- rowSums(wide_data[, species_columns_sorted], na.rm = TRUE)

# Sample one observation per location
filtered_data <- wide_data %>%
  dplyr::filter(!is.na(elevation)) %>%  # Assuming 'strata' means elevation strata
  group_by(location) %>%
  slice_sample(n = 1)

# 5. Model: Run NMDS using Jaccard distance ----
set.seed(999)  # For reproducibility
nmds_result <- metaMDS(filtered_data[, species_columns_sorted], 
                       distance = "jaccard", 
                       k = 2, 
                       trymax = 100)  # Increase trymax if convergence issues occur

# 6. Save Results ----
# Extract NMDS coordinates
coordinates <- scores(nmds_result, display = "sites")
filtered_data$mdsA <- coordinates[, 1]
filtered_data$mdsB <- coordinates[, 2]

# Save the NMDS results to a CSV
write.csv(filtered_data[, c("location", "elevation", "doy", "year", "mdsA", "mdsB")], 
          "Output/waterton_nmds_2axis_EMDA.csv", row.names = FALSE)

# 7. Plot Results ----
# Calculate means for each strata (elevation bins)
means <- aggregate(cbind(mdsA, mdsB) ~ elevation, data = filtered_data, mean)

# Prepare NMDS plot
nmdsplot <- ggplot(filtered_data, aes(x = mdsA, y = mdsB, color = factor(elevation))) +
  geom_point(aes(shape = factor(elevation))) +  # Different shapes for each strata
  scale_color_manual(values = c("green", "green", "green", "orange", "orange", "orange")) +  # Color strata
  geom_text(data = means, aes(label = factor(elevation)), vjust = 2, color = "black") +  # Label means
  # Add ellipses for each strata group
  stat_ellipse(aes(fill = factor(elevation)), geom = "polygon", level = 0.95, alpha = 0.2) +
  labs(color = "Elevation Strata", shape = "Elevation Strata", fill = "Confidence") +
  theme_minimal() +
  theme(legend.position = "right") +
  xlim(c(-3, 0))

# Save the plot
ggsave("Output/waterton_nmds_plot.png", plot = nmdsplot, width = 8, height = 6)

# Show the plot in RStudio
print(nmdsplot)
