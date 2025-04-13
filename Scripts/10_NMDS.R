# NMDS and PERMANOVA using single visit type A dataset and bray-curtis 
# Updated March 25 2025

# 1. Load Required Packages ----
library(dplyr)
library(vegan)
library(ggplot2)
library(readr)

# 2. Load Filtered Dataset ----
emda4 <- read_csv("Output/07_cleaned_single_visit_filtered.csv")

# 3. Add Elevation × Burn Strata ----
if (!"strata" %in% colnames(emda4)) {
  location_data <- emda4 |>
    select(location, latitude, longitude, elevation, grid.code) |>
    distinct()
  
  location_data$altbin <- cut(location_data$elevation,
                              breaks = quantile(location_data$elevation, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                              labels = FALSE, include.lowest = TRUE)
  
  location_data$strata <- ifelse(location_data$grid.code > 0,
                                 location_data$altbin + 3,
                                 location_data$altbin)
  
  emda4 <- left_join(emda4, location_data[, c("location", "strata")], by = "location")
}

# 4. Prepare Data for NMDS ----
species_columns <- grep("^[A-Z]{4}$", colnames(emda4), value = TRUE)

nmds_data <- emda4 |>
  filter(!is.na(strata)) |>
  select(location, strata, latitude, longitude, elevation, all_of(species_columns))

species_matrix <- nmds_data[, species_columns]

# Remove rows with no species present
non_zero_rows <- rowSums(species_matrix) > 0
species_matrix <- species_matrix[non_zero_rows, , drop = FALSE]
nmds_data <- nmds_data[non_zero_rows, ]

# 5. Run NMDS ----
set.seed(999)
# Use Bray-Curtis for abundance data (instead of Jaccard which is for presence/absence)
nmds_result <- metaMDS(species_matrix, distance = "bray", k = 2, trymax = 100)

# Extract scores
nmds_coords <- scores(nmds_result, display = "sites") |> as.data.frame()
nmds_data$mdsA <- nmds_coords$NMDS1
nmds_data$mdsB <- nmds_coords$NMDS2

# Save scores
write_csv(nmds_data[, c("location", "strata", "mdsA", "mdsB")], "Output/10_NMDS/waterton_nmds_2axis_EMDA.csv")

# 6. NMDS Plot (Elly Style) ----
dat3 <- read_csv("Output/10_NMDS/waterton_nmds_2axis_EMDA.csv") |>
  mutate(
    treatment = ifelse(strata %in% c(1:3), "Unburned", "Burned"),
    elevation = case_when(
      strata %in% c(1, 4) ~ "Low elevation",
      strata %in% c(2, 5) ~ "Mid elevation",
      strata %in% c(3, 6) ~ "High elevation"
    ),
    elevation = factor(elevation, levels = c("Low elevation", "Mid elevation", "High elevation"))
  )

plot3 <- ggplot(dat3, aes(x = mdsA, y = mdsB)) +
  stat_ellipse(aes(colour = treatment, linetype = elevation), level = 0.68, linewidth = 1) +
  scale_color_manual(values = c("Burned" = "tomato3", "Unburned" = "steelblue3"), name = "") +
  scale_linetype_manual(values = c("solid", "dotted", "dashed"), name = "") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "NMDS axis 1", y = "NMDS axis 2")

ggsave("Output/10_NMDS/10_nmds_plot_final.png", plot = plot3, width = 8, height = 6, dpi = 300, bg = "white")

# 7. PERMANOVA ----
env_data <- nmds_data |>
  mutate(
    treatment = ifelse(strata %in% c(1:3), "Unburned", "Burned"),
    elevation = case_when(
      strata %in% c(1, 4) ~ "Low elevation",
      strata %in% c(2, 5) ~ "Mid elevation",
      strata %in% c(3, 6) ~ "High elevation"
    )
  ) |>
  select(treatment, elevation)

set.seed(999)
# Again, use Bray-Curtis here to match the NMDS method
permanova_result <- adonis2(species_matrix ~ treatment * elevation,
                            data = env_data,
                            method = "bray",
                            permutations = 999)

write_csv(as.data.frame(permanova_result), "Output/10_NMDS/permanova_results.csv")
print(permanova_result)


##get delta AIC values 

# Full model with interaction
full_model <- adonis2(species_matrix ~ treatment * elevation, data = env_data, method = "bray", permutations = 999)

# Additive model
additive_model <- adonis2(species_matrix ~ treatment + elevation, data = env_data, method = "bray", permutations = 999)

# Elevation only
elevation_model <- adonis2(species_matrix ~ elevation, data = env_data, method = "bray", permutations = 999)

# Treatment only
treatment_model <- adonis2(species_matrix ~ treatment, data = env_data, method = "bray", permutations = 999)

# Extract AIC values and calculate ΔAIC manually (pseudo-AIC based on residuals and df)
aics <- c(
  full_model = full_model$aic,
  additive_model = additive_model$aic,
  elevation_model = elevation_model$aic,
  treatment_model = treatment_model$aic
)
delta_aic <- aics - min(aics)

# Save delta AIC results
write.csv(data.frame(Model = names(delta_aic), AIC = aics, Delta_AIC = delta_aic), 
          "Output/10_NMDS/permanova_model_comparisons.csv", row.names = FALSE)

















####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
#old code for nmds and permanova refering to both the multi and single visit dfs 

# NMDS Analysis Script - Step 10

#re-run without removing the zeros but use Bray-Curtis instead (Feb 5 2025 note)

# Set working directory
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# Load required packages
packs <- c("dplyr", "vegan", "ggplot2")
for (q in 1:length(packs)) {
  if (!require(packs[q], character.only = TRUE)) {
    install.packages(packs[q])
    require(packs[q])
  }
}

# Read in the filtered single visit dataset
emda4 <- read.csv("Output/07_cleaned_single_visit_filtered.csv")

# Add strata column if not already present
if (!"strata" %in% colnames(emda4)) {
  location_data <- data.frame(location = unique(emda4$location))
  
  for (q in 1:nrow(location_data)) {
    location_data$latitude[q] <- emda4$latitude[emda4$location == location_data$location[q]][1]
    location_data$longitude[q] <- emda4$longitude[emda4$location == location_data$location[q]][1]
    location_data$altitude[q] <- emda4$elevation[emda4$location == location_data$location[q]][1]
    location_data$severity[q] <- emda4$grid.code[emda4$location == location_data$location[q]][1]
  }
  
  location_data$altbin <- cut(location_data$altitude, 
                              breaks = quantile(location_data$altitude, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                              labels = FALSE, include.lowest = TRUE)
  location_data$strata <- ifelse(location_data$severity > 0, 
                                 location_data$altbin + 3, 
                                 location_data$altbin)
  
  emda4 <- left_join(emda4, location_data[, c("location", "strata")], by = "location")
}

# Define species columns
species_columns <- grep("^[A-Z]{4}$", colnames(emda4), value = TRUE)

# Create nmds_data and species_matrix
nmds_data <- emda4 %>%
  filter(!is.na(strata)) %>%
  select(location, strata, latitude, longitude, elevation, all_of(species_columns))

species_matrix <- nmds_data[, species_columns]

# Remove rows with no species present (zeros) from both datasets
non_zero_rows <- rowSums(species_matrix) > 0
species_matrix <- species_matrix[non_zero_rows, , drop = FALSE]  # Ensure matrix format
nmds_data <- nmds_data[non_zero_rows, ]

# Verify dimensions match
if (nrow(species_matrix) != nrow(nmds_data)) {
  stop("Mismatch between species_matrix and nmds_data rows after filtering.")
}

# Run NMDS
set.seed(999)
nmds_result <- metaMDS(species_matrix, distance = "jaccard", k = 2, trymax = 100)

# Extract NMDS scores and add to nmds_data
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_data$mdsA <- nmds_coords[, 1]
nmds_data$mdsB <- nmds_coords[, 2]

# Save the NMDS result for plotting
write.csv(nmds_data[, c("location", "strata", "mdsA", "mdsB")], 
          "Output/10_NMDS/waterton_nmds_2axis_EMDA.csv", row.names = FALSE)

#######plot

# Read in the NMDS result CSV
dat3 <- read.csv("Output/10_NMDS/waterton_nmds_2axis_EMDA.csv") |> 
  mutate(treatment = ifelse(strata %in% c(1:3), "Unburned", "Burned"),
         elevation = case_when(strata %in% c(1, 4) ~ "Low elevation", 
                               strata %in% c(2, 5) ~ "Mid elevation",
                               strata %in% c(3, 6) ~ "High elevation"),
         elevation = factor(elevation, levels=c("Low elevation", "Mid elevation", "High elevation")))

# Filter out the outlier AFTER reading the CSV
dat3 <- dat3 %>% filter(location != "WLNP-12-1")

# Check summary to confirm outlier is removed
summary(dat3$mdsA)

# Create the legend data for plot
linelegend <- data.frame(expand.grid(x = 1, y = 1, elevation = unique(dat3$elevation)))

plot3 <- ggplot(dat3, aes(x = mdsA, y = mdsB)) +
  stat_ellipse(aes(colour = treatment, linetype = elevation), level = 0.68, linewidth = 1) +
  scale_color_manual(values = c("tomato3", "steelblue3"), name = "") +  # Burned = Red, Unburned = Blue
  scale_linetype_manual(values = c("solid", "dotted", "dashed"), name = "") +  # More distinct line types
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "NMDS axis 1",
       y = "NMDS axis 2")
plot3
# Save the plot
ggsave(plot3, filename = "Output/10_NMDS/NMDS_EMDA.jpeg", width = 8, height = 6)


#PERMANOVA

# Ensure species_matrix is a matrix
species_matrix <- as.matrix(species_matrix)

# Prepare environmental variables for PERMANOVA
env_data <- nmds_data %>%
  mutate(
    treatment = ifelse(strata %in% c(1, 2, 3), "Unburned", "Burned"),
    elevation = case_when(
      strata %in% c(1, 4) ~ "Low elevation",
      strata %in% c(2, 5) ~ "Mid elevation",
      strata %in% c(3, 6) ~ "High elevation"
    )
  ) %>%
  select(treatment, elevation)  # Select only relevant columns for PERMANOVA

# Check dimensions match before PERMANOVA
if (nrow(species_matrix) != nrow(env_data)) {
  stop("Mismatch between species_matrix and env_data rows. Check filtering steps.")
}

# Run PERMANOVA
set.seed(999)
permanova_result <- adonis2(
  species_matrix ~ treatment * elevation, 
  data = env_data, 
  method = "jaccard", 
  permutations = 999
)

# Display PERMANOVA results
print(permanova_result)

# Save PERMANOVA results to file
write.csv(as.data.frame(permanova_result), "Output/10_NMDS/permanova_results.csv", row.names = TRUE)
