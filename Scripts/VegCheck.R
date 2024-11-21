# Step 1: Set working directory and load libraries
setwd("~/Documents/local-git/Waterton_FireElevation")  
library(sf)  # For spatial data manipulation
library(dplyr)  # For data manipulation
library(units)  # For handling distance units

# Step 2: Load input data
elc <- st_read("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Input/Veg_Plots_Corrected_2013.shp")
locations <- read.csv("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Input/Locations.csv")

# Step 3: Convert locations to spatial data and match CRS
locations_sf <- st_as_sf(locations, coords = c("longitude", "latitude"), crs = 4326)  # WGS84 (EPSG:4326)
locations_sf <- st_transform(locations_sf, crs = st_crs(elc))  # Match the CRS of `elc`

# Step 4: Filter ELC data
elc_filtered <- elc %>% filter(!is.na(Classified))  # Remove rows where Classified is NA

# Step 5: Create buffers for ELC and locations
elc_filtered$buffer <- st_buffer(elc_filtered$geometry, dist = 1000)  # 1000 m buffer for ELC points
locations_sf$buffer <- st_buffer(locations_sf$geometry, dist = 300)  # 300 m buffer for locations

# Step 6: Find intersections
intersections <- st_intersects(locations_sf$buffer, elc_filtered$buffer)

# Step 7: Assign values based on intersections
locations_sf$Classified <- NA
locations_sf$nearest_elevation <- NA

for (i in seq_along(intersections)) {
  if (length(intersections[[i]]) > 0) {
    # Get all overlapping ELC indices
    elc_indices <- intersections[[i]]
    
    # Extract `Classified` values for these indices
    classified_values <- elc_filtered$Classified[elc_indices]
    
    # Assign the most frequent Classified value to the location
    locations_sf$Classified[i] <- names(sort(table(classified_values), decreasing = TRUE))[1]
    
    # Assign the corresponding elevation value (mean of overlapping elevations)
    locations_sf$nearest_elevation[i] <- mean(elc_filtered$Elevation[elc_indices], na.rm = TRUE)
  }
}

# Step 8: Handle unmatched locations with nearest neighbors
unmatched <- which(is.na(locations_sf$Classified))
cat("Number of unmatched rows before loop:", length(unmatched), "\n")

# First pass: Assign values within 1000 m
if (length(unmatched) > 0) {
  for (i in unmatched) {
    # Find the nearest ELC feature
    nearest_idx <- st_nearest_feature(locations_sf[i, ], elc_filtered)
    nearest_distance <- st_distance(locations_sf[i, ], elc_filtered[nearest_idx, ], by_element = TRUE)
    threshold <- units::set_units(1000, "m")
    
    # Debugging: Print nearest index and distance
    cat("Processing unmatched row:", i, "\n")
    cat("Nearest index:", nearest_idx, "\n")
    cat("Nearest distance:", nearest_distance, "\n")
    
    # Assign values if within the threshold
    if (nearest_distance <= threshold) {
      locations_sf$Classified[i] <- elc_filtered$Classified[nearest_idx]
      locations_sf$nearest_elevation[i] <- elc_filtered$Elevation[nearest_idx]
    }
  }
}

# Second pass: Relax threshold for remaining unmatched rows
unmatched <- which(is.na(locations_sf$Classified))
cat("Number of unmatched rows after first pass:", length(unmatched), "\n")

if (length(unmatched) > 0) {
  for (i in unmatched) {
    # Find the nearest ELC feature
    nearest_idx <- st_nearest_feature(locations_sf[i, ], elc_filtered)
    nearest_distance <- st_distance(locations_sf[i, ], elc_filtered[nearest_idx, ], by_element = TRUE)
    threshold <- units::set_units(3000, "m")  # Increased threshold for second pass
    
    # Debugging: Print nearest index and distance
    cat("Processing unmatched row (second pass):", i, "\n")
    cat("Nearest index:", nearest_idx, "\n")
    cat("Nearest distance:", nearest_distance, "\n")
    
    # Assign values if within the relaxed threshold
    if (nearest_distance <= threshold) {
      locations_sf$Classified[i] <- elc_filtered$Classified[nearest_idx]
      locations_sf$nearest_elevation[i] <- elc_filtered$Elevation[nearest_idx]
    }
  }
}

# Step 9: Add latitude and longitude columns
locations_sf$latitude <- st_coordinates(locations_sf)[, 2]
locations_sf$longitude <- st_coordinates(locations_sf)[, 1]

# Step 10: Resolve and validate final results
final_results <- locations_sf %>%
  st_drop_geometry() %>%
  select(location, latitude, longitude, elevation, Classified, nearest_elevation)

# Save to CSV
write.csv(final_results, "/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Output/VegCheckResults.csv", row.names = FALSE)

# Debugging: Check number of rows with NA in final results
cat("Number of rows with NA in Classified:", sum(is.na(final_results$Classified)), "\n")


####################collinearity test 


library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Read in the final results CSV
final_results <- read.csv("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Output/VegCheckResults.csv")

# Step 2: Bin elevation into intervals
bin_width <- 100  # Define the bin size (100 meters in this example)
final_results <- final_results %>%
  mutate(elevation_bin = cut(nearest_elevation, 
                             breaks = seq(floor(min(nearest_elevation, na.rm = TRUE)),
                                          ceiling(max(nearest_elevation, na.rm = TRUE)),
                                          by = bin_width),
                             include.lowest = TRUE,
                             labels = FALSE))  # Assign numeric labels to bins

# Step 3: Calculate percentages of vegetation types in each bin
veg_proportions <- final_results %>%
  filter(!is.na(Classified)) %>%  # Remove rows with NA vegetation types
  group_by(elevation_bin, Classified) %>% 
  summarise(count = n(), .groups = "drop") %>%  # Count occurrences of each vegetation type per bin
  group_by(elevation_bin) %>%
  mutate(percentage = count / sum(count) * 100)  # Calculate percentages

# Step 4: Create a wide format dataset for correlation analysis
veg_matrix <- veg_proportions %>%
  select(elevation_bin, Classified, percentage) %>%
  pivot_wider(names_from = Classified, values_from = percentage, values_fill = 0)  # Fill missing combinations with 0

# Step 5: Calculate correlations with elevation_bin
# Add numeric elevation bin values for correlation
veg_matrix <- veg_matrix %>%
  mutate(elevation_bin_numeric = as.numeric(elevation_bin))  # Ensure elevation bin is numeric

# Remove the elevation_bin column for correlation
correlation_data <- veg_matrix %>%
  select(-elevation_bin)

# Calculate correlation matrix
correlation_matrix <- cor(correlation_data, use = "complete.obs")

# Extract correlations with elevation_bin
elevation_correlations <- correlation_matrix["elevation_bin_numeric", -1]  # Skip the diagonal element

# Step 6: Check for strong collinearity
strong_correlations <- elevation_correlations[abs(elevation_correlations) > 0.7]

cat("Correlations between elevation bins and vegetation types:\n")
print(elevation_correlations)

if (length(strong_correlations) > 0) {
  cat("\nStrong collinearity detected (R > 0.7):\n")
  print(strong_correlations)
} else {
  cat("\nNo strong collinearity detected (R > 0.7).\n")
}

# Step 7: Optional: Visualize the data
ggplot(veg_proportions, aes(x = factor(elevation_bin), y = percentage, fill = Classified)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Elevation Bin", y = "Percentage", fill = "Vegetation Type") +
  theme_minimal() +
  ggtitle("Vegetation Percentages Across Elevation Bins")
















############################ extra checks

# Check if all geometries are valid
all(st_is_valid(locations_sf))  # Should return TRUE
all(st_is_valid(elc_filtered))  # Should return TRUE

# Plot all locations
plot(st_geometry(locations_sf), col = "red", main = "Unmatched Locations and ELC Features")
plot(st_geometry(elc_filtered), col = "blue", add = TRUE)

# Highlight unmatched locations
unmatched_locations <- locations_sf[remaining_na, ]
plot(st_geometry(unmatched_locations), col = "yellow", add = TRUE)



# Repair invalid geometries if needed
locations_sf <- st_make_valid(locations_sf)
elc_filtered <- st_make_valid(elc_filtered)

