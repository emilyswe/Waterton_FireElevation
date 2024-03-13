
library(sf)
library(dplyr)
library(readr)

# Load point count locations
point_counts <- read_csv("Input/02_Standardized.csv")

# Convert to an sf object
point_counts_sf <- st_as_sf(point_counts, coords = c("longitude", "latitude"), crs = 4326)

# Load the Kenow Wildfire severity shapefile
kenow_severity <- st_read("Input/Kenow 2017 Burn severity/Kenow_severity_classes.shp")


# Check the CRS of both datasets
crs_points <- st_crs(point_counts_sf)
crs_severity <- st_crs(kenow_severity)

print(crs_points)
print(crs_severity)

# Transform the CRS of the point dataset to match the shapefile
point_counts_sf <- st_transform(point_counts_sf, crs = st_crs(kenow_severity))

# Now try the spatial join
point_severity <- st_join(point_counts_sf, kenow_severity, join = st_within)

# Proceed with integrating the severity data into your original dataset...

# Spatial join to find fire severity for each point
point_severity <- st_join(point_counts_sf, kenow_severity, join = st_within)

# Convert back to a data frame and select relevant columns, including the new fire severity column
final_data <- as.data.frame(point_severity) %>%
  select(-geometry, everything()) # Adjust to include only necessary columns

#change NA values for gridcode to 0 
final_data <- final_data %>%
  mutate(`gridcode` = ifelse(is.na(`gridcode`), 0, `gridcode`))

# Export the updated dataset
write_csv(final_data, "Input/03_Standardized_Severity.csv")
