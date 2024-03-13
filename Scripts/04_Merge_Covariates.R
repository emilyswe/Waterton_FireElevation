library(dplyr)
library(readr)

# Load datasets
severity_data <- read_csv("Input/03_Standardized_Severity.csv")
terrain_data <- read_csv("Input/03_Standardized_Terrain_GEE.csv")
landcover_data <- read_csv("Input/03_Standardized_LandCover_2016_2019_GEE.csv")
coordinates_data <- read_csv("Input/02_Standardized.csv") # Load the coordinates data

# Ensure the 'location' column is consistently named across datasets
terrain_data <- rename(terrain_data, location = locatin)
landcover_data <- rename(landcover_data, location = locatin)

# Keep only the first row of each duplicate entry based on 'location' and 'year'
severity_data_unique <- severity_data %>%
  distinct(location, year, .keep_all = TRUE)

terrain_data_unique <- terrain_data %>%
  distinct(location, year, .keep_all = TRUE)

landcover_data_unique <- landcover_data %>%
  distinct(location, year, .keep_all = TRUE)

# Now join the datasets on 'location' and 'year'
final_data <- severity_data_unique %>%
  inner_join(terrain_data_unique, by = c("location", "year")) %>%
  inner_join(landcover_data_unique, by = c("location", "year"))

# Prepare the coordinates_data by ensuring it only contains unique rows for each 'location'
coordinates_data_unique <- coordinates_data %>%
  distinct(location, .keep_all = TRUE)

# Join the coordinates from the coordinates_data_unique
final_data_with_coordinates <- final_data %>%
  inner_join(coordinates_data_unique %>% select(location, latitude, longitude), by = "location")

# Write the final dataset with coordinates to a new CSV file
write_csv(final_data_with_coordinates, "Input/04_Cleaned_Waterton.csv")
