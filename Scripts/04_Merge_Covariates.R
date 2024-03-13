library(dplyr)
library(readr)

# Load datasets
severity_data <- read_csv("Input/03_Standardized_Severity.csv")
terrain_data <- read_csv("Input/03_Standardized_Terrain_GEE.csv")
landcover_data <- read_csv("Input/03_Standardized_LandCover_2016_2019_GEE.csv")

# Ensure the 'location' column is consistently named across datasets.
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
  inner_join(landcover_data_unique, by = c("location", "year")) %>%
  select(everything()) # Adjust the select statement as per your specific column requirements

# Write the final dataset to a new CSV file
write_csv(final_data, "Input/04_Cleaned_Waterton.csv")




#######checked for duplicates 
# Check for duplicate 'location' and 'year' combinations for severity
duplicates_severity <- severity_data %>%
  group_by(location, year) %>%
  filter(n() > 1) %>%
  ungroup()

# If duplicates exist, they will be shown here
print(duplicates_severity)

# Check for duplicate 'location' and 'year' combinations for landcover
duplicates_landcover <- landcover_data %>%
  group_by(location, year) %>%
  filter(n() > 1) %>%
  ungroup()

# If duplicates exist, they will be shown here
print(duplicates_landcover)


# Check for duplicate 'location' and 'year' combinations in the terrain_data
duplicates_terrain <- terrain_data %>%
  group_by(location, year) %>%
  filter(n() > 1) %>%
  ungroup()

# If duplicates exist, they will be shown here
print(duplicates_terrain)

