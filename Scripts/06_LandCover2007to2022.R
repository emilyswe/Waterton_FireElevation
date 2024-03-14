# Load necessary libraries
library(dplyr)
library(readr)

# Load the datasets
waterton_data <- read_csv("Input/04_Cleaned_Waterton.csv")
land_cover_data <- read_csv("Input/LandCover_2018.csv")

# Correct the column name from 'locatin' to 'location' in the land cover dataset
land_cover_data <- rename(land_cover_data, location = locatin)

# Ensure there are no duplicate entries for each location in the land cover dataset
unique_land_cover_data <- distinct(land_cover_data, location, .keep_all = TRUE)

# Merge the waterton data with the unique land cover data
merged_data <- left_join(waterton_data, unique_land_cover_data[, c("location", "landCover_2018")], by = "location")

# You can then write this merged data to a new CSV file if desired
write_csv(merged_data, "Input/04_Cleaned_Waterton.csv")

###############summarize landcover data per site per year 

# Load your dataset
data <- read_csv('Input/04_Cleaned_Waterton.csv')

# Create a lookup table for land cover codes and their corresponding category names
landcover_lookup <- data.frame(
  code = c(0, 20, 31, 32, 33, 40, 50, 80, 81, 100, 210, 220, 230),
  category = c("Unclassified", "Water", "Snow/Ice", "Rock/Rubble", "Exposed/Barren Land", 
               "Bryoids", "Shrubs", "Wetland", "Wetland Treed", "Herbs", 
               "Coniferous", "Broadleaf", "Mixedwood")
)

# Function to retrieve land cover code and category name based on the year
get_landcover_info <- function(row, lookup_table) {
  if (row$year >= 2020 && row$year <= 2022) {
    return(data.frame(code = NA, category = NA))
  } else {
    year_col <- paste0('landCover_', row$year)
    if (!year_col %in% names(row)) {
      return(data.frame(code = NA, category = NA))
    }
    code <- row[[year_col]]
    category <- lookup_table$category[lookup_table$code == code]
    
    if (length(category) == 0) {
      category <- NA
    }
    
    return(data.frame(code = code, category = category))
  }
}

# Apply the function to each row and bind the results as new columns
data_with_landcover <- data %>%
  rowwise() %>%
  mutate(landcover_info = list(get_landcover_info(cur_data(), landcover_lookup))) %>%
  unnest_wider(landcover_info)

# Check the structure of the new data
str(data_with_landcover)

# Write the new dataset to a CSV file
write_csv(data_with_landcover, 'Input/04_Cleaned_Waterton.csv')