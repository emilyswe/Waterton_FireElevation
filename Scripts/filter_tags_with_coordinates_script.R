
# Load required library
library(dplyr)

# Function to filter tag data based on site locations and add location coordinates
filter_tags_with_coordinates <- function(site_file, tag_file, output_file) {
  # Read the CSV files
  sites_df <- read.csv(site_file, stringsAsFactors = FALSE)
  tags_df <- read.csv(tag_file, stringsAsFactors = FALSE)

  # Extract unique location IDs from the sites dataframe
  unique_location_ids <- unique(sites_df$location)

  # Filter the tags dataframe
  filtered_tags_df <- tags_df %>% filter(location %in% unique_location_ids)

  # Add location coordinates
  final_df <- merge(filtered_tags_df, sites_df[, c('location', 'latitude', 'longitude')], by = 'location', all.x = TRUE)

  # Write the final dataframe to a new CSV file
  write.csv(final_df, output_file, row.names = FALSE)
}

# Example usage of the function
# filter_tags_with_coordinates("path_to_site_file.csv", "path_to_tag_file.csv", "path_to_output_file.csv")
