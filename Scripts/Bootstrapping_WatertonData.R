library(dplyr)
library(lubridate)

# Read the data from the merged CSV file
merged_df <- read.csv("~/Documents/local-git/Waterton/Output/MasterSpeciesMatrix_Waterton.csv")

# Convert recording_date_time to Date-Time object
# It's important to specify the format to match your data
merged_df$recording_date_time <- dmy(merged_df$recording_date_time, tz = "UTC")

# Define good detectability criteria (e.g., within 2 hours of sunrise, in the month of June)
# Since sunrise is approximately at 5:30 AM, we define a time window from 3:30 AM to 7:30 AM
good_detectability <- merged_df %>%
  filter(hour(recording_date_time) >= 3 & hour(recording_date_time) < 8, # 2 hours of sunrise
         month(recording_date_time) == 6) # Month of June

# Randomly select one visit for each unique ARU location for each year
set.seed(123) # for reproducibility
unique_visits <- good_detectability %>%
  group_by(location, year = year(recording_date_time)) %>%
  do(sample_n(., size = 1))

# If you want to include locations that are new in the most recent year and not present in the previous year,
# you would first identify the unique locations in the most recent year and then perform the same
# random selection for those locations.
locations_recent_year <- unique(good_detectability %>% filter(year(recording_date_time) == max(year(recording_date_time)))$location)
locations_previous_years <- unique(good_detectability %>% filter(year(recording_date_time) < max(year(recording_date_time)))$location)
new_locations <- setdiff(locations_recent_year, locations_previous_years)

new_visits_recent_year <- good_detectability %>%
  filter(location %in% new_locations) %>%
  group_by(location) %>%
  do(sample_n(., size = 1))

# Combine the new visits of the most recent year with the previous unique visits
final_dataset <- bind_rows(unique_visits, new_visits_recent_year)

# Write the final dataset to a new CSV file
write.csv(final_dataset, "~/Documents/local-git/Waterton/Output/BootstrappedDataset_Waterton.csv", row.names = FALSE)

