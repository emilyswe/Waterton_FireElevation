#setwd
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# 1. Load Packages ----
# This section loads all necessary packages for the analysis.
packs <- c("haven", "ggplot2", "dplyr", "broom", "lubridate", "vegan", "lme4", 
           "tidyr", "mgcv", "GGally", "tidyverse", "MuMIn", "data.table", "ggeffects")

for (q in 1:length(packs)) {
  if (!require(packs[q], character.only = TRUE)) {
    install.packages(packs[q])  # Install the package if not already installed
    require(packs[q])           # Load the package
  }
}

# 2. Read in Erin's Data ----
# Erin's data is used for column consistency checks later on.
final_data <- read.csv("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Input/Raw/waterton_alpha_richness_perlocation_peryear.csv")

# 3. Read in Emily's Data ----
# Multi-visit dataset: Read and filter the dataset containing multiple visits per location.
emda1.raw <- read.csv("Input/01_Wrangled.csv") |> 
  mutate(location = toupper(location),  # Standardize location names to uppercase
         year = year(date_time),        # Extract the year from the date
         doy = yday(date_time),         # Extract the day of the year (doy)
         hod = hour(date_time),         # Extract the hour of the day (hod)
         task_duration = as.numeric(str_sub(task_duration, -100, -2))) |>  # Extract numeric task duration
  dplyr::filter(hod < 10,               # Only include records between 4 am and 10 am
                hod > 4,
                doy >= 152,             # Only include records between June 1st (doy 152) and July 1st (doy 182)
                doy <= 182,
                longitude > -114.2)     # Filter out locations based on longitude

# Single-visit dataset: Read and filter the dataset containing one visit per location (median count across visits).
emda4 <- read.csv("Input/04_Cleaned_Waterton.csv") |> 
  dplyr::filter(type != "PC_historic")  # Remove records labeled as 'PC_historic' which differ from other data types

# 4. Wrangle ----
# Remove rare species from the multi-visit data, make it wide, and add environmental variables like elevation.
# This block ensures species that are rare (less than 5 observations) are excluded.

spp.count <- emda1.raw |> 
  group_by(species) |> 
  summarize(counts = n()) |> 
  ungroup() |> 
  dplyr::filter(counts >= 5)  # Only keep species with 5 or more observations

# Filter the multi-visit dataset for common species, pivot it to wide format, and join with environmental data from the single-visit dataset
emda1 <- dplyr::filter(emda1.raw, species %in% spp.count$species) |> 
  pivot_wider(names_from = "species", values_from = "count", values_fill = 0) |> 
  left_join(emda4 |> 
              dplyr::select(location, latitude, longitude, elevation, slope, aspect, TPI, northness, grid.code, type) |> 
              unique(),  # Ensure the environmental data is unique per location
            by = "location") |> 
  dplyr::filter(type != "PC_historic")

# Remove duplicate columns for latitude and longitude ----
# Keep one set of latitude/longitude (we'll keep the .y columns for now)
emda1 <- emda1 %>% 
  dplyr::select(-latitude.x, -longitude.x)

# Rename latitude.y and longitude.y to latitude and longitude ----
emda1 <- emda1 %>% 
  rename(latitude = latitude.y, 
         longitude = longitude.y)

# 5. Check for Duplicates ----
# Ensure no duplicate rows for the same location and timestamp.
emda1 <- emda1 %>%
  distinct(location, date_time, .keep_all = TRUE)  # Keep only unique rows based on location and timestamp

# 6. Add Burn Variables ----
# Add a binary burn variable and before/after burn variable for both datasets.
# The burnYN condition accounts for spatial (grid.code) and temporal (year) rules.
emda1$burnYN <- as.factor(ifelse(((emda1$grid.code == 0) | (emda1$year < 2018)), 0, 1))
emda4$burnYN <- as.factor(ifelse(((emda4$grid.code == 0) | (emda4$year < 2018)), 0, 1))

# Add before/after burn variable
emda1$burnBA <- as.factor(ifelse(emda1$year < 2018, 0, 1))
emda4$burnBA <- as.factor(ifelse(emda4$year < 2018, 0, 1))

# 7. Explore ----
# Exploratory Plot 1: Histogram of Elevation by burnBA
HistogramElev <- ggplot(emda1) +
  geom_histogram(aes(x = elevation)) +
  facet_grid(. ~ burnBA)

# Save Plot 1
ggsave(filename = "Output/07_HistogramElev.png", plot = HistogramElev, width = 8, height = 6)

# Exploratory Plot 2: Spatial Distribution of Elevation by burnBA
EleSpatialDistribution_byburnBA <- ggplot(emda1) + 
  geom_point(aes(x = longitude, y = latitude, colour = elevation)) +
  facet_grid(. ~ burnBA) +
  scale_colour_viridis_c()

# Save Plot 2
ggsave(filename = "Output/07_EleSpatialDistribution_byburnBA.png", plot = EleSpatialDistribution_byburnBA, width = 8, height = 6)

# 8. Output ----
# Save the cleaned dataframes for use in the next scripts
write.csv(emda1, "Output/07_cleaned_multivisit_data.csv", row.names = FALSE)  # Multi-visit dataset
write.csv(emda4, "Output/07_cleaned_single_visit_data.csv", row.names = FALSE)  # Single-visit dataset
