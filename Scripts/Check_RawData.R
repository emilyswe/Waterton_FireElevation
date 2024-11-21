# Load necessary packages
library(dplyr)
library(tidyr)
library(lubridate)

#setwd
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# Load glmer_data2
glmer_data2 <- readRDS("Output/glmer_data2.RDS")

# Load predictions2
predictions2 <- readRDS("Output/predictions2.RDS")

# Step 1: Filter `glmer_data2` directly for high elevation data for the problematic species
problematic_species <- c("AMCR", "CANG", "GRCA", "SACR", "TRES", "VGSW", "WEWP")

# Extract entries from glmer_data2 for problematic species and elevations > 1800m
high_elevation_data_points <- glmer_data2 |> 
  dplyr::filter(altitude > 1800) |> 
  dplyr::select(location, hod, doy, year, altitude, all_of(problematic_species))

# Step 2: Identify which species were detected (presence = 1)
high_elevation_species_detections <- high_elevation_data_points |> 
  tidyr::pivot_longer(cols = all_of(problematic_species), 
                      names_to = "species", 
                      values_to = "presence") |> 
  dplyr::filter(presence == 1) # Only keep rows where species were detected

# Step 3: Convert `doy` to a standard date format using `lubridate`
high_elevation_species_detections <- high_elevation_species_detections |> 
  mutate(date = ymd(paste(year, doy, sep = "-")))

# Step 4: Save the filtered data for reference in WildTrax
write.csv(high_elevation_species_detections, "Output/High_Elevation_Species_Detections.csv", row.names = FALSE)

###################### keep investigating species not in high_elevation_species_detections##################################

# List of species to double-check
species_to_check <- c("AMCR", "CANG", "GRCA", "SACR","WEWP")

# Step 1: Directly inspect `glmer_data2` for these species at high elevations
missing_species_data_points <- glmer_data2 |> 
  dplyr::filter(altitude > 1800) |> 
  dplyr::select(location, hod, doy, altitude, all_of(species_to_check))

# Step 2: Pivot longer to see detections
missing_species_detections <- missing_species_data_points |> 
  tidyr::pivot_longer(cols = all_of(species_to_check), 
                      names_to = "species", 
                      values_to = "presence") |> 
  dplyr::filter(presence == 1) # Keep where species were detected

# Output this to see what we find
write.csv(missing_species_detections, "Output/Missing_Species_High_Elevation.csv", row.names = FALSE)

##################
# Directly check predictions for these species at high elevations, regardless of confirmed presence
u_shaped_high_elevations <- predictions2 |> 
  dplyr::filter(species %in% species_to_check, x > 1800)

# Output to examine what's there
View(u_shaped_high_elevations)
write.csv(u_shaped_high_elevations, "Output/U_Shaped_High_Elevation_Examination.csv", row.names = FALSE)

# -------------------------------------
# Check Raw Data for Missing Species
# -------------------------------------

# List of species to double-check
species_to_check <- c("AMCR", "CANG", "GRCA", "SACR", "WEWP")

# Step 1: Directly inspect `glmer_data2` for these species at high elevations
missing_species_data_points <- glmer_data2 |> 
  dplyr::filter(altitude > 1800) |> 
  dplyr::select(location, hod, doy, year, altitude, all_of(species_to_check))

# Step 2: Pivot longer to see detections
missing_species_detections <- missing_species_data_points |> 
  tidyr::pivot_longer(cols = all_of(species_to_check), 
                      names_to = "species", 
                      values_to = "presence") |> 
  dplyr::filter(presence == 1) # Keep where species were detected

# Step 3: Save this to recheck WildTrax
write.csv(missing_species_detections, "Output/Missing_Species_High_Elevation.csv", row.names = FALSE)


#########save objects

# Save glmer_data2
saveRDS(glmer_data2, "Output/glmer_data2.RDS")

# Save predictions2
saveRDS(predictions2, "Output/predictions2.RDS")


###########################################################
# Specifically check predictions for AMCR at any elevation
amcr_predictions <- predictions2 |> 
  dplyr::filter(species == "AMCR")

# Save output to review
write.csv(amcr_predictions, "Output/AMCR_Predictions_Check.csv", row.names = FALSE)


####################### glmer2 detection checks ###########################
# Load glmer_data2
glmer_data2 <- readRDS("Output/glmer_data2.RDS")

# Step 1: Check presence records in `glmer2` for the problematic species
problematic_species <- c("AMCR", "CANG", "GRCA", "SACR", "WEWP")

# Extract all entries for the problematic species across all elevations
all_elevations_glmer2 <- glmer_data2 |> 
  dplyr::select(location, hod, doy, year, altitude, all_of(problematic_species))

# Pivot to see detections
glmer2_detection_check <- all_elevations_glmer2 |> 
  tidyr::pivot_longer(cols = all_of(problematic_species), 
                      names_to = "species", 
                      values_to = "presence") |> 
  dplyr::filter(presence == 1) # Keep where species were detected

# Check detections above 1800m 
glmer2_above1800m_detections <- glmer2_detection_check |> 
  dplyr::filter(altitude > 1800)

# Output to view or save
View(glmer2_above1800m_detections)
write.csv(glmer2_above1700m_detections, "Output/glmer2_above1700m_detections.csv", row.names = FALSE)


##################### emda4 detection checks #############################
# Step 2: Check `emda4` for abundance of these species to validate raw data
emda4 <- read.csv("Output/07_cleaned_single_visit_data.csv")
all_elevations_emda4 <- emda4 |> 
  dplyr::select(location, year, elevation, all_of(problematic_species))

# Pivot to see abundances
emda4_detection_check <- all_elevations_emda4 |> 
  tidyr::pivot_longer(cols = all_of(problematic_species), 
                      names_to = "species", 
                      values_to = "abundance") |> 
  dplyr::filter(abundance > 0) # Keep where species were detected

#check detections above 1800m
emda4_above1800m_detections <- emda4_detection_check |> 
  dplyr::filter(elevation > 1800)

# Output to view or save
View(emda4_above1800m_detections)
write.csv(emda4_above1800m_detections, "emda4_above1700m_detections", row.names = FALSE)


# Save the outputs if needed
write.csv(all_elevation_species_detections, "Output/All_Elevation_Species_Detections.csv", row.names = FALSE)
write.csv(all_elevation_species_abundance, "Output/All_Elevation_Species_Abundance.csv", row.names = FALSE)

# Validation Checks: Findings were consistent across these scripts:
# 1. `glmer_data2` indicated no high-elevation detections for problematic species.
# 2. Raw data `emda4` still had 14 rows — Suggesting possible data loss during conversion.
# 3. Model trends still visible, suggesting overfit or irregular predictions.



