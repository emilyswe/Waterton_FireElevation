#setwd
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# 1. Load Packages ----
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)  
library(lubridate)


# 2. Read in Single-Visit Dataset ----
emda4 <- read.csv("Input/04_Cleaned_Waterton.csv") |> 
  dplyr::filter(type != "PC_historic")  # Remove 'PC_historic' records

# 3. Add Burn Variables ----
emda4$burnYN <- as.factor(ifelse(((emda4$grid.code == 0) | (emda4$year < 2018)), 0, 1))
emda4$burnBA <- as.factor(ifelse(emda4$year < 2018, 0, 1))

# 4. Filter Out Non-Territorial and Unreliable Species ----
species_to_remove <- c(
  "RWBL", "BHCO", "YHBL", "BRBL",  # Icterids
  "AMCR", "GRAJ", "BLJA", "STJA", "CORA", "CAJA",  # Corvids
  "BANS", "CLSW", "TRES", "VGSW", "NRWS", # Hirundinids
  "SACR", "CANG", "PBGR", "MALL", "AMWI", "AMBI",  # Waterfowl/Non-Territorial
  "KILL", "SPSA", "COLO", "SORA",  # Shorebirds, Marsh Birds, Loons
  "PUFI", "MODO", "WISN", "AMCO",  # PUFI false ID, Mourning Dove, American Coot
  "GHOW", "RTHA", "BBMA", "DUGR", "RUGR", "CLNU", "STGR"  # Owls, Raptors, Grouse, etc.
)

emda4_filtered <- emda4 %>%
  select(-all_of(species_to_remove))

# 5. Save the Final Dataset ----
write.csv(emda4_filtered, "Output/07_singlevisit_typeA.csv", row.names = FALSE)

# 6. Verify Final Species List ----
final_species_columns <- grep("^[A-Z]{4}$", colnames(emda4_filtered), value = TRUE)
final_species_list <- data.frame(species = final_species_columns)
print(final_species_list)
