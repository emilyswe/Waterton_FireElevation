# 1. Load libraries ----
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(stringr)

# 2. Load and filter raw detections ----
detections <- read_csv("Input/01_Wrangled.csv") %>%
  mutate(
    date_time = ymd_hms(date_time),
    year = year(date_time),
    doy = yday(date_time),
    hour = hour(date_time)
  ) %>%
  filter(
    hour >= 4 & hour < 10,
    doy >= 152 & doy <= 182,
    longitude > -114.2
  )

# 3. Remove rare species (fewer than 5 detections) ----
spp_counts <- detections %>%
  count(species, name = "n") %>%
  filter(n >= 5)

filtered <- detections %>%
  filter(species %in% spp_counts$species)

# 4. Aggregate to mean count per species × location × year ----
agg_counts <- filtered %>%
  group_by(location, year, species) %>%
  summarise(count = ceiling(mean(count, na.rm = TRUE)), .groups = "drop")

# 5. Pivot to wide format ----
wide_counts <- agg_counts %>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0)

# 6. Load covariates, remove old species columns ----
covariates <- read_csv("Input/04_Cleaned_Waterton.csv") %>%
  filter(type != "PC_historic") %>%
  distinct(location, year, .keep_all = TRUE)

# Remove any AOU species code columns (4-letter uppercase codes)
spp_cols_to_remove <- grep("^[A-Z]{4}$", colnames(covariates), value = TRUE)
covariates_clean <- covariates %>%
  select(-all_of(spp_cols_to_remove))

# 7. Join with wide-format counts ----
joined <- left_join(covariates_clean, wide_counts, by = c("location", "year"))

# 8. Create burn variables ----
joined <- joined %>%
  mutate(
    `grid code` = ifelse(is.na(`grid code`), 0, `grid code`),
    burnYN = as.factor(ifelse(`grid code` == 0 | year < 2018, 0, 1)),
    burnBA = as.factor(ifelse(year < 2018, 0, 1))
  )

# 9. Filter out non-Territory A species ----
species_to_remove <- c(
  "RWBL", "BHCO", "YHBL", "BRBL",      # Icterids
  "AMCR", "GRAJ", "BLJA", "STJA", "CORA", "CAJA",  # Corvids
  "BANS", "CLSW", "TRES", "VGSW", "NRWS",          # Hirundinids
  "SACR", "CANG", "PBGR", "MALL", "AMWI", "AMBI",  # Waterfowl
  "KILL", "SPSA", "COLO", "SORA",                  # Shorebirds, Loons
  "PUFI", "MODO", "WISN", "WIWR","AMCO",                  # False ID, Doves
  "GHOW", "RTHA", "BBMA", "DUGR", "RUGR", "CLNU", "STGR"  # Raptors, Owls, Grouse
)

species_columns <- grep("^[A-Z]{4}$", colnames(joined), value = TRUE)
species_keep <- setdiff(species_columns, species_to_remove)

# 9.5 Filter out species with <5 detections after all filtering
species_check <- joined %>%
  select(all_of(species_keep)) %>%
  summarise_all(~sum(. > 0, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "species", values_to = "detections")

species_keep_final <- species_check %>%
  filter(detections >= 5) %>%
  pull(species)

# 10. Final dataset ----
final_data <- joined %>%
  select(location, year, elevation, slope, aspect, TPI, northness, type,
         `grid code`, burnYN, burnBA, all_of(species_keep_final)) %>%
  filter(!is.na(elevation)) %>%
  distinct()

# Save it
write_csv(final_data, "Output/07_singlevisit_typeA.csv")
cat("✅ Done! Final dataset has", nrow(final_data), "rows and", length(species_keep_final), "species.\n")
