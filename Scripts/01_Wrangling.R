#NOTES

#PC_WLNP projects have an abundance cap of 3



#PREAMBLE####

#1. Load library----
library(tidyverse)
library(wildRtrax)
library(lubridate)

#2. Log in----
source("login.R")
wt_auth()

#3. Read in functional wildRtrax code---
source("Scripts/00_wildRtrax_functions.R")

#READ IN WILDTRAX DATA#####

#1. Get list of files----
files.wt <- list.files("Input", pattern="*main_report.csv", full.names = TRUE)

#2. Read in files----
dat.wt <- purrr::map(.x=files.wt, .f = ~ read.csv(.x)) %>% 
  do.call(rbind, .)

#3. Clean up the species list----
dat.tidy <- wt_tidy_species(dat.wt, remove=c("mammal", "amphibian", "abiotic", "insect", "unknown"))

#4. Replace the tmtts----
table(dat.tidy$individual_count)
dat.tmtt <- wt_replace_tmtt(dat.tidy) %>% 
  mutate(individual_count = as.numeric(ifelse(individual_count=="CI 1", 1, individual_count)))
table(dat.tmtt$individual_count)

#5. Sum to abundance per species per survey
dat.sum <- wt_make_wide(dat.tmtt) %>% 
  pivot_longer(ALFL:YRWA, names_to="species", values_to="count") %>% 
  dplyr::filter(count!=0)

#READ IN HISTORICAL DATA####

#1. Get the data----
dat.hist <- read.csv("Input/PC_WLNP_EI Monitoring_2007to2019.csv")

#2. Wrangle----
dat.pc <- dat.hist %>% 
  mutate(location = paste0(Transect, "_", PointID),
         date_time = ymd_hms(paste0(Year, "-", Month, "-", Day, " 06:00:00")),
         organization = "PC",
         project_id = 999,
         task_duration = "180s",
         latitude = NA,
         longitude = NA) %>% 
  rename(count = X.detected_0.3m20s,
         species = AOU_Code) %>% 
  dplyr::filter(Time.1st.Detected <= 3,
                species %in% dat.sum$species)

#PUT TOGETHER####

#1. List the columns we want----
cols <- c("organization", "project_id", "location", "latitude", "longitude", "date_time", "task_duration", "species", "count")

#2. Put it together----
dat.out <- dat.sum %>% 
  rename(date_time = recording_date_time) %>% 
  dplyr::select(all_of(cols)) %>% 
  rbind(dat.pc %>% 
          dplyr::select(all_of(cols)))

#3. Write it out----
write.csv(dat.out, "Input/01_Wrangled.csv", row.names = FALSE)
