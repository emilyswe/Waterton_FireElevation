#NOTES

#PC_WLNP projects have an abundance cap of 3



#PREAMBLE####

#1. Load library----
library(tidyverse)
library(wildRtrax)
library(lubridate)
library(sf)

#2. Log in----
source("login.R")
wt_auth()

#3. Read in functional wildRtrax code---
source("Scripts/00_wildRtrax_functions.R")

#READ IN WILDTRAX DATA#####

#1. Get list of files----
files.wt <- list.files("Input/Raw", pattern="*main_report.csv", full.names = TRUE)

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
dat.hist <- read.csv("Input/Raw/PC_WLNP_EI Monitoring_2007to2019.csv")

#2. Get the locations----
loc.raw <- read.csv("Input/Raw/WLNP_tbl_point_locations.csv") %>% 
  dplyr::filter(!is.na(UTM.Northing),
                UTM.Northing > 5000000)
  
#3. Reproject to lat long----
loc.11 <- loc.raw %>% 
  dplyr::filter(UTM.Zone==11) %>% 
  st_as_sf(coords=c("UTM.Easting", "UTM.Northing"), crs=26911) %>% 
  st_transform(crs=4326)

loc.12 <- loc.raw %>% 
  dplyr::filter(UTM.Zone==12) %>% 
  st_as_sf(coords=c("UTM.Easting", "UTM.Northing"), crs=26912) %>% 
  st_transform(crs=4326)

loc.use <- rbind(loc.11, loc.12) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(longitude = X, latitude = Y) %>% 
  cbind(loc.raw) %>% 
  dplyr::select(-Comments)

#4. Get the time data----
time.raw <- read.csv("Input/Raw/WLNP_tbl_Field_information.csv") %>% 
  mutate(date_time = ymd_hm(paste0(Year, "-", Month, "-", Day, " ", Time))) %>% 
  dplyr::filter(!is.na(date_time)) %>% 
  dplyr::select(-TransectID)

#3. Wrangle----
dat.pc <- dat.hist %>% 
  inner_join(loc.use) %>% 
  inner_join(time.raw) %>% 
  mutate(location = paste0(Transect, "_", PointID),
         organization = "PC",
         project_id = 999,
         task_duration = "180s") %>% 
  rename(count = X.detected_0.3m20s,
         species = AOU_Code) %>% 
  dplyr::filter(Time.1st.Detected <= 3,
                species %in% dat.sum$species)

#PUT TOGETHER####

#1. List the columns we want----
cols <- c("organization", "project_id", "location", "latitude", "longitude", "date_time", "task_duration", "species", "count")

#2. Put it together----
dat.out <- dat.sum %>% 
  mutate(date_time = ymd_hms(recording_date_time)) %>% 
  dplyr::select(all_of(cols)) %>% 
  rbind(dat.pc %>% 
          dplyr::select(all_of(cols))) %>% 
  mutate(hour = hour(date_time))

#3. Sanity checks----

#year
ggplot(dat.out) +
  geom_histogram(aes(x=date_time, fill=factor(project_id)))

#month
ggplot(dat.out) +
  geom_histogram(aes(x=yday(date_time), fill=factor(project_id)))

#hour
ggplot(dat.out) +
  geom_histogram(aes(x=hour(date_time), fill=factor(project_id)))

#location
ggplot(dat.out) +
  geom_point(aes(x=longitude, y=latitude, colour=factor(project_id)))

#3. Write it out----
write.csv(dat.out, "Input/01_Wrangled.csv", row.names = FALSE)
