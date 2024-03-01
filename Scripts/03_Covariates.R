#PREAMBLE####

#1. Load library----
library(tidyverse)
library(sf)
library(terra)

#2. Load the data & project----
raw <- read.csv("Input/02_Standardized.csv") %>% 
  st_as_sf(coords=c("longitude", "latitude"), crs=4326)

#FIRE SEVERITY######

#1. Read in the layer----
fire <- read_sf("Input/Kenow 2017 Burn severity/Kenow_severity_classes.shp") %>% 
  st_transform(crs=crs(raw)) %>% 
  st_make_valid() %>% 
  dplyr::select(gridcode) %>% 
  rename(severity = gridcode)

#2. Extract----
dat.fire <- raw %>% 
  st_intersection(fire) %>% 
  full_join(raw) %>% 
  mutate(severity = ifelse(is.na(severity), 0, severity))

#TOPOGRAPHY#####

#1. Read in the layer----
dem <- rast("File.tiff")

#2. Extract----
dat.dem <- dat.fire %>% 
  terra::extract(dem)

#LANDCOVER#######


#FINALIZE######

#1. Sanity checks----
ggplot(dat.out) + 
  geom_point(aes(x=longitude, y=latitude, colour=severity))



#X. Save----
write.csv(dat.out, "Input/03_Covariates.csv", row.names=FALSE)