#PREAMBLE####

#1. Load library----
library(tidyverse)
library(lubridate)

#2. Read in data----
raw <- read.csv("Input/01_Wrangled.csv") %>% 
  mutate(type = case_when(organization=="BU" ~ "BU",
                          organization=="PC" & project_id %in% c(598, 353, 140, 10) ~ "BU",
                          organization=="PC" & project_id %in% c(925, 1113) ~ "PC_truncated",
                          project_id==999 & year(date_time) <= 2017 ~ "PC_historic",
                          project_id==999 & year(date_time) > 2017 ~ "PC_postfire"),
         doy = yday(date_time),
         hour = hour(date_time))

#FILTER####

#1. Remove things we don't want----
dat1 <- raw %>% 
  dplyr::filter(hour < 10,
                doy >= 152,
                doy <= 182,
                longitude > -114.2)

#2. Remove rare species----
spp.count <- dat1 %>% 
  group_by(species) %>% 
  summarize(counts = n()) %>% 
  ungroup() %>% 
  dplyr::filter(counts >= 5)

dat2 <- dplyr::filter(dat1, species %in% spp.count$species)

#ACCOUNT FOR SURVEY EFFORT####

#1. Write a function to find the mode----
mode = function(){
  return(sort(-table(dat))[1])
}

#1. Take median count per species per location per year----
med <- dat2 %>% 
  mutate(year = year(date_time)) %>% 
  group_by(organization, type, project_id, location, latitude, longitude, year, species) %>% 
  summarize(count = ceiling(mean(count))) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(count))

#POKE AROUND#####

#1. Visualize----
ggplot(med) +
  geom_histogram(aes(x=count)) +
  facet_wrap(~type, scales="free_y")

table(med$type, med$count)

#2. Test difference in counts----
m1 <- glm(count ~ type + species, family="poisson", data=med)
summary(m1)

#FINAL FILTERING & MAKE WIDE#####

#1. Do the thing----
use <- med %>% 
  dplyr::filter(!type %in% c("postfire", "historic")) %>% 
  pivot_wider(names_from="species", values_from="count", values_fill = 0)

#2. Save the thing----
write.csv(use, "Input/02_Standardized.csv", row.names = FALSE)
