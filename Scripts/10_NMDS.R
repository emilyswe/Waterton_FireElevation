#STEP 10: NMDS

#setwd
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# 1. Load required packages ----
 packs <- c("vegan", "dplyr", "ggplot2", "tidyverse")
 for (q in 1:length(packs)) {
    if (!require(packs[q], character.only = TRUE)) {
       install.packages(packs[q])
        require(packs[q])
      }
   }

# 2. Read in Cleaned Multi-visit and Single-visit Datasets ----
emda1 <- read.csv("Output/07_cleaned_multivisit_data.csv")  # Multi-visit dataset
emda4 <- read.csv("Output/07_cleaned_single_visit_data.csv")  # Single-visit dataset

# 3. Merge Necessary Columns ----
# Merge latitude, longitude, elevation, etc. from emda4 into emda1 (multi-visit)
superdata <- left_join(emda1, 
                       emda4 %>% dplyr::select(location, latitude, longitude, elevation, slope, aspect, TPI, northness, grid.code, burnYN, burnBA) %>%
                         unique(), 
                       by = "location")

# 4. Fix Column Names After Merge ----
# Remove .x and .y suffixes
colnames(superdata) <- gsub("\\.x$", "", colnames(superdata))  # Remove .x suffixes
colnames(superdata) <- gsub("\\.y$", "", colnames(superdata))  # Remove .y suffixes

# 5. Check if Elevation Exists and Is Numeric ----
# Validate that 'elevation' column was successfully merged and is numeric
if (!"elevation" %in% colnames(superdata)) {
  stop("Error: Elevation column is missing from superdata.")
}

# Check for missing elevation values
if (any(is.na(superdata$elevation))) {
  stop("Error: Missing values detected in the elevation column after merging.")
}

# Ensure 'elevation' is numeric
if (!is.numeric(superdata$elevation)) {
  superdata$elevation <- as.numeric(superdata$elevation)
}

# 5. Create Strata Based on Altitude and Burn Severity ----
# This replicates Elly's logic of using elevation and grid code to define strata
superdata$altbin <- cut(superdata$elevation, 
                        breaks = quantile(superdata$elevation, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                        labels = FALSE, include.lowest = TRUE)

# Adjust strata based on burn severity
superdata$strata <- ifelse(superdata$grid.code > 0, superdata$altbin + 3, superdata$altbin)

# 6. Identify Species Columns ----
# Use a regular expression to identify the 4-letter uppercase bird species codes
species_columns <- colnames(superdata)[grep("^[A-Z]{4}$", colnames(superdata))]

# 7. Prepare Data for NMDS ----
# Remove non-species columns, ensuring species columns are sorted alphabetically
non_species_columns <- setdiff(colnames(superdata), species_columns)
species_columns_sorted <- sort(species_columns)

# Create wide data format needed for NMDS
wide_data <- superdata[c(non_species_columns, species_columns_sorted)]

# Filter to ensure no missing strata and sample one observation per location
filtered_data <- wide_data %>%
  dplyr::filter(!is.na(strata)) %>%
  group_by(location) %>%
  slice_sample(n = 1)

# 8. Perform NMDS ----
# Run NMDS using species data, distance = jaccard, and k = 2
set.seed(999)
nmds_result <- metaMDS(filtered_data[, species_columns_sorted], distance = "jaccard", k = 2, trymax = 100)

# 9. Save NMDS Results ----
# Extract NMDS coordinates for each site
coordinates <- scores(nmds_result, display = "sites")
filtered_data$mdsA <- coordinates[, 1]
filtered_data$mdsB <- coordinates[, 2]

# Save the results as a CSV
write.csv(filtered_data[, c("location", "strata", "doy", "year", "mdsA", "mdsB")],  "Output/10_nmds_2axis.csv", row.names = FALSE)

# 10. NMDS Plot ----
# Calculate mean NMDS scores for each strata to add group labels
means <- aggregate(cbind(mdsA, mdsB) ~ strata, data = filtered_data, mean)

# Create the NMDS plot
nmds_plot <- ggplot(filtered_data, aes(x = mdsA, y = mdsB, color = factor(strata))) +
  geom_point(aes(shape = factor(strata))) +
  scale_color_manual(values = c("green", "green", "green", "orange", "orange", "orange")) +
  geom_text(data = means, aes(label = factor(strata)), vjust = 2, color = "black") +
  stat_ellipse(aes(fill = factor(strata)), geom = "polygon", level = 0.95, alpha = 0.2) +
  labs(color = "Strata", shape = "Strata", fill = "Confidence") +
  theme_bw() +
  theme(legend.position = "right")

# Save the plot
ggsave("Output/10_nmds_plot.png", plot = nmds_plot, width = 8, height = 6)

#cleaned plot ----

dat3 <- read.csv("Output/10_waterton_mds_2axis.csv") |> 
  mutate(treatment = ifelse(strata %in% c(1:3), "Unburned", "Burned"),
         elevation = case_when(strata %in% c(1, 4) ~ "Low elevation", 
                               strata %in% c(2, 5) ~ "Mid elevation",
                               strata %in% c(3, 6) ~ "High elevation"),
         elevation = factor(elevation, levels=c("Low elevation", "Mid elevation", "High elevation")))

linelegend <- data.frame(expand.grid(x=1, y=1, elevation=unique(dat3$elevation)))

plot3 <- ggplot(dat3, aes(x = mdsA, y = mdsB)) +
  stat_ellipse(aes(colour=treatment, linetype=elevation), geom = "polygon", level = 0.5, fill="white", alpha = 0, linewidth=1)
  #  geom_point(aes(shape = elevation, colour=treatment)) + 
  geom_line(aes(x=x, y=y, linetype=elevation), data=linelegend) +
  scale_color_manual(values=c("steelblue3", "tomato3"),name="") +
  scale_linetype_manual(values=c("solid", "dotted", "dashed"), name="") +
  theme_minimal() +
  theme(legend.position = "right") +
  xlim(c(-0.003, 0.006)) +
  ylim(c(-0.006, 0.005)) +
  labs(x = "NMDS axis 1",
       y = "NMDS axis 2")
plot3

ggsave(plot3, filename="Output/10_NMDS_Clean.jpeg", width=8, height=6)



#gpt check
sum(is.na(dat3$mdsA) | is.infinite(dat3$mdsA))
sum(is.na(dat3$mdsB) | is.infinite(dat3$mdsB))

unique(dat3$strata)

table(dat3$treatment, dat3$elevation)

dat3 %>%
  group_by(treatment, elevation) %>%
  summarize(count = n())


linelegend

print(unique(dat3$elevation))
print(linelegend)


dat3 %>%
  group_by(treatment, elevation) %>%
  summarize(count = n()) %>%
  print()

# Check for rows where mdsA or mdsB is non-finite (NA, NaN, Inf)
dat3 %>%
  filter(!is.finite(mdsA) | !is.finite(mdsB)) %>%
  print()

dat3 %>%
  group_by(strata) %>%
  summarize(count = n()) %>%
  print()

plot3 <- ggplot(dat3, aes(x = mdsA, y = mdsB, color = treatment)) +
  geom_point(aes(shape = elevation)) + 
  scale_color_manual(values = c("steelblue3", "tomato3"), name = "") +
  theme_minimal() +
  theme(legend.position = "right") +
  xlim(c(-0.003, 0.006)) +
  ylim(c(-0.006, 0.005)) +
  labs(x = "NMDS axis 1", y = "NMDS axis 2")

# Save the plot
ggsave(plot3, filename = "Output/10_NMDS_PointsOnly.jpeg", width = 8, height = 6)
print(plot3)

# Check for missing or non-finite values in other relevant columns
sum(is.na(dat3$strata) | is.infinite(dat3$strata))
sum(is.na(dat3$treatment) | is.infinite(dat3$treatment))
sum(is.na(dat3$elevation) | is.infinite(dat3$elevation))










