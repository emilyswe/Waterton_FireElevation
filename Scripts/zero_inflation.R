# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

###########single visit 

# Load the single-visit dataset
single_visit_typeA <- read.csv("Output/07_cleaned_single_visit_filtered.csv")

# Identify species columns (assume these represent the median counts for each species at each site)
species_columns <- names(single_visit_typeA)[grep("^[A-Z]{4}$", names(single_visit_typeA))]

# Step 1: Count occurrences of zeros, ones, twos, etc., across all species columns
value_counts <- single_visit_typeA %>%
  select(all_of(species_columns)) %>%
  pivot_longer(cols = everything(), names_to = "species", values_to = "median_count") %>%
  group_by(median_count) %>%
  summarise(count = n(), .groups = "drop")

# Print value counts
print(value_counts)

# Step 2: Scatterplot of Median Count (y-axis) vs. Number of Visits (x-axis)
# For single-visit data, assume the number of visits is constant (1 for all rows)
single_visit_typeA <- single_visit_typeA %>%
  mutate(number_of_visits = 1)

# Flatten the dataset into species and median counts for scatterplot
species_median_counts <- single_visit_typeA %>%
  select(location, number_of_visits, all_of(species_columns)) %>%
  pivot_longer(cols = all_of(species_columns), names_to = "species", values_to = "median_count")

# Scatterplot
plot1 <- ggplot(species_median_counts, aes(x = number_of_visits, y = median_count)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(
    title = "Median Count vs. Number of Visits (Single Visit Type A)",
    x = "Number of Visits",
    y = "Median Count (per Species)"
  ) +
  theme_bw()

# Save Scatterplot
ggsave("Output/12_Effort/single_visit_median_vs_visits_fixed.png", plot = plot1, width = 8, height = 6)

# Step 3: Histogram of Species-Level Median Counts
# Use the long-format species_median_counts dataframe
plot2 <- ggplot(species_median_counts, aes(x = median_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Median Counts (Single Visit Type A)",
    x = "Median Count",
    y = "Frequency"
  ) +
  theme_bw()

# Save Histogram
ggsave("Output/12_Effort/single_visit_histogram_median_counts_fixed.png", plot = plot2, width = 8, height = 6)

# Display plots (optional)
print(plot1)
print(plot2)

###############multivist
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the dataset
multi_visit_typeA <- read.csv("Output/07_multivist_typeA.csv")

# Identify species columns (assume these represent the median counts for each species at each site)
species_columns <- names(multi_visit_typeA)[grep("^[A-Z]{4}$", names(multi_visit_typeA))]

# Step 1: Count occurrences of zeros, ones, twos, etc., across all species columns
value_counts <- multi_visit_typeA %>%
  select(all_of(species_columns)) %>%
  pivot_longer(cols = everything(), names_to = "species", values_to = "median_count") %>%
  group_by(median_count) %>%
  summarise(count = n(), .groups = "drop")

# Print value counts
print(value_counts)

# Step 2: Scatterplot of Median Count (y-axis) vs. Number of Visits (x-axis)
# For multi-visit data, calculate the number of visits per site
multi_visit_typeA <- multi_visit_typeA %>%
  group_by(location) %>%
  mutate(number_of_visits = n()) %>%  # Count the number of visits for each site
  ungroup()

# Flatten the dataset into species and median counts for scatterplot
species_median_counts <- multi_visit_typeA %>%
  select(location, number_of_visits, all_of(species_columns)) %>%
  pivot_longer(cols = all_of(species_columns), names_to = "species", values_to = "median_count")

# Scatterplot
plot1 <- ggplot(species_median_counts, aes(x = number_of_visits, y = median_count)) +
  geom_jitter(alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue") +
  labs(
    title = "Median Count vs. Number of Visits (Multi-Visit Type A)",
    x = "Number of Visits",
    y = "Median Count (per Species)"
  ) +
  theme_bw()

###uopdated code
# Scatterplot with LOESS smoothing and simplified computation
plot1 <- ggplot(species_median_counts, aes(x = number_of_visits, y = median_count)) +
  geom_jitter(alpha = 0.6, width = 0.2) +  # Add some jitter to spread out points
  geom_smooth(method = "loess", color = "blue", se = FALSE) +  # Disable confidence interval to save memory
  labs(
    title = "Median Count vs. Number of Visits (Multi-Visit Type A)",
    x = "Number of Visits",
    y = "Median Count (per Species)"
  ) +
  theme_bw()

# Save the corrected scatterplot
ggsave("Output/12_Effort/multi_visit_median_vs_visits_corrected.png", plot = plot1, width = 8, height = 6)


# Save Scatterplot
ggsave("Output/12_Effort/multi_visit_median_vs_visits_fixed.png", plot = plot1, width = 8, height = 6)

# Step 3: Histogram of Species-Level Median Counts
# Use the long-format species_median_counts dataframe
plot2 <- ggplot(species_median_counts, aes(x = median_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Median Counts (Multi-Visit Type A)",
    x = "Median Count",
    y = "Frequency"
  ) +
  theme_bw()

# Save Histogram
ggsave("Output/12_Effort/multi_visit_histogram_median_counts_fixed.png", plot = plot2, width = 8, height = 6)

# Display plots (optional)
print(plot1)
print(plot2)

##########raw, cleaned all visits
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the raw dataset (file path provided)
raw_dataset <- read.csv("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Input/04_Cleaned_Waterton.csv")

# Identify species columns (assume these represent the median counts for each species at each site)
species_columns <- names(raw_dataset)[grep("^[A-Z]{4}$", names(raw_dataset))]

# Step 1: Calculate the number of visits per site
raw_dataset <- raw_dataset %>%
  group_by(location) %>%
  mutate(number_of_visits = n()) %>%  # Count the number of visits for each site
  ungroup()

# Step 2: Flatten the dataset for plotting
species_median_counts <- raw_dataset %>%
  select(location, number_of_visits, all_of(species_columns)) %>%
  pivot_longer(cols = all_of(species_columns), names_to = "species", values_to = "median_count")

# Step 3: Scatterplot of Median Count (y-axis) vs. Number of Visits (x-axis)
plot1 <- ggplot(species_median_counts, aes(x = number_of_visits, y = median_count)) +
  geom_jitter(alpha = 0.6, width = 0.2) +  # Add jitter for better visibility
  geom_smooth(method = "loess", color = "blue", se = FALSE) +  # LOESS trend line without confidence interval
  labs(
    title = "Median Count vs. Number of Visits (Raw Dataset)",
    x = "Number of Visits",
    y = "Median Count (per Species)"
  ) +
  theme_bw()

# Save Scatterplot
ggsave("Output/12_Effort/07_allvisits_scatterplot.png", plot = plot1, width = 8, height = 6)

# Step 4: Histogram of Species-Level Median Counts
plot2 <- ggplot(species_median_counts, aes(x = median_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Median Counts (Raw Dataset)",
    x = "Median Count",
    y = "Frequency"
  ) +
  theme_bw()

# Save Histogram
ggsave("Output/12_Effort/07_allvisits_histogram.png", plot = plot2, width = 8, height = 6)

# Display plots (optional)
print(plot1)
print(plot2)


##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
#REMOVE ZEROS and check relationships again -

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# === SINGLE-VISIT DATASET ===
single_visit_typeA <- read.csv("Output/07_cleaned_single_visit_filtered.csv")
species_columns <- names(single_visit_typeA)[grep("^[A-Z]{4}$", names(single_visit_typeA))]

# Filter out zero counts
non_zero_single_visit <- single_visit_typeA %>%
  pivot_longer(cols = all_of(species_columns), names_to = "species", values_to = "median_count") %>%
  filter(median_count > 0)

# Scatterplot (non-zero data) - Add jitter to simulate variation
plot1_single <- ggplot(non_zero_single_visit, aes(x = factor(1), y = median_count)) +
  geom_jitter(width = 0.2, alpha = 0.6) +  # Add jitter for visibility
  labs(
    title = "Median Count for Non-Zero Single Visit Data",
    x = "Single Visit (Constant)",
    y = "Median Count (per Species)"
  ) +
  theme_bw()

# Save the scatterplot
ggsave("Output/12_Effort/single_visit_non_zero_scatterplot.png", plot = plot1_single, width = 8, height = 6)

# Histogram (non-zero data)
plot2_single <- ggplot(non_zero_single_visit, aes(x = median_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Median Counts (Non-Zero Single Visit)",
    x = "Median Count",
    y = "Frequency"
  ) +
  theme_bw()

# Save the histogram
ggsave("Output/12_Effort/single_visit_non_zero_histogram.png", plot = plot2_single, width = 8, height = 6)

# Print plots (optional)
print(plot1_single)
print(plot2_single)


# === MULTI-VISIT DATASET ===
multi_visit_typeA <- read.csv("Output/07_multivist_typeA.csv")
species_columns <- names(multi_visit_typeA)[grep("^[A-Z]{4}$", names(multi_visit_typeA))]

# Calculate number of visits and filter out zeros
multi_visit_typeA <- multi_visit_typeA %>%
  group_by(location) %>%
  mutate(number_of_visits = n()) %>%
  ungroup()

non_zero_multi_visit <- multi_visit_typeA %>%
  pivot_longer(cols = all_of(species_columns), names_to = "species", values_to = "median_count") %>%
  filter(median_count > 0)

# Scatterplot
plot1_multi <- ggplot(non_zero_multi_visit, aes(x = number_of_visits, y = median_count)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(
    title = "Median Count vs. Number of Visits (Non-Zero Multi-Visit)",
    x = "Number of Visits",
    y = "Median Count (per Species)"
  ) +
  theme_bw()

ggsave("Output/12_Effort/non_zero_multi_visit_scatterplot.png", plot = plot1_multi, width = 8, height = 6)

# Histogram
plot2_multi <- ggplot(non_zero_multi_visit, aes(x = median_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Median Counts (Non-Zero Multi-Visit)",
    x = "Median Count",
    y = "Frequency"
  ) +
  theme_bw()

ggsave("Output/12_Effort/non_zero_multi_visit_histogram.png", plot = plot2_multi, width = 8, height = 6)

# === RAW ALL-VISIT DATASET ===
raw_dataset <- read.csv("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Input/04_Cleaned_Waterton.csv")
species_columns <- names(raw_dataset)[grep("^[A-Z]{4}$", names(raw_dataset))]

# Calculate number of visits and filter out zeros
raw_dataset <- raw_dataset %>%
  group_by(location) %>%
  mutate(number_of_visits = n()) %>%
  ungroup()

non_zero_raw <- raw_dataset %>%
  pivot_longer(cols = all_of(species_columns), names_to = "species", values_to = "median_count") %>%
  filter(median_count > 0)

# Scatterplot
plot1_raw <- ggplot(non_zero_raw, aes(x = number_of_visits, y = median_count)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(
    title = "Median Count vs. Number of Visits (Non-Zero Raw Dataset)",
    x = "Number of Visits",
    y = "Median Count (per Species)"
  ) +
  theme_bw()

ggsave("Output/12_Effort/non_zero_raw_scatterplot.png", plot = plot1_raw, width = 8, height = 6)

# Histogram
plot2_raw <- ggplot(non_zero_raw, aes(x = median_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Median Counts (Non-Zero Raw Dataset)",
    x = "Median Count",
    y = "Frequency"
  ) +
  theme_bw()

ggsave("Output/12_Effort/non_zero_raw_histogram.png", plot = plot2_raw, width = 8, height = 6)

###############
###############
#rename files for clarity
# Set the directory containing the files
file_directory <- "/Users/Bronwyn/Documents/local-git/Waterton_FireElevation/Output/12_Effort"

# Define the current and new file names
file_renames <- list(
  "non_zero_raw_histogram.png" = "Raw_NonZero_Histogram.png",
  "non_zero_raw_scatterplot.png" = "Raw_NonZero_Scatterplot.png",
  "non_zero_multi_visit_histogram.png" = "MultiVisit_NonZero_Histogram.png",
  "non_zero_multi_visit_scatterplot.png" = "MultiVisit_NonZero_Scatterplot.png",
  "non_zero_single_visit_histogram.png" = "SingleVisit_NonZero_Histogram.png",
  "non_zero_single_visit_scatterplot.png" = "SingleVisit_NonZero_Scatterplot.png",
  "07_allvisits_histogram.png" = "Raw_AllVisits_Histogram.png",
  "07_allvisits_scatterplot.png" = "Raw_AllVisits_Scatterplot.png",
  "multivisit_median_vs_visits.png" = "MultiVisit_Scatterplot.png",
  "singlevisit_histogram_median_counts.png" = "SingleVisit_Histogram.png",
  "single_visit_median_vs_visits.png" = "SingleVisit_Scatterplot.png",
  "multivisit_histogram_median_counts.png" = "MultiVisit_Histogram.png"
)

# Rename the files
for (old_name in names(file_renames)) {
  old_path <- file.path(file_directory, old_name)
  new_path <- file.path(file_directory, file_renames[[old_name]])
  
  # Check if the old file exists before renaming
  if (file.exists(old_path)) {
    file.rename(old_path, new_path)
  } else {
    cat("File not found:", old_name, "\n")
  }
}





