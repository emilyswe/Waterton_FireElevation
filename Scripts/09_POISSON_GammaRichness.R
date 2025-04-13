# Gamma Richness Script - Step 09
#Uses the single visit dataset

setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# 1. Load required packages ----
packs <- c("dplyr", "ggplot2", "tidyverse")
for (q in 1:length(packs)) {
  if (!require(packs[q], character.only = TRUE)) {
    install.packages(packs[q])
    require(packs[q])
  }
}

# 2. Read in the required data from previous steps ----
emda4 <- read.csv("Output/07_cleaned_single_visit_filtered.csv")

#emda1 <- read.csv("Output/07_cleaned_multivisit_data.csv")

# 3. Wrangle data ----

# Create elevation strata (bins based on altitude for gamma richness analysis)
location_data <- data.frame(location = unique(emda4$location))

# How Many Locations
hml <- length(location_data$location)

# Clunky but working solution to populate location_data and get strata
for (q in 1:hml) {
  location_data$latitude[q] <- emda4$latitude[which(emda4$location == location_data$location[q])[1]]
  location_data$longitude[q] <- emda4$longitude[which(emda4$location == location_data$location[q])[1]]
  location_data$altitude[q] <- emda4$elevation[which(emda4$location == location_data$location[q])[1]]
  location_data$severity[q] <- emda4$grid.code[which(emda4$location == location_data$location[q])[1]]
}

# Compute tertiles for altitude
location_data$altbin <- cut(location_data$altitude, breaks = quantile(location_data$altitude, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), labels = FALSE, include.lowest = TRUE)

# Adjust bins based on severity
location_data$strata <- ifelse(location_data$severity > 0, location_data$altbin + 3, location_data$altbin)

# Ensure all columns are correctly handled
superdata <- left_join(emda4, location_data, by = "location")

# Ensure the species columns are correctly identified
species_columns <- colnames(emda4)[grep("^[A-Z]{4}$", colnames(emda4))] # Regex to capture species column names using uppercase 4-letter bird codes

# 4. Calculate gamma richness ----

# Define the number of replicates and sampling points
n_replicates <- 100
n_sampling_points <- seq(10, 100, by = 10)

# Function to calculate richness for a given number of samples within a specific stratum
calculate_richness <- function(data, n) {
  if (nrow(data) < n) {  # Check if the sample size is larger than available data
    n <- nrow(data)  # Adjust n to the maximum possible
  }
  # Handle species columns safely to avoid missing values
  sample_data <- data[sample(nrow(data), n, replace = TRUE), species_columns, drop = FALSE]
  unique_species_count <- sum(colSums(sample_data > 0, na.rm = TRUE) > 0)  # Count columns with any presence
  return(unique_species_count)
}

# Prepare a dataframe to store the results
gari <- data.frame(
  strata = integer(),
  replicate = integer(),
  numpc = integer(),
  richness = integer()
)

# Loop over each stratum
strata <- sort(unique(superdata$strata))
for (stratum in strata) {
  stratum_data <- superdata[superdata$strata == stratum, ]
  for (n in n_sampling_points) {
    for (rep in 1:n_replicates) {
      if (nrow(stratum_data) > 0) {  # Ensure there is data in this stratum
        richness_value <- calculate_richness(stratum_data, n)
        gari <- rbind(gari, data.frame(strata = stratum, replicate = rep, numpc = n, richness = richness_value))
      }
    }
  }
}

# Rename columns appropriately
colnames(gari) <- c("strata", "replicate", "numpc", "richness")

# Calculate summary statistics with Poisson CIs
gari_stats <- gari %>%
  group_by(strata, numpc) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness = sd(richness, na.rm = TRUE),
    
    # Poisson-based lower and upper 95% CIs
    l95 = qpois(0.025, mean_richness),  # Lower 95% Poisson CI
    u95 = qpois(0.975, mean_richness),  # Upper 95% Poisson CI
    
    .groups = 'drop'
  )


# Label definitions for strata
strata_labels <- c("Low Elev. Unburned", "Mid. Elev. Unburned", "High Elev. Unburned",
                   "Low Elev. Burned", "Mid. Elev. Burned", "High Elev. Burned")
names(strata_labels) <- 1:6

# Add readable strata labels to the data
gari_stats$strata_label <- strata_labels[as.character(gari_stats$strata)]

# Save the summary results
write.csv(gari_stats, "Output/09_Gamma/filtered_poisson_gamma_richness_6strata.csv", row.names = FALSE)

# 5. Plot ----

# Function to generate plot
generate_plot <- function(data, strata_id) {
  ggplot(data[data$strata == strata_id, ], aes(x = numpc, y = mean_richness)) +
    geom_ribbon(aes(ymin = l95, ymax = u95), fill = ifelse(strata_id <= 3, "green", "red"), alpha = 0.6) +
    geom_line(color = ifelse(strata_id <= 3, "green", "red")) +
    labs(title = paste("Cumulative number of species -", strata_labels[strata_id]),
         x = "Number of surveys", y = "Cumulative number of species") +
    theme_bw() +
    theme(legend.position = "none")
}

# Generate and store each plot
plot_low_unburned <- generate_plot(gari_stats, 1)
plot_mid_unburned <- generate_plot(gari_stats, 2)
plot_high_unburned <- generate_plot(gari_stats, 3)
plot_low_burned <- generate_plot(gari_stats, 4)
plot_mid_burned <- generate_plot(gari_stats, 5)
plot_high_burned <- generate_plot(gari_stats, 6)

# Save the plots
ggsave("Output/09_Gamma/POISSON_gamma_richness_low_unburned.png", plot = plot_low_unburned, width = 8, height = 6)
ggsave("Output/09_Gamma/POISSON_gamma_richness_mid_unburned.png", plot = plot_mid_unburned, width = 8, height = 6)
ggsave("Output/09_Gamma/POISSON_gamma_richness_high_unburned.png", plot = plot_high_unburned, width = 8, height = 6)
ggsave("Output/09_Gamma/POISSON_gamma_richness_low_burned.png", plot = plot_low_burned, width = 8, height = 6)
ggsave("Output/09_Gamma/POISSON_gamma_richness_mid_burned.png", plot = plot_mid_burned, width = 8, height = 6)
ggsave("Output/09_Gamma/POISSON_gamma_richness_high_burned.png", plot = plot_high_burned, width = 8, height = 6)

## If there are warnings, investigate 
summary(gari_stats)
head(gari)
summary(gari)
summary(superdata)
head(species_columns)


# 6. Comparative Plots ----

generate_comparative_plot <- function(data, unburned_strata, burned_strata) {
  data$strata_label <- ifelse(data$strata == unburned_strata, "Unburned", ifelse(data$strata == burned_strata, "Burned", NA))
  
  ggplot() +
    geom_ribbon(data = data[data$strata %in% c(unburned_strata, burned_strata), ],
                aes(x = numpc, ymin = l95, ymax = u95, fill = strata_label),
                alpha = 0.5) +
    geom_line(data = data[data$strata %in% c(unburned_strata, burned_strata), ],
              aes(x = numpc, y = mean_richness, color = strata_label)) +
    scale_fill_manual(values = c("Unburned" = "green", "Burned" = "red")) +
    scale_color_manual(values = c("Unburned" = "green", "Burned" = "red")) +
    labs(title = "Comparison of Cumulative Species Richness - Unburned vs. Burned",
         x = "Number of Surveys",
         y = "Cumulative Number of Species") +
    theme_bw() +
    theme(legend.position = "bottom")
}

# Generate and store comparative plots
plot_low_comparison <- generate_comparative_plot(gari_stats, 1, 4)
plot_mid_comparison <- generate_comparative_plot(gari_stats, 2, 5)
plot_high_comparison <- generate_comparative_plot(gari_stats, 3, 6)

# Save the comparative plots
ggsave("Output/09_Gamma/POISSON_gamma_richness_low_comparison.png", plot = plot_low_comparison, width = 8, height = 6)
ggsave("Output/09_Gamma/POISSON_gamma_richness_mid_comparison.png", plot = plot_mid_comparison, width = 8, height = 6)
ggsave("Output/09_Gamma/POISSON_gamma_richness_high_comparison.png", plot = plot_high_comparison, width = 8, height = 6)

#summarize key results
key_results <- gari_stats %>%
  filter(numpc == max(numpc)) %>%
  select(strata_label, mean_richness, l95, u95)
print(key_results)






