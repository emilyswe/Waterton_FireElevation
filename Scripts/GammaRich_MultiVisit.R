setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# 1. Load required packages ----
packs <- c("dplyr", "ggplot2", "tidyverse")
for (q in 1:length(packs)) {
  if (!require(packs[q], character.only = TRUE)) {
    install.packages(packs[q])
    require(packs[q])
  }
}

# 2. Load data 
# Load the multi-visit dataset
emda1 <- read.csv("Output/07_cleaned_multivisit_data.csv")

# 3. Wrangle data ----
# Create location_data dataframe to store unique locations and related altitude
location_data <- data.frame(location = unique(emda1$location))

# Populate location_data with necessary columns from emda1
hml <- length(location_data$location)
for (q in 1:hml) {
  location_data$latitude[q] <- emda1$latitude[which(emda1$location == location_data$location[q])[1]]
  location_data$longitude[q] <- emda1$longitude[which(emda1$location == location_data$location[q])[1]]
  location_data$altitude[q] <- emda1$elevation[which(emda1$location == location_data$location[q])[1]]
}

# Compute tertiles for altitude
location_data$altbin <- cut(location_data$altitude, 
                            breaks = quantile(location_data$altitude, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                            labels = FALSE, include.lowest = TRUE)

# Merge strata and burn status data with emda1
superdata <- left_join(emda1, location_data, by = "location")

# Ensure strata assignment by altitude bin and burn status
superdata$strata <- superdata$altbin

# Verify strata and burnYN assignment
print(table(superdata$strata, superdata$burnYN))

# Identify species columns in the dataset
species_columns_emda1 <- colnames(emda1)[grep("^[A-Z]{4}$", colnames(emda1))]

# Verify that species columns are properly identified
print(species_columns_emda1)

# 4. Gamma Richness ----

# Function to calculate richness for a given number of samples within a specific stratum
calculate_richness <- function(data, n) {
  if (nrow(data) < n) {  # Check if the sample size is larger than available data
    n <- nrow(data)  # Adjust n to the maximum possible
  }
  sample_data <- data[sample(nrow(data), n, replace = TRUE), species_columns_emda1, drop = FALSE]
  unique_species_count <- sum(colSums(sample_data > 0, na.rm = TRUE) > 0)  # Count columns with any presence
  return(unique_species_count)
}

# Prepare a dataframe to store the results
gari <- data.frame(
  strata = integer(),
  replicate = integer(),
  numpc = integer(),
  richness = integer(),
  burnYN = integer()
)

# Set the number of sampling points and replicates
n_sampling_points <- c(25, 50, 75, 100)
n_replicates <- 100

# Loop over each stratum (elevation and burn status)
for (stratum in sort(unique(superdata$strata))) {
  for (burn_status in unique(superdata$burnYN)) {
    stratum_data <- superdata[superdata$strata == stratum & superdata$burnYN == burn_status, ]
    for (n in n_sampling_points) {
      for (rep in 1:n_replicates) {
        if (nrow(stratum_data) > 0) {  # Ensure there is data in this stratum
          richness_value <- calculate_richness(stratum_data, n)
          gari <- rbind(gari, data.frame(strata = stratum, replicate = rep, numpc = n, richness = richness_value, burnYN = burn_status))
        }
      }
    }
  }
}

# Summarize the richness statistics
gari_stats <- gari %>%
  group_by(strata, numpc, burnYN) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness = sd(richness, na.rm = TRUE),
    l95 = mean_richness - 1.96 * sd_richness,
    u95 = mean_richness + 1.96 * sd_richness,
    .groups = 'drop'
  )

# Save gamma richness results as a CSV
write.csv(gari_stats, "Output/09_gamma_richness_multivisit.csv", row.names = FALSE)

# Check the result
print(head(gari_stats))

# 5. Plot ----
# 5. Plot ----

# Function to generate plot for each stratum
generate_plot <- function(data, strata_id, burnYN_status) {
  ggplot(data[data$strata == strata_id & data$burnYN == burnYN_status, ], aes(x = numpc, y = mean_richness)) +
    geom_ribbon(aes(ymin = l95, ymax = u95), fill = ifelse(burnYN_status == 0, "green", "red"), alpha = 0.6) +
    geom_line(color = ifelse(burnYN_status == 0, "green", "red")) +
    labs(title = paste("Cumulative number of species -", ifelse(burnYN_status == 0, "Unburned", "Burned")),
         x = "Number of surveys", y = "Cumulative number of species") +
    theme_bw() +
    theme(legend.position = "none")
}

# Generate and store each plot (for burned and unburned separately)
plot_low_unburned <- generate_plot(gari_stats, 1, 0)
plot_mid_unburned <- generate_plot(gari_stats, 2, 0)
plot_high_unburned <- generate_plot(gari_stats, 3, 0)
plot_low_burned <- generate_plot(gari_stats, 1, 1)
plot_mid_burned <- generate_plot(gari_stats, 2, 1)
plot_high_burned <- generate_plot(gari_stats, 3, 1)

# Save the plots
ggsave("Output/09_gamma_richness_low_unburned_multivisit.png", plot = plot_low_unburned, width = 8, height = 6)
ggsave("Output/09_gamma_richness_mid_unburned_multivisit.png", plot = plot_mid_unburned, width = 8, height = 6)
ggsave("Output/09_gamma_richness_high_unburned_multivisit.png", plot = plot_high_unburned, width = 8, height = 6)
ggsave("Output/09_gamma_richness_low_burned_multivisit.png", plot = plot_low_burned, width = 8, height = 6)
ggsave("Output/09_gamma_richness_mid_burned_multivisit.png", plot = plot_mid_burned, width = 8, height = 6)
ggsave("Output/09_gamma_richness_high_burned_multivisit.png", plot = plot_high_burned, width = 8, height = 6)

# 6. Comparative Plots ----

# Function for comparison between unburned and burned
generate_comparative_plot <- function(data, unburned_strata, burned_strata) {
  ggplot() +
    geom_ribbon(data = data[data$strata == unburned_strata & data$burnYN == 0, ],
                aes(x = numpc, ymin = l95, ymax = u95, fill = "Unburned"), alpha = 0.5) +
    geom_line(data = data[data$strata == unburned_strata & data$burnYN == 0, ],
              aes(x = numpc, y = mean_richness, color = "Unburned")) +
    geom_ribbon(data = data[data$strata == burned_strata & data$burnYN == 1, ],
                aes(x = numpc, ymin = l95, ymax = u95, fill = "Burned"), alpha = 0.5) +
    geom_line(data = data[data$strata == burned_strata & data$burnYN == 1, ],
              aes(x = numpc, y = mean_richness, color = "Burned")) +
    scale_fill_manual(values = c("Unburned" = "green", "Burned" = "red")) +
    scale_color_manual(values = c("Unburned" = "green", "Burned" = "red")) +
    labs(title = "Comparison of Cumulative Species Richness - Unburned vs. Burned",
         x = "Number of Surveys",
         y = "Cumulative Number of Species") +
    theme_bw() +
    theme(legend.position = "bottom")
}

# Generate and store comparative plots
plot_low_comparison <- generate_comparative_plot(gari_stats, 1, 1)
plot_mid_comparison <- generate_comparative_plot(gari_stats, 2, 2)
plot_high_comparison <- generate_comparative_plot(gari_stats, 3, 3)

# Save the comparative plots
ggsave("Output/09_gamma_richness_low_comparison_multivisit.png", plot = plot_low_comparison, width = 8, height = 6)
ggsave("Output/09_gamma_richness_mid_comparison_multivisit.png", plot = plot_mid_comparison, width = 8, height = 6)
ggsave("Output/09_gamma_richness_high_comparison_multivisit.png", plot = plot_high_comparison, width = 8, height = 6)

