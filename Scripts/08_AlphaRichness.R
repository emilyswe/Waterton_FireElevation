# ALPHA RICHNESS ####

# 1. Load Required Packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(mgcv)
library(MuMIn)
library(ggeffects)

# 2. Read in Cleaned Data ----
# Load the cleaned single-visit and multi-visit datasets
emda1 <- read.csv("Output/07_cleaned_multivisit_data.csv")
emda4 <- read.csv("Output/07_cleaned_single_visit_data.csv")

# 3. Calculate Alpha Richness ----
# Alpha richness represents the species diversity at each site
# Initialize alpha richness column with a default value
emda1$alpharich <- 1
emda4$alpharich <- 1

# Identify the species columns using regex (exactly four uppercase letters)
species_columns_emda1 <- colnames(emda1)[grep("^[A-Z]{4}$", colnames(emda1))]
species_columns_emda4 <- colnames(emda4)[grep("^[A-Z]{4}$", colnames(emda4))]

# Calculate alpha richness for multi-visit data (emda1)
for (q in 1:nrow(emda1)) {
  emda1$alpharich[q] <- length(which(emda1[q, species_columns_emda1] > 0))
}

# Calculate alpha richness for single-visit data (emda4)
for (q in 1:nrow(emda4)) {
  emda4$alpharich[q] <- length(which(emda4[q, species_columns_emda4] > 0))
}

# 4. Explore Alpha Richness ----
# Visualize alpha richness against elevation and burn status
# Pairwise plots for single-visit data
pairwise_plot <- ggpairs(emda4 |> 
                           dplyr::select(alpharich, elevation, burnYN, burnBA) |> 
                           mutate(burnYN = as.factor(burnYN),
                                  burnBA = as.factor(burnBA)))

# Save the plot
ggsave(filename = "Output/08_pairwise_plots_single_visit_data.png", plot = pairwise_plot, width = 10, height = 8)


# Scatter plot with smoothed lines for alpha richness by elevation and burn status
smoothed_scatterplot <- ggplot(emda4) +
  geom_smooth(aes(x = elevation, y = alpharich, colour = factor(burnYN))) +
  geom_point(aes(x = elevation, y = alpharich, colour = factor(burnYN))) +
  facet_wrap(~type)

# Save the plot
ggsave(filename = "Output/08_smoothed_scatterplot_alpharich.png", plot = smoothed_scatterplot, width = 10, height = 8)

# Histogram of alpha richness for single-visit data
histogram_singlevisit_alpharich <- ggplot(emda4) +
  geom_histogram(aes(x = alpharich)) +
  facet_wrap(~type)

ggsave(filename = "Output/08_histogram_singlevist_alpharich.png", plot = histogram_singlevisit_alpharich, width = 10, height = 8)

# 5. Linear Model for Multi-Visit Data ----
# Model the effect of elevation and burn status on alpha richness using a linear mixed-effects model
model <- lmer(alpharich ~ elevation * burnYN + burnBA + doy + hod + (1 | location), 
              data = emda1, na.action = "na.fail")

# Output the summary of the model to assess fit and coefficients
summary(model)

# Use model dredging to explore all possible models and rank them by AICc
dredge(model)

# Create a new dataset for predictions
# This dataset covers a range of elevations and burn statuses for prediction purposes
new_data <- expand.grid(
  elevation = seq(1200, 2200, by = 50),  # Range of elevations
  burnYN = as.factor(c(0, 1)),           # Burn status: 0 = Unburned, 1 = Burned
  burnBA = as.factor(c(0, 1)),           # Include both pre-burn (0) and post-burn (1) statuses
  year = mean(emda1$year, na.rm = TRUE), # Mean year for prediction
  hod = mean(emda1$hod, na.rm = TRUE),   # Mean hour of day for prediction
  doy = mean(emda1$doy, na.rm = TRUE)    # Mean day of year for prediction
)

# Predict alpha richness using the linear model
new_data$predicted <- predict(model, newdata = new_data, re.form = NA)

# Plot the predicted effect of burn status on alpha diversity by elevation
plot_preds <- ggplot(new_data, aes(x = elevation, y = predicted, color = factor(burnYN))) +
  geom_line() +
  scale_color_manual(values = c("green", "red"), labels = c("Unburned", "Burned")) +
  labs(title = "Effect of Burn on Alpha Diversity by Elevation", 
       y = "Mean number of species", 
       x = "Elevation - metres") +
  theme_classic() +  # Use a classic white background
  theme(legend.title = element_blank())  # Optional: further customize legend

# Save the plot
ggsave(filename = "Output/08_predburnstatus_alpharich_byelev.png", plot = plot_preds, width = 8, height = 6)

# 6. Generalized Additive Models (GAMs) for Single-Visit Data ----
# Fit GAMs to model the non-linear relationship between elevation and alpha richness
# GAM for unburned sites
gam_NOBURN <- gam(alpharich ~ s(elevation, k = 4), data = emda4[emda4$burnYN == 0, ], method = "REML")
summary(gam_NOBURN)

# GAM for burned sites
gam_BURN <- gam(alpharich ~ s(elevation, k = 4), data = emda4[emda4$burnYN == 1, ], method = "REML")
summary(gam_BURN)

# Predict alpha richness for GAMs
# For unburned sites
pred_NOBURN <- predict(gam_NOBURN, newdata = new_data[new_data$burnYN == 0, ], se.fit = TRUE)
new_data$s_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit
new_data$l95_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit - 1.96 * pred_NOBURN$se.fit
new_data$u95_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit + 1.96 * pred_NOBURN$se.fit

# For burned sites
pred_BURN <- predict(gam_BURN, newdata = new_data[new_data$burnYN == 1, ], se.fit = TRUE)
new_data$s_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit
new_data$l95_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit - 1.96 * pred_BURN$se.fit
new_data$u95_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit + 1.96 * pred_BURN$se.fit

# Save the predictions to a CSV for further analysis
write.csv(new_data, "Output/08_alpha_richness_predictions.csv", row.names = FALSE)

# Plot the GAM predictions
alpharichvsburnele <- ggplot() +
  geom_ribbon(data = new_data[new_data$burnYN == 0, ], aes(x = elevation, ymin = l95_NOBURNalt, ymax = u95_NOBURNalt, fill = "Unburned"), alpha = 0.5) +
  geom_line(data = new_data[new_data$burnYN == 0, ], aes(x = elevation, y = s_NOBURNalt, color = "Unburned")) +
  geom_ribbon(data = new_data[new_data$burnYN == 1, ], aes(x = elevation, ymin = l95_BURNalt, ymax = u95_BURNalt, fill = "Burned"), alpha = 0.5) +
  geom_line(data = new_data[new_data$burnYN == 1, ], aes(x = elevation, y = s_BURNalt, color = "Burned")) +
  labs(y = "Mean number of species", x = "Elevation - metres", title = "") +
  scale_fill_manual(values = c("Unburned" = "green", "Burned" = "red")) +
  scale_color_manual(values = c("Unburned" = "green", "Burned" = "red")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))
alpharichvsburnele

# Save the final GAM plot
ggsave(filename = "Output/08_alpharich_vs_elevation_burn_status.png", plot = alpharichvsburnele, width = 8, height = 6)
