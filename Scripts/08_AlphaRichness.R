#Alpha richness modeling using single visit type A dataset and gams with additive model structure ####
#Updated March 25 2025 by Emily Swerdfager

# 1. Load Packages ----
library(ggplot2)
library(dplyr)
library(mgcv)
library(readr)

# 2. Read and Prep Data ----
emda4 <- read_csv("Output/07_singlevisit_typeA.csv") |>
  mutate(burnYN = as.factor(burnYN))

# Identify species columns (exactly four uppercase letters)
species_columns <- colnames(emda4)[grepl("^[A-Z]{4}$", colnames(emda4))]

# Calculate alpha richness
emda4$alpharich <- apply(emda4[, species_columns], 1, function(x) sum(x > 0, na.rm = TRUE))

# 3. Fit GAM Model ----
gam_model <- gam(alpharich ~ s(elevation, by = burnYN) + burnYN,
                 data = emda4, method = "REML", family = gaussian())
summary(gam_model)
saveRDS(gam_model, "Output/08_Alpha/08_gam_model_alpharich.rds")

# 4. Generate Predictions (manual) ----
# Create new data frame for predictions
new_data <- expand.grid(
  elevation = seq(min(emda4$elevation, na.rm = TRUE),
                  max(emda4$elevation, na.rm = TRUE),
                  by = 25),
  burnYN = factor(c(0, 1))
)

# Predict for unburned (burnYN == 0)
pred_NOBURN <- predict(gam_model, newdata = new_data[new_data$burnYN == 0, ], se.fit = TRUE)
new_data$s_NOBURNalt <- NA
new_data$l95_NOBURNalt <- NA
new_data$u95_NOBURNalt <- NA
new_data$s_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit
new_data$l95_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit - 1.96 * pred_NOBURN$se.fit
new_data$u95_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit + 1.96 * pred_NOBURN$se.fit

# Predict for burned (burnYN == 1)
pred_BURN <- predict(gam_model, newdata = new_data[new_data$burnYN == 1, ], se.fit = TRUE)
new_data$s_BURNalt <- NA
new_data$l95_BURNalt <- NA
new_data$u95_BURNalt <- NA
new_data$s_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit
new_data$l95_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit - 1.96 * pred_BURN$se.fit
new_data$u95_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit + 1.96 * pred_BURN$se.fit

# Save predictions
write_csv(new_data, "Output/08_Alpha/08_alpha_richness_predictions.csv")

# 5. Elly-Style Plot ----
plot1 <- ggplot() +
  geom_ribbon(data = new_data[new_data$burnYN == 0, ],
              aes(x = elevation, ymin = l95_NOBURNalt, ymax = u95_NOBURNalt, colour = "Unburned"),
              alpha = 0.2, linetype = "dashed", linewidth = 0.5) +
  geom_line(data = new_data[new_data$burnYN == 0, ],
            aes(x = elevation, y = s_NOBURNalt, color = "Unburned"),
            linewidth = 1) +
  geom_ribbon(data = new_data[new_data$burnYN == 1, ],
              aes(x = elevation, ymin = l95_BURNalt, ymax = u95_BURNalt, colour = "Burned"),
              alpha = 0.2, linetype = "dashed", linewidth = 0.5) +
  geom_line(data = new_data[new_data$burnYN == 1, ],
            aes(x = elevation, y = s_BURNalt, color = "Burned"),
            linewidth = 0.5) +
  labs(y = "Mean number of species", x = "Elevation (metres)", title = "") +
  scale_color_manual(values = c("Unburned" = "steelblue3", "Burned" = "tomato3")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))

# Save final plot with white background + high resolution
ggsave(
  filename = "Output/08_Alpha/08_alpha_richness_EllyPlot.png",
  plot = plot1,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"  # This fixes the dark background issue
)

# Display plot in RStudio
plot1










####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
# old code for ALPHA RICHNESS referring to both the multi (lme) and single visit (gam) datasets ####

# 1. Load Required Packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(mgcv)
library(MuMIn)
library(ggeffects)
library(GGally)


# 2. Read in filtered cleaned single visit data ----
emda4 <- read.csv("Output/07_singlevisit_typeA.csv")

# Convert burnYN and burnBA to factors
emda4 <- emda4 |> 
  mutate(
    burnYN = as.factor(burnYN),
    burnBA = as.factor(burnBA)
  )

str(emda4$burnYN)
str(emda4$burnBA)

# 3. Calculate Alpha Richness ----
# Alpha richness represents the species diversity at each site
# Initialize alpha richness column with a default value
emda4$alpharich <- 1

# Identify the species columns using regex (exactly four uppercase letters)
species_columns_emda4 <- colnames(emda4)[grep("^[A-Z]{4}$", colnames(emda4))]

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
ggsave(filename = "Output/08_Alpha/08_filtered_pairwise_plots_single_visit_data.png", plot = pairwise_plot, width = 10, height = 8)


# Scatter plot with smoothed lines for alpha richness by elevation and burn status
smoothed_scatterplot <- ggplot(emda4) +
  geom_smooth(aes(x = elevation, y = alpharich, colour = factor(burnYN))) +
  geom_point(aes(x = elevation, y = alpharich, colour = factor(burnYN))) +
  facet_wrap(~type)

# Save the plot
ggsave(filename = "Output/08_Alpha/08_filtered_smoothed_scatterplot_alpharich.png", plot = smoothed_scatterplot, width = 10, height = 8)

# Histogram of alpha richness for single-visit data
histogram_singlevisit_alpharich <- ggplot(emda4) +
  geom_histogram(aes(x = alpharich)) +
  facet_wrap(~type)

ggsave(filename = "Output/08_Alpha/08_filtered_histogram_singlevist_alpharich.png", plot = histogram_singlevisit_alpharich, width = 10, height = 8)

# 5. Linear Model for Multi-Visit Data ----
# Model the effect of elevation and burn status on alpha richness using a linear mixed-effects model
model <- lmer(alpharich ~ elevation * burnYN + (1 | location), 
              data = emda4, na.action = "na.fail")

# Output the summary of the model to assess fit and coefficients
summary(model)

saveRDS(model, "Output/08_Alpha/08_filtered_singlevisit_lmer_alpharich.rds")

#save dredge glm single visit results
dredge_results_df <- as.data.frame(dredge_results)
write.csv(dredge_results_df, "Output/08_Alpha/08_filtered_dredge_results.csv", row.names = FALSE)


# Create new data for predictions
# Generate a range of elevations and include burn status for prediction purposes
new_data <- expand.grid(
  elevation = seq(min(emda4$elevation, na.rm = TRUE), max(emda4$elevation, na.rm = TRUE), by = 50),
  burnYN = as.factor(c(0, 1))  # Burn status: 0 = Unburned, 1 = Burned
)

# Predict alpha richness using the linear mixed-effects model
new_data$predicted <- predict(model, newdata = new_data, re.form = NA)

# Plot the predicted effect of burn status on alpha diversity by elevation
plot_preds <- ggplot(new_data, aes(x = elevation, y = predicted, color = factor(burnYN))) +
  geom_line(size = 1) +
  scale_color_manual(values = c("green", "red"), labels = c("Unburned", "Burned")) +
  labs(
    title = "Effect of Burn on Alpha Diversity by Elevation",
    y = "Predicted Alpha Richness",
    x = "Elevation (metres)"
  ) +
  theme_classic() +  # Consistent theme for a clean look
  theme(legend.title = element_blank(), legend.position = "bottom")

# Save the plot
ggsave(
  filename = "Output/08_Alpha/08_filtered_predburnstatus_alpharich_byelev.png",
  plot = plot_preds,
  width = 8,
  height = 6
)

# Display the plot
plot_preds

# 6. Generalized Additive Models (GAMs) for Single-Visit Data ----

# Fit GAMs to model the non-linear relationship between elevation and alpha richness

# GAM for unburned sites
gam_NOBURN <- gam(alpharich ~ s(elevation, k = 4), data = emda4[emda4$burnYN == 0, ], method = "REML")
summary(gam_NOBURN)

# Save the summary of the GAM for unburned sites
capture.output(summary(gam_NOBURN), file = "Output/08_Alpha/08_gam_summary_unburned.txt")

# GAM for burned sites
gam_BURN <- gam(alpharich ~ s(elevation, k = 4), data = emda4[emda4$burnYN == 1, ], method = "REML")
summary(gam_BURN)

# Save the summary of the GAM for burned sites
capture.output(summary(gam_BURN), file = "Output/08_Alpha/08_gam_summary_burned.txt")

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
write.csv(new_data, "Output/08_Alpha/08_filtered_alpha_richness_predictions.csv", row.names = FALSE)

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

# Save the final GAM plot
ggsave(filename = "Output/08_Alpha/08_filtered_alpharich_vs_elevation_burn_status.png", plot = alpharichvsburnele, width = 8, height = 6)


