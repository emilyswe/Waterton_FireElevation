# STEP 11: Single Species Models – Type A
#March 24 2025

# Set working directory
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

# Load packages
library(ggplot2)
library(dplyr)
library(data.table)
library(ggeffects)

# Load cleaned dataset (already filtered to Type A species)
data <- read.csv("Output/07_singlevisit_typeA.csv")

# Rename and standardize
data <- data |> 
  rename(altitude = elevation) |> 
  mutate(burnYN = factor(burnYN, levels = c(0, 1)),
         year = factor(year))

# Identify species columns (4-letter bird codes)
species_columns <- grep("^[A-Z]{4}$", colnames(data), value = TRUE)

# Convert to presence/absence
data[species_columns] <- lapply(data[species_columns], function(x) ifelse(x > 0, 1, 0))

# --- Step 1: Convergence check ----
check_convergence <- function(species_name, df) {
  f <- as.formula(paste(species_name, "~ altitude + burnYN + year"))
  model <- tryCatch(glm(f, data = df, family = binomial()), error = function(e) NULL)
  if (!is.null(model) && model$converged) {
    cat("✅", species_name, "kept\n")
    return(species_name)
  } else {
    cat("❌", species_name, "dropped\n")
    return(NULL)
  }
}

converged_species <- unlist(lapply(species_columns, check_convergence, df = data))

# --- Step 2: Model fitting and prediction function ----
run_models <- function(species, df) {
  models <- list(
    quadINT   = try(glm(as.formula(paste(species, "~ poly(altitude, 2) * burnYN + year")), df, family = binomial())),
    linearINT = try(glm(as.formula(paste(species, "~ altitude * burnYN + year")), df, family = binomial())),
    quadMAIN  = try(glm(as.formula(paste(species, "~ poly(altitude, 2) + burnYN + year")), df, family = binomial())),
    linearMAIN= try(glm(as.formula(paste(species, "~ altitude + burnYN + year")), df, family = binomial())),
    quadALT   = try(glm(as.formula(paste(species, "~ poly(altitude, 2) + year")), df, family = binomial())),
    linearALT = try(glm(as.formula(paste(species, "~ altitude + year")), df, family = binomial())),
    burnONLY  = try(glm(as.formula(paste(species, "~ burnYN + year")), df, family = binomial())),
    base      = try(glm(as.formula(paste(species, "~ year")), df, family = binomial()))
  )
  
  # Model selection using BIC
  results <- lapply(models, function(m) {
    if (inherits(m, "try-error")) return(list(AIC = NA, BIC = NA, df = NA))
    list(AIC = AIC(m), BIC = BIC(m), df = m$df.residual)
  })
  
  aicbic <- rbindlist(results) |> 
    mutate(model = names(models), species = species) |> 
    arrange(BIC) |> 
    mutate(delta = BIC - min(BIC, na.rm = TRUE))
  
  best_model_name <- aicbic$model[1]
  best_model <- models[[best_model_name]]
  
  # Predict using ggeffects
  if (best_model_name %in% c("quadINT", "linearINT", "quadMAIN", "linearMAIN")) {
    preds <- ggeffects::predict_response(best_model, terms = c("altitude [all]", "burnYN")) |> 
      as.data.frame() |> 
      rename(elevation = x, predicted = predicted, conf.low = conf.low, conf.high = conf.high, burnYN = group) |> 
      mutate(species = species, model = best_model_name)
  } else if (best_model_name %in% c("quadALT", "linearALT")) {
    preds <- ggeffects::predict_response(best_model, terms = c("altitude [all]")) |> 
      as.data.frame() |> 
      rename(elevation = x, predicted = predicted, conf.low = conf.low, conf.high = conf.high) |> 
      mutate(burnYN = NA, species = species, model = best_model_name)
  } else if (best_model_name == "burnONLY") {
    preds <- ggeffects::predict_response(best_model, terms = c("burnYN")) |> 
      as.data.frame() |> 
      rename(elevation = x, predicted = predicted, conf.low = conf.low, conf.high = conf.high, burnYN = group) |> 
      mutate(species = species, model = best_model_name)
  } else {
    preds <- data.frame(elevation = NA, predicted = NA, conf.low = NA, conf.high = NA, burnYN = NA, species = species, model = best_model_name)
  }
  
  return(list(aicbic = aicbic, predictions = preds))
}

# --- Step 3: Run model loop ----
results <- lapply(converged_species, function(species) run_models(species, df = data))
names(results) <- converged_species

# Save model selection results
bic_all <- rbindlist(lapply(results, function(x) x$aicbic), fill = TRUE)
write.csv(bic_all, "Output/11_SingleSpecies/11_BestModelSelection.csv", row.names = FALSE)

# Save predictions
pred_all <- rbindlist(lapply(results, function(x) x$predictions), fill = TRUE)
write.csv(pred_all, "Output/11_SingleSpecies/11_Predictions.csv", row.names = FALSE)

cat("✅ Models complete. Files saved to Output/11_SingleSpecies/\n")






# Load packages
library(ggplot2)
library(dplyr)
library(readr)

# Load predictions
dat4 <- read_csv("Output/11_SingleSpecies/11_Predictions.csv") |>
  mutate(
    burnYN = factor(burnYN, levels = c(1, 0), labels = c("Burned", "Unburned")),
    model = as.character(model)
  ) |>
  filter(!is.na(species), !model %in% c("burnONLY", "linearALT")) |>
  arrange(species, elevation, burnYN)

# Custom color palette
burn_colors <- c("Burned" = "tomato3", "Unburned" = "steelblue3")

# --- Plot 1: Burn Interaction (quadINT / linearINT) ---
plot_burnint <- ggplot(dat4 |> filter(model %in% c("linearINT", "quadINT"))) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high, color = burnYN), 
              linetype = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_line(aes(x = elevation, y = predicted, color = burnYN), linewidth = 1) +
  scale_color_manual(values = burn_colors, name = "") +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom")

ggsave("Output/11_SingleSpecies/11_SingleSpecies_BurnInteraction_Fixed.jpeg", 
       plot_burnint, width = 10, height = 6)

# --- Plot 2: Burn Additive (quadMAIN / linearMAIN) ---
plot_burnmain <- ggplot(dat4 |> filter(model %in% c("linearMAIN", "quadMAIN"))) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high, color = burnYN), 
              linetype = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_line(aes(x = elevation, y = predicted, color = burnYN), linewidth = 1) +
  scale_color_manual(values = burn_colors, name = "") +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom")

ggsave("Output/11_SingleSpecies/11_SingleSpecies_BurnAdditive_Fixed.jpeg", 
       plot_burnmain, width = 10, height = 6)

# --- Plot 3: Altitude Only (quadALT) ---
plot_alt <- ggplot(dat4 |> filter(model == "quadALT")) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high), 
              fill = "grey80", alpha = 0.3, linetype = "dashed") +
  geom_line(aes(x = elevation, y = predicted), color = "black", linewidth = 1) +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "none")

ggsave("Output/11_SingleSpecies/11_SingleSpecies_Altitude_Fixed.jpeg", 
       plot_alt, width = 10, height = 12)

cat("✅ All plots saved to Output/11_SingleSpecies/\n")


###############
###############
#Updated April 8 2025
#Combined plots code:

# Load packages
library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)

# Load prediction data
dat4 <- read_csv("Output/11_SingleSpecies/11_Predictions.csv") |>
  mutate(
    burnYN = factor(burnYN, levels = c(1, 0), labels = c("Burned", "Unburned")),
    model = as.character(model)
  ) |>
  filter(!is.na(species), !model %in% c("burnONLY", "linearALT")) |>
  arrange(species, elevation, burnYN)

burn_colors <- c("Burned" = "tomato3", "Unburned" = "steelblue3")

# Panel A: Altitude-only
plot_alt <- ggplot(dat4 |> filter(model == "quadALT")) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high), 
              fill = "grey80", alpha = 0.3, linetype = "dashed") +
  geom_line(aes(x = elevation, y = predicted), color = "black", linewidth = 1) +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "none") +
  ggtitle("A. Altitude Model")

# Panel B: Burn Additive
plot_burnmain <- ggplot(dat4 |> filter(model %in% c("linearMAIN", "quadMAIN"))) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high, color = burnYN), 
              linetype = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_line(aes(x = elevation, y = predicted, color = burnYN), linewidth = 1) +
  scale_color_manual(values = burn_colors, name = "") +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom") +
  ggtitle("B. Burn Additive Model")

# Panel C: Burn x Elevation Interaction
plot_burnint <- ggplot(dat4 |> filter(model %in% c("linearINT", "quadINT"))) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high, color = burnYN), 
              linetype = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_line(aes(x = elevation, y = predicted, color = burnYN), linewidth = 1) +
  scale_color_manual(values = burn_colors, name = "") +
  facet_wrap(~species, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom") +
  ggtitle("C. Interaction Model")

# Combine into single figure
final_combined <- plot_alt / (plot_burnmain | plot_burnint) + 
  plot_layout(heights = c(1.3, 1))  # Make Panel A a bit taller for readability

# Save to file
ggsave("Output/11_SingleSpecies/11_SingleSpecies_FinalCombinedPlot.jpeg", 
       final_combined, width = 12, height = 14)