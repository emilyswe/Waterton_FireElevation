################################# Updated March 21 2025
# Run models and collect results
results_list <- lapply(species_columns, function(species) compare_models(species, model_data))
names(results_list) <- species_columns

# Combine all model comparisons
bic_all <- rbindlist(lapply(results_list, function(x) x[[1]]), use.names = TRUE, fill = TRUE)
write.csv(bic_all, "Output/11_SingleSpecies/11_BestModelSelection.csv", row.names = FALSE)

# Combine all predictions
pred_all <- rbindlist(lapply(results_list, function(x) x[[2]]), use.names = TRUE, fill = TRUE)
write.csv(pred_all, "Output/11_SingleSpecies/11_Predictions.csv", row.names = FALSE)

# Done!
cat("✅ Model fitting complete. BIC and prediction outputs saved.")



###########plotting code
library(ggplot2)
library(dplyr)

# Load predictions
dat4 <- read.csv("Output/11_SingleSpecies/11_Predictions.csv") |> 
  rename(elevation = x, burnYN = group) |> 
  mutate(burnYN = factor(burnYN, levels = c("1", "0"), labels = c("Burned", "Unburned"))) |> 
  filter(!is.na(species), !model %in% c("burnONLY", "linearALT")) |> 
  arrange(species, elevation, burnYN)

# Define correct colors: Burned = red, Unburned = blue
correct_colors <- c("Burned" = "tomato3", "Unburned" = "steelblue3")

# Burn interaction models
plot_burn_int <- ggplot(dat4 |> filter(model %in% c("linearINT", "quadINT"))) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high, fill = burnYN), 
              linetype = "dashed", linewidth = 0.5, alpha = 0.2) +
  geom_line(aes(x = elevation, y = predicted, color = burnYN), linewidth = 1) +
  scale_color_manual(values = correct_colors, name = "") +
  scale_fill_manual(values = correct_colors, name = "") +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (meters)") +
  theme(legend.position = "bottom")

ggsave("Output/11_SingleSpecies/11_SingleSpecies_BurnInteraction.jpeg", plot = plot_burn_int, width = 10, height = 6)

# Burn additive models
plot_burn_main <- ggplot(dat4 |> filter(model %in% c("linearMAIN", "quadMAIN"),
                                        !species %in% c("CAVI"))) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high, fill = burnYN), 
              linetype = "dashed", linewidth = 0.5, alpha = 0.2) +
  geom_line(aes(x = elevation, y = predicted, color = burnYN), linewidth = 1) +
  scale_color_manual(values = correct_colors, name = "") +
  scale_fill_manual(values = correct_colors, name = "") +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (meters)") +
  theme(legend.position = "bottom")

ggsave("Output/11_SingleSpecies/11_SingleSpecies_BurnAdditive.jpeg", plot = plot_burn_main, width = 10, height = 6)

# Altitude-only models
plot_alt <- ggplot(dat4 |> filter(model %in% c("quadALT"))) +
  geom_ribbon(aes(x = elevation, ymin = conf.low, ymax = conf.high), 
              linetype = "dashed", linewidth = 0.5, alpha = 0.2) +
  geom_line(aes(x = elevation, y = predicted), linewidth = 1) +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (meters)")

ggsave("Output/11_SingleSpecies/11_SingleSpecies_Altitude.jpeg", plot = plot_alt, width = 10, height = 12)




























#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################



#Updated March 20 2025 

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

# Load dataset
single_species_data <- read.csv("Output/07_singlevisit_typeA.csv")

# Ensure burnYN is a factor
single_species_data$burnYN <- factor(single_species_data$burnYN, levels = c("0", "1"))

# Automatically detect species columns using regex (4-letter uppercase codes)
species_columns <- grep("^[A-Z]{4}$", colnames(single_species_data), value = TRUE)

# Convert abundance to presence/absence
single_species_data[species_columns] <- lapply(single_species_data[species_columns], function(x) ifelse(x > 0, 1, 0))

# Function to compare all models and select the best one
compare_models <- function(species, model_data) {
  if (!(species %in% names(model_data))) {
    return(NULL)
  }
  
  # Define model formulas
  formulas <- list(
    quadratic_interaction = as.formula(paste0(species, " ~ poly(elevation, 2) * burnYN")),
    linear_interaction = as.formula(paste0(species, " ~ elevation * burnYN")),
    quadratic_additive = as.formula(paste0(species, " ~ poly(elevation, 2) + burnYN")),
    linear_additive = as.formula(paste0(species, " ~ elevation + burnYN"))
  )
  
  # Fit models
  models <- lapply(formulas, function(f) {
    tryCatch(glm(f, data = model_data, family = binomial), error = function(e) NULL)
  })
  
  # Compute BIC for each model
  bics <- sapply(models, function(m) if (!is.null(m)) BIC(m) else NA)
  
  # Select the best model (lowest BIC)
  best_model_name <- names(which.min(bics))
  best_model <- models[[best_model_name]]
  
  # Prediction grid
  pred_data <- expand.grid(
    elevation = seq(min(model_data$elevation), max(model_data$elevation), length.out = 100),
    burnYN = levels(model_data$burnYN)
  )
  pred_data$predicted <- predict(best_model, newdata = pred_data, type = "response")
  
  return(list(
    best_model_name = best_model_name,
    bic_values = bics,
    pred_data = pred_data,
    model = best_model
  ))
}

# Run models for all species
species_results <- list()

for (species in species_columns) {
  cat("Running model for:", species, "\n")
  species_results[[species]] <- compare_models(species, single_species_data)
}

# Extract best model results
best_model_info <- data.frame(
  species = names(species_results),
  best_model = sapply(species_results, function(x) if (!is.null(x)) x$best_model_name else NA),
  BIC_values = sapply(species_results, function(x) if (!is.null(x)) min(x$bic_values, na.rm = TRUE) else NA)
)

# Save results
write.csv(best_model_info, "Output/11_SingleSpecies/11_BestModelSelection.csv", row.names = FALSE)

# Save predictions for plotting
all_predictions <- do.call(rbind, lapply(names(species_results), function(species) {
  if (!is.null(species_results[[species]]$pred_data)) {
    data.frame(species = species, species_results[[species]]$pred_data)
  }
}))

write.csv(all_predictions, "Output/11_SingleSpecies/11_Predictions.csv", row.names = FALSE)

# Print summary outputs
print(head(best_model_info))
print(head(all_predictions))






#####################slightly working plot code
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Load predictions
all_predictions <- read.csv("Output/11_SingleSpecies/11_Predictions.csv")

# Load best model selection info
best_models <- read.csv("Output/11_SingleSpecies/11_BestModelSelection.csv")

# Merge model info into predictions
all_predictions <- all_predictions %>%
  left_join(best_models, by = "species") %>%
  rename(model = best_model) %>%
  mutate(burnYN = factor(burnYN, levels = c("1", "0"), labels = c("Burned", "Unburned")))


#plot 1 burn x elevation interaction 
plot_burn_interaction <- ggplot(all_predictions %>%
                                  filter(model %in% c("quadratic_interaction", "linear_interaction"))) +
  geom_ribbon(aes(x = elevation, ymin = predicted - 0.05, ymax = predicted + 0.05, colour = burnYN), 
              linetype = "dashed", linewidth = 0.5, alpha = 0.2) +
  geom_line(aes(x = elevation, y = predicted, colour = burnYN), linewidth = 1) +
  scale_color_manual(values = c("steelblue3", "tomato3"), name = "") +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom")

ggsave(plot_burn_interaction, filename = "Output/11_SingleSpecies/11_SingleSpecies_TypeA_BurnInteraction.jpeg", width = 5, height = 3)

#plot 2 burn additive
plot_burn_additive <- ggplot(all_predictions %>%
                               filter(model %in% c("quadratic_additive", "linear_additive"),
                                      !species %in% c("CAVI"))) +
  geom_ribbon(aes(x = elevation, ymin = predicted - 0.05, ymax = predicted + 0.05, colour = burnYN), 
              linetype = "dashed", linewidth = 0.5, alpha = 0.2) +
  geom_line(aes(x = elevation, y = predicted, colour = burnYN), linewidth = 1) +
  scale_color_manual(values = c("steelblue3", "tomato3"), name = "") +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom")

ggsave(plot_burn_additive, filename = "Output/11_SingleSpecies/11_SingleSpecies_TypeA_BurnAdditive.jpeg", width = 10, height = 6)

#plot 3 altitude only
plot_altitude <- ggplot(all_predictions %>%
                          filter(model == "quadratic_additive")) +
  geom_ribbon(aes(x = elevation, ymin = predicted - 0.05, ymax = predicted + 0.05), 
              linetype = "dashed", linewidth = 0.5, alpha = 0.2) +
  geom_line(aes(x = elevation, y = predicted), linewidth = 1) +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)")

ggsave(plot_altitude, filename = "Output/11_SingleSpecies/11_SingleSpecies_TypeA_Altitude.jpeg", width = 10, height = 12)


































#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
##################### OLD CODE: 
#STEP 11_SINGLE SPECIES MODELS

#setwd
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

#load packages
library(ggplot2)
library(dplyr)

# Read in Cleaned Data ----
# Load the cleaned single-visit and multi-visit datasets
emda1 <- read.csv("Output/07_cleaned_multivisit_data.csv")
emda4 <- read.csv("Output/07_cleaned_single_visit_data.csv")

# Ensure the species columns are correctly identified using regex for 4-letter bird codes
species_columns <- colnames(emda4)[grep("^[A-Z]{4}$", colnames(emda4))] # Applies to the `emda4` dataset

#1. Wrangle altitude----
emda1$stdaltitude <- scale(emda1$elevation, center = TRUE, scale = TRUE)
emda1$stddoy <- emda1$doy - min(emda1$doy)
emda1$stdhod <- emda1$hod - min(emda1$hod)

#2. Function to run logistic regression and check convergence----
check_convergence <- function(species_name, data) {
  formula <- as.formula(paste(species_name, "~ altitude + factor(burnYN)  + factor(year)"))
  model <- tryCatch({
    glm(formula, data = data, family = binomial())
  }, error = function(e) {
    cat(paste(species_name, "- Error: ", e$message, "\n"))
    return(NULL)
  })
  
  # Check convergence
  if (!is.null(model) && !model$converged) {
    cat(paste(species_name, "- Dropped due to non-convergence\n"))
    return(NULL)
  } else if (!is.null(model)) {
    cat(paste(species_name, "- Kept\n"))
    return(species_name)
  }
}

#3. Remove species that don't converge on basic model----
glmer_data <- emda1 |> 
  rename(altitude = elevation)
glmer_data2 <- emda4 |> 
  rename(altitude = elevation)

# convert abundance to presence/absence using `species_columns`
glmer_data[species_columns] <- lapply(glmer_data[species_columns], function(x) ifelse(x > 0, 1, 0))
glmer_data2[species_columns] <- lapply(glmer_data2[species_columns], function(x) ifelse(x > 0, 1, 0))

# apply functions
converged_species <- sapply(species_columns, check_convergence, data = glmer_data)
converged_species2 <- sapply(species_columns, check_convergence, data = glmer_data2)

# filter out non-converged species
converged_species <- names(Filter(Negate(is.null), converged_species))
converged_species2 <- names(Filter(Negate(is.null), converged_species2))

glmer_data <- glmer_data[c("location", "hod", "doy", "year", "latitude", "longitude", "altitude", "burnYN", "burnBA", converged_species)] |> 
  data.frame()
glmer_data2 <- glmer_data[c("location", "hod", "doy", "year", "latitude", "longitude", "altitude", "burnYN", converged_species2)] |> 
  data.frame()

#4. Function to run mixed effect models----
run_models_re <- function(species, model_data) {
  
  # Models configuration 
  models <- list(
    quadINT = try(glmer(as.formula(paste(species, "~ poly(altitude, 2) * burnYN + burnBA + hod + doy + (1|location)")), data = model_data, family = binomial())),
    linearINT = try(glmer(as.formula(paste(species, "~ altitude * burnYN + burnBA + hod + doy + (1|location)")), data = model_data, family = binomial())),
    quadMAIN = try(glmer(as.formula(paste(species, "~ poly(altitude, 2) + burnBA + burnYN + hod + doy + (1|location)")), data = model_data, family = binomial())),
    linearMAIN = try(glmer(as.formula(paste(species, "~ altitude + burnBA + burnYN + hod + doy + (1|location)")), data = model_data, family = binomial())),
    quadALT = try(glmer(as.formula(paste(species, "~ poly(altitude, 2) + burnBA + hod + doy + (1|location)")), data = model_data, family = binomial())),
    linearALT = try(glmer(as.formula(paste(species, "~ altitude + burnBA + hod + doy + (1|location)")), data = model_data, family = binomial())),
    burnONLY = try(glmer(as.formula(paste(species, "~ burnYN + burnBA + hod + doy + (1|location)")), data = model_data, family = binomial())),
    base = try(glmer(as.formula(paste(species, "~ burnBA + hod + doy + (1|location)")), data = model_data, family = binomial()))
  )
  
  # Output AIC and BIC for each model, with error handling
  results <- lapply(models, function(model) {
    if (class(model)[1] == "try-error") {
      cat("Error in model fitting for", species, "\n")
      return(list(AIC = NA, BIC = NA))
    } else {
      return(list(AIC = AIC(model), BIC = BIC(model)))
    }
  })
  
  # Identify model with lowest BIC
  aicbic <- data.table::rbindlist(results) |> 
    mutate(model = names(models),
           species = species) |> 
    arrange(BIC)
  
  best <- models[aicbic$model[1]][[1]]
  
  # Make predictions from selected model
  if(aicbic$model[1] %in% names(models)[1:4]){
    preds <- predict_response(best, terms=c("altitude[all]", "burnYN")) |> 
      data.frame() |> 
      mutate(species = species,
             model = aicbic$model[1])
  }
  if(aicbic$model[1] %in% names(models)[5:6]){
    preds <- predict_response(best, terms=c("altitude[all]")) |> 
      data.frame() |> 
      mutate(species = species,
             model = aicbic$model[1])
  }
  if(aicbic$model[1] %in% names(models)[7]){
    preds <- predict_response(best, terms=c("burnYN")) |>
      data.frame() |> 
      mutate(species = species,
             model = aicbic$model[1])
    
  }
  if(aicbic$model[1] %in% names(models)[8]){
    preds <- data.frame(x=NA, predicted=NA, conf.low=NA, conf.high=NA, group=NA, species=NA, model=NA)
  }
  
  return(list(aicbic, preds))
}

#5. Function to run fixed effect models----
run_models <- function(species, model_data) {
  
  # Models configuration 
  models <- list(
    quadINT = try(glm(as.formula(paste(species, "~ poly(altitude, 2) * burnYN + factor(year)")), data = model_data, family = binomial())),
    linearINT = try(glm(as.formula(paste(species, "~ altitude * burnYN + factor(year)")), data = model_data, family = binomial())),
    quadMAIN = try(glm(as.formula(paste(species, "~ poly(altitude, 2) + burnYN + factor(year)")), data = model_data, family = binomial())),
    linearMAIN = try(glm(as.formula(paste(species, "~ altitude + burnYN + factor(year)")), data = model_data, family = binomial())),
    quadALT = try(glm(as.formula(paste(species, "~ poly(altitude, 2) + factor(year)")), data = model_data, family = binomial())),
    linearALT = try(glm(as.formula(paste(species, "~ altitude + factor(year)")), data = model_data, family = binomial())),
    burnONLY = try(glm(as.formula(paste(species, "~ burnYN + factor(year)")), data = model_data, family = binomial())),
    base = try(glm(as.formula(paste(species, "~ factor(year)")), data = model_data, family = binomial()))
  )
  
  # Output AIC and BIC for each model, with error handling
  results <- lapply(models, function(model) {
    if (class(model)[1] == "try-error") {
      cat("Error in model fitting for", species, "\n")
      return(list(AIC = NA, BIC = NA))
    } else {
      return(list(AIC = AIC(model), BIC = BIC(model), df = model$df.residual))
    }
  })
  
  # Identify model with lowest BIC
  aicbic <- data.table::rbindlist(results) |> 
    mutate(model = names(models),
           species = species) |> 
    arrange(BIC) |> 
    mutate(delta = BIC - lag(BIC),
           delta = ifelse(is.na(delta), 0, delta))
  
  select <- aicbic |> 
    dplyr::filter(delta < 2) |> 
    dplyr::filter(df==max(df)) |> 
    arrange(delta) |> 
    head(1)
  
  best <- models[select$model]
  
  if(length(best)==1){
    # Make predictions from selected model
    if(select$model[1] %in% names(models)[1:4]){
      preds <- predict_response(best[[1]], terms=c("altitude[all]", "burnYN")) |> 
        data.frame() |> 
        mutate(species = species,
               model = select$model[1])
    }
    if(select$model[1] %in% names(models)[5:6]){
      preds <- predict_response(best[[1]], terms=c("altitude[all]")) |> 
        data.frame() |> 
        mutate(species = species,
               model = select$model[1])
    }
    if(select$model[1] %in% names(models)[7]){
      preds <- predict_response(best[[1]], terms=c("burnYN")) |>
        data.frame() |> 
        mutate(species = species,
               model = select$model[1])
      
    }
    if(select$model[1] %in% names(models)[8]){
      preds <- data.frame(x=NA, predicted=NA, conf.low=NA, conf.high=NA, group=NA, species=NA, model=select$model[1])
    }
    
    return(list(aicbic, preds))
  }
  
}


#6. Run models----
# Running the combined models function for all species
#results <- lapply(converged_species, run_models_re, model_data = glmer_data)
#names(results) <- converged_species

results2 <- lapply(converged_species2, run_models, model_data = glmer_data2)
names(results2) <- converged_species2

#7. Package results for plotting----
#bic_results <- do.call(rbind, lapply(results, function(x) x[[1]]))
bic_results2 <- do.call(rbind, lapply(results2, function(x) x[[1]]))

#predictions <- lapply(results, function(x) x[[2]]) |> 
#  data.table::rbindlist(fill=TRUE, use.names=TRUE)

predictions2 <- lapply(results2, function(x) x[[2]]) |> 
  data.table::rbindlist(fill=TRUE, use.names=TRUE)

#Ended up using linear models with a year effect from the single visit dataset instead of the mixed effects ones on the multivisit dataset because the latter weren't making sense and didn't have enough power

write.csv(bic_results2, "Output/11_singlespecies_BIC_EMDA.csv", row.names = FALSE)

write.csv(predictions2, "Output/11_singlespecies_predictions_EMDA.csv", row.names = FALSE)

#8. Plot----

#Single species predictions EMDA
dat4 <- read.csv("Output/11_singlespecies_predictions_EMDA.csv") |> 
  rename(elevation = x, burnYN = group) |> 
  mutate(burnYN = factor(burnYN, levels=c("1", "0"), labels=c("Burned", "Unburned"))) |> 
  dplyr::filter(!is.na(species),
                !model %in% c("burnONLY", "linearALT")) |> 
  arrange(species, elevation, burnYN)

mod <- dat4 |> 
  dplyr::select(species, model) |> 
  unique()
table(mod$model)

#Single Species Burn Interaction EMDA
plot4.burnint <- ggplot(dat4 |> dplyr::filter(model %in% c("linearINT", "quadINT"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high, colour=burnYN), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted, colour=burnYN), linewidth=1) +
  scale_color_manual(values=c("steelblue3", "tomato3"),name="") +
  facet_wrap(~species, scales="free", ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom")
plot4.burnint

ggsave(plot4.burnint, filename="Output/11_SingleSpecies_BurnInteraction_EMDA.jpeg", width=5, height=3)

#Single Species Burn Additive EMDA
plot4.burnmain <- ggplot(dat4 |> dplyr::filter(model %in% c("linearMAIN", "quadMAIN"),
                                               !species %in% c("CAVI"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high, colour=burnYN), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted, colour=burnYN), linewidth=1) +
  scale_color_manual(values=c("steelblue3", "tomato3"),name="") +
  facet_wrap(~species, scales="free",  ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") + 
  theme(legend.position="bottom")
plot4.burnmain

ggsave(plot4.burnmain, filename="Output/11_SingleSpecies_BurnAdditive_EMDA.jpeg", width=10, height=6)

#Single Species Altitude EMDA
plot4.alt <- ggplot(dat4 |> dplyr::filter(model %in% c("quadALT"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted), linewidth=1) +
  facet_wrap(~species, scales="free",  ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)")
plot4.alt

ggsave(plot4.alt, filename="Output/11_SingleSpecies_Altitude_EMDA.jpeg", width=10, height=12)


#################################################################### changed colours of burn and unburned 

# 8. Plot ---- 

# Single species predictions EMDA
dat4 <- read.csv("Output/11_singlespecies_predictions_EMDA.csv") |> 
  rename(elevation = x, burnYN = group) |> 
  mutate(burnYN = factor(burnYN, levels=c("1", "0"), labels=c("Burned", "Unburned"))) |> 
  dplyr::filter(!is.na(species),
                !model %in% c("burnONLY", "linearALT")) |> 
  arrange(species, elevation, burnYN)

mod <- dat4 |> 
  dplyr::select(species, model) |> 
  unique()
table(mod$model)

# Single Species Burn Interaction EMDA
plot4.burnint <- ggplot(dat4 |> dplyr::filter(model %in% c("linearINT", "quadINT"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high, colour=burnYN), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted, colour=burnYN), linewidth=1) +
  scale_color_manual(values=c("tomato3", "steelblue3"),name="") +  # Switched colors here
  facet_wrap(~species, scales="free", ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom")
plot4.burnint

ggsave(plot4.burnint, filename="Output/11_SingleSpecies_BurnInteraction_EMDA.jpeg", width=5, height=3)

# Single Species Burn Additive EMDA
plot4.burnmain <- ggplot(dat4 |> dplyr::filter(model %in% c("linearMAIN", "quadMAIN"),
                                               !species %in% c("CAVI"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high, colour=burnYN), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted, colour=burnYN), linewidth=1) +
  scale_color_manual(values=c("tomato3", "steelblue3"),name="") +  # Switched colors here
  facet_wrap(~species, scales="free",  ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") + 
  theme(legend.position="bottom")
plot4.burnmain

ggsave(plot4.burnmain, filename="Output/11_SingleSpecies_BurnAdditive_EMDA.jpeg", width=10, height=6)

# Single Species Altitude EMDA
plot4.alt <- ggplot(dat4 |> dplyr::filter(model %in% c("quadALT"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted), linewidth=1) +
  facet_wrap(~species, scales="free",  ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)")
plot4.alt

ggsave(plot4.alt, filename="Output/11_SingleSpecies_Altitude_EMDA.jpeg", width=10, height=12)

####################################################################
####################################################################
