#STEP 11_SINGLE SPECIES MODELS

#setwd
setwd("/Users/Bronwyn/Documents/local-git/Waterton_FireElevation")

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

#4. Function to run mixed effect models with third-order and square-root polynomials----
run_models_re <- function(species, model_data) {
  
  # Models configuration 
  models <- list(
    third_order = try(glmer(as.formula(paste(species, "~ poly(altitude, 3) * burnYN + burnBA + hod + doy + (1|location)")), data = model_data, family = binomial())),
    sqrt_poly = try(glmer(as.formula(paste(species, "~ sqrt(altitude) * burnYN + burnBA + hod + doy + (1|location)")), data = model_data, family = binomial())),
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

#5. Function to run fixed effect models with third-order and square-root polynomials----
run_models <- function(species, model_data) {
  
  # Models configuration f
  models <- list(
    third_order = try(glm(as.formula(paste(species, "~ poly(altitude, 3) * burnYN + factor(year)")), data = model_data, family = binomial())),
    sqrt_poly = try(glm(as.formula(paste(species, "~ sqrt(altitude) * burnYN + factor(year)")), data = model_data, family = binomial())),
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

write.csv(bic_results2, "Output/11_singlespecies_BIC_POLYS.csv", row.names = FALSE)

write.csv(predictions2, "Output/11_singlespecies_predictions_POLYS.csv", row.names = FALSE)

#8. Plot----

#Single species predictions EMDA
dat4 <- read.csv("Output/11_singlespecies_predictions_POLYS.csv") |> 
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

ggsave(plot4.burnint, filename="Output/11_SingleSpecies_BurnInteraction_POLYS.jpeg", width=5, height=3)

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

ggsave(plot4.burnmain, filename="Output/11_SingleSpecies_BurnAdditive_POLYS.jpeg", width=10, height=6)

#Single Species Altitude EMDA
plot4.alt <- ggplot(dat4 |> dplyr::filter(model %in% c("quadALT"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted), linewidth=1) +
  facet_wrap(~species, scales="free",  ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)")
plot4.alt

ggsave(plot4.alt, filename="Output/11_SingleSpecies_Altitude_POLYS.jpeg", width=10, height=12)