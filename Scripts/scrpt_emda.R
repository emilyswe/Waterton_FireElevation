#PREAMBLE####

#1. Load packages----
#load('../current.RData')
packs<-c("haven",
         #"emmeans",
         "ggplot2","dplyr","broom","lubridate","vegan","lme4","tidyr","mgcv", "GGally", "tidyverse", "MuMIn", "data.table", "ggeffects")
for(q in 1:length(packs)) {
	if (!require(packs[q], character.only = TRUE)) 
		{
			install.packages(packs[q])  
			require(packs[q])
	}
}

#2. read in Erin's data----
final_data <- read.csv("G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_alpha_richness_perlocation_peryear.csv")

#3. read in Emily's data----

setwd("~/ABMI/Waterton_FireElevation")

#Take out PC historic, it's different from the other datasets it seems

#This is the data object with multiple visits per location
emda1.raw <- read.csv("Input/01_Wrangled.csv") |> 
  mutate(location = toupper(location),
         year = year(date_time),
         doy=yday(date_time),
         hod=hour(date_time),
         task_duration = as.numeric(str_sub(task_duration, -100, -2))) |> 
  dplyr::filter(hod < 10,
                hod > 4,
                doy >= 152,
                doy <= 182,
                longitude > -114.2)

#This is the data object with one visit per location (median count across visits)
emda4 <- read.csv("Input/04_Cleaned_Waterton.csv") |> 
  dplyr::filter(type != "PC_historic")

#4. Wrangle----

#Remove rare species from multi data, make wide, add elevation
spp.count <- emda1.raw |> 
  group_by(species) |> 
  summarize(counts = n()) |> 
  ungroup() |> 
  dplyr::filter(counts >= 5)

emda1 <- dplyr::filter(emda1.raw, species %in% spp.count$species) |> 
  pivot_wider(names_from="species", values_from="count", values_fill = 0) |> 
  left_join(emda4 |> 
              dplyr::select(latitude, longitude, elevation, slope, aspect, TPI, northness, grid.code, type) |>
              unique(),
            multiple="all") |> 
  dplyr::filter(type != "PC_historic")

#Check columns against Erin's data
colnames(final_data)[which(!colnames(final_data)%in%colnames(emda4))]
colnames(final_data)[which(!colnames(final_data)%in%colnames(emda1))]

#added the condition of "prior to 2017", else we'd assume places has burned prior to the actual fire |needed because Erin's data was all post-fire so it didn't matter then
emda1$burnYN <- as.factor(ifelse(((emda1$grid.code == 0)|(emda1$year<2018)),0,1))
emda4$burnYN <- as.factor(ifelse(((emda4$grid.code == 0)|(emda4$year<2018)),0,1))

#add before/afer burn variable
emda1$burnBA <- as.factor(ifelse(emda1$year < 2018, 0, 1))
emda4$burnBA <- as.factor(ifelse(emda4$year < 2018, 0, 1))

#5. Explore----
ggplot(emda1) +
  geom_histogram(aes(x=elevation)) +
  facet_grid(. ~ burnBA)

ggplot(emda1) + 
  geom_point(aes(x=longitude, y=latitude, colour=elevation)) +
  facet_grid(. ~ burnBA) +
  scale_colour_viridis_c()

#ALPHA RICHNESS####

#Uses single visit df for the GAMS and multi df for glms

#1. Calculate alpha richness----
emda1$alpharich<-1
emda4$alpharich<-1

#first and last species' columns
fsp4<-which(colnames(emda1)=="COLO")
lsp4<-which(colnames(emda1)=="STGR")
fsp1<-which(colnames(emda1)=="COLO")
lsp1<-which(colnames(emda1)=="CLSW")
emda1$location<-toupper(emda1$location)
emda4$location<-toupper(emda4$location)

for(q in 1:nrow(emda1)) emda1$alpharich[q]<-length(which(emda1[q,fsp1:lsp1]>0))

for(q in 1:nrow(emda4)) emda4$alpharich[q]<-length(which(emda4[q,fsp4:lsp4]>0))

#Not sure what this is doing
aggregated_data1 <- aggregate(alpharich ~ location + year + doy + hod, emda1, sum)
aggregated_data4 <- aggregate(alpharich ~ location + year, emda4, sum)

#2. Explore----
ggpairs(emda4 |> 
          dplyr::select(alpharich, elevation, burnYN, burnBA) |> 
          mutate(burnYN = as.factor(burnYN),
                 burnBA = as.factor(burnBA)))

ggplot(emda4) +
  geom_smooth(aes(x=elevation, y=alpharich, colour=factor(burnYN))) +
  geom_point(aes(x=elevation, y=alpharich, colour=factor(burnYN))) +
  facet_wrap(~type)

ggplot(emda4) +
  geom_histogram(aes(x=alpharich)) +
  facet_wrap(~type)

#3. Linear model----

#Model
model <- lmer(alpharich ~ elevation * burnYN + burnBA + doy + hod + (1 | location), data = emda1, na.action="na.fail")
summary(model)
dredge(model)

#New data
new_data <- expand.grid(
  elevation = seq(1200, 2200, by = 50),
  burnYN = as.factor(c(0, 1)),
  burnBA = as.factor(c(1)),
  year = mean(emda1$year, na.rm = TRUE),
  hod = mean(emda1$hod, na.rm = TRUE),
  doy = mean(emda1$doy, na.rm = TRUE)
)

#Predict
new_data$predicted <- predict(model, newdata = new_data, re.form = NA)

#Plot predictions
plot_preds <- ggplot(new_data, aes(x = elevation, y = predicted, color = factor(burnYN))) +
  geom_line() +
  scale_color_manual(values = c("green", "red"), labels = c("Unburned", "Burned")) +
  labs(title = "Effect of Burn on Alpha Diversity by Elevation", y = "Mean number of species", x = "Elevation - metres") +
  theme_minimal() +
  theme(legend.title = element_blank())
plot_preds
  
#4. GAM----

#Model
gam_NOBURN <- gam(alpharich ~ s(elevation, k = 4), data = emda4[emda4$burnYN == 0, ], method = "REML")
summary(gam_NOBURN)

gam_BURN <- gam(alpharich ~ s(elevation, k = 4), data = emda4[emda4$burnYN == 1, ], method = "REML")
summary(gam_BURN)

#Predict
pred_NOBURN <- predict(gam_NOBURN, newdata = new_data[new_data$burnYN == 0, ], se.fit = TRUE)
new_data$s_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit
new_data$l95_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit - 1.96 * pred_NOBURN$se.fit
new_data$u95_NOBURNalt[new_data$burnYN == 0] <- pred_NOBURN$fit + 1.96 * pred_NOBURN$se.fit

pred_BURN <- predict(gam_BURN, newdata = new_data[new_data$burnYN == 1, ], se.fit = TRUE)
new_data$s_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit
new_data$l95_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit - 1.96 * pred_BURN$se.fit
new_data$u95_BURNalt[new_data$burnYN == 1] <- pred_BURN$fit + 1.96 * pred_BURN$se.fit

#Save results
write.csv(new_data, "G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_alpha_richness_predictions_EMDA.csv", row.names = FALSE)

#Plot
alpharichvsburnele<-ggplot() +
    geom_ribbon(data = new_data[new_data$burnYN == 0, ], aes(x = elevation, ymin = l95_NOBURNalt, ymax = u95_NOBURNalt, fill = "Unburned"), alpha = 0.5) +
    geom_line(data = new_data[new_data$burnYN == 0, ], aes(x = elevation, y = s_NOBURNalt, color = "Unburned")) +
    geom_ribbon(data = new_data[new_data$burnYN == 1, ], aes(x = elevation, ymin = l95_BURNalt, ymax = u95_BURNalt, fill = "Burned"), alpha = 0.5) +
    geom_line(data = new_data[new_data$burnYN == 1, ], aes(x = elevation, y = s_BURNalt, color = "Burned")) +
    labs(y = "Mean number of species", x = "Elevation - metres", title = "") +
    scale_fill_manual(values = c("Unburned" = "green", "Burned" = "red")) +
    scale_color_manual(values = c("Unburned" = "green", "Burned" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))
alpharichvsburnele

#GAMMA RICHNESS########

#Uses the single visit dataset

#1. Create elevation strata----

#location   latitude longitude altitude severity
location_data<-data.frame(location=unique(emda4$location))
# How Many Locations
hml<-length(location_data$location)
#clunky but working solution to populate location_data and getting strata
for(q in 1:hml)
	{
	location_data$latitude[q]<-emda4$latitude[which(emda4$location==location_data$location[q])[1]]
	location_data$longitude[q]<-emda4$longitude[which(emda4$location==location_data$location[q])[1]]
	location_data$altitude[q]<-emda4$elevation[which(emda4$location==location_data$location[q])[1]]
	location_data$severity[q]<-emda4$grid.code[which(emda4$location==location_data$location[q])[1]]
	}

location_data$altbin <- cut(location_data$altitude, breaks = quantile(location_data$altitude, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), labels = FALSE, include.lowest = TRUE)

# Compute tertiles for altitude
location_data$altbin <- cut(location_data$altitude, breaks = quantile(location_data$altitude, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), labels = FALSE, include.lowest = TRUE)

# Adjust bins based on severity
location_data$strata <- ifelse(location_data$severity > 0, location_data$altbin + 3, location_data$altbin)

superdata<-left_join(emda1, location_data, by = "location")

species_columns <- colnames(emda1)[fsp1:lsp1]

#2. Calculate gamma richness----

# Define the number of replicates and sampling points
n_replicates <- 100
n_sampling_points <- seq(10, 100, by = 10)

# Function to calculate richness for given number of samples within a specific stratum
calculate_richness <- function(data, n) {
  if (nrow(data) < n) {  # Check if the sample size is larger than available data
    n <- nrow(data)  # Adjust n to the maximum possible
  }
  sample_data <- data[sample(nrow(data), n, replace = TRUE), species_columns]
  unique_species_count <- sum(colSums(sample_data > 0) > 0)  # Count columns with any presence
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
      richness_value <- calculate_richness(stratum_data, n)
      gari <- rbind(gari, c(stratum, rep, n, richness_value))
    }
  }
}

# Rename columns appropriately
colnames(gari) <- c("strata", "replicate", "numpc", "richness")

gari_stats <- gari %>%
  group_by(strata, numpc) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness = sd(richness, na.rm = TRUE),
    l95 = mean_richness - 1.96 * sd_richness,
    u95 = mean_richness + 1.96 * sd_richness,
    .groups = 'drop'
  )
  
# Label definitions for strata
strata_labels <- c("Low Elev. Unburned", "Mid. Elev. Unburned", "High Elev. Unburned",
                   "Low Elev. Burned", "Mid. Elev. Burned", "High Elev. Burned")
names(strata_labels) <- 1:6

# Add readable strata labels to the data
gari_stats$strata_label <- strata_labels[as.character(gari_stats$strata)]

write.csv(gari_stats, "G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_gamma_richness_6strata_EMDA.csv", row.names=FALSE)

#3. Plot----
generate_plot <- function(data, strata_id) {
  ggplot(data[data$strata == strata_id, ], aes(x = numpc, y = mean_richness)) +
    geom_ribbon(aes(ymin = l95, ymax = u95), fill = ifelse(strata_id <= 3, "green", "red"), alpha = 0.6) +
    geom_line(color = ifelse(strata_id <= 3, "green", "red")) +
    labs(title = paste("Cumulative number of species -", strata_labels[strata_id]),
         x = "Number of surveys", y = "Cumulative number of species") +
    theme_minimal() +
    theme(legend.position = "none")
}


# Create and store each plot
plot_low_unburned <- generate_plot(gari_stats, 1)
plot_mid_unburned <- generate_plot(gari_stats, 2)
plot_high_unburned <- generate_plot(gari_stats, 3)
plot_low_burned <- generate_plot(gari_stats, 4)
plot_mid_burned <- generate_plot(gari_stats, 5)
plot_high_burned <- generate_plot(gari_stats, 6)


# Function to generate comparative plots with correct colors
generate_comparative_plot <- function(data, unburned_strata, burned_strata) {
  # Ensure strata_label is added for clarity in plots
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
    theme_minimal() +
    theme(legend.position = "bottom")
}
# Generate and store plots for each elevation comparison
plot_low_comparison <- generate_comparative_plot(gari_stats, 1, 4)
plot_mid_comparison <- generate_comparative_plot(gari_stats, 2, 5)
plot_high_comparison <- generate_comparative_plot(gari_stats, 3, 6)

#NMDS#############

#1. Wrangle----

# Identify species columns by excluding non-species columns
non_species_columns <- setdiff(colnames(superdata), species_columns)

superdata$richness<-1

# Sort these species columns alphabetically
species_columns_sorted <- sort(species_columns)

# Arrange the DataFrame so species columns are sorted
wide_data <- superdata[c(non_species_columns, species_columns_sorted)]

# update first and last species' columns
fsp<-which(colnames(wide_data)==species_columns_sorted[1])
lsp<-which(colnames(wide_data)==species_columns_sorted[length(species_columns_sorted)])

# Calculate total richness across all species
wide_data$richness <- rowSums(wide_data[, species_columns_sorted], na.rm = TRUE)

# Sample one observation per location if needed
filtered_data <- wide_data %>%
  dplyr::filter(!is.na(strata)) |> 
  group_by(location) %>%
  slice_sample(n = 1)

#2. Model----
# Run NMDS using all species data
set.seed(999)
nmds_result <- metaMDS(filtered_data[, species_columns_sorted], 
                       distance = "jaccard", 
                       k = 2, 
                       trymax = 100)  # Increase trymax if convergence issues occur

#3. Save----
# Extract NMDS coordinates
coordinates <- scores(nmds_result, display = "sites")
filtered_data$mdsA <- coordinates[, 1]
filtered_data$mdsB <- coordinates[, 2]

write.csv(filtered_data[, c("location", "strata", "doy", "year", "mdsA", "mdsB")], "G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_mds_2axis_EMDA.csv", row.names = FALSE)

#4. Plot----
means <- aggregate(cbind(mdsA, mdsB) ~ strata, data = filtered_data, mean)

# Prepare plot 
nmdsplot<- ggplot(filtered_data, aes(x = mdsA, y = mdsB, color = factor(strata))) +
  geom_point(aes(shape = factor(strata))) +  # Different shapes for each strata
  scale_color_manual(values = c("green", "green", "green", "orange", "orange", "orange")) +  # Colors for unburned and burned
  geom_text(data = means, aes(label = factor(strata)), vjust = 2, color = "black") +  # Label means
  # Add ellipses for each strata group
  stat_ellipse(aes(fill = factor(strata)), geom = "polygon", level = 0.95, alpha = 0.2) +
  labs(color = "Strata", shape = "Strata", fill = "Confidence") +
  theme_minimal() +
  theme(legend.position = "right") +
  xlim(c(-3, 0))

#plot
plot(nmdsplot)

#SINGLE SPECIES MODELS####

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

# convert abundance to presence/absence
glmer_data[species_columns] <- lapply(glmer_data[species_columns], function(x) ifelse(x > 0, 1, 0))
glmer_data2[species_columns] <- lapply(glmer_data2[species_columns], function(x) ifelse(x > 0, 1, 0))
# apply functions
converged_species <- sapply(species_columns, check_convergence, data = glmer_data)
converged_species2 <- sapply(species_columns, check_convergence, data = glmer_data2)
#filter out non-converged species
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

#4. Function to run fixed effect models----
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


#5. Run models----
# Running the combined models function for all species
#results <- lapply(converged_species, run_models_re, model_data = glmer_data)
#names(results) <- converged_species

results2 <- lapply(converged_species2, run_models, model_data = glmer_data2)
names(results2) <- converged_species2

#6. Package results for plotting----
#bic_results <- do.call(rbind, lapply(results, function(x) x[[1]]))
bic_results2 <- do.call(rbind, lapply(results2, function(x) x[[1]]))

#predictions <- lapply(results, function(x) x[[2]]) |> 
#  data.table::rbindlist(fill=TRUE, use.names=TRUE)

predictions2 <- lapply(results2, function(x) x[[2]]) |> 
  data.table::rbindlist(fill=TRUE, use.names=TRUE)

#Ended up using linear models with a year effect from the single visit dataset instead of the mixed effects ones on the multivisit dataset because the latter weren't making sense and didn't have enough power

write.csv(bic_results2, "G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_singlespecies_BIC_EMDA.csv", row.names = FALSE)

write.csv(predictions2, "G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_singlespecies_predictions_EMDA.csv", row.names = FALSE)
