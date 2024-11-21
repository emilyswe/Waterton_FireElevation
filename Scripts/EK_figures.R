library(tidyverse)

setwd("G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/")

#1. Alpha richness----

dat1 <- read.csv("G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_alpha_richness_predictions_EMDA.csv")

plot1 <- ggplot() +
  geom_ribbon(data = dat1[dat1$burnYN == 0, ], aes(x = elevation, ymin = l95_NOBURNalt, ymax = u95_NOBURNalt, colour = "Unburned"), alpha = 0.2, linetype="dashed", linewidth=0.5) +
  geom_line(data = dat1[dat1$burnYN == 0, ], aes(x = elevation, y = s_NOBURNalt, color = "Unburned"), linewidth=1) +
  geom_ribbon(data = dat1[dat1$burnYN == 1, ], aes(x = elevation, ymin = l95_BURNalt, ymax = u95_BURNalt, colour = "Burned"), alpha = 0.2, linetype="dashed", linewidth=0.5) +
  geom_line(data = dat1[dat1$burnYN == 1, ], aes(x = elevation, y = s_BURNalt, color = "Burned"), linewidth=0.5) +
  labs(y = "Mean number of species", x = "Elevation (metres)", title = "") +
  scale_color_manual(values = c("Unburned" = "steelblue3", "Burned" = "tomato3")) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))
plot1

ggsave(plot1, filename="FIGURES/AlphaRichness_EMDA.jpeg", width=6, height=5)

#2. Gamma richness----

dat2 <-read.csv("G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_gamma_richness_6strata_EMDA.csv") |> 
  mutate(treatment = ifelse(strata %in% c(1:3), "Unburned", "Burned"),
         elevation = case_when(strata %in% c(1, 4) ~ "Low elevation", 
                               strata %in% c(2, 5) ~ "Mid elevation",
                               strata %in% c(3, 6) ~ "High elevation"),
         elevation = factor(elevation, levels=c("Low elevation", "Mid elevation", "High elevation")))

plot2 <- ggplot(dat2) +
  geom_ribbon(aes(x=numpc, ymin=l95, ymax=u95, colour=treatment), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=numpc, y=mean_richness, colour=treatment), linewidth=1) +
  scale_color_manual(values=c("steelblue3", "tomato3"),name="") +
  facet_wrap(~elevation) +
  theme_minimal() +
  labs(x = "Number of Surveys",
       y = "Cumulative Number of Species") +
  theme(legend.position = "bottom")
plot2

ggsave(plot2, filename="G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/FIGURES/GammaRichness_EMDA.jpeg", width=10, height=5)

#3. NMDS----

dat3 <- read.csv("G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_mds_2axis_EMDA.csv") |> 
  mutate(treatment = ifelse(strata %in% c(1:3), "Unburned", "Burned"),
         elevation = case_when(strata %in% c(1, 4) ~ "Low elevation", 
                               strata %in% c(2, 5) ~ "Mid elevation",
                               strata %in% c(3, 6) ~ "High elevation"),
         elevation = factor(elevation, levels=c("Low elevation", "Mid elevation", "High elevation")))

linelegend <- data.frame(expand.grid(x=1, y=1, elevation=unique(dat3$elevation)))

plot3 <- ggplot(dat3, aes(x = mdsA, y = mdsB)) +
  stat_ellipse(aes(colour=treatment, linetype=elevation), geom = "polygon", level = 0.68, fill="white", alpha = 0, linewidth=1) +
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

ggsave(plot3, filename="G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/FIGURES/NMDS_EMDA.jpeg", width=8, height=6)

#4. Single species----

dat4 <- read.csv("G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/SUMMARIZED/waterton_singlespecies_predictions_EMDA.csv") |> 
  rename(elevation = x, burnYN = group) |> 
  mutate(burnYN = factor(burnYN, levels=c("1", "0"), labels=c("Burned", "Unburned"))) |> 
  dplyr::filter(!is.na(species),
                !model %in% c("burnONLY", "linearALT")) |> 
  arrange(species, elevation, burnYN)

mod <- dat4 |> 
  dplyr::select(species, model) |> 
  unique()
table(mod$model)

plot4.burnint <- ggplot(dat4 |> dplyr::filter(model %in% c("linearINT", "quadINT"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high, colour=burnYN), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted, colour=burnYN), linewidth=1) +
  scale_color_manual(values=c("steelblue3", "tomato3"),name="") +
  facet_wrap(~species, scales="free", ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)") +
  theme(legend.position = "bottom")
plot4.burnint

ggsave(plot4.burnint, filename="G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/FIGURES/SingleSpecies_BurnInteraction_EMDA.jpeg", width=5, height=3)

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

ggsave(plot4.burnmain, filename="G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/FIGURES/SingleSpecies_BurnAdditive_EMDA.jpeg", width=10, height=6)

plot4.alt <- ggplot(dat4 |> dplyr::filter(model %in% c("quadALT"))) +
  geom_ribbon(aes(x=elevation, ymin=conf.low, ymax=conf.high), linetype="dashed", linewidth=0.5, alpha=0.2) +
  geom_line(aes(x=elevation, y=predicted), linewidth=1) +
  facet_wrap(~species, scales="free",  ncol=4) +
  theme_minimal() +
  labs(y = "Probability of occurrence", x = "Elevation (metres)")
plot4.alt

ggsave(plot4.alt, filename="G:/.shortcut-targets-by-id/1Fe3s33yF_NClLt-5Q52EVLrXNN4nT2uZ/WATERTON_FIRE/2022_TALK/DATA/FIGURES/SingleSpecies_Altitude_EMDA.jpeg", width=10, height=12)
