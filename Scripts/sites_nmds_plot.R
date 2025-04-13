# Load NMDS results
library(tidyverse)

dat3 <- read_csv("Output/10_NMDS/waterton_nmds_2axis_EMDA.csv", show_col_types = FALSE) |>
  filter(!location %in% c("WLNP-12-1", "WLNP-22-5")) |>  # remove visual outliers
  mutate(
    treatment = ifelse(strata %in% c(1:3), "Unburned", "Burned"),
    elevation = case_when(
      strata %in% c(1, 4) ~ "Low elevation",
      strata %in% c(2, 5) ~ "Mid elevation",
      strata %in% c(3, 6) ~ "High elevation"
    ),
    elevation = factor(elevation, levels = c("Low elevation", "Mid elevation", "High elevation"))
  )

# Make plot with site points and ellipses
plot3_sites <- ggplot(dat3, aes(x = mdsA, y = mdsB)) +
  geom_point(aes(color = treatment), alpha = 0.5, size = 2) +
  stat_ellipse(aes(colour = treatment, linetype = elevation), level = 0.68, linewidth = 1) +
  scale_color_manual(values = c("Burned" = "tomato3", "Unburned" = "steelblue3"), name = "") +
  scale_linetype_manual(values = c("solid", "dotted", "dashed"), name = "") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "NMDS axis 1", y = "NMDS axis 2")

# Save it
ggsave("Output/10_NMDS/10_nmds_plot_with_sites_FINAL_FIXED.png", plot = plot3_sites,
       width = 8, height = 6, dpi = 300, bg = "white")

