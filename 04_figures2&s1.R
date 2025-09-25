#### generate figures 2 and suppl. figure 1  ####

rm(list = ls())

# load libraries
library(tidyverse)
# Run git lfs pull to fetch large files before loading
system("git lfs pull")
# load data
load("tmp/tolerance_specific_pumas.rda")
load("tmp/ind_coefs_and_weights.rda")
# Define a function to get coefficient values from weighted_df
get_coefficients <- function(weighted_df, var, quantile, sex) {
  unname(weighted_df[which(weighted_df$var == var & weighted_df$sex == sex), 
                     quantile])
}
#### investigate specific pumas with high tolerance ####
high_males<- coef_db %>% 
  filter(grepl("M", pumaID)) %>%
  mutate(tolerance_category = factor(case_when(
    # high (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "m") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "f") ~ "high",
    
    # low (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "m") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "f") ~ "low",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("high", "mixed", "low"))) %>%
  filter(tolerance_category == "high") %>%
  dplyr::select(pumaID)
high_females<- coef_db %>% 
  filter(grepl("F", pumaID)) %>%
  mutate(tolerance_category = factor(case_when(
    
    # high (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "f") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "f") ~ "high",
    
    # low (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "f") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "f") ~ "low",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("high", "mixed", "low")))  %>%
  filter(tolerance_category == "high") %>%
  dplyr::select(pumaID)

print(high_males)
print(high_females)

#### make figure 2 ####
# Add sex identifier and combine datasets
real_m_pumas_combined <- real_m_pumas_combined_forplot %>% mutate(Sex = "Male")
real_f_pumas_combined <- real_f_pumas_combined_forplot %>% mutate(Sex = "Female")
combined_long <- bind_rows(real_m_pumas_combined, real_f_pumas_combined) %>%
  dplyr::select(-c(pumas:sl_scaled)) %>%
  rename(`Housing Density_mean` = hd_150_mean,
         `Housing Density_sd` = hd_150_sd,
         `Dist. to Urban Edge_mean` = prox_urbanedge_mean,
         `Dist. to Urban Edge_sd` = prox_urbanedge_sd,
         `Human Landcover_mean` = landcover_simp2,
         `Human Landcover_sd` = landcover_simphuman_sd,
         `Veg. Cover_mean` = cover_mean,
         `Veg. Cover_sd` = cover_sd,
         `Anthropogenic Tolerance` = tolerance_category,
         Slope_mean = slope_mean,
         Slope_sd = slope_sd) %>%
  pivot_longer(
    cols = -c(Sex, `Anthropogenic Tolerance`),              # pivot all covariate columns
    names_to = c("covariate", ".value"),           # split names into covariate and a value type
    names_sep = "_"                                # underscore separates covariate name and type (mean or sd)
  ) %>%
  mutate(`Anthropogenic Tolerance` = factor(case_when(`Anthropogenic Tolerance` == "low" ~ "low",
                                                      `Anthropogenic Tolerance` == "mixed" ~ "medium",
                                                      `Anthropogenic Tolerance` == "high" ~ "high"), levels = 
                                              c("low", "medium", "high")),
         `Covariate Type` = factor(ifelse(covariate %in% c("Veg. Cover", "Slope"),
                                          "Natural", "Anthropogenic")),
         mean = ifelse(`Covariate Type` == "Natural",
                       mean*10, mean),
         sd = ifelse(`Covariate Type` == "Natural",
                     sd*10, sd))
pos_dodge <- position_dodge(width = 0.6, preserve = "total")
p <- ggplot(combined_long, aes(
  x = mean,
  y = covariate,
  color = `Anthropogenic Tolerance`,
  group = `Anthropogenic Tolerance`  # <- key fix
)) +
  geom_errorbarh(
    aes(xmin = mean - 1.96 * sd, xmax = mean + 1.96 * sd),
    height = 0.2,
    size = 0.6,
    position = pos_dodge
  ) +
  geom_point(size = 3, position = pos_dodge) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  facet_grid(`Covariate Type`~ Sex, scales = "free") +
  labs(
    x = "Relative Selection", y = "Covariate",
    subtitle = "Weighted Coefficient Estimates by Sex and Tolerance"
  ) +
  theme_bw(base_size = 24) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 24),      # Increase legend text size
    legend.title = element_text(size = 24)
  ) +
  scale_color_manual(
    name = "Anthropogenic Tolerance",
    values = RColorBrewer::brewer.pal(3, "Dark2"),
    breaks = c("high", "medium", "low")
  ) +
  scale_x_continuous(
    breaks = c(-10, 0, 10),
    labels = c("avoid", "0", "select")
  ) +
  ggtitle("Figure 2")
p
ggsave("output/figure2.png", plot = p, width = 16, height = 8, dpi = 300)

#### create supplementary figures figure ####
categorized_males<- coef_db %>% 
  filter(grepl("M", pumaID)) %>%
  mutate(tolerance_category = factor(case_when(
    # high (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "m") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "m") ~ "high",
    
    # low (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "m") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "m") ~ "low",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("high", "mixed", "low"))) %>%
  mutate(sex = "male")
categorized_females<- coef_db %>% 
  filter(grepl("F", pumaID)) %>%
  mutate(tolerance_category = factor(case_when(
    
    # high (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "f") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "f") ~ "high",
    
    # low (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "f") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "f") ~ "low",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("high", "mixed", "low")))  %>%
  mutate(sex = "female")
theory_cats<- rbind(real_f_pumas_combined, real_m_pumas_combined) %>%
  rename(sex = Sex,
         prox_urbanedge = prox_urbanedge_mean,
         landcover_simphuman = landcover_simp2) %>%
  dplyr::select(tolerance_category, sex, prox_urbanedge, landcover_simphuman) %>%
  mutate(cluster = "theoretical average",
         sex = tolower(sex))
cats<- rbind(categorized_females, categorized_males) %>%
  dplyr::select(tolerance_category, sex, prox_urbanedge, landcover_simphuman) %>%
  mutate(cluster = "actual")
all_cats<- rbind(theory_cats, cats) %>%
  mutate(
    # Recode tolerance_category as a factor with new labels
    disturbance = fct_recode(tolerance_category,
                             "low" = "low",
                             "medium" = "mixed",
                             "high" = "high"),
    # Capitalize sex for facet labels
    sex = recode(sex,
                 "male" = "Males",
                 "female" = "Females")
  )
thresholds <- tibble::tibble(
  sex = rep(c("Females", "Males"), each = 3 * 2),  # 3 percentiles x 2 vars
  var = rep(c("prox_urbanedge", "landcover_simphuman"), times = 3, each = 2),
  label = rep(c("p40", "p50", "p60"), each = 2, times = 2),
  value = c(
    # Female prox_urbanedge
    get_coefficients(weighted_df, "prox_urbanedge", "p40", "f"),
    get_coefficients(weighted_df, "landcover_simphuman", "p60", "f"),
    get_coefficients(weighted_df, "prox_urbanedge", "p50", "f"),
    get_coefficients(weighted_df, "landcover_simphuman", "p50", "f"),
    get_coefficients(weighted_df, "prox_urbanedge", "p60", "f"),
    get_coefficients(weighted_df, "landcover_simphuman", "p40", "f"),
    
    # Male prox_urbanedge
    get_coefficients(weighted_df, "prox_urbanedge", "p40", "m"),
    get_coefficients(weighted_df, "landcover_simphuman", "p60", "m"),
    get_coefficients(weighted_df, "prox_urbanedge", "p50", "m"),
    get_coefficients(weighted_df, "landcover_simphuman", "p50", "m"),
    get_coefficients(weighted_df, "prox_urbanedge", "p60", "m"),
    get_coefficients(weighted_df, "landcover_simphuman", "p40", "m")
  )
)
# create line (MAKE THIS SEX SPECIFIC)
# Filter
theoretical_cats <- filter(all_cats, cluster == "theoretical average")
# Fit linear model
mod_m <- lm(landcover_simphuman ~ prox_urbanedge, data = theoretical_cats %>% 
              filter(sex == "Males"))
summary(mod_m)
# Fit linear model females
mod_f <- lm(landcover_simphuman ~ prox_urbanedge, data = theoretical_cats %>% 
              filter(sex == "Females"))
summary(mod_m)
summary(mod_f)
# Create prediction data over full x-range of plot
newdata <- data.frame(
  prox_urbanedge = seq(-10, 
                       max(all_cats$prox_urbanedge), 
                       length.out = 100)
)
# Add predictions and confidence intervals
preds_m <- predict(mod_m, newdata, interval = "confidence")
preds_f <- predict(mod_f, newdata, interval = "confidence")
newdata_m <- cbind(newdata, as.data.frame(preds_m)) %>%
  mutate(sex = "Males")
newdata_f <- cbind(newdata, as.data.frame(preds_f)) %>%
  mutate(sex = "Females")
newdata<- rbind(newdata_f, newdata_m)

# make plot
p2<- all_cats %>%
  ggplot(aes(x = prox_urbanedge, y = landcover_simphuman)) +
  geom_line(data = newdata, 
            aes(x = prox_urbanedge, y = fit), 
            inherit.aes = FALSE,
            color = "grey40", linetype = "dashed", size = 1) +
  # Shaded high tolerance zone ( for females)
  geom_rect(
    data = data.frame(sex = "Females"),
    aes(xmin = -Inf, xmax = get_coefficients(weighted_df, "prox_urbanedge", "p50", "f"), 
        ymin = get_coefficients(weighted_df, "landcover_simphuman", "p50", "f"), ymax = Inf),
    fill = "red", alpha = 0.1, inherit.aes = FALSE
  ) +
  # Shaded low tolerance zone ( for females)
  geom_rect(
    data = data.frame(sex = "Females"),
    aes(xmax = Inf, xmin = get_coefficients(weighted_df, "prox_urbanedge", "p60", "f"), 
        ymax = get_coefficients(weighted_df, "landcover_simphuman", "p40", "f"), ymin = -Inf),
    fill = "blue", alpha = 0.1, inherit.aes = FALSE
  ) +
  # Shaded high tolerance zone ( for females)
  geom_rect(
    data = data.frame(sex = "Males"),
    aes(xmin = -Inf, xmax = get_coefficients(weighted_df, "prox_urbanedge", "p50", "m"), 
        ymin = get_coefficients(weighted_df, "landcover_simphuman", "p50", "m"), ymax = Inf),
    fill = "red", alpha = 0.1, inherit.aes = FALSE
  ) +
  # Shaded low tolerance zone ( for males)
  geom_rect(
    data = data.frame(sex = "Males"),
    aes(xmax = Inf, xmin = get_coefficients(weighted_df, "prox_urbanedge", "p60", "m"), 
        ymax = get_coefficients(weighted_df, "landcover_simphuman", "p40", "m"), ymin = -Inf),
    fill = "blue", alpha = 0.1, inherit.aes = FALSE
  ) +
  # Hollow circles around theoretical points for emphasis
  geom_point(
    data = filter(all_cats, cluster == "theoretical average"),
    aes(color = disturbance),
    size = 6, shape = 1, stroke = 1.2, alpha = 0.7
  ) +
  
  # Main points for all pumas (theoretical and real)
  geom_point(
    aes(color = disturbance, 
        shape = cluster, 
        alpha = cluster),
    size = 3
  ) +
  
  # Adjust alpha manually
  scale_alpha_manual(values = c(actual = 0.5, `theoretical average` = 1)) +
  
  # Optional: make shapes more informative
  scale_shape_manual(values = c(actual = 16, `theoretical average` = 17)) +
  
  facet_wrap(~ sex, scales = "free") +
  theme_minimal() +
  labs(
    x = "Distance to Urban Edge",
    y = "Simplified Human Landcover",
    color = "Anthropogenic Tolerance",
    shape = "Puma Type",
    alpha = "Puma Type",
    subtitle = "Coefficient Estimates Relative to Tolerance"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  ) +
  ggtitle("Figure S1") +
  theme_bw(base_size = 24) +
  theme(
    strip.text = element_text(face = "bold"),      # high facet titles
    panel.grid = element_blank()                   # remove background grid for clarity
  ) +
  xlim(c(-10,5))
p2
ggsave("output/supp_fig1.png", plot = p2, width = 16, height = 8, dpi = 300)

