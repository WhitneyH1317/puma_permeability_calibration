
# load libraries
rm(list = ls())
library(sf)
library(lubridate)
library(mapview)
library(tidyverse)
library(raster)
# useful function
`%!in%`<- Negate('%in%')
# load data
load("tmp/data_for_projections.rda")
load("tmp/ind_coefs_and_weights.rda")
# Create a logical mask of non-NA cells
valid_cells <- !is.na(getValues(stk[[1]])) & !is.na(getValues(lulc_current_reclassed))

### generate personality-specific predictions ### 
# sex-specific predictions
fem<- weighted_df %>%
  filter(sex == "f")
mal<- weighted_df %>%
  filter(sex == "m")
# Define a function to get coefficient values from weighted_df
get_coefficients <- function(weighted_df, var, quantile, sex) {
  unname(weighted_df[which(weighted_df$var == var & weighted_df$sex == sex), 
                     quantile])
}

#####################################################################
# investigate female/male covariate spread
test<- weighted_df %>% 
  pivot_longer(., cols = c("min":"max"), names_to = "percentile") %>% 
  filter(var %in% c("hd_150", "prox_urbanedge", "slope", "cover", "landcover_simphuman")) 
test$percentile<- factor(test$percentile, levels = c("min", "p10", "p20", "p30", "p40", "p50", "avg", "p60", "p70", "p80", "p90", "max"))
test %>% 
  filter(percentile %!in% c("avg", "min", "max", "p10", "p90")) %>% 
  ggplot(., aes(x = percentile, y = value, group = sex, color = sex)) + 
  geom_point() + 
  facet_wrap(~ Covariate, scales = "free")


#### FEMALE REAL COEFFS ####
### combination ###
(real_f_pumas_combined <- coef_db %>% 
   filter(grepl("F", pumaID)) %>%
   mutate(shyness_category = factor(case_when(
     
     # Bold (selects human landcover and is close to urban edge)
     prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "f") & 
       landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "f") ~ "bold",
     
     # shy (avoids landcover and urban edge)
     prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "f") &
       landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "f") ~ "shy",
     
     # Neutral
     TRUE ~ "mixed"
   ), levels = c("bold", "mixed", "shy")))  %>%
   group_by(shyness_category) %>%
   na.omit() %>%
   na.omit() %>%
   summarize(
     across(c(hd_150:landcover_simphuman),
            list(mean = ~mean(.), sd = ~sd(.)),
            .names = "{.col}_{.fn}"),
     pumas = n(),
     .groups = "drop"
   ) %>%
   mutate(log_sl = weighted_df$avg[which(weighted_df$var == "log_sl" & weighted_df$sex == "f")],
          ta_scaled = 0,
          sl_scaled = 0) %>%
   arrange(shyness_category)) 
#### MALE REAL COEFFS ####
# Create a combined dataset
(real_m_pumas_combined <- coef_db %>% 
   filter(grepl("M", pumaID)) %>%
   mutate(shyness_category = factor(case_when(
     # Bold (selects human landcover and is close to urban edge)
     prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "m") & 
       landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "m") ~ "bold",
     
     # shy (avoids landcover and urban edge)
     prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "m") &
       landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "m") ~ "shy",
     
     # Neutral
     TRUE ~ "mixed"
   ), levels = c("bold", "mixed", "shy"))) %>%
   group_by(shyness_category) %>%
   na.omit() %>%
   summarize(
     across(c(hd_150:landcover_simphuman),
            list(mean = ~mean(.), sd = ~sd(.)),
            .names = "{.col}_{.fn}"),
     pumas = n(),
     .groups = "drop"
   ) %>%
   mutate(log_sl = weighted_df$avg[which(weighted_df$var == "log_sl" & weighted_df$sex == "f")],
          ta_scaled = 0,
          sl_scaled = 0) %>%
   arrange(shyness_category)) 

# investigate specific identifies
bold_males<- coef_db %>% 
  filter(grepl("M", pumaID)) %>%
  mutate(shyness_category = factor(case_when(
    # Bold (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "m") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "f") ~ "bold",
    
    # shy (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "m") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "f") ~ "shy",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("bold", "mixed", "shy"))) %>%
  filter(shyness_category == "bold") %>%
  dplyr::select(pumaID)
bold_females<- coef_db %>% 
  filter(grepl("F", pumaID)) %>%
  mutate(shyness_category = factor(case_when(
    
    # Bold (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "f") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "f") ~ "bold",
    
    # shy (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "f") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "f") ~ "shy",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("bold", "mixed", "shy")))  %>%
  filter(shyness_category == "bold") %>%
  dplyr::select(pumaID)

print(bold_males)
print(bold_females)

# make figure 2
# Add sex identifier and combine datasets
real_m_pumas_combined <- real_m_pumas_combined %>% mutate(Sex = "Male")
real_f_pumas_combined <- real_f_pumas_combined %>% mutate(Sex = "Female")
combined_long <- bind_rows(real_m_pumas_combined, real_f_pumas_combined) %>%
  dplyr::select(-c(pumas:sl_scaled)) %>%
  rename(`Housing Density_mean` = hd_150_mean,
         `Housing Density_sd` = hd_150_sd,
         `Dist. to Urban Edge_mean` = prox_urbanedge_mean,
         `Dist. to Urban Edge_sd` = prox_urbanedge_sd,
         `Human Landcover_mean` = landcover_simphuman_mean,
         `Human Landcover_sd` = landcover_simphuman_sd,
         `Veg. Cover_mean` = cover_mean,
         `Veg. Cover_sd` = cover_sd,
         `Anthropogenic Tolerance` = shyness_category,
         Slope_mean = slope_mean,
         Slope_sd = slope_sd) %>%
  pivot_longer(
    cols = -c(Sex, `Anthropogenic Tolerance`),              # pivot all covariate columns
    names_to = c("covariate", ".value"),           # split names into covariate and a value type
    names_sep = "_"                                # underscore separates covariate name and type (mean or sd)
  ) %>%
  mutate(`Anthropogenic Tolerance` = factor(case_when(`Anthropogenic Tolerance` == "shy" ~ "low",
                                                      `Anthropogenic Tolerance` == "mixed" ~ "medium",
                                                      `Anthropogenic Tolerance` == "bold" ~ "high"), levels = 
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
ggsave("~/output/coef_var_tolerance_sex.png", plot = p, width = 16, height = 8, dpi = 300)

# make different figure
categorized_males<- coef_db %>% 
  filter(grepl("M", pumaID)) %>%
  mutate(shyness_category = factor(case_when(
    # Bold (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "m") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "m") ~ "bold",
    
    # shy (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "m") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "m") ~ "shy",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("bold", "mixed", "shy"))) %>%
  mutate(sex = "male")
categorized_females<- coef_db %>% 
  filter(grepl("F", pumaID)) %>%
  mutate(shyness_category = factor(case_when(
    
    # Bold (selects human landcover and is close to urban edge)
    prox_urbanedge < get_coefficients(weighted_df, "prox_urbanedge", "p50", "f") & 
      landcover_simphuman > get_coefficients(weighted_df, "landcover_simphuman", "p50", "f") ~ "bold",
    
    # shy (avoids landcover and urban edge)
    prox_urbanedge > get_coefficients(weighted_df, "prox_urbanedge", "p60", "f") &
      landcover_simphuman < get_coefficients(weighted_df, "landcover_simphuman", "p40", "f") ~ "shy",
    
    # Neutral
    TRUE ~ "mixed"
  ), levels = c("bold", "mixed", "shy")))  %>%
  mutate(sex = "female")
theory_cats<- rbind(real_f_pumas_combined, real_m_pumas_combined) %>%
  rename(sex = Sex,
         prox_urbanedge = prox_urbanedge_mean,
         landcover_simphuman = landcover_simphuman_mean) %>%
  dplyr::select(shyness_category, sex, prox_urbanedge, landcover_simphuman) %>%
  mutate(cluster = "theoretical average",
         sex = tolower(sex))
cats<- rbind(categorized_females, categorized_males) %>%
  dplyr::select(shyness_category, sex, prox_urbanedge, landcover_simphuman) %>%
  mutate(cluster = "actual")
all_cats<- rbind(theory_cats, cats) %>%
  mutate(
    # Recode shyness_category as a factor with new labels
    disturbance = fct_recode(shyness_category,
                             "low" = "shy",
                             "medium" = "mixed",
                             "high" = "bold"),
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
  # geom_ribbon(data = newdata, 
  #             aes(x = prox_urbanedge, ymin = lwr, ymax = upr), 
  #             inherit.aes = FALSE,
  #             fill = "grey70", alpha = 0.5) +
  geom_line(data = newdata, 
            aes(x = prox_urbanedge, y = fit), 
            inherit.aes = FALSE,
            color = "grey40", linetype = "dashed", size = 1) +
  # geom_smooth(
  #   data = all_cats,
  #   aes(x = prox_urbanedge, y = landcover_simphuman),
  #   method = "lm",
  #   linetype = "dashed",        # Dashed line
  #   color = "black",           # Grey line
  #   size = 1
  # ) +
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
    strip.text = element_text(face = "bold"),      # bold facet titles
    panel.grid = element_blank()                   # remove background grid for clarity
  ) +
  xlim(c(-10,5))
p2
ggsave("~/puma_proj/puma_data/figs/landcover_proxurbanedge.png", plot = p2, width = 16, height = 8, dpi = 300)


# remove extraneous column and rename
real_f_pumas_combined<- real_f_pumas_combined %>% dplyr::select(-c(pumas)) %>% 
  rename( cover = cover_mean,
          landcover_simp = landcover_simphuman_mean,
          prox_urbanedge = prox_urbanedge_mean,
          slope = slope_mean,
          hd_150 = hd_150_mean)
real_m_pumas_combined<- real_m_pumas_combined %>% dplyr::select(-c(pumas)) %>% 
  rename( cover = cover_mean,
          landcover_simp = landcover_simphuman_mean,
          prox_urbanedge = prox_urbanedge_mean,
          slope = slope_mean,
          hd_150 = hd_150_mean)

# Print dimensions and names to verify
print(dim(design_matrix))
print(colnames(real_m_pumas_combined))
print(colnames(real_f_pumas_combined))

#investigate values
print(real_m_pumas_combined)
print(real_f_pumas_combined)
# Compute linear predictors and apply transformations

######## MALES #######
# Split the dataframe into a list based on 'percentile_lc'
real_m_pumas_combined <-  split(as.data.frame(real_m_pumas_combined[, -1]), real_m_pumas_combined$shyness_category)
# Convert the list of dataframes to a list of named numeric vectors
real_m_pumas_combined <- lapply(real_m_pumas_combined, function(df) {
  setNames(as.numeric(df[1, ]), colnames(df))  # Convert first row to numeric vector with names
})
# m_puma_coeffs_list
predictions <- lapply(real_m_pumas_combined, function(x) { 
  # change coef list...
  linear_predictors <- design_matrix %*% x
  plogis(linear_predictors)
})
# Create a template raster filled with NA values
predict_rasters_m <- lapply(predictions, function(predicted_values) {
  r <- setValues(stk[[1]], rep(NA, ncell(stk[[1]])))
})
# Assign predicted values to the corresponding cells in the raster
for (i in seq_along(predictions)) {
  predict_rasters_m[[i]][valid_cells] <- predictions[[i]] 
  predict_rasters_m[[i]] <-overlay(predict_rasters_m[[i]], weight_mask, fun = function(con, weight) {
    ifelse(is.na(weight), con, con * weight) # using scaled traffic values to multiple (and shrink) connectivity values
  })  
  predict_rasters_m[[i]][!is.na(city_mask)] <- predict_rasters_m[[i]][!is.na(city_mask)]/4 # Set permeability to value divided by 4 
  predict_rasters_m[[i]][r_final_masked == 1] <- 1  # Set permeability to 1 where source pop is 
}
# Check results
lapply(predictions, summary)
par(mfrow = c(1,length(predict_rasters_m)))
for (i in 1:length(real_m_pumas_combined)) {
  plot(predict_rasters_m[[i]], main = paste(names(real_m_pumas_combined)[i], "male", sep = " "), col = color_palette)
}

#### FEMALES f_puma_coeffs_list ######

# Split the dataframe into a list based on 'percentile_lc'
real_f_pumas_combined <- split(as.data.frame(real_f_pumas_combined[, -1]), real_f_pumas_combined$shyness_category)

# Convert the list of dataframes to a list of named numeric vectors
real_f_pumas_combined <- lapply(real_f_pumas_combined, function(df) {
  setNames(as.numeric(df[1, ]), colnames(df))  # Convert first row to numeric vector with names
})

# generate predictions
predictions <- lapply(real_f_pumas_combined, function(x) {
  linear_predictors <- design_matrix %*% x
  plogis(linear_predictors)
})
# Create a template raster filled with NA values
predict_rasters_f <- lapply(predictions, function(predicted_values) {
  r <- setValues(stk[[1]], rep(NA, ncell(stk[[1]])))
  r
})
# Assign predicted values to the corresponding cells in the raster
for (i in seq_along(predictions)) {
  predict_rasters_f[[i]][valid_cells] <- predictions[[i]] 
  predict_rasters_f[[i]] <-overlay(predict_rasters_f[[i]], weight_mask, fun = function(con, weight) {
    ifelse(is.na(weight), con, con * weight) # using scaled traffic values to multiple (and shrink) connectivity values
  })  
  predict_rasters_f[[i]][!is.na(city_mask)] <- predict_rasters_f[[i]][!is.na(city_mask)]/2  # Set permeability to value divided by 4 (arbitrary??)
  predict_rasters_f[[i]][r_final_masked == 1] <- 1  # Set permeability to 1 where source pop is 
}
# plot
par(mfrow = c(1,length(predict_rasters_f)))
for (i in seq_along(predict_rasters_f)) {
  plot(predict_rasters_f[[i]], main = paste(names(predict_rasters_f)[i], "female", sep = " "), col = color_palette
  )
}
# Check results
lapply(predict_rasters_f, summary)


################### Raster investigation #############
# find low values:
low_value_cells <- which(predict_rasters_f[[3]][] < 0.15, arr.ind = TRUE)
low_value_indices <- valid_cells[low_value_cells]  # Get design matrix row indices
low_value_design_matrix <- design_matrix[low_value_indices, ]  # Extract corresponding design matrix rows
# Create an empty raster with the same dimensions as the original raster
low_value_raster <- predict_rasters_f[[1]]
low_value_raster[] <- NA  # Set all values to NA
# Assign 1 to the low-value cells
low_value_raster[low_value_cells] <- 1
# Plot the raster to visualize low values
par(mfrow = c(1,1))
plot(low_value_raster, main = "Low Permeability Pixels (<0.15)")

# scaling occurs in collab

# calculate FCH
f_quart<- max(sapply(predict_rasters_f, function(x) max(values(x), na.rm = TRUE)))*0.25
m_quart<- max(sapply(predict_rasters_m, function(x) max(values(x), na.rm = TRUE)))*0.25

female_FCH<- lapply(predict_rasters_f, function(x) { ((length(which(values(x) > f_quart)) * 300^2)/(1000000))/( (length(stk[[1]]) * 300^2)/1000000) }) # get number of pixels with certain connectivity value, multiple by pixel size, then divide by conversion to kilometers squared (1e^6). Lastly, divide by total area of region
male_FCH<- lapply(predict_rasters_m, function(x) { ((length(which(values(x) > m_quart)) * 300^2)/(1000000))/( (length(stk[[1]]) * 300^2)/1000000) }) # get number of pixels with certain connectivity value, multiple by pixel size, then divide by conversion to kilometers squared (1e^6). Lastly, divide by total area of region
FCHs<- data.frame("tolerance" = rep(c("very_high", "average", "very_low"), 2), 
                  "FCH" = c(unlist(female_FCH), unlist(male_FCH)), 
                  "sex" = c(rep("Females", each = 3), rep("Males", each = 3))) %>%
  mutate(tolerance = factor(tolerance, levels = c("very_high", "average", "very_low")))
FCHs
# plot difference in FCH by tolerance and sex
par(mfrow = c(1,1))
ggplot(FCHs, aes(x = tolerance, y = FCH*100, group = sex)) + 
  geom_line(aes(color = sex), size = 2) + 
  facet_wrap(~ sex) +
  ylab("Percent of Functionally Connected Habitat in Region") +
  xlab("Degree of Anthropogenic Tolerance") +
  ggtitle("Functionally Connected Habitat by sex and tolerance") +
  theme_bw()

# Compute a difference raster
diff_raster <- predict_rasters_m[[3]] - predict_rasters_m[[2]]
par(mfrow = c(1,1))
plot(diff_raster, main = "Resistance Improvement from Bold/Neutral to Shy Male Puma")  
summary(diff_raster)
# Compute the difference raster
diff_raster <- predict_rasters_f[[3]] - predict_rasters_f[[1]]
par(mfrow = c(1,1))
plot(diff_raster, main = "Resistance Improvement from Shy to Bold Female Puma")  
summary(diff_raster)

# save rasters
setwd("~/rasters_for_collab")

writeRaster(predict_rasters_f[[1]], filename = "perm_rasters_female/permeability_bold_fem_raw.tif", format = "GTiff", 
            overwrite = TRUE)
# writeRaster(predict_rasters_f[[2]], filename = "permeability_medbold_fem_raw.tif", format = "GTiff", 
#             overwrite = TRUE)
writeRaster(predict_rasters_f[[2]], filename = "perm_rasters_female/permeability_mix_fem_raw.tif", format = "GTiff", 
            overwrite = TRUE)
# writeRaster(predict_rasters_f[[7]], filename = "permeability_medshy_fem_raw.tif", format = "GTiff", 
#             overwrite = TRUE)
writeRaster(predict_rasters_f[[3]], filename = "perm_rasters_female/permeability_shy_fem_raw.tif", format = "GTiff", 
            overwrite = TRUE)

# males
writeRaster(predict_rasters_m[[1]], filename = "perm_rasters_male/permeability_bold_mal_raw.tif", format = "GTiff", 
            overwrite = TRUE)
writeRaster(predict_rasters_m[[2]], filename = "perm_rasters_male/permeability_mix_mal_raw.tif", format = "GTiff", 
            overwrite = TRUE)
writeRaster(predict_rasters_m[[3]], filename = "perm_rasters_male/permeability_shy_mal_raw.tif", format = "GTiff", 
            overwrite = TRUE)
# generate plot
perm_stack<- stack(predict_rasters_f[[1]], predict_rasters_f[[2]], predict_rasters_f[[3]],
                   predict_rasters_m[[1]], predict_rasters_m[[2]], predict_rasters_m[[3]])
names(perm_stack)<- c("High - Female", "Medium - Female", "Low - Female",
                      "High - Male", "Medium - Male", "Low - Male")
# Create a data.frame with facet labels
facet_df <- data.frame(
  name = names(perm_stack),
  Sex = rep(c("Female", "Male"), each = 3),
  Tolerance = rep(c("High", "Medium", "Low"), times = 2)
)
# Assign names that match raster layer names
names(perm_stack) <- paste(facet_df$Sex, facet_df$Tolerance, sep = "_")

# Convert rasters to data frame with facet labels
df_list <- list()
sexes <- c("Female", "Male")
tols <- c("High", "Medium", "Low")
for (i in 1:nlayers(perm_stack)) {
  df <- as.data.frame(perm_stack[[i]], xy = TRUE)
  colnames(df)[3] <- "value"
  df$Sex <- rep(sexes, each = 3)[i]
  df$Tolerance <- rep(tols, times = 2)[i]
  df_list[[i]] <- df
}
df_all <- do.call(rbind, df_list) %>%
  mutate(Sex = factor(Sex, levels = c("Female", "Male")),
         Tolerance = factor(Tolerance, levels = c("Low", "Medium", "High"))) %>%
  rename(Permeability = value)
p<- ggplot(df_all, aes(x = x, y = y, fill = Permeability)) +
  geom_raster() +
  facet_grid(Sex ~ Tolerance, labeller = label_value, switch = "y") +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal() +
  theme(
    strip.text.x = element_text(face = "bold", size = 14),  # Column labels: Low, Medium, High
    strip.text.y.left = element_text(face = "bold", angle = 0, size = 14),  # Row labels: Male, Female
    strip.placement = "outside",  # Pushes strip labels to left
    axis.ticks = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    #axis.text = element_text(size = 14),  # Increase size of axis tick labels
    #axis.title,
    title = element_text(size = 20)
  ) +
  labs(title = "Figure 3", 
       subtitle = "Sex and Tolerance Specific Permeability",
       x = "Tolerance", y = "Sex")
print(p)
ggsave("~/output/permeabilities.png", plot = p, width = 20, height = 8, dpi = 1000)


