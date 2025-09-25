#### generate permeability layers per sex per tolerance level ####

# load libraries
rm(list = ls())
library(sf)
library(lubridate)
library(mapview)
library(tidyverse)
library(raster)
# useful function
`%!in%`<- Negate('%in%')
color_palette <- colorRampPalette(c("blue", "darkgreen", "yellow", "orange", "red"))(100)
# Run git lfs pull to fetch large files before loading
system("git lfs pull")
# load data
load("tmp/data_for_projections.rda")
load("tmp/ind_coefs_and_weights.rda")
load("data/env.rda")
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


#### FEMALE REAL COEFFS ####
### combination ###
(real_f_pumas_combined_forplot <- coef_db %>% 
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
   group_by(tolerance_category) %>%
   na.omit() %>%
   summarize(
     across(c(landcover_simphuman, cover, 
              prox_urbanedge, slope, hd_150),
            list(mean = ~mean(.), sd = ~sd(.)),
            .names = "{.col}_{.fn}"),
     pumas = n(),
     .groups = "drop"
   ) %>%
   mutate(log_sl = weighted_df$avg[which(weighted_df$var == "log_sl" & weighted_df$sex == "f")],
          ta_scaled = 0,
          sl_scaled = 0) %>%
   rename(landcover_simp2 = landcover_simphuman_mean) %>%
   arrange(tolerance_category)) 
#### MALE REAL COEFFS ####
# Create a combined dataset
(real_m_pumas_combined_forplot <- coef_db %>% 
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
   group_by(tolerance_category) %>%
   na.omit() %>%
   summarize(
     across(c(landcover_simphuman, cover, 
              prox_urbanedge, slope, hd_150),
            list(mean = ~mean(.), sd = ~sd(.)),
            .names = "{.col}_{.fn}"),
     pumas = n(),
     .groups = "drop"
   ) %>%
   mutate(log_sl = weighted_df$avg[which(weighted_df$var == "log_sl" & weighted_df$sex == "f")],
          ta_scaled = 0,
          sl_scaled = 0) %>%
   rename(landcover_simp2 = landcover_simphuman_mean) %>%
   arrange(tolerance_category)) 

#### remove extraneous column and rename ####
real_f_pumas_combined<- real_f_pumas_combined_forplot %>% 
  dplyr::select(-c(pumas)) %>% 
  dplyr::select(-matches("_sd")) %>%
  rename( cover = cover_mean,
          prox_urbanedge = prox_urbanedge_mean,
          slope = slope_mean,
          hd_150 = hd_150_mean)
real_m_pumas_combined<- real_m_pumas_combined_forplot %>% 
  dplyr::select(-matches("_sd")) %>%
  dplyr::select(-c(pumas)) %>% 
  rename( cover = cover_mean,
          prox_urbanedge = prox_urbanedge_mean,
          slope = slope_mean,
          hd_150 = hd_150_mean)

# Print dimensions and names to verify
print(dim(design_matrix))
print(colnames(real_m_pumas_combined))
print(colnames(real_f_pumas_combined))

# Compute linear predictors and apply transformations
#### MALES ####
# Split the dataframe into a list based on 'tolerance'
real_m_pumas_combined <-  split(as.data.frame(real_m_pumas_combined[, -1]), real_m_pumas_combined$tolerance_category)
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

#### FEMALES ####

# Split the dataframe into a list based on 'tolerance'
real_f_pumas_combined <- split(as.data.frame(real_f_pumas_combined[, -1]), real_f_pumas_combined$tolerance_category)

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
# Check results
lapply(predict_rasters_f, summary)

#### generate figure 3 ####
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
ggsave("output/figure3.png", plot = p, width = 20, height = 8, dpi = 1000)


#### save rasters for inputting into EcoScape ####
setwd("rasters_for_collab")

writeRaster(predict_rasters_f[[1]], filename = "perm_rasters_female/permeability_high_fem_raw.tif", format = "GTiff", 
            overwrite = TRUE)
writeRaster(predict_rasters_f[[2]], filename = "perm_rasters_female/permeability_mix_fem_raw.tif", format = "GTiff", 
            overwrite = TRUE)
writeRaster(predict_rasters_f[[3]], filename = "perm_rasters_female/permeability_low_fem_raw.tif", format = "GTiff", 
            overwrite = TRUE)

# males
writeRaster(predict_rasters_m[[1]], filename = "perm_rasters_male/permeability_high_mal_raw.tif", format = "GTiff", 
            overwrite = TRUE)
writeRaster(predict_rasters_m[[2]], filename = "perm_rasters_male/permeability_mix_mal_raw.tif", format = "GTiff", 
            overwrite = TRUE)
writeRaster(predict_rasters_m[[3]], filename = "perm_rasters_male/permeability_low_mal_raw.tif", format = "GTiff", 
            overwrite = TRUE)

#### save objects for plotting ####
save(real_f_pumas_combined_forplot, real_m_pumas_combined_forplot, 
     file = "tmp/tolerance_specific_pumas.rda")

