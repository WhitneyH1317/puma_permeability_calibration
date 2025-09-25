#### prepare spatial data for projections ####

rm(list =ls())
# load libraries
library(amt)
library(tidyverse)
library(sf)
library(lubridate)
library(mapview)
library(raster)
library(lwgeom)

# useful function
`%!in%`<- Negate('%in%')
# Run git lfs pull to fetch large files before loading
system("git lfs pull")
# load data
load("data/env.rda")
load("data/implicit_data.rda")
# generate mean movement metrics
means<- implicit_used_control %>%
  summarize(mean_step_dist = mean(step_dist),
            sd_step_dist = sd(step_dist),
            mean_turning_angle = mean(turning_angle),
            sd_turning_angle = mean(turning_angle),
            mean_log_sl = mean(log_sl),
            sd_log_sl = sd(log_sl)) %>%
  pivot_longer(cols = everything(),
               names_to = c(".value", "cov"),
               names_pattern = "(mean|sd)_(.*)") %>%
  rename(mean = mean, sd = sd)

#### generate landcover simp raster ####

# read in new projections
lulc_reprojected <- projectRaster(lulc_current, crs = crs(stk[[4]]), method = "ngb")
lulc_current <- resample(lulc_reprojected, stk[[4]], method = "ngb") # Use "bilinear" or "ngb" (nearest neighbor) depending on your data
# reclassify according to key
reclass_matrix <- as.matrix(iclus_key[, c("code", "landcover_simp")])
lulc_current_reclassed<- reclassify(lulc_current, reclass_matrix) 

#### create roads mask weighted with traffic values ####
# Crop roads shapefile
ext<- terra::ext(stk[[1]])
raster_bbox <- terra::as.polygons(ext, crs = as.character(crs(stk[[1]])))
roadc<- st_crop(road_shp, raster_bbox)
trafficc<- st_crop(traffic, raster_bbox)
# snap points to roads
traffic_snapped<- st_snap(trafficc, roadc, 10)
# Compute road segmentation based on point locations
# Split roads at traffic points
road_segments <- roadc %>%
  st_cast("LINESTRING") %>%
  lwgeom::st_split(., traffic_snapped) %>%
  st_collection_extract("LINESTRING")
# Join traffic points to segments (nearest snapped point)
road_segments <- st_join(road_segments, traffic_snapped, join = st_nearest_feature)
# Replace geometry with centroids for distance-based imputation
seg_centroids <- st_centroid(road_segments)
# Find NA and non-NA rows
na_rows <- which(is.na(seg_centroids$BACK_PEAK_MADT))
non_na_rows <- which(!is.na(seg_centroids$BACK_PEAK_MADT))

# Fill in missing traffic values using nearest from same ROUTE
for (i in na_rows) {
  this_route <- seg_centroids$ROUTE[i]
  candidates <- non_na_rows[seg_centroids$ROUTE[non_na_rows] == this_route]
  
  if (length(candidates) > 0) {
    dists <- st_distance(seg_centroids[i, ], seg_centroids[candidates, ])
    nearest <- candidates[which.min(dists)]
    seg_centroids$BACK_PEAK_MADT[i] <- seg_centroids$BACK_PEAK_MADT[nearest]
  }
}

# Assign filled values back to road_segments
road_segments$traffic_value <- seg_centroids$BACK_PEAK_MADT
road_segments_buff<- road_segments %>%
  st_buffer(., 300) %>%
  rename(traffic = BACK_PEAK_MADT)
# Create a blank raster with the extent and resolution of the connectivity raster
blank_raster <- setValues(stk[[1]], rep(NA, ncell(stk[[1]])))
# Rasterize the road traffic sf object (assuming a "traffic" field in sf)
road_mask_traffic_raster <- rasterize(road_segments_buff, blank_raster, field = "traffic")
# Define a scaling function (adjust as needed)
scaling_function <- function(traffic) {
    1 / (1 + log1p(traffic))
  }  # Example: higher traffic reduces connectivity more
 # and apply
weight_mask <- calc(road_mask_traffic_raster, fun = scaling_function)

#### make city mask ####
# Create a logical mask of non-NA cells
valid_cells <- !is.na(getValues(stk[[1]])) & !is.na(getValues(lulc_current_reclassed))
# Create a binary mask of the city
city_mask<- stk[[6]]
city_mask[city_mask > 0] <- NA  # 

#### Create a binary mask of the raster where values == 1 ###
r<- terra::rast(r) # mcp of all puma collar data (same locations in implicit location data)
r_binary <- terra::classify(r, cbind(-Inf, 0.9999, NA))  # Set everything not 1 to NA
r_binary[!is.na(r_binary)] <- 1                   # Set 1s explicitly
sc_protect# publicly available protected areas in Santa Cruz county
# Mask using shapefile (only keeps values within polygons)
r_masked <- mask(r_binary, sc_protect, updatevalue = 0, inverse = TRUE)
# Create a binary mask of where roads exist
road_presence <- !is.na(weight_mask)
# Set r_final values to NA where road pixels are present
r_final_masked <- mask(raster(r_masked), road_presence, maskvalue = 1)
# now write it for use in googlecolab
writeRaster(r_final_masked, "output/puma_pop_seed_raster.tif",
            overwrite = T)

#### generate new data to be used in projections ####
new_data<- data.frame(hd_150 = getValues(stk[[4]]),
                      prox_urbanedge = getValues(stk[[6]]),
                      slope = getValues(stk[[1]]), 
                      cover = getValues(stk[[2]]), 
                      landcover_simp = as.factor(getValues(lulc_current_reclassed)),
                      step_dist = means[means$cov == "step_dist",]$mean,
                      turning_angle = means[means$cov == "turning_angle",]$mean,
                      log_sl =means[means$cov == "log_sl",]$mean,
                      end_gpsID = 1) %>%
  mutate(hd_150 = ifelse(is.na(hd_150), 0, hd_150),
  ) %>%
  # scaling data according to same averages we used before
  mutate(hd_150 = (hd_150 - unname(scaling_averages[1,1]))/unname(scaling_sds[1,1]), 
         slope = (slope - unname(scaling_averages[4,1]))/unname(scaling_sds[4,1]),
         cover = (cover - unname(scaling_averages[5,1]))/unname(scaling_sds[5,1]),
         prox_urbanedge = (prox_urbanedge - unname(scaling_averages[2,1]))/unname(scaling_sds[2,1]),
         sl_scaled = 0,
         ta_scaled = 0) %>%
  filter(!is.na(slope) & !is.na(landcover_simp))
# Create the design matrix using model.matrix
design_matrix <- model.matrix(~ hd_150 + #hd_150_sq + 
                                prox_urbanedge  +
                                slope + cover  + 
                                landcover_simp +
                                sl_scaled + ta_scaled +
                                log_sl,
                              data = new_data)
# Ensure the design matrix does not include an intercept if the model does not have one
design_matrix <- design_matrix[, -1]

# save it all
save(design_matrix, 
     lulc_current_reclassed,
     weight_mask, 
     city_mask,
     r_final_masked,
     file = "~/tmp/data_for_projections.rda")


