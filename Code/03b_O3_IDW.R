#-------------------------------------------------------------------------------
# title: "03b O3 Inverse Distance Weighting"  
# author: E Lunsford  
# date: 2025-08-06 
#  
# This code is to spatially model O3 daily measurements over a 5 year period
# (2017 - 2021) to both the neighborhood level and census tract level.
#
# Method 1: IDW point estimate model to neighborhood centroid.
# Method 2: IDW point estimate model to census tract centroid.
# Method 3: IDW as a raster over Denver area with extraction to both neighborhood and census tract centroid
#
#
# Last Run: 08/06/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------


#################################################################################
# Load Libraries.                                                               #
#################################################################################
library(broom)
library(ggmap)
library(ggplot2)
library(ggspatial)
library(gstat)
library(ggthemes)
library(knitr)
library(lubridate)
library(purrr)
library(RAQSAPI)
library(raster)
library(sf)
library(skimr)
library(sp)
library(spatialEco)
library(stars)
library(stringr)
library(terra)
library(tidycensus)
library(tidyverse)
library(tigris)
library(units)

#################################################################################
# Load parameters of interest from "01" file. 
#################################################################################

# Years of interest
years <- c(2017:2021)

# Counties
load(file = "Outputs/denver_metro_counties_wgs84.RData") 

# Census tract boundaries
load(file = "Outputs/metro_tract_wgs84.RData") ## Denver metropolitan census tracts
load(file = "Outputs/den_tract_wgs84.RData") ## CCoD census tracts
load(file = "Outputs/den_nbhd_wgs84.RData") ## CCoD neighborhoods

# Load combined census tract and nbhd sf
load(file = "Outputs/den_tract_nbhd.RData") # Many-to-one 
load(file = "Outputs/den_tract_nbhd2.RData") # w/o DIA
load(file = "Outputs/den_ct_nbhd.RData") # One-to-one
load(file = "Outputs/den_ct_nbhd2.RData") #w/o DIA

# Set bboxes
metro_bbox <- st_bbox(metro_tract_wgs84) + c(-0.005,-0.005,0.005,0.005)
den_bbox <- st_bbox(den_tract_wgs84) + c(-0.005,-0.005,0.005,0.005)

#################################################################################
# Load EPA O3 data from "02" file.
#################################################################################

load(file = "Outputs/denver_o3_2017_2021.RData") ## Loads as ozone_aqs2

# Rename for consistency
den_o3 <- ozone_aqs2

#################################################################################
# Preview Data.                                                                 #
# Use the glimpse function from dplyr to view data.                             #
## Alternatively, use the skim function from the skimr package to see the       #
## distribution of data.                                                        #
#################################################################################

dplyr::glimpse(den_o3)

#################################################################################
# Filter Data.                                                                  #
## Filter to the pollutant standard of choice.                                  #
## Note that the PM2.5 standard was updated from 12ug/m3 to 9ug/m3 in 2024.     #
## For PM2.5, stick with the 2012 standard as the 2024 standard                 #
## was released after our data.                                                 #
## For ozone, use the Ozone 8-hour 2015.                                        #
## For PM10 there is no pollutant standard listed.                              #
## Note for NO2 filter to NO2 1-hour 2010.                                      #
## Note that for SO2, use the SO2 1-hour 2010.                                  #
## Note that for carbon monoxide, use the CO 8-hour 1971.                       #
#################################################################################

den_o3 <- den_o3 %>%
  filter(pollutant_standard == "Ozone 8-hour 2015")

#################################################################################
# Calculate summary statistics by county, site number, and year.                #
#################################################################################

sum_stats <- den_o3 %>%
  st_drop_geometry() %>%
  group_by(site_number) %>%
  summarize(count = n(),
            min = min(o3_max_ppm),
            max = max(o3_max_ppm),
            average = mean(o3_max_ppm),
            med = median(o3_max_ppm))

knitr::kable(sum_stats, digits = 2)

#################################################################################
# Plot the full area, with CCOD highlighted and monitor locations.              #
#################################################################################

ggplot() +
  # geom_sf(data = metro_tract) +
  geom_sf(data = den_ct_nbhd2, fill = "green", col = "black") +
  geom_sf(data = den_o3, aes(col = o3_max_ppm))


#################################################################################
## Grids, buffers, and points.                                                  #
## Create grid and buffer for Denver metropolitan area.                         #
#################################################################################

# Create grid for metro area.
metro_grid <- st_make_grid(metro_tract_wgs84)

# Create grid for CCoD.
den_grid <- st_make_grid(den_tract_nbhd2)

# Create 50 km buffer boundary box for monitor locations.
o3_bb <- st_bbox(den_o3 %>% 
                   st_union() %>%
                   st_buffer(dist = 50000))

# Create 50 km buffer for CCoD.
den_buffer <- st_buffer(den_grid, dist = 50000)

#################################################################################
## Calculate the centroids.                                                     #
#################################################################################

#co_centroid <- sf::st_centroid(colorado_sf_wgs84)
metro_centroid <- sf::st_centroid(metro_tract_wgs84)
den_tract_centroid <- sf::st_centroid(den_tract_wgs84)
den_nbhd_centroid <- sf::st_centroid(den_nbhd_wgs84)

den_ct_nbhd2_centroid <- sf::st_centroid(den_ct_nbhd2)

# Check with plot.
plot(st_geometry(den_ct_nbhd2)) 
plot(st_geometry(den_ct_nbhd2_centroid), col = "black", add = T)
plot(st_geometry(den_buffer), add = T, border = "green")

#################################################################################
# Start the IDW Process.                                                        #
#################################################################################

#################################################################################
# Method 1. IDW as points from neighborhood centroids.                          #
#################################################################################

# 1. Create new sf object with points of neighborhood centroids. 
## Pull the x and y coordinates using sf::st_coordinates.
sf_nbhd_centroids <- den_nbhd_centroid %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

# 2. Use IDW formula to predict concentrations (ppm) at neighborhood centroids.
o3_mean_idw_ppm_nbhd <- idw(formula = o3_mean_ppm ~1,
                            locations = den_o3,
                            newdata = sf_nbhd_centroids,
                            idp = 2.0)

o3_max_idw_ppm_nbhd <- idw(formula = o3_max_ppm ~1,
                           locations = den_o3,
                           newdata = sf_nbhd_centroids,
                           idp = 2)

predicted_values <- o3_mean_idw_ppm_nbhd$var1.pred
observed_values <- den_o3$o3_mean_ppm

predicted_observed_comparison <- data.frame(
  observed = observed_values[1:78],
  predicted = predicted_values[1:78]
)

rmse2 <- sqrt(mean((predicted_observed_comparison$observed -
                      predicted_observed_comparison$predicted)^2))
print(paste("RMSE: ", rmse2))
# RMSE: 0.0140


# 3. Spatially join neighborhood names with IDW prediction values.
IDW_o3_nbhd_level_ppm <- sf_nbhd_centroids %>%
  st_join(o3_mean_idw_ppm_nbhd) %>%
  dplyr::rename(o3_mean_pred_ppm = var1.pred) %>%
  select(-var1.var)

IDW_o3_nbhd_level_ppm <- IDW_o3_nbhd_level_ppm %>%
  st_join(o3_max_idw_ppm_nbhd) %>%
  dplyr::rename(o3_max_pred_ppm = var1.pred) %>%
  select(-var1.var)

# 4. Map it.
ggplot() +
  geom_sf(data = den_nbhd_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_o3_nbhd_level_ppm,
             mapping = aes(color = o3_mean_pred_ppm,
                           #size = 1.0,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")

#################################################################################
#  Method 2. IDW as points from census tract centroids.                         #
#################################################################################

# 1. Create new sf object with points of census tract centroids. 
## Pull the x and y coordinates using sf::st_coordinates.
sf_census_tract_centroids <- den_tract_centroid %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

# 2. Use IDW formula to predict concentrations (ppm) at census tract centroids.
o3_mean_census_tract_idw_ppm <- idw(formula = o3_mean_ppm ~ 1,
                                    locations = den_o3,
                                    newdata = sf_census_tract_centroids,
                                    idp = 2.0)

o3_max_census_tract_idw_ppm <- gstat::idw(formula = o3_max_ppm ~ 1,
                                          locations = den_o3,
                                          newdata = sf_census_tract_centroids,
                                          idp = 2.0)

# 3. Spatially join the census tract information with the IDW prediction values.
IDW_o3_census_tract_level_ppm <- sf_census_tract_centroids %>%
  st_join(o3_mean_census_tract_idw_ppm) %>%
  dplyr::rename(o3_mean_pred_ppm = var1.pred) %>%
  select(-var1.var)

IDW_o3_census_tract_level_ppm <- IDW_o3_census_tract_level_ppm %>%
  sf::st_join(o3_max_census_tract_idw_ppm) %>%
  dplyr::rename(o3_max_pred_ppm = var1.pred) %>%
  select(-var1.var)

# 4. Map it.
# 8-hr mean
ggplot() +
  geom_sf(data = den_tract_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_o3_census_tract_level_ppm,
             mapping = aes(color = o3_mean_pred_ppm,
                           #size = 1.0,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")

# 8-hr max
ggplot() +
  geom_sf(data = den_tract_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_o3_census_tract_level_ppm,
             mapping = aes(color = o3_max_pred_ppm,
                           #size = 1.0,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")



#################################################################################
#  Method 3. Raster over Denver area and extract points.                        #
#################################################################################

#################################################################################
## A. Prep data                                                                 #
## The gstat package uses sp objects instead of sf objects.                     #
## Convert sf objects to sp objects.                                            #
#################################################################################

# A1. Convert sf object to Spatial Points Data Frame.
den_o3_sp <- as_Spatial(den_o3)

# A2. Convert Metro Denver tracts to Spatial Polygons Data Frame.
metro_tract_sp <- as_Spatial(metro_tract_wgs84)

# A3. Convert CCoD tracts to Spatial Polygons Data Frame.
den_tract_sp <- as_Spatial(den_ct_nbhd)

#################################################################################
## B. Create grid, points, and pixels.                                          #
## The gstat package uses sp objects instead of sf objects.                     #
## Convert sf objects to sp objects.                                            #
# Code follows: https://github.com/adeel1997/Spatial-Interpolation-on-Air-Quality-data/blob/master/IDW_Interpolation.R
#################################################################################

# B1. Create grid for Denver metro.
grdpts <- makegrid(metro_tract_sp)

# B2. Create spatial points.
spgrd <- SpatialPoints(grdpts, proj4string = CRS(proj4string(metro_tract_sp)))

# B3. Create SpatialPixels.
spgrdWithin <- SpatialPixels(spgrd[metro_tract_sp,])

#################################################################################
## C. Interpolation.                                                            #
## The gstat package uses sp objects instead of sf objects.                     #
## Convert sf objects to sp objects.                                            #
#################################################################################

# C1. Interpolate the grid cells using a power value of 2 (idp=2.0).
## IDW = function for inverse distance weighted interpolation. 
## Idp = numeric; specify the inverse distance weighting power.
o3_mean_idw_ppm <- idw(formula = o3_mean_ppm ~ 1, 
                       locations = den_o3_sp,
                       newdata = spgrdWithin,
                       idp=2.0)

o3_max_idw_ppm <- gstat::idw(formula = o3_max_ppm ~ 1,
                             locations = den_o3_sp,
                             newdata = spgrdWithin,
                             idp = 2.0)

# C2. Set IDW as raster.
raster_o3_mean_idw <- raster(o3_mean_idw_ppm)

raster_o3_max_idw <- raster(o3_max_idw_ppm)

# C3. Get raster values.
val1_mean = as.numeric(c(
  minValue(raster_o3_mean_idw):maxValue(raster_o3_mean_idw)
))

val1_max = as.numeric(c(
  minValue(raster_o3_max_idw):maxValue(raster_o3_max_idw)
))

# C4. Clip to CCoD.
## Crop function
crop_raster_o3_mean_idw <- terra::crop(x = raster_o3_mean_idw,
                                       y = den_tract_sp)
crop_raster_o3_max_idw <- terra::crop(x = raster_o3_max_idw,
                                      y = den_tract_sp)

## Mask function
mask_raster_o3_mean_idw <- terra::mask(x = raster_o3_mean_idw,
                                       mask = den_tract_sp)
mask_raster_o3_max_idw <- terra::mask(x = raster_o3_max_idw,
                                      mask = den_tract_sp)

#################################################################################
## D. Mapping.                                                                  #
#################################################################################

# D1. Set cropped raster and masked raster as data frames for mapping.
## cropped
crop_raster_o3_mean_idw_df <- as.data.frame(crop_raster_o3_mean_idw, xy=TRUE)

## masked
mask_raster_o3_mean_idw_df <- as.data.frame(mask_raster_o3_mean_idw, xy = TRUE) %>%
  na.omit()

# D2. Map IDW raster with CCoD census tracts.
## cropped
ggplot() +
  geom_raster(data = crop_raster_o3_mean_idw_df,
              aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis_c()+
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 1) 

## masked
ggplot() +
  geom_raster(data = mask_raster_o3_mean_idw_df,
              aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis_c()+
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 1) 

#################################################################################
## E. Extraction                                                                #
# E1. Extract raster IDW to neighborhood centroids.                             #
#################################################################################

# E1 Step 1. Change sf neighborhood centroids (sf_nbhd_centroid to spatial points).
sp_nbhd_centroid <- as_Spatial(sf_nbhd_centroids)

# E1 Step 2. Extract raster values from the pm IDW raster model.
rasValue_mean_nbhd <- terra::extract(raster_o3_mean_idw, 
                                     sp_nbhd_centroid)

rasValue_max_nbhd <- terra::extract(raster_o3_max_idw,
                                    sp_nbhd_centroid)

# E1 Step 3. Combine raster values to Method 1.
IDW_o3_nbhd_centroids <- cbind(IDW_o3_nbhd_level_ppm, 
                               rasValue_mean_nbhd,
                               rasValue_max_nbhd)

IDW_o3_nbhd_centroids <- IDW_o3_nbhd_centroids %>%
  dplyr::rename(o3_mean_ras_ppm = rasValue_mean_nbhd,
                o3_max_ras_ppm = rasValue_max_nbhd)

# E1 Step 4. Map it.
ggplot() +
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_o3_nbhd_centroids,
             mapping = aes(color = o3_mean_ras_ppm,
                           #size = 1.0,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")

#################################################################################
## E. Extraction.                                                               #
# E2. Extract raster to census tract centroids.                                 #
#################################################################################

# E2 Step 1. Change sf census tract centroids (df_census_tracts) to spatial points.
sp_den_tract_centroid <- as_Spatial(sf_census_tract_centroids)

# E2 Step 2. Extract raster values from the IDW raster model.
rasValue_ct_mean <- terra::extract(raster_o3_mean_idw, 
                                   sp_den_tract_centroid)

rasValue_ct_max <- terra::extract(raster_o3_max_idw,
                                  sp_den_tract_centroid)

# E2 Step 3. Combine raster values to Method 2.
IDW_o3_census_tract_centroids <- cbind(IDW_o3_census_tract_level_ppm, 
                                       rasValue_ct_mean,
                                       rasValue_ct_max)

IDW_o3_census_tract_centroids <- IDW_o3_census_tract_centroids %>%
  dplyr::rename(o3_mean_ras_ct_ppm = rasValue_ct_mean,
                o3_max_ras_ct_ppm = rasValue_ct_max)

# E2 Step 4. Map it.
ggplot() +
  geom_sf(data = den_tract_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_o3_census_tract_centroids,
             mapping = aes(color = o3_mean_ras_ct_ppm,
                           #size = 1.0,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")


#################################################################################
# Finalize both IDW datasets for future consistency and usability.              #
#################################################################################

# Save both datasets.
save(file = "Outputs/IDW_O3_nbhd_centroids_08_2025.RData",
     x = IDW_o3_nbhd_centroids)

save(file = "Outputs/IDW_O3_ct_centroids.RData",
     x = IDW_o3_census_tract_centroids)








