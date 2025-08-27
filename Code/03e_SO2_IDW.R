#-------------------------------------------------------------------------------
# title: "03e SO2 Inverse Distance Weighting"  
# author: E Lunsford  
# date: 2025-08-07 
#  
# This code is to spatially model SO2 daily measurements over a 5 year period
# (2017 - 2021) to both the neighborhood level and census tract level.
#
# Method 1: IDW point estimate model to neighborhood centroid.
# Method 2: IDW point estimate model to census tract centroid.
# Method 3: IDW as a raster over Denver area with extraction to both neighborhood and census tract centroid
#
#
# Last Run: 08/07/2025 and code was in working order 
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
# Load EPA SO2 data from "02" file.
#################################################################################

load(file = "Outputs/denver_so2_2017_2021.RData") # Loads as so2_aqs2

# Rename for consistency
den_so2 <- so2_aqs2

#################################################################################
# Preview Data.                                                                 #
# Use the glimpse function from dplyr to view data.                             #
## Alternatively, use the skim function from the skimr package to see the       #
## distribution of data.                                                        #
#################################################################################

dplyr::glimpse(den_so2)

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

den_so2 <- den_so2 %>%
  filter(pollutant_standard == "SO2 1-hour 2010")

#################################################################################
# Calculate summary statistics by county, site number, and year.                #
#################################################################################

sum_stats <- den_so2 %>%
  st_drop_geometry() %>%
  group_by(site_number) %>%
  summarize(count = n(),
            min = min(so2_ppb),
            max = max(so2_ppb),
            average = mean(so2_ppb),
            med = median(so2_ppb))

knitr::kable(sum_stats, digits = 2)

#################################################################################
# Plot the full area, with CCOD highlighted and monitor locations.              #
#################################################################################

ggplot() +
  # geom_sf(data = metro_tract) +
  geom_sf(data = den_nbhd_wgs84, fill = "green", col = "black") +
  geom_sf(data = den_so2, aes(col = so2_ppb))

#################################################################################
## Grids, buffers, and points.                                                  #
## Create grid and buffer for Denver metropolitan area.                         #
#################################################################################

# Create grid for metro area.
metro_grid <- st_make_grid(metro_tract_wgs84)

# Create grid for CCoD.
den_grid <- st_make_grid(den_tract_nbhd2)

# Create 50 km buffer boundary box for monitor locations.
so2_bb <- st_bbox(den_so2 %>% 
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
plot(st_geometry(den_tract_nbhd2))
plot(st_geometry(den_ct_nbhd2)) 
plot(st_geometry(den_ct_nbhd2_centroid), col = "black", add = T)
#plot(st_geometry(den_buffer), add = T, border = "green")

plot(st_geometry(den_nbhd_wgs84))
plot(st_geometry(den_nbhd_centroid), col = "black", add = T)

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

# 2. Use IDW formula to predict concentrations (ppb) at neighborhood centroids.
so2_nbhd_idw_ppb <- idw(formula = so2_ppb ~ 1,
                        locations = den_so2,
                        newdata = sf_nbhd_centroids,
                        idp = 2.0)

# 3. Spatially join the neighborhood names with the IDW prediction values.
IDW_so2_nbhd_level_ppb <- sf_nbhd_centroids %>%
  st_join(so2_nbhd_idw_ppb) %>%
  dplyr::rename(so2_pred_ppb = var1.pred)

# 4. Map it.
ggplot() +
  geom_sf(data = den_nbhd_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_so2_nbhd_level_ppb,
             mapping = aes(color = so2_pred_ppb,
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

# 2. Use IDW to predict concentrations (ppb) at census tract centroids.
so2_census_tract_idw_ppb <- idw(formula = so2_ppb ~ 1,
                                locations = den_so2,
                                newdata = sf_census_tract_centroids,
                                idp = 2.0)

# 3. Spatially join the census tract information with the IDW prediction values.
IDW_so2_census_tract_level_ppb <- sf_census_tract_centroids %>%
  st_join(so2_census_tract_idw_ppb) %>%
  dplyr::rename(so2_pred_ppb = var1.pred)

# 4. Map it.
ggplot() +
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_so2_census_tract_level_ppb,
             mapping = aes(color = so2_pred_ppb,
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
den_so2_sp <- as_Spatial(den_so2)

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
so2_idw <- idw(formula = so2_ppb ~ 1, 
               locations = den_so2_sp,
               newdata = spgrdWithin,
               idp = 2.0)

# C2. Set IDW as raster.
raster_so2_idw <- raster(so2_idw)

# C3. Get raster values.
val1 = as.numeric(c(minValue(raster_so2_idw):maxValue(raster_so2_idw)))

# C4. Clip to CCoD.
## Crop function
crop_raster_so2_idw <- crop(raster_so2_idw, den_tract_sp)

## Mask function
mask_raster_so2_idw <- mask(raster_so2_idw, den_tract_sp)

#################################################################################
## D. Mapping.                                                                  #
#################################################################################

# D1. Set cropped raster and masked raster as data frames for mapping.
## cropped
crop_raster_so2_idw_df <- as.data.frame(crop_raster_so2_idw, xy=TRUE)

## masked
mask_raster_so2_idw_df <- as.data.frame(mask_raster_so2_idw, xy = TRUE) %>%
  na.omit()

# D2. Map IDW raster with CCoD census tracts.
## cropped
ggplot() +
  geom_raster(data = crop_raster_so2_idw_df,
              aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis_c()+
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 1) 

## masked
ggplot() +
  geom_raster(data = mask_raster_so2_idw_df,
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

# E1 Step 1. Change sf neighborhood centroids to spatial points.
sp_nbhd_centroid <- as_Spatial(sf_nbhd_centroids)

# E1 Step 2. Extract raster values from the pm IDW raster model.
rasValue_nbhd <- terra::extract(raster_so2_idw, sp_nbhd_centroid)

# E1 Step 3. Combine raster values to Method 1.
IDW_so2_nbhd_centroids <- cbind(IDW_so2_nbhd_level_ppb, rasValue_nbhd)

IDW_so2_nbhd_centroids <- IDW_so2_nbhd_centroids %>%
  dplyr::rename(so2_ras_ppb = rasValue_nbhd)

# E1 Step 4. Map it.
ggplot() +
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_so2_nbhd_centroids,
             mapping = aes(color = so2_ras_ppb,
                           #size = 1.0,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")

#################################################################################
## E. Extraction.                                                               #
# E2. Extract raster to census tract centroids.                                 #
#################################################################################

# E2 Step 1. Change sf census tract centroids to spatial points.
sp_den_tract_centroid <- as_Spatial(sf_census_tract_centroids)

# E2 Step 2. Extract raster values from the IDW raster model.
rasValue_census_tract <- terra::extract(raster_so2_idw, sp_den_tract_centroid)

# E2 Step 3. Combine raster values to Method 2.
IDW_so2_census_tract_centroids <- cbind(IDW_so2_census_tract_level_ppb, 
                                        rasValue_census_tract)

IDW_so2_census_tract_centroids <- IDW_so2_census_tract_centroids %>%
  dplyr::rename(so2_ras_ppb = rasValue_census_tract)

# E2 Step 4. Map it.
ggplot() +
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_point(data = IDW_so2_census_tract_centroids,
             mapping = aes(color = so2_ras_ppb,
                           #size = 1.0,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")

#################################################################################
# Finalize both IDW data sets for future consistency and usability.             #
#################################################################################

# Save both datasets.
save(file = "Outputs/IDW_SO2_nbhd_centroids_08_2025.RData",
     x = IDW_so2_nbhd_centroids)

save(file = "Outputs/IDW_SO2_ct_centroids_08_2025.RData", 
     x = IDW_so2_census_tract_centroids)


