#-------------------------------------------------------------------------------
# title: "04 Exposure Difference"  
# author: E Lunsford  
# date: 2025-08-08 
#  
# This code is to finalize all IDW data for PM2.5, O3, PM10, NO2, SO2, & CO,
# and find the differences between the lowest and highest exposure.
#
# Last Run: 08/08/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------


#################################################################################
# Load Libraries.
#################################################################################
library(broom)
library(ggmap)
library(ggplot2)
library(ggspatial)
library(gstat)
library(ggthemes)
library(knitr)
library(purrr)
library(raster)
library(sf)
library(skimr)
library(sp)
library(stars)
library(stringr)
library(tidyverse)

#################################################################################
# Load data from "01" file.
#################################################################################

# Census tract boundaries
load(file = "Outputs/metro_tract_wgs84.RData") ## Denver metropolitan census tracts
load(file = "Outputs/den_tract_wgs84.RData") ## CCoD census tracts
load(file = "Outputs/den_nbhd_wgs84.RData") ## CCoD neighborhoods

# Load combined census tract and nbhd sf
load(file = "Outputs/den_tract_nbhd.RData") # Many-to-one 
load(file = "Outputs/den_tract_nbhd2.RData") # w/o DIA
load(file = "Outputs/den_ct_nbhd.RData") # One-to-one
load(file = "Outputs/den_ct_nbhd2.RData") #w/o DIA

#################################################################################
# Load data from O3a through 03f
#################################################################################

# PM2.5
load(file = "Outputs/IDW_PM25_nbhd_centroids_08_2025.RData")
load(file = "Outputs/IDW_PM25_ct_centroids_08_2025.RData")

# O3
load(file = "Outputs/IDW_O3_nbhd_centroids_08_2025.RData")
load(file = "Outputs/IDW_O3_ct_centroids.RData")

# PM10
load(file = "Outputs/IDW_PM10_nbhd_centroids_08_2025.RData")
load(file = "Outputs/IDW_PM10_ct_centroids_08_2025.RData")

# NO2
load(file = "Outputs/IDW_NO2_nbhd_centroids_08_2025.RData")
load(file = "Outputs/IDW_NO2_ct_centroids_08_2025.RData")

# SO2
load(file = "Outputs/IDW_SO2_nbhd_centroids_08_2025.RData")
load(file = "Outputs/IDW_SO2_ct_centroids_08_2025.RData")

# CO
load(file = "Outputs/IDW_CO_nbhd_centroids_08_2025.RData")
load(file = "Outputs/IDW_CO_ct_centroids_08_2025.RData")

#-------------------------------------------------------------------------------
# For neighborhoods
#-------------------------------------------------------------------------------

#################################################################################
# 1. Combine all neighborhood centroid IDW measurements
#################################################################################

# Combine all nbhd centroid IDW measurements
IDW_nbhd <- cbind(IDW_PM25_nbhd_centroids,
                  IDW_o3_nbhd_centroids$o3_mean_pred_ppm,
                  IDW_o3_nbhd_centroids$o3_mean_ras_ppm,
                  IDW_o3_nbhd_centroids$o3_max_pred_ppm,
                  IDW_o3_nbhd_centroids$o3_max_ras_ppm,
                  IDW_pm10_nbhd_centroids$pm10_pred_ugm3,
                  IDW_pm10_nbhd_centroids$pm10_ras_ugm3,
                  IDW_no2_nbhd_centroids$no2_pred_ppb,
                  IDW_no2_nbhd_centroids$no2_ras_ppb,
                  IDW_so2_nbhd_centroids$so2_pred_ppb,
                  IDW_so2_nbhd_centroids$so2_ras_ppb,
                  IDW_co_nbhd_centroids$co_pred_ppm,
                  IDW_co_nbhd_centroids$co_ras_ppm)

# Rename columns for consistency
IDW_nbhd2 <- IDW_nbhd %>%
  dplyr::select(-GLOBALID) %>%
  dplyr::rename(o3_mean_pred_ppm = IDW_o3_nbhd_centroids.o3_mean_pred_ppm,
                o3_mean_ras_ppm = IDW_o3_nbhd_centroids.o3_mean_ras_ppm,
                o3_max_pred_ppm = IDW_o3_nbhd_centroids.o3_max_pred_ppm,
                o3_max_ras_ppm = IDW_o3_nbhd_centroids.o3_max_ras_ppm,
                pm10_pred_ugm3 = IDW_pm10_nbhd_centroids.pm10_pred_ugm3,
                pm10_ras_ugm3 = IDW_pm10_nbhd_centroids.pm10_ras_ugm3,
                no2_pred_ppb = IDW_no2_nbhd_centroids.no2_pred_ppb,
                no2_ras_ppb = IDW_no2_nbhd_centroids.no2_ras_ppb,
                so2_pred_ppb = IDW_so2_nbhd_centroids.so2_pred_ppb,
                so2_ras_ppb = IDW_so2_nbhd_centroids.so2_ras_ppb,
                co_pred_ppm = IDW_co_nbhd_centroids.co_pred_ppm,
                co_ras_ppm = IDW_co_nbhd_centroids.co_ras_ppm)

# Map for check
ggplot() +
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_sf(data = den_nbhd_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 0.75) +
  geom_point(data = IDW_nbhd2,
             mapping = aes(color = pm25_idp2_ugm3,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")

# Save data
save(file = "Outputs/IDW_nbhd.RData",
     x = IDW_nbhd2)

write_csv(x = IDW_nbhd2,
          file = "Outputs/IDW_nbhd.csv")

#################################################################################
# 2. Find the min air pollutant values for neighborhoods
#################################################################################

# Store first 4 columns (unchanged)
df_fixed <- IDW_nbhd2[, 1:4]

# Initialize empty list
filtered_rows_min <- list()

# Loop through each air pollutant to find minimum estimated level
for (i in 5:19) {
  # Set column name
  col_name <- colnames(IDW_nbhd2)[i]
  
  # Find minimum value
  col_min <- min(IDW_nbhd2[[i]], na.rm = TRUE)
  
  # Filter rows where column i has the minimum
  filtered_row_min <- IDW_nbhd2[IDW_nbhd2[[i]] == col_min, ]
  
  #Add column to show which pollutant
  filtered_row_min$pollutant <- col_name
  
  # Add column to show min value
  filtered_row_min$min_value <- col_min
  
  # Keep first 4 columns + pollutant + min value
  filtered_rows_min[[col_name]] <- filtered_row_min[, c(1:4, 
                                                        ncol(filtered_row_min)-1, 
                                                        ncol(filtered_row_min))]
  }

# Combine into a single data frame
nbhd_result_min <- do.call(rbind, filtered_rows_min)

nbhd_result_min2 <- nbhd_result_min %>%
  dplyr::select(-NBHD_ID, -lon, -lat) %>%
  sf::st_drop_geometry() 

#################################################################################
# 3. Find the max air pollutant values for neighborhoods
#################################################################################

# Store first 4 columns (unchanged)
df_fixed <- IDW_nbhd2[, 1:4]

# Initialize empty list
filtered_rows_max <- list()

# Loop through each air pollutant to find minimum estimated level
for (i in 5:19) {
  # Set column name
  col_name <- colnames(IDW_nbhd2)[i]
  
  # Find minimum value
  col_max <- max(IDW_nbhd2[[i]], na.rm = TRUE)
  
  # Filter rows where column i has the minimum
  filtered_row_max <- IDW_nbhd2[IDW_nbhd2[[i]] == col_max, ]
  
  #Add column to show which pollutant
  filtered_row_max$pollutant <- col_name
  
  # Add column to show min value
  filtered_row_max$max_value <- col_max
  
  # Keep first 4 columns + pollutant + min value
  filtered_rows_max[[col_name]] <- filtered_row_max[, c(1:4, 
                                                        ncol(filtered_row_max)-1, 
                                                        ncol(filtered_row_max))]
  }

# Combine into a single data frame
nbhd_result_max <- do.call(rbind, filtered_rows_max)  

#################################################################################
# 4. Find GES measurements from neighborhood IDW results
#################################################################################

IDW_GES_nbhd <- IDW_nbhd2 %>%
  dplyr::filter(NBHD_NAME == NBHD_NAME[str_detect(string = NBHD_NAME, 
                                                  pattern = "Globeville")] |
                  NBHD_NAME == NBHD_NAME[str_detect(string = NBHD_NAME,
                                                    pattern = "Elyria Swansea")]
                ) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(-NBHD_ID, -lat, -lon)

# Format data for merge
IDW_GES_nbhd2 <- IDW_GES_nbhd %>%
  # Pivot longer to match nbhd results
  tidyr::pivot_longer(cols = pm25_idp2_ugm3:co_ras_ppm,
                      names_to = "pollutant",
                      values_to = "value") %>%
  # Pivot wider to "collapse" data from NA's that occur from pivot longer
  pivot_wider(
    names_from = NBHD_NAME,
    values_from = value,
    names_glue = "{NBHD_NAME}_{.value}"
  ) %>%
  # Rename columns for consistency
  dplyr::rename(glob_exp = `Globeville_value`,
                elyr_exp = `Elyria Swansea_value`)

#################################################################################
# 5. Create DF that contains minimum levels & GES from nbhd data
#################################################################################

nbhd_exp_diff <- dplyr::left_join(x = nbhd_result_min2,
                                  y = IDW_GES_nbhd2,
                                  by = "pollutant") %>%
  # Calculate difference in exposure
  mutate(glob_diff = glob_exp - min_value,
         elyr_diff = elyr_exp - min_value)

# Save data
save(file = "Outputs/nbhd_exp_diff.RData",
     x = nbhd_exp_diff)

#-------------------------------------------------------------------------------
# For census tracts
#-------------------------------------------------------------------------------

#################################################################################
# 1. Combine all census tract IDW measurements
#################################################################################

# Combine all census tract IDW measurements
IDW_ct <- cbind(IDW_PM25_census_tract_centroids,
                IDW_o3_census_tract_centroids$o3_mean_pred_ppm,
                IDW_o3_census_tract_centroids$o3_mean_ras_ct_ppm,
                IDW_o3_census_tract_centroids$o3_max_pred_ppm,
                IDW_o3_census_tract_centroids$o3_max_ras_ct_ppm,
                IDW_pm10_census_tract_centroids$pm10_pred_ugm3,
                IDW_pm10_census_tract_centroids$pm10_ras_ugm3,
                IDW_no2_census_tract_centroids$no2_pred_ppb,
                IDW_no2_census_tract_centroids$no2_ras_ppb,
                IDW_so2_census_tract_centroids$so2_pred_ppb,
                IDW_so2_census_tract_centroids$so2_ras_ppb,
                IDW_co_census_tract_centroids$co_pred_ppm,
                IDW_co_census_tract_centroids$co_ras_ppm)

# Rename columns for consistency
IDW_ct2 <- IDW_ct %>%
  dplyr::rename(o3_mean_pred_ppm = IDW_o3_census_tract_centroids.o3_mean_pred_ppm,
                o3_mean_ras_ppm = IDW_o3_census_tract_centroids.o3_mean_ras_ct_ppm,
                o3_max_pred_ppm = IDW_o3_census_tract_centroids.o3_max_pred_ppm,
                o3_max_ras_ppm = IDW_o3_census_tract_centroids.o3_max_ras_ct_ppm,
                pm10_pred_ugm3 = IDW_pm10_census_tract_centroids.pm10_pred_ugm3,
                pm10_ras_ugm3 = IDW_pm10_census_tract_centroids.pm10_ras_ugm3,
                no2_pred_ppb = IDW_no2_census_tract_centroids.no2_pred_ppb,
                no2_ras_ppb = IDW_no2_census_tract_centroids.no2_ras_ppb,
                so2_pred_ppb = IDW_so2_census_tract_centroids.so2_pred_ppb,
                so2_ras_ppb = IDW_so2_census_tract_centroids.so2_ras_ppb,
                co_pred_ppm = IDW_co_census_tract_centroids.co_pred_ppm,
                co_ras_ppm = IDW_co_census_tract_centroids.co_ras_ppm)

# Map for check
ggplot() +
  geom_sf(data = den_ct_nbhd,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 0.5) +
  geom_sf(data = den_nbhd_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 0.75) +
  geom_point(data = IDW_ct2,
             mapping = aes(color = o3_mean_pred_ppm,
                           geometry = geometry),
             stat = "sf_coordinates") +
  scale_color_viridis_c(option = "C")

# Save data
save(file = "Outputs/IDW_ct.RData",
     x = IDW_ct2)

write_csv(x = IDW_ct2,
          file = "Outputs/IDW_ct.csv")

#################################################################################
# 2. Find the min air pollutant values for census tracts
#################################################################################

# Store first 15 columns (unchanged)
df_fixed <- IDW_ct2[, 1:15]

# Initialize empty list
filtered_rows_min <- list()

# Loop through each air pollutant to find minimum estimated level
for (i in 16:29) {
  # Set column name
  col_name <- colnames(IDW_ct2)[i]
  
  # Find minimum value
  col_min <- min(IDW_ct2[[i]], na.rm = TRUE)
  
  # Filter rows where column i has the minimum
  filtered_row_min <- IDW_ct2[IDW_ct2[[i]] == col_min, ]
  
  #Add column to show which pollutant
  filtered_row_min$pollutant <- col_name
  
  # Add column to show min value
  filtered_row_min$min_value <- col_min
  
  # Keep first 4 columns + pollutant + min value
  filtered_rows_min[[col_name]] <- filtered_row_min[, c(5:7, 
                                                        ncol(filtered_row_min)-1, 
                                                        ncol(filtered_row_min))]
}

# Combine into a single data frame
ct_result_min <- do.call(rbind, filtered_rows_min)

#################################################################################
# 3. Find the min & max air pollutant values for census tracts
#################################################################################

# Store first 15 columns (unchanged)
df_fixed <- IDW_ct2[, 1:15]

# Initialize empty list
filtered_rows_max <- list()

# Loop through each air pollutant to find minimum estimated level
for (i in 16:29) {
  # Set column name
  col_name <- colnames(IDW_ct2)[i]
  
  # Find minimum value
  col_max <- max(IDW_ct2[[i]], na.rm = TRUE)
  
  # Filter rows where column i has the minimum
  filtered_row_max <- IDW_ct2[IDW_ct2[[i]] == col_max, ]
  
  #Add column to show which pollutant
  filtered_row_max$pollutant <- col_name
  
  # Add column to show min value
  filtered_row_max$max_value <- col_max
  
  # Keep first 4 columns + pollutant + min value
  filtered_rows_max[[col_name]] <- filtered_row_max[, c(1:4, 
                                                        ncol(filtered_row_max)-1, 
                                                        ncol(filtered_row_max))]
}

# Combine into a single data frame
ct_result_max <- do.call(rbind, filtered_rows_max)  

#################################################################################
# 4. Find GES measurements from census tract IDW results
#################################################################################

GES_geoid <- c("08031001500", "08031003501", "08031003502")

IDW_GES_ct <- IDW_ct2 %>%
  # Remove geometry
  sf::st_drop_geometry() %>%
  # Remove columns
  dplyr::select(-AFFGEOID,-STATEFP, -COUNTYFP, -TRACTCE, -STUSPS, -NAMELSADCO,
                -STATE_NAME, -LSAD, -ALAND, -AWATER, -lon, -lat, -NAME) %>%
  # Filter to select GES
  dplyr::filter(GEOID %in% GES_geoid)

# Format for data merge
IDW_GES_ct2 <- IDW_GES_ct %>%
  # Pivot longer to match ct results
  tidyr::pivot_longer(cols = pm25_pred_ugm3:co_ras_ppm,
                      names_to = "pollutant",
                      values_to = "value") %>%
  # Pivot wider to collapse NA's
  dplyr::select(-GEOID) %>%
  pivot_wider(names_from = NAMELSAD,
              values_from = value,
              names_glue = "{NAMELSAD}_{.value}"
              ) %>%
  # Rename columns
  dplyr::rename(globe_exp = `Census Tract 15_value`,
                elyr1_exp = `Census Tract 35.01_value`,
                elyr2_exp = `Census Tract 35.02_value`)

#################################################################################
# 5. Create DF that contains minimum levels & GES from ct data
#################################################################################

ct_result_min2 <- ct_result_min %>%
  sf::st_drop_geometry()

ct_exp_diff <- dplyr::left_join(x = ct_result_min2,
                                y = IDW_GES_ct2,
                                by = "pollutant")

ct_exp_diff2 <- ct_exp_diff %>%
  mutate(elyr1_diff = elyr1_exp - min_value,
         elyr2_diff = elyr2_exp - min_value,
         glob_diff = globe_exp - min_value)

# Save data
save(file = "Outputs/ct_exp_diff.RData",
     x = ct_exp_diff2)

