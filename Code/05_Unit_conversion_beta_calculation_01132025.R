#################################################################################
# title: "05 Unit Conversions and beta calculations"                            #
# author: "Beth Lunsford"                                                       #
# last updated: "2025-01-13"                                                    #
#                                                                               #
# Last Run: 01/13/2025 and code was in working order using R 4.4.2 and          #
# RStudio 2024.12.0+467.                                                        #
#                                                                               #
#################################################################################


#################################################################################
# Neighborhood Level Analysis.                                                  #
#                                                                               #
# Burden calculation follows the formula:                                       #
#                                                                               #
# Delta-y = pop * (Y_0/365)*(1-exp(-beta * delta-c))                            #
# Where:                                                                        #
# Delta-y is the estimated number of excess cases of disease or mortality,      #
# pop is the population exposed to air pollution                                #
# (i.e., Globeville population, Elyria-Swansea population),                     #
# Y_0 is the annual baseline rate (i.e., incidence) of deaths or illness,       #
# Delta-c is the exposure difference.                                           #
#                                                                               #
#                                                                               #
# This code is to convert literature obtained effect estimates unit of change   #
# from ug/m3 to ppb and ppm units dependent on the air pollutant of             #
# interest.                                                                     #
#                                                                               #
# Calculation of Beta                                                           #
# 1. Load exposure response functions                                           #
# 2. Calculate per X ppb or ppm based on EPA and unit conversions               #
# 2. Use formula:                                                               #
#       Beta = ln(effect size) / delta-X                                        #
#       delta-X refers to the per unit change.                                  #
#                                                                               #
#################################################################################



#################################################################################
# Set working directory to external hard drive CEBA folder.                     #
#################################################################################

setwd("D:/CEBA/CEBA_Project")

#################################################################################
# Load Libraries.                                                               #
#################################################################################
library(broom)
library(ggmap)
library(ggplot2)
library(ggspatial)
library(gstat)
library(ggthemes)
library(keyring)
library(knitr)
library(leaflet)
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
# Load exposure response functions.                                             #
#################################################################################

# Load Exposure response functions
exp_fctn <- as.data.frame(read_csv("exp_resp_fnct.csv",
                                   col_names = TRUE))

# Create new column for air pollutant ID for future merge of datasets.
exp_fctn2 <- exp_fctn %>%
  mutate(air_var = tolower(air_var)) %>%
  mutate(unit = str_split_fixed(string = `Exposure unit or Comparator`,
                                pattern = " ",
                                n = 3)) %>%
  rename(low_ci = `Low CI`,
         upper_ci = `Upper CI`)

# Rename split "unit" columns for delta-x and unit text.
exp_fctn3 <- exp_fctn2 %>%
  mutate(delta_x = as.numeric(unit[,2]),
         unit_text = unit[,3]) %>%
  select(-unit)


#################################################################################
# Load unit conversion csv file.                                                #
#################################################################################

# Load unit conversions
unit_conv <- as.data.frame(read_csv("Unit_conversions.csv",
                                    col_names = TRUE)) %>%
  na.omit

# Add columns for molecular weight and molecular volume.
unit_conv2 <- unit_conv %>%
  mutate(mol_weight = `mole_weight_g/mol`,
         mol_volume = 22.41 * (298/273) * (1013.25/840)) %>%
  rename(air_var = air_pollutant) %>%
  select(air_var,
         mol_weight,
         mol_volume)

# Join exposure functions and unit conversions into 1 dataframe.
exp_fctn4 <- left_join(x = exp_fctn3,
                       y = unit_conv2)

#################################################################################
# Convert effect size unit of change to EPA units and calculate betas by air    #
# pollutant.                                                                    #
#                                                                               #
# Where:                                                                        #
# β = ln(effect size)/ΔX.                                                       #
# ΔX refers to the per unit change.                                             #
#################################################################################



#################################################################################
# PM2.5                                                                         #
# PM2.5 does not need any unit conversion.                                      #
#################################################################################

pm25_exp_fctn <- exp_fctn4 %>%
  dplyr::filter(air_var == "pm25") %>%
  mutate(delta_x_ugm3 = delta_x,
         delta_x_mgm3 = delta_x_ugm3 / 1000,
         delta_x_ppb = NA,
         delta_x_ppm = NA,
         beta = log(Effect_size)/delta_x_ugm3,
         beta_lc = log(low_ci)/delta_x_ugm3,
         beta_uc = log(upper_ci)/delta_x_ugm3)

#################################################################################
# PM10                                                                          #
# PM10 does not need any unit conversion.                                       #
#################################################################################

pm10_exp_fctn <- exp_fctn4 %>%
  dplyr::filter(air_var == "pm10") %>%
  mutate(delta_x_ugm3 = delta_x,
         delta_x_mgm3 = delta_x_ugm3 / 1000,
         delta_x_ppb = NA,
         delta_x_ppm = NA,
         beta = log(Effect_size)/delta_x_ugm3,
         beta_lc = log(low_ci)/delta_x_ugm3,
         beta_uc = log(upper_ci)/delta_x_ugm3)

#################################################################################
# O3                                                                            #
# O3 needs to be converted to ppm units.                                        #
# Effect estimates obtained from literature are in ppb, and ug/m3.              #
#################################################################################

# Filter out O3
o3_exp_fctn <- exp_fctn4 %>%
  dplyr::filter(air_var == "o3") 

# Filter effect estimates of ppm units
o3_exp_fctn_ppb <- o3_exp_fctn[grep("ppb", o3_exp_fctn$unit_text), ]

# Calculate beta
o3_exp_fctn2_ppb <- o3_exp_fctn_ppb %>%
  mutate(delta_x_ppb = delta_x,
         delta_x_ppm = 1000 * delta_x_ppb,
         delta_x_ugm3 = NA,
         delta_x_mgm3 = NA,
         beta = log(Effect_size)/delta_x_ppm,
         beta_lc = log(low_ci)/delta_x_ppm,
         beta_uc = log(upper_ci)/delta_x_ppm)

# Filter effect estimates of ug/m3 units
o3_exp_fctn_ugm3 <- o3_exp_fctn[grep("ug/m3", o3_exp_fctn$unit_text), ]

# Calculate beta
o3_exp_fctn2_ugm3 <- o3_exp_fctn_ugm3 %>%
  mutate(delta_x_ugm3 = delta_x,
         delta_x_mgm3 = delta_x_ugm3 / 1000,
         delta_x_ppb = delta_x_ugm3 * (mol_volume/mol_weight),
         delta_x_ppm = delta_x_ppb / 1000,
         beta = log(Effect_size) / delta_x_ppm,
         beta_lc = log(low_ci) / delta_x_ppm,
         beta_uc = log(upper_ci) / delta_x_ppm)

# Combine both O3 df into 1
o3_exp_fctn_units <- rbind(o3_exp_fctn2_ppb,
                           o3_exp_fctn2_ugm3)

#################################################################################
# CO                                                                            #
# CO needs to be converted to ppm units.                                        #
# Effect estimates obtained from literature are in ppb, mg/m3 and ug/m3.        #
#################################################################################

# Filter out CO
co_exp_fctn <- exp_fctn4 %>%
  filter(air_var == "co")

# Filter effect estimates of ppm units
co_exp_fctn_ppb <- co_exp_fctn[grep("ppb", co_exp_fctn$unit_text), ]

# Calculate beta
co_exp_fctn2_ppb <- co_exp_fctn_ppb %>%
  mutate(delta_x_ppb = delta_x,
         delta_x_ppm = 1000 * delta_x_ppb,
         delta_x_ugm3 = NA,
         delta_x_mgm3 = NA,
         beta = log(Effect_size) / delta_x_ppm,
         beta_lc = log(low_ci) / delta_x_ppm,
         beta_uc = log(upper_ci) / delta_x_ppm)

# Filter effect estimates of ug/m3 units.
co_exp_fctn_ugm3 <- co_exp_fctn[grep("ug/m3", co_exp_fctn$unit_text), ]

# Calculate beta
co_exp_fctn2_ugm3 <- co_exp_fctn_ugm3 %>%
  mutate(delta_x_ugm3 = delta_x,
         delta_x_mgm3 = delta_x_ugm3 / 1000,
         delta_x_ppb = delta_x_ugm3 * (mol_volume/mol_weight),
         delta_x_ppm = delta_x_ppb / 1000,
         beta = log(Effect_size) / delta_x_ppm,
         beta_lc = log(low_ci) / delta_x_ppm,
         beta_uc = log(upper_ci) / delta_x_ppm)

# Filter effect estimates of mg/m3 units.
co_exp_fctn_mgm3 <- co_exp_fctn[grep("mg/m3", co_exp_fctn$unit_text), ]

# Calculate beta
co_exp_fctn2_mgm3 <- co_exp_fctn_mgm3 %>%
  mutate(delta_x_mgm3 = delta_x,
         delta_x_ugm3 = delta_x_mgm3 * 1000,
         delta_x_ppm = delta_x_mgm3 * (mol_volume/mol_weight),
         delta_x_ppb = delta_x_ppm * 1000,
         beta = log(Effect_size) / delta_x_ppm,
         beta_lc = log(low_ci) / delta_x_ppm,
         beta_uc = log(upper_ci) / delta_x_ppm)

# Combine all CO beta df into 1.
co_exp_fctn_units <- rbind(co_exp_fctn2_mgm3,
                           co_exp_fctn2_ppb,
                           co_exp_fctn2_ugm3)

#################################################################################
# SO2                                                                           #
# SO2 needs to be converted to ppb units.                                       #
# Literature effect estimates obtained are in ppb and ug/m3 units.              #
#################################################################################

# Filter to SO2 effect estimates.
so2_exp_fctn <- exp_fctn4 %>%
  filter(air_var == "so2")

# Filter effect estimates of ppb units.
so2_exp_fctn_ppb <- so2_exp_fctn[grep("ppb", so2_exp_fctn$unit_text), ]

# Calculate beta
so2_exp_fctn2_ppb <- so2_exp_fctn_ppb %>%
  mutate(delta_x_ppb = delta_x,
         delta_x_ppm = 1000 * delta_x_ppb,
         delta_x_ugm3 = NA,
         delta_x_mgm3 = NA,
         beta = log(Effect_size) / delta_x_ppb,
         beta_lc = log(low_ci) / delta_x_ppb,
         beta_uc = log(upper_ci) / delta_x_ppb)

# Filter effect estimates of ug/m3 units,
so2_exp_fctn_ugm3 <- so2_exp_fctn[grep("ug/m3", so2_exp_fctn$unit_text), ]

# Calculate beta.
so2_exp_fctn2_ugm3 <- so2_exp_fctn_ugm3 %>%
  mutate(delta_x_ugm3 = delta_x,
         delta_x_mgm3 = delta_x_ugm3 / 1000,
         delta_x_ppb = delta_x_ugm3 * (mol_volume/mol_weight),
         delta_x_ppm = delta_x_ppb / 1000,
         beta = log(Effect_size) / delta_x_ppb,
         beta_lc = log(low_ci) / delta_x_ppb,
         beta_uc = log(upper_ci) / delta_x_ppb)

# Combine SO2 dfs into 1.
so2_exp_fctn_units <- rbind(so2_exp_fctn2_ppb,
                            so2_exp_fctn2_ugm3)

#################################################################################
# NO2                                                                           #
# NO2 effect estimates need to be converted to ppb units.                       #
# Literature obtained effect estimates are in ppb and ug/m3 units.              #
#################################################################################

# Filter to NO2 effect estimates. 
no2_exp_fctn <- exp_fctn4 %>%
  filter(air_var == "no2")

# Filter to ppb units.
no2_exp_fctn_ppb <- no2_exp_fctn[grep("ppb", no2_exp_fctn$unit_text), ]

# Calculate beta.
no2_exp_fctn2_ppb <- no2_exp_fctn_ppb %>%
  mutate(delta_x_ppb = delta_x,
         delta_x_ppm = 1000 * delta_x_ppb,
         delta_x_ugm3 = NA,
         delta_x_mgm3 = NA,
         beta = log(Effect_size) / delta_x_ppb,
         beta_lc = log(low_ci) / delta_x_ppb,
         beta_uc = log(upper_ci) / delta_x_ppb)

# Filter to ug/m3 units.
no2_exp_fctn_ugm3 <- no2_exp_fctn[grep("ug/m3", no2_exp_fctn$unit_text), ]

# Calculate beta.
no2_exp_fctn2_ugm3 <- no2_exp_fctn_ugm3 %>%
  mutate(delta_x_ugm3 = delta_x,
         delta_x_mgm3 = delta_x_ugm3 / 1000,
         delta_x_ppb = delta_x_ugm3 * (mol_volume/mol_weight),
         delta_x_ppm = delta_x_ppb / 1000,
         beta = log(Effect_size) / delta_x_ppb,
         beta_lc = log(low_ci) / delta_x_ppb,
         beta_uc = log(upper_ci) / delta_x_ppb)

# Combine NO2 dfs into 1.
no2_exp_fctn_units <- rbind(no2_exp_fctn2_ppb,
                            no2_exp_fctn2_ugm3)

#################################################################################
# Finalize literature obtained effect estimates for all pollutants into 1 DF.   #
#################################################################################

# Combine into 1 DF.
exp_fctn_units <- rbind(pm25_exp_fctn,
                        pm10_exp_fctn,
                        o3_exp_fctn_units,
                        co_exp_fctn_units,
                        so2_exp_fctn_units,
                        no2_exp_fctn_units)

# Export file as .RData file and .csv file.

save(file = "exp_fctn_units.RData", exp_fctn_units)

write_csv(file = "exp_fctn_units.csv",
          x = exp_fctn_units)


load(file = "exp_fctn_units.RData")






