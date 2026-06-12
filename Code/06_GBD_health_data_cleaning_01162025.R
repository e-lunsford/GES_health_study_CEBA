#################################################################################
# title: "06 GBD Health Data Cleaning"                                          #
# author: "Beth Lunsford"                                                       #
# last updated: "2025-01-16"                                                    #
#                                                                               #
# Last Run: 01/16/2025 and code was in working order using R 4.4.2 and          #
# RStudio 2024.12.0+467.                                                        #
#                                                                               #
#################################################################################

#################################################################################
#                                                                               #
# This code is to obtain and clean health data obtained from GBD and BRFFS,     #
# and birth and death records from CHED data request.                           #
#                                                                               #
# GBD: https://vizhub.healthdata.org/gbd-results/
# BRFFS: https://data.cdc.gov/Nutrition-Physical-Activity-and-Obesity/Nutrition-Physical-Activity-and-Obesity-Behavioral/hn4x-zwk7/about_data
#                                                                               #
#                                                                               #
# Health variables include:                                                     #
# Preterm birth (CHED)                                                          #
# Low birth weight (CHED)                                                       #
# All-cause mortality (GBD)                                                     #
# All-cause cancer mortality (GBD)                                              #
# CVD Incidence and mortality (GBD)                                             #
# Hypertension (GBD)                                                            #
# Stroke Incidence and mortality (GBD)                                          #
# Obesity incidence (BRFFS)                                                     #
# T2D incidence (GBD)                                                           #
# All-cause respiratory mortality (GBD)                                         #
# COPD mortality (GBD)                                                          #
# childhood asthma incidence (GBD)                                              #
#                                                                               #
#################################################################################




#################################################################################
# Set working directory.                                                        #
#################################################################################

setwd("D:/CEBA/CEBA_Project")

#################################################################################
# Load Libraries.                                                               #
#################################################################################
library(broom)
library(ggmap)
library(ggplot2)
library(gstat)
library(ggthemes)
library(knitr)
library(lubridate)
library(purrr)
library(readr)
library(skimr)
library(stringr)
library(tidycensus)
library(tidyverse)
library(units)


#################################################################################
# GBD Data                                                                      #
#                                                                               #
# Load GBD obtained health data.                                                #
# Add column of health outcome to match exp fctns from code 05.                 #
# Clean each health data DF to match exp_fctn_units for future combining.       #
#                                                                               #
#################################################################################

# All cause mortality, Respiratory infections & TB, Chronic respiratory diseases
GBD1_health <- as.data.frame(read_csv(
  file = "health_data/IHME-GBD_2021_co_respiratory2.csv",
  col_names = TRUE))

# CVD incidence and mortality, Stroke, T2DM, hypertension, COPD, Cancer
GBD2_health <- 
  as.data.frame(read_csv(file = "health_data/IHME-GBD_2021_GBD2.csv",
                         col_names = TRUE))

#################################################################################
# Mortality (All-cause)                                                         #
#                                                                               #
# Calculate overall mortality rate                                              #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD1 to all-cause mortality rate
mortality1 <- GBD1_health %>%
  filter(measure_id == 1, # Deaths
         cause_id == 294, # All causes,
         metric_id == 3)  %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "all cause mortality")

# Create an 18+ age category
mortality_15_19 <- mortality1 %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

mortality_20plus <- mortality1 %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = val,
         upper2 = upper,
         lower2 = lower)

mortality_all_ages <- mortality1 %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = val,
         upper2 = upper,
         lower2 = lower)

mortality_18plus <- rbind(mortality_15_19,
                          mortality_20plus)

mortality2_18plus <- mortality_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "Deaths",
         Health_outcome = "all cause mortality")

mortality2_20plus <- mortality_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "Deaths",
         Health_outcome = "all cause mortality")

mortality2_all_ages <- mortality_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "Deaths",
         Health_outcome = "all cause mortality")

mortality_final <- rbind(mortality2_18plus,
                         mortality2_20plus,
                         mortality2_all_ages)


#################################################################################
# Cardiovascular diseases (CVD)                                                 #
#                                                                               #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD2 to CVD mortality
cvd_mortality <- GBD2_health %>%
  dplyr::filter(measure_id == 1, # Deaths
                cause_id == 491, # Cardiovascular diseases
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "cvd mortality")

# Create an 18+ age category
cvd_mortality_15_19 <- cvd_mortality %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

cvd_mortality_20plus <- cvd_mortality %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

cvd_mortality_all_ages <- cvd_mortality %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

cvd_mortality_18plus <- rbind(cvd_mortality_15_19,
                              cvd_mortality_20plus)

cvd_mortality2_18plus <- cvd_mortality_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "Deaths",
         Health_outcome = "cvd mortality")

cvd_mortality2_20plus <- cvd_mortality_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "Deaths",
         Health_outcome = "cvd mortality")

cvd_mortality2_all_ages <- cvd_mortality_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "Deaths",
         Health_outcome = "cvd mortality")

cvd_mortality_final <- rbind(cvd_mortality2_18plus,
                             cvd_mortality2_20plus,
                             cvd_mortality2_all_ages)

# Filter GBD2 to CVD incidence
cvd_incidence <- GBD2_health %>%
  dplyr::filter(measure_id == 6, # Incidence
                cause_id == 491, # Cardiovascular diseases
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "cvd incidence")

# Create an 18+ age category
cvd_incidence_15_19 <- cvd_incidence %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

cvd_incidence_20plus <- cvd_incidence %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

cvd_incidence_all_ages <- cvd_incidence %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

cvd_incidence_18plus <- rbind(cvd_incidence_15_19,
                              cvd_incidence_20plus)

cvd_incidence2_18plus <- cvd_incidence_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "incidence",
         Health_outcome = "cvd incidence")

cvd_incidence2_20plus <- cvd_incidence_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "incidence",
         Health_outcome = "cvd incidence")

cvd_incidence2_all_ages <- cvd_incidence_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "incidence",
         Health_outcome = "cvd incidence")

cvd_incidence_final <- rbind(cvd_incidence2_18plus,
                             cvd_incidence2_20plus,
                             cvd_incidence2_all_ages)

cvd_final <- rbind(cvd_mortality_final,
                   cvd_incidence_final)

#################################################################################
# COPD mortality                                                                #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD2 to COPD mortality
copd_mortality <- GBD2_health %>%
  dplyr::filter(measure_id == 1, # Deaths
                cause_id == 509, # Chronic obstructive pulmonary disease
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "copd mortality")

# Create an 18+ age category
copd_mortality_15_19 <- copd_mortality %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

copd_mortality_20plus <- copd_mortality %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

copd_mortality_all_ages <- copd_mortality %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

copd_mortality_18plus <- rbind(copd_mortality_15_19,
                               copd_mortality_20plus)

copd_mortality2_18plus <- copd_mortality_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "Deaths",
         Health_outcome = "copd mortality")

copd_mortality2_20plus <- copd_mortality_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "Deaths",
         Health_outcome = "copd mortality")

copd_mortality2_all_ages <- copd_mortality_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "Deaths",
         Health_outcome = "copd mortality")

copd_mortality_final <- rbind(copd_mortality2_18plus,
                              copd_mortality2_20plus,
                              copd_mortality2_all_ages)

#################################################################################
# Diabetes (T2DM) incidence                                                     #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD2 to T2DM incidence
t2dm_incidence <- GBD2_health %>%
  dplyr::filter(measure_id == 6, # Incidence
                cause_id == 976, # Diabetes mellitus type 2
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "t2dm")

# Create an 18+ age category
t2dm_incidence_15_19 <- t2dm_incidence %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

t2dm_incidence_20plus <- t2dm_incidence %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

t2dm_incidence_all_ages <- t2dm_incidence %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

t2dm_incidence_18plus <- rbind(t2dm_incidence_15_19,
                               t2dm_incidence_20plus)

t2dm_incidence2_18plus <- t2dm_incidence_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "incidence",
         Health_outcome = "t2dm")

t2dm_incidence2_20plus <- t2dm_incidence_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "incidence",
         Health_outcome = "t2dm")

t2dm_incidence2_all_ages <- t2dm_incidence_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "incidence",
         Health_outcome = "t2dm")

t2dm_incidence_final <- rbind(t2dm_incidence2_18plus,
                              t2dm_incidence2_20plus,
                              t2dm_incidence2_all_ages)

#################################################################################
# Hypertension incidence                                                        #
#                                                                               #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD2 to hypertension incidence
hypertension_incidence <- GBD2_health %>%
  dplyr::filter(measure_id == 6, # Incidence
                cause_id == 1004, # Pulmonary Arterial Hypertension
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "hypertension")

# Create an 18+ age category
hypertension_incidence_15_19 <- hypertension_incidence %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

hypertension_incidence_20plus <- hypertension_incidence %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

hypertension_incidence_all_ages <- hypertension_incidence %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

# mean(hypertension_incidence_all_ages$val)

hypertension_incidence_18plus <- rbind(hypertension_incidence_15_19,
                                       hypertension_incidence_20plus)

hypertension_incidence2_18plus <- hypertension_incidence_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "incidence",
         Health_outcome = "hypertension")

hypertension_incidence2_20plus <- hypertension_incidence_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "incidence",
         Health_outcome = "hypertension")

hypertension_incidence2_all_ages <- hypertension_incidence_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "incidence",
         Health_outcome = "hypertension")

hypertension_incidence_final <- rbind(hypertension_incidence2_18plus,
                                      hypertension_incidence2_20plus,
                                      hypertension_incidence2_all_ages)


#################################################################################
# Stroke mortality and incidence                                                #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD2 to stroke mortality
stroke_mortality <- GBD2_health %>%
  dplyr::filter(measure_id == 1, # Deaths
                cause_id == 494, # Stroke
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "stroke mortality")

# Create an 18+ age category
stroke_mortality_15_19 <- stroke_mortality %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

stroke_mortality_20plus <- stroke_mortality %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

stroke_mortality_all_ages <- stroke_mortality %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

stroke_mortality_18plus <- rbind(stroke_mortality_15_19,
                                 stroke_mortality_20plus)

stroke_mortality2_18plus <- stroke_mortality_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "Deaths",
         Health_outcome = "stroke mortality")

stroke_mortality2_20plus <- stroke_mortality_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "Deaths",
         Health_outcome = "stroke mortality")

stroke_mortality2_all_ages <- stroke_mortality_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "Deaths",
         Health_outcome = "stroke mortality")

stroke_mortality_final <- rbind(stroke_mortality2_18plus,
                                stroke_mortality2_20plus,
                                stroke_mortality2_all_ages)

# Filter GBD2 to stroke incidence
stroke_incidence <- GBD2_health %>%
  dplyr::filter(measure_id == 6, # Incidence
                cause_id == 494, # stroke
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "stroke incidence")

# Create an 18+ age category
stroke_incidence_15_19 <- stroke_incidence %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

stroke_incidence_20plus <- stroke_incidence %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

stroke_incidence_all_ages <- stroke_incidence %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

stroke_incidence_18plus <- rbind(stroke_incidence_15_19,
                                 stroke_incidence_20plus)

stroke_incidence2_18plus <- stroke_incidence_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "incidence",
         Health_outcome = "stroke incidence")

stroke_incidence2_20plus <- stroke_incidence_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "incidence",
         Health_outcome = "stroke incidence")

stroke_incidence2_all_ages <- stroke_incidence_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "incidence",
         Health_outcome = "stroke incidence")

stroke_incidence_final <- rbind(stroke_incidence2_18plus,
                                stroke_incidence2_20plus,
                                stroke_incidence2_all_ages)

stroke_final <- rbind(stroke_mortality_final,
                      stroke_incidence_final)



#################################################################################
# Cancer mortality (Total cancers)                                              #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD2 to total cancers mortality rate
cancer_mortality1 <- GBD2_health %>%
  filter(measure_id == 1, # Deaths
         cause_id == 1029, # Total cancers,
         metric_id == 3)  %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "cancer mortality")

# Create an 18+ age category
cancer_mortality_15_19 <- cancer_mortality1 %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

cancer_mortality_20plus <- cancer_mortality1 %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = val,
         upper2 = upper,
         lower2 = lower)

cancer_mortality_all_ages <- cancer_mortality1 %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = val,
         upper2 = upper,
         lower2 = lower)

cancer_mortality_18plus <- rbind(cancer_mortality_15_19,
                                 cancer_mortality_20plus)

cancer_mortality2_18plus <- cancer_mortality_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "Deaths",
         Health_outcome = "cancer mortality")

cancer_mortality2_20plus <- cancer_mortality_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "Deaths",
         Health_outcome = "cancer mortality")

cancer_mortality2_all_ages <- cancer_mortality_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "Deaths",
         Health_outcome = "cancer mortality")

cancer_mortality_final <- rbind(cancer_mortality2_18plus,
                                cancer_mortality2_20plus,
                                cancer_mortality2_all_ages)

#################################################################################
# Respiratory disease mortality (All-cause)                                     #
# Note: GBD divides respiratory into chronic respiratory diseases and           #
# Respiratory infections and tuberculosis                                       #
# Calculate a combined mortality rate                                           #
#                                                                               #
#                                                                               #
# Calculate overall rates                                                       #
# Need age range all ages and 18+                                               #
# Have 15-19 and 20+                                                            #
#                                                                               #
# 5 ages covered in 15-19; need 2 (18, 19)                                      #
# resulting in 0.4 (2/5) of 15-19 age range                                     #
# need 100% of 20+                                                              #
#################################################################################

# Filter GBD1 to respiratory mortality
respiratory_mortality <- GBD1_health %>%
  dplyr::filter(measure_id == 1, # Deaths
                !cause_id == 294, # All causes
              #  cause_id == 956,  # Respiratory infections and tuberculosis
              #  cause_id == 508, # Chronic respiratory diseases
                metric_id == 3) %>% # Rate
  # Add column for Health outcome
  mutate(Health_outcome = "respiratory mortality")

# Create an 18+ age category
respiratory_mortality_15_19 <- respiratory_mortality %>%
  dplyr::filter(age_name == "15-19 years") %>%
  mutate(val2 = (val * 0.4),
         upper2 = (upper * 0.4),
         lower2 = (lower * 0.4))

respiratory_mortality_20plus <- respiratory_mortality %>%
  dplyr::filter(age_name == "20+ years") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

respiratory_mortality_all_ages <- respiratory_mortality %>%
  dplyr::filter(age_name == "All ages") %>%
  mutate(val2 = (val),
         upper2 = (upper),
         lower2 = (lower))

respiratory_mortality_18plus <- rbind(respiratory_mortality_15_19,
                                      respiratory_mortality_20plus)

respiratory_mortality2_18plus <- respiratory_mortality_18plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "18+ years",
         measure_name = "Deaths",
         Health_outcome = "respiratory mortality")

respiratory_mortality2_20plus <- respiratory_mortality_20plus %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "20+ years",
         measure_name = "Deaths",
         Health_outcome = "respiratory mortality")

respiratory_mortality2_all_ages <- respiratory_mortality_all_ages %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "All ages",
         measure_name = "Deaths",
         Health_outcome = "respiratory mortality")

respiratory_mortality_final <- rbind(respiratory_mortality2_18plus,
                                     respiratory_mortality2_20plus,
                                     respiratory_mortality2_all_ages)

#################################################################################
# Childhood asthma incidence                                                    #
#                                                                               #
# Calculate overall incidence rate                                              #
# Need age range as 1-18 years                                                  #
# Have <20, 0-14, 10-14, 10-19, 5-14, 5-9                                       #
#                                                                               #
# 15 ages covered in 0-14; need 14 (1-14)                                       #
# resulting in 0.9 of 0-14 range                                                #
# 20 ages covered in 10-19; need 4 ages (15, 16, 17, 18)                        #
# resulting in 0.2 of 10-19 range                                               #
#################################################################################

# Load asthma GBD data
asthma <- 
  as.data.frame(readr::read_csv(
    file = "health_data/IHME-GBD_2021_co_child_ashtma.csv",
    col_names = TRUE))

# Filter to incidence rates
asthma2 <- asthma %>%
  dplyr::filter(measure_name == "Incidence",
                metric_name == "Rate") %>%
  # Add column for health outcome
  mutate(Health_outcome = "childhood asthma incidence")

asthma_1_14 <- asthma2 %>%
  dplyr::filter(age_name == "0-14 years") %>%
  mutate(val2 = (val * 0.9),
         upper2 = (upper * 0.9),
         lower2 = (lower * 0.9))
  
asthma_0_14 <- asthma2 %>%
  dplyr::filter(age_name == "0-14 years") %>%
  mutate(val2 = (val),
         upper2 = upper,
         lower2 = lower)

asthma_10_19 <- asthma2 %>%
  dplyr::filter(age_name == "10-19 years") %>%
  mutate(val2 = (val * 0.2),
         upper2 = (upper * 0.2),
         lower2 = (lower * 0.2))

asthma_1_18 <- rbind(asthma_1_14,
                     asthma_10_19)

asthma2_1_18 <- asthma_1_18 %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "1-18 years",
         measure_name = "Incidence",
         Health_outcome = "childhood asthma incidence")

asthma2_0_18 <- rbind(asthma_0_14,
                      asthma_10_19)

asthma3_0_18 <- asthma2_0_18 %>%
  reframe(sum_y0 = sum(val2),
          sum_y0_upper = sum(upper2),
          sum_y0_lower = sum(lower2),
          avg_y0 = mean(val2),
          avg_y0_upper = mean(upper2),
          avg_y0_lower = mean(lower2)) %>%
  mutate(age_name = "0-18 years",
         measure_name = "Incidence",
         Health_outcome = "childhood asthma incidence")

asthma_final <- rbind(asthma2_1_18,
                      asthma3_0_18)

#################################################################################
# Join all GBD data together                                                    #
#                                                                               #
#################################################################################

GBD_health_combo <- rbind(asthma_final,
                          cancer_mortality_final,
                          copd_mortality_final,
                          cvd_incidence_final,
                          cvd_mortality_final,
                          hypertension_incidence_final,
                          mortality_final,
                          respiratory_mortality_final,
                          stroke_incidence_final,
                          stroke_mortality_final,
                          t2dm_incidence_final)

GBD_health_0116 <- GBD_health_combo %>%
  mutate(source = "GBD")

save(file = "GBD_health_0116.RData",
     GBD_health_0116)

readr::write_csv(x = GBD_health_0116,
                 file = "GBD_health_0116.csv")

