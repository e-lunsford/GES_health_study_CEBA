#################################################################################
# title: "07 GBD Delta y"                                                       #
# author: "Beth Lunsford"                                                       #
# last updated: "2025-02-14"                                                    #
#                                                                               #
# Last Run: 02/14/2025 and code was in working order using R 4.4.2 and          #
# RStudio 2024.12.0+467.                                                        #
#                                                                               #
#################################################################################

#################################################################################
#                                                                               #
# This code is to calculate excess cases of disease and mortality by air        #
# pollutant and health outcome.                                                 #
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
library(knitr)
library(lubridate)
library(purrr)
library(raster)
library(sf)
library(sp)
library(stringr)
library(tidycensus)
library(tidyverse)
library(units)

# Function: Delta-Y = pop * Y0 * 1 - exp(-beta * delta-C)

# Δy = Pop * (Y_0/365)*(1-exp(-β*∆C))
# Where Δy is the estimated number of excess cases of disease or mortality, 
# pop is the population exposed to air pollution 
# (i.e., Globeville population, Elyria-Swansea population), 
# Y0 is the annual baseline rate (i.e., incidence) of deaths or illness, 
# ΔC is the exposure difference. 



#################################################################################
# Population
#################################################################################


# Load ACS populations

ges_population <- as.data.frame(read_csv(file = "ges_population.csv",
                                         col_names = TRUE))



#################################################################################
# Y0
#################################################################################

# load health data
# OLD HEALTH
GBD_health5 <- as.data.frame(read_csv(file = "GBD_health5.csv",
                                      col_names = TRUE))
# CURRENT HEALTH
# Has ages 18+, 20+, and all ages
GBD_health1 <- as.data.frame(read_csv(file = "GBD_health_0116.csv",
                                      col_names = TRUE))

#################################################################################
# Beta
#################################################################################

# Load exposure functions
load(file = "exp_fctn_units.RData")

exp_fctn_units2 <- exp_fctn_units %>%
  select(air_var,
         `Exposure unit or Comparator`,
         Exposure_temporality,
         Population,
         Health_outcome,
         Reference,
         Effect_size,
         low_ci, upper_ci,
         beta, beta_lc, beta_uc) %>%
  mutate(Health_outcome = tolower(Health_outcome))

#################################################################################
# Delta C: Load air exposure estimates and exposure differences                 #
#################################################################################

# Load neighborhood level exposure differences
load(file = "nbhds_exp_dif_01_14_2025.RData")

nbhds_exp_dif2 <- nbhds_exp_dif_01_14_2025 %>%
  rename(parameter = id) %>%
  dplyr::mutate(air_var = stringr::str_extract(parameter, "[^_]+"),
                method = stringr::str_extract(parameter, "^[^_]+_[^_]+"))

nbhds_exp_dif3 <- nbhds_exp_dif2 %>%
  rename(g_delta_c = dif_exp_G,
         es_delta_c = dif_exp_ES)



#################################################################################
#                                                                               #
# join data together
# Filter out asthma, birth outcomes, and obesity (These will be calculated
# later). All other health outcomes retained are to be applied to total 
# population counts
#                                                                               #
#################################################################################

# Δy = Pop * (Y_0/100000/365)*(1-exp(-β*∆C))
GBD_health2 <- GBD_health1 %>%
  filter(!Health_outcome == "childhood asthma incidence")


exp_fctn_units3 <- exp_fctn_units2 %>%
  filter(!Health_outcome == "childhood asthma incidence",
         !Health_outcome == "low birth weight" ,
         !Health_outcome == "preterm birth" ,
         !Health_outcome == "obesity")

health_exp <- left_join(x = GBD_health2,
                        y = exp_fctn_units3,
                        relationship = "many-to-many")




health_exp2 <- left_join(health_exp,
                         y = nbhds_exp_dif3,
                         relationship = "many-to-many")

# Rename columns for delta-y calculations
health_exp3 <- health_exp2 %>%
  rename(beta_lower = beta_lc,
         beta_upper = beta_uc,
         glbvle_delta_c = g_delta_c,
         elsw_delta_c = es_delta_c)

# Add population to data frame

####################################
# 2/17/25

health_exp4 <- health_exp3 %>%
  dplyr::mutate(
    glbvle_pop = case_when(age_name == "18+ years" ~ ges_population$adult18_pop[ges_population$GES == "glbvle"],
                           age_name == "20+ years" ~ ges_population$plus20_pop[ges_population$GES == "glbvle"],
                           age_name == "All ages" ~ ges_population$total_pop[ges_population$GES == "glbvle"]),
    glbvle_pop_upper = case_when(age_name == "18+ years" ~ ges_population$adult18_pop_upper[ges_population$GES == "glbvle"],
                             age_name == "20+ years" ~ ges_population$plus20_pop_upper[ges_population$GES == "glbvle"],
                             age_name == "All ages" ~ ges_population$total_pop_upper[ges_population$GES == "glbvle"]),
    glbvle_pop_lower = case_when(age_name == "18+ years" ~ ges_population$adult18_pop_lower[ges_population$GES == "glbvle"],
                             age_name == "20+ years" ~ ges_population$plus20_pop_lower[ges_population$GES == "glbvle"],
                             age_name == "All ages" ~ ges_population$total_pop_lower[ges_population$GES == "glbvle"]),
    elsw_pop = case_when(age_name == "18+ years" ~ ges_population$adult18_pop[ges_population$GES == "elsw"],
                         age_name == "20+ years" ~ ges_population$plus20_pop[ges_population$GES == "elsw"],
                         age_name == "All ages" ~ ges_population$total_pop[ges_population$GES == "elsw"]),
    elsw_pop_upper = case_when(age_name == "18+ years" ~ ges_population$adult18_pop_upper[ges_population$GES == "elsw"],
                           age_name == "20+ years" ~ ges_population$plus20_pop_upper[ges_population$GES == "elsw"],
                           age_name == "All ages" ~ ges_population$total_pop_upper[ges_population$GES == "elsw"]),
    elsw_pop_lower = case_when(age_name == "18+ years" ~ ges_population$adult18_pop_lower[ges_population$GES == "elsw"],
                           age_name == "20+ years" ~ ges_population$plus20_pop_lower[ges_population$GES == "elsw"],
                           age_name == "All ages" ~ ges_population$total_pop_lower[ges_population$GES == "elsw"])
    )
  
GBD_df <- health_exp4 %>%
  select(-c(sum_y0,
            sum_y0_upper,
            sum_y0_lower,
            Exposure_temporality,
            Population
            ))

#################################################################################
# Calculate Delta y

# Delta-y = Pop * (Y_0/1000/5) * (1-exp(-beta * delta-C))

# Delta-y is the estimated number of excess cases of disease or mortality, 

# pop is the population exposed to air pollution 
# From CHED data cleaning, we found that the overall birth rate in CO 
# between 2017 and 2021 is ~36% -> Multiply pop by 0.36

# Y0 is the annual baseline rate (i.e., incidence) of deaths or illness, 
# For birth outcomes rate is per 1,000

# Delta-C is the exposure difference. 
#
#################################################################################


## Using Monte Carlo ##


## Using Monte Carlo ##

# Monte Carlo function to calculate delta_y and 95% CI
monte_carlo_delta_y <- function(avg_y0_lower, avg_y0_upper, 
                                beta_lower, beta_upper,
                                pop_lower, pop_upper, 
                                delta_c, n_iter = 10000) {
  
  # Generate random samples for avg_y0, beta, and pop
  avg_y0_samples <- runif(n_iter, avg_y0_lower, avg_y0_upper)
  beta_samples <- runif(n_iter, beta_lower, beta_upper)
  pop_samples <- runif(n_iter, pop_lower, pop_upper)
  
  # Compute delta_y for each sample
  delta_y_samples <- (pop_samples) * ((avg_y0_samples/100000)*5) * (1 - exp(-beta_samples * delta_c))
  
  # Compute mean and 95% CIs (2.5th and 97.5th percentiles)
  delta_y_mean <- mean(delta_y_samples)
  ci_lower <- quantile(delta_y_samples, 0.025)
  ci_upper <- quantile(delta_y_samples, 0.975)
  
  return(c(delta_y_mean, ci_lower, ci_upper))
}

# Apply Monte Carlo simulation row-wise to compute delta_y and CIs 
# for both Globeville and Elyria-Swansea
mc_GBD_df <- GBD_df %>%
  rowwise() %>%
  mutate(
    glbvle_results = list(monte_carlo_delta_y(avg_y0_lower, avg_y0_upper, 
                                              beta_lower, beta_upper,
                                              glbvle_pop_lower, glbvle_pop_upper, 
                                              glbvle_delta_c)),
    elsw_results = list(monte_carlo_delta_y(avg_y0_lower, avg_y0_upper, 
                                            beta_lower, beta_upper,
                                            elsw_pop_lower, elsw_pop_upper, 
                                            elsw_delta_c))
  ) %>%
  ungroup()

# Extract computed delta_y and confidence intervals into separate columns
mc_GBD_df2 <- as.data.frame(mc_GBD_df) %>%
  mutate(
    glbvle_delta_y = map_dbl(glbvle_results, 1),
    glbvle_delta_y_ci_lower = map_dbl(glbvle_results, 2),
    glbvle_delta_y_ci_upper = map_dbl(glbvle_results, 3),
    elsw_delta_y = map_dbl(elsw_results, 1),
    elsw_delta_y_ci_lower = map_dbl(elsw_results, 2),
    elsw_delta_y_ci_upper = map_dbl(elsw_results, 3)
  ) %>%
  select(-glbvle_results, -elsw_results,
         -nbhd_min_name)  # Remove lists for clarity


write_csv(x = mc_GBD_df2,
          file = "delta_y_GBD_MC.csv")

#################################################################################
# OLD CODE:

health_exp4 <- health_exp3 %>%
  mutate(delta_y_g = g_pop * ((avg_y0/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_low = g_pop * ((avg_y0_lower/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_high = g_pop * ((avg_y0_upper/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_es = es_pop * ((avg_y0/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_low = es_pop * ((avg_y0_lower/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_high = es_pop * ((avg_y0_upper/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_g_20 = g_pop_20plus * ((avg_y0/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_20_low = g_pop_20plus * ((avg_y0_lower/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_20_high = g_pop_20plus * ((avg_y0_upper/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_es_20 = es_pop_20plus * ((avg_y0/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_20_low = es_pop_20plus * ((avg_y0_lower/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_20_high = es_pop_20plus * ((avg_y0_upper/100000)*5) * (1-exp(-beta*es_delta_c))
  )

write_csv(x = health_exp4,
          file = "delta_y_GBD1_nbhds_01202025.csv")

write_csv(x = health_exp4,
          file = "delta_y_GBD1_nbhds_01272025.csv")


