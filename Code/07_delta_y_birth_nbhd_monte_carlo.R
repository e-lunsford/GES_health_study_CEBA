#################################################################################
# title: "06 delta y - Health at Birth"                                         #
# author: "Beth Lunsford"                                                       #
# last updated: "2025-02-14"                                                    #
#                                                                               #
# Last Run: 02/14/2025 and code was in working order using R 4.4.2 and          #
# RStudio 2024.12.0+467.                                                        #
#                                                                               #
#################################################################################

#################################################################################
#                                                                               #
# This code is to calculate excess cases of low birth weight and preterm        #
# birth.                                                                        #

# Function: Delta-Y = pop * Y0 * 1 - exp(-beta * delta-C)

# Δy = Pop * (Y_0/365)*(1-exp(-β*∆C))
# Where Δy is the estimated number of excess cases of disease or mortality, 
# pop is the population exposed to air pollution 
# Y0 is the annual baseline rate (i.e., incidence) of deaths or illness, 
# ΔC is the exposure difference. 
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
library(dplyr)
library(ggmap)
library(ggplot2)
library(gstat)
library(ggthemes)
library(knitr)
library(purrr)
library(raster)
library(sf)
library(sp)
library(stringr)
library(tidycensus)
library(tidyverse)
library(units)





#################################################################################
# Population
#################################################################################

# Load ACS populations

ges_population <- as.data.frame(read_csv(file = "ges_population.csv",
                                         col_names = TRUE))

#################################################################################
# Y0
#################################################################################

# load birth data
load(file = "ched_birth5.RData")

#################################################################################
# Beta
#################################################################################

# Load exposure functions
load(file = "exp_fctn_units.RData")

exp_fctn_units2 <- exp_fctn_units %>%
  dplyr::select(air_var,
         `Exposure unit or Comparator`,
         Exposure_temporality,
         Population,
         Health_outcome,
         Reference,
         Effect_size,
         low_ci, upper_ci,
         beta, beta_lc, beta_uc) %>%
  mutate(Health_outcome = tolower(Health_outcome))

exp_fctn_birth <- exp_fctn_units2 %>%
  dplyr::filter(Health_outcome == "low birth weight" |
           Health_outcome == "preterm birth")

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
# join data together                                                            #
#################################################################################

birth_exp <- left_join(x = ched_birth5,
                       y = exp_fctn_birth,
                       relationship = "many-to-many")

birth_env_exp <- left_join(x = birth_exp,
                           y = nbhds_exp_dif3,
                           relationship = "many-to-many")

# Add population
# For birth outcomes, the population at risk is females aged 15-44

birth_env_exp2 <- birth_env_exp %>%
  dplyr::mutate(
    glbvle_pop = ges_population$female_15_44_pop[ges_population$GES == "glbvle"],
    glbvle_pop_lower = ges_population$female_15_44_pop_lower[ges_population$GES == "glbvle"],
    glbvle_pop_upper = ges_population$female_15_44_pop_upper[ges_population$GES == "glbvle"],
    elsw_pop = ges_population$female_15_44_pop[ges_population$GES == "elsw"],
    elsw_pop_lower = ges_population$female_15_44_pop_lower[ges_population$GES == "elsw"],
    elsw_pop_upper = ges_population$female_15_44_pop_upper[ges_population$GES == "elsw"]
  ) %>%
  dplyr::select(-nbhd_min_name)

# Finalize data for delta-y calculations
birth_df <- birth_env_exp2 %>%
  rename(avg_y0 = avg_inc,
         avg_y0_lower = avg_inc_low,
         avg_y0_upper = avg_inc_high,
         beta_lower = beta_lc,
         beta_upper = beta_uc,
         glbvle_delta_c = g_delta_c,
         elsw_delta_c = es_delta_c)

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
  delta_y_samples <- (pop_samples*0.36) * ((avg_y0_samples/1000)*5) * (1 - exp(-beta_samples * delta_c))
  
  # Compute mean and 95% CIs (2.5th and 97.5th percentiles)
  delta_y_mean <- mean(delta_y_samples)
  ci_lower <- quantile(delta_y_samples, 0.025)
  ci_upper <- quantile(delta_y_samples, 0.975)
  
  return(c(delta_y_mean, ci_lower, ci_upper))
}

# Apply Monte Carlo simulation row-wise to compute delta_y and CIs 
# for both Globeville and Elyria-Swansea
mc_birth_df <- birth_df %>%
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
mc_birth_df2 <- mc_birth_df %>%
  mutate(
    glbvle_delta_y = map_dbl(glbvle_results, 1),
    glbvle_delta_y_ci_lower = map_dbl(glbvle_results, 2),
    glbvle_delta_y_ci_upper = map_dbl(glbvle_results, 3),
    elsw_delta_y = map_dbl(elsw_results, 1),
    elsw_delta_y_ci_lower = map_dbl(elsw_results, 2),
    elsw_delta_y_ci_upper = map_dbl(elsw_results, 3)
  ) %>%
  select(-glbvle_results, -elsw_results)  # Remove lists for clarity

write_csv(x = mc_birth_df2,
          file = "delta_y_birth_MC.csv")

######################################
## Using bootstrap ##

# Define function
calculate_delta_y <- function(pop, avg_y0, beta, delta_c) {
  (pop*0.36) * ((avg_y0/1000)*5) * (1 - exp(-beta * delta_c))
}

# Test function
test_birth_df <- birth_df %>%
  mutate(
    glbvle_delta_y = calculate_delta_y(pop = glbvle_pop,
                                       avg_y0 = avg_y0,
                                       beta = beta,
                                       delta_c = glbvle_delta_c),
    elsw_delta_y = calculate_delta_y(pop = elsw_pop,
                                     avg_y0 = avg_y0,
                                     beta = beta,
                                     delta_c = elsw_delta_c)
  )


# Define bootstrap function
bootstrap_delta_y <- function(data, n_iter = 1000) {
  # Use boot to resample data and calculate delta_y
  boot_results <- boot(data, statistic = function(data, indices) {
    sampled_data <- data[indices, ]
    
    # Calculate delta_y for glbvle and elsw using sampled data
    glbvle_delta_y <- calculate_delta_y(
      pop = sampled_data$glbvle_pop,
      avg_y0 = sampled_data$avg_y0,
      beta = sampled_data$beta,
      delta_c = sampled_data$glbvle_delta_c
    )
    
    elsw_delta_y <- calculate_delta_y(
      pop = sampled_data$elsw_pop,
      avg_y0 = sampled_data$avg_y0,
      beta = sampled_data$beta,
      delta_c = sampled_data$elsw_delta_c
    )
    
    return(c(glbvle_delta_y, elsw_delta_y))
  }, R = n_iter)
  
  return(boot_results)
}

# Apply bootstrapping to each row of asthma DF
bootstrap_results <- birth_df %>%
  rowwise() %>%
  mutate(
    bootstrap_result = list(bootstrap_delta_y(cur_data()))  # Store as a list
  )

# Extract the bootstrap results and calculate the 95% CI

##### WARNING: extraction step does not work ######
boot_CI_results <- bootstrap_results %>%
  mutate(
    # Extracting the 95% CI for both glbvle and elsw using safe extraction
    glbvle_delta_y_ci = map(bootstrap_result, ~ {
      ci <- boot.ci(.x, type = "bca")
      # Extract BCa lower and upper confidence interval values (4:5 for BCa)
      c(lower = ci$bca[4], upper = ci$bca[5])
    }),
    elsw_delta_y_ci = map(bootstrap_result, ~ {
      ci <- boot.ci(.x, type = "bca")
      # Extract BCa lower and upper confidence interval values (4:5 for BCa)
      c(lower = ci$bca[4], upper = ci$bca[5])
    })
  )




#################################################################################
# OLD CODE:

birth_env_exp3 <- birth_env_exp2 %>%
  dplyr::mutate(# Delta Y Globeville
    delta_y_g = g_female_pop * ((avg_inc/1000)*5) * (1-exp(-beta*g_delta_c)),
    # Delta Y Globeville lower CI
    delta_y_g_low = g_female_pop * ((avg_inc_low/1000)*5) * (1-exp(-beta*g_delta_c)),
    # Delta Y Globeville upper CI
    delta_y_g_high = g_female_pop * ((avg_inc_high/1000)*5) * (1-exp(-beta*g_delta_c)),
    # Delta Y ES 
    delta_y_es = es_female_pop * ((avg_inc/1000)*5) * (1-exp(-beta*es_delta_c)),
    # Delta Y ES lower CI
    delta_y_es_low = es_female_pop * ((avg_inc_low/1000)*5) * (1-exp(-beta*es_delta_c)),
    # Delta Y ES upper CI
    delta_y_es_high = es_female_pop * ((avg_inc_high/1000)*5) * (1-exp(-beta*es_delta_c)),
    delta_y_test = g_female_pop * avg_prop_test1 * 5
  )

write_csv(x = birth_env_exp3,
          file = "delta_y_birth1_nbhds_01272025.csv")

birth_env_exp4 <- birth_env_exp3 %>%
  dplyr::select(Health_outcome,air_var, method, 
         Reference, Effect_size, low_ci, upper_ci,
         delta_y_g,delta_y_g_low,delta_y_g_high,
         delta_y_es, delta_y_es_low, delta_y_es_high
  )


write_csv(x = birth_env_exp4,
          file = "delta_y_birth2_nbhds_01202025.csv")


#################################################################################
# Coffman et. al. 2024
# https://ehp.niehs.nih.gov/doi/full/10.1289/EHP12969#sec-2

# Δy = Pop * Y_0 * (1-(sum(exp(-β_i*∆C_i))

# Δy = estimated joint health impact
# pop is the population exposed to air pollution 
# (i.e., Globeville population, Elyria-Swansea population), 
# Y_0 is baseline incidence
# β_i and ΔC_i are the effect estimate and change in air pollutant given
# air pollutant (i) in a set of pollutants (n).

birth_multi <- birth_env_exp2 %>%
  filter(Health_outcome == "low birth weight")

birth_multi2 <- birth_multi %>%
  slice(c(2,5,10,14,15,16))

birth_multi3 <- birth_multi2 %>%
  dplyr::mutate(delta_y_g = g_female_pop * ((avg_inc/1000)*5) * (1-exp(-sum(beta*g_delta_c)))
                )

birth_multi4 <- birth_env_exp2 %>%
  filter(Health_outcome == "preterm birth") %>%
  slice(-c(3,4))

birth_multi5 <- birth_multi4 %>%
  dplyr::mutate(delta_y_g = g_female_pop * ((avg_inc/1000)*5) * (1-exp(-sum(beta*g_delta_c)))
  )       


birth_multi6 <- birth_env_exp2 %>%
  group_by(Health_outcome) %>%
  dplyr::mutate(delta_y_g = g_female_pop * ((avg_inc/1000)*5) * (1-exp(-sum(beta*g_delta_c)))
  )



