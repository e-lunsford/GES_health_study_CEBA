#################################################################################
# title: "06 delta y - Asthma"                                                  #
# author: "Beth Lunsford"                                                       #
# last updated: "2025-02-13"                                                    #
#                                                                               #
# Last Run: 02/13/2025 and code was in working order using R 4.4.2 and          #
# RStudio 2024.12.0+467.                                                        #
#                                                                               #
#################################################################################

#################################################################################
#                                                                               #
# This code is to calculate excess cases of childhood asthma.                   #
#                                                                               #
#                                                                               #
# Function: Delta-Y = pop * Y0 * 1 - exp(-beta * delta-C)
#                                                                               #
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
library(boot)
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
# Population
#################################################################################

ges_population <- as.data.frame(read_csv(file = "ges_population.csv",
                                         col_names = TRUE))


#################################################################################
# Y0
#################################################################################

# load health data
GBD_health1 <- as.data.frame(read_csv(file = "GBD_health_0116.csv",
                                      col_names = TRUE))

GBD_health_asthma <- GBD_health1 %>%
  dplyr::filter(Health_outcome == "childhood asthma incidence")

GBD_health_asthma2 <- GBD_health_asthma %>%
  filter(age_name == "1-18 years")

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

exp_fctn_asthma <- exp_fctn_units2 %>%
  filter(Health_outcome == "childhood asthma incidence")

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

asthma_exp <- left_join(GBD_health_asthma2,
                        exp_fctn_asthma,
                        relationship = "many-to-many")


asthma_env_exp <- left_join(x = asthma_exp,
                            y = nbhds_exp_dif3,
                            relationship = "many-to-many")

# Add population
asthma_env_exp2 <- asthma_env_exp %>%
  dplyr::mutate(glbvle_pop = ges_population$under18_pop[ges_population$GES == "glbvle"],
                glbvle_pop_lower = ges_population$under18_pop_lower[ges_population$GES == "glbvle"],
                glbvle_pop_upper = ges_population$under18_pop_upper[ges_population$GES == "glbvle"],
                elsw_pop = ges_population$under18_pop[ges_population$GES == "elsw"],
                elsw_pop_lower = ges_population$under18_pop_lower[ges_population$GES == "elsw"],
                elsw_pop_upper = ges_population$under18_pop_upper[ges_population$GES == "elsw"]
  ) %>%
  dplyr::select(-nbhd_min_name)

# Finalize data frame for delta-y calculations
asthma_df <- asthma_env_exp2 %>%
  select(-c(sum_y0, sum_y0_upper, sum_y0_lower,
            Exposure_temporality)
        ) %>%
  rename(beta_lower = beta_lc,
         beta_upper = beta_uc,
         glbvle_delta_c = g_delta_c,
         elsw_delta_c = es_delta_c)

#################################################################################

# Calculate Delta-y
# Delta-y = Pop * (Y_0/100000/5)*(1-exp(-beta * delta-C))

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
  delta_y_samples <- pop_samples * ((avg_y0_samples/100000)*5) * (1 - exp(-beta_samples * delta_c))
  
  # Compute mean and 95% CIs (2.5th and 97.5th percentiles)
  delta_y_mean <- mean(delta_y_samples)
  ci_lower <- quantile(delta_y_samples, 0.025)
  ci_upper <- quantile(delta_y_samples, 0.975)
  
  return(c(delta_y_mean, ci_lower, ci_upper))
}

# Apply Monte Carlo simulation row-wise to compute delta_y and CIs 
# for both Globeville and Elyria-Swansea
mc_asthma_df <- asthma_df %>%
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
mc_asthma_df2 <- mc_asthma_df %>%
  mutate(
    glbvle_delta_y = map_dbl(glbvle_results, 1),
    glbvle_delta_y_ci_lower = map_dbl(glbvle_results, 2),
    glbvle_delta_y_ci_upper = map_dbl(glbvle_results, 3),
    elsw_delta_y = map_dbl(elsw_results, 1),
    elsw_delta_y_ci_lower = map_dbl(elsw_results, 2),
    elsw_delta_y_ci_upper = map_dbl(elsw_results, 3)
  ) %>%
  select(-glbvle_results, -elsw_results)  # Remove lists for clarity

write_csv(x = mc_asthma_df2,
          file = "delta_y_asthma_MC.csv")

######################################
## Using bootstrap ##

# Define function
calculate_delta_y <- function(pop, avg_y0, beta, delta_c) {
  pop * ((avg_y0/100000)*5) * (1 - exp(-beta * delta_c))
}

# Test function
test_asthma_df <- asthma_df %>%
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
bootstrap_results <- asthma_df %>%
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


################################################################################
## OLD CODE ##

asthma_env_exp3 <- asthma_env_exp2 %>%
  mutate(delta_y_g = g_pop_under18 * ((avg_y0/100000) * 5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_low = g_pop_under18 * ((avg_y0_lower/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_high = g_pop_under18 * ((avg_y0_upper/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_es = es_pop_under18 * ((avg_y0/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_low = es_pop_under18 * ((avg_y0_lower/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_high = es_pop_under18 * ((avg_y0_upper/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_g_0_17 = g_pop_0_17 * ((avg_y0/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_0_17_low = g_pop_0_17 * ((avg_y0_lower/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_g_0_17_high = g_pop_0_17 * ((avg_y0_upper/100000)*5) * (1-exp(-beta*g_delta_c)),
         delta_y_es_0_17 = es_pop_0_17 * ((avg_y0/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_0_17_low = es_pop_0_17 * ((avg_y0_lower/100000)*5) * (1-exp(-beta*es_delta_c)),
         delta_y_es_0_17_high = es_pop_0_17 * ((avg_y0_upper/100000)*5) * (1-exp(-beta*es_delta_c))
  )

write_csv(x = asthma_env_exp3,
          file = "delta_y_asthma1_nbhd_01202025.csv")

write_csv(x = asthma_env_exp3,
          file = "delta_y_asthma1_nbhd_01272025.csv")

asthma <- as.data.frame(read_csv("delta_y_asthma1_nbhd_01272025.csv",
                                 col_names = TRUE))





