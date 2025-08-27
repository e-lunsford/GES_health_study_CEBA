#-------------------------------------------------------------------------------
# title: "Obtain ACS data"  
# author: E Lunsford  
# date: 2025-08-04 
#  
# This code is to obtain population data from ACS.
# ACS 5 year estimates for 2017-2021 for total population and under 18.                                                             
#
# Last Run: 08/04/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
# Methods
# https://www2.census.gov/programs-surveys/acs/tech_docs/statistical_testing/2023_Instructions_for_Stat_Testing_ACS.pdf
# 
# SE for each population age group is derived from MOE with Z = 1.645.
# Standard Error = Margin of Error / 1.645
# MOE(A+B) = SQRT(MOE_A + MOE_B)
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
library(keyring)
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
#                                                                               #
# Population: Get census counts for GES for age categories.                     #
#                                                                               #
#################################################################################

# Load keyring to obtain saved ACS key
library(keyring)
key_list()

# Load 2021 ACS5 variable list
acs5_var_list <- load_variables(2021, "acs5", cache = TRUE)

# Create new list with obtained codes
acs5_vars <- c(total_male = "B01001_002",
               male_under5 = "B01001_003",
               male_5_9 = "B01001_004",
               male_10_14 = "B01001_005",
               male_15_17 = "B01001_006",
               male_18_19 = "B01001_007",
               male_20 = "B01001_008",
               male_21 = "B01001_009",
               male_22_24 = "B01001_010",
               male_25_29 = "B01001_011",
               male_30_34 = "B01001_012",
               male_35_39 = "B01001_013",
               male_40_44 = "B01001_014",
               male_45_49 = "B01001_015",
               male_50_54 = "B01001_016",
               male_55_59 = "B01001_017",
               male_60_61 = "B01001_018",
               male_62_64 = "B01001_019",
               male_65_66 = "B01001_020",
               male_67_69 = "B01001_021",
               male_70_74 = "B01001_022",
               male_75_79 = "B01001_023",
               male_80_84 = "B01001_024",
               male_85plus = "B01001_025",
               total_female = "B01001_026",
               female_under5 = "B01001_027",
               female_5_9 = "B01001_028",
               female_10_14 = "B01001_029",
               female_15_17 = "B01001_030",
               female_18_19 = "B01001_031",
               female_20 = "B01001_032",
               female_21 = "B01001_033",
               female_22_24 = "B01001_034",
               female_25_29 = "B01001_035",
               female_30_34 = "B01001_036",
               female_35_39 = "B01001_037",
               female_40_44 = "B01001_038",
               female_45_49 = "B01001_039",
               female_50_54 = "B01001_040",
               female_55_59 = "B01001_041",
               female_60_61 = "B01001_042",
               female_62_64 = "B01001_043",
               female_65_66 = "B01001_044",
               female_67_69 = "B01001_045",
               female_70_74 = "B01001_046",
               female_75_79 = "B01001_047",
               female_80_84 = "B01001_048",
               female_85plus = "B01001_049"
)

# Get Denver ACS at the census 
den_acs5_2021_nad <- get_acs(geography = "tract", #tract = census tracts
                         variables = acs5_vars, # variable list
                         geometry = TRUE,
                         key = keyring::key_get(service = "get_acs"),
                         state = "Colorado",
                         county = "Denver",
                         year = 2021, # 2017-2021
                         survey = "acs5") #default of get acs is 5

# Save ACS data to outputs folder
readr::write_csv(file = "Outputs/den_acs5_2021_nad.csv", 
                 x = den_acs5_2021_nad)


#################################################################################
# Match coordinate reference systems between ACS and air pollutant estimates    #
## EPA air data CRS is WGS84, therefore the exposure estimates are in WGS84.    #
#################################################################################

# Use st_transform to change ACS from NAD83 to WGS84
den_acs5_2021_wgs <- den_acs5_2021_nad %>%
  st_transform(crs = 4326)

#################################################################################
# Pivot wider for age summaries.                                                #
#################################################################################

den_acs_wgs <- den_acs5_2021_wgs %>%
  pivot_wider(id_cols = c(GEOID, NAME, geometry),
              names_from = variable,
              values_from = c(estimate, moe))

#################################################################################
# Calculate total count and SE (male and female) by age.                        #
#################################################################################

Z = 1.645

den_acs_pop <- den_acs_wgs %>%
  mutate(# Total population
    total_pop = estimate_total_male + estimate_total_female,
    total_pop_se = sqrt(moe_total_male + moe_total_female) / Z,
    # Under 5 population
    under5 = estimate_male_under5 + estimate_female_under5,
    under5_se = sqrt(moe_male_under5 + moe_female_under5) / Z,
    # ages 5 to 9 population
    age_5_9 = estimate_male_5_9 + estimate_female_5_9,
    age_5_9_se = sqrt(moe_male_5_9 + moe_female_5_9) / Z,
    # ages 10 to 14 population
    age_10_14 = estimate_male_10_14 + estimate_female_10_14,
    age_10_14_se = sqrt(moe_male_10_14 + moe_female_10_14) / Z,
    # ages 15 to 17 population
    age_15_17 = estimate_male_15_17 + estimate_female_15_17,
    age_15_17_se = sqrt(moe_male_15_17 + moe_female_15_17) / Z,
    # age 18 and 19 population
    age_18_19 = estimate_male_18_19 + estimate_female_18_19,
    age_18_19_se = sqrt(moe_male_18_19 + moe_female_18_19) / Z,
    # age 20 population
    age_20 = estimate_male_20 + estimate_female_20,
    age_20_se = sqrt(moe_male_20 + moe_female_20) / Z,
    # age 18 plus
    age_18plus = (estimate_male_18_19 * 0.5) + (estimate_female_18_19 * 0.5) +
      estimate_male_20 + estimate_male_21 + estimate_male_22_24 + 
      estimate_male_25_29 + estimate_male_30_34 + estimate_male_35_39 + 
      estimate_male_40_44 + estimate_male_45_49 + estimate_male_50_54 + 
      estimate_male_55_59 + estimate_male_60_61 + estimate_male_62_64 +
      estimate_male_65_66 + estimate_male_70_74 + estimate_male_75_79 + 
      estimate_male_80_84 + estimate_male_85plus +
      estimate_female_20 + estimate_female_21 + estimate_female_22_24 + 
      estimate_female_25_29 + estimate_female_30_34 + 
      estimate_female_35_39 + estimate_female_40_44 +
      estimate_female_45_49 + estimate_female_50_54 + 
      estimate_female_55_59 + estimate_female_60_61 + 
      estimate_female_62_64 + estimate_female_65_66 +
      estimate_female_67_69 + estimate_female_70_74 + 
      estimate_female_75_79 + estimate_female_80_84 + 
      estimate_female_85plus,
    age_18plus_se = sqrt((moe_male_18_19  *0.5) + (moe_female_18_19 * 0.5) +
                           moe_male_20 + moe_male_21 + moe_male_22_24 + 
                           moe_male_25_29 + moe_male_30_34 + moe_male_35_39 + 
                           moe_male_40_44 + moe_male_45_49 + moe_male_50_54 + 
                           moe_male_55_59 + moe_male_60_61 + moe_male_62_64 +
                           moe_male_65_66 + moe_male_70_74 + moe_male_75_79 + 
                           moe_male_80_84 + moe_male_85plus +
                           moe_female_20 + moe_female_21 + moe_female_22_24 + 
                           moe_female_25_29 + moe_female_30_34 + 
                           moe_female_35_39 + moe_female_40_44 +
                           moe_female_45_49 + moe_female_50_54 + 
                           moe_female_55_59 + moe_female_60_61 + 
                           moe_female_62_64 + moe_female_65_66 +
                           moe_female_67_69 + moe_female_70_74 + 
                           moe_female_75_79 + moe_female_80_84 + 
                           moe_female_85plus) / Z,
    # age 20 plus
    age_20plus = estimate_male_20 + estimate_male_21 + estimate_male_22_24 + 
      estimate_male_25_29 + estimate_male_30_34 + estimate_male_35_39 + 
      estimate_male_40_44 + estimate_male_45_49 + estimate_male_50_54 + 
      estimate_male_55_59 + estimate_male_60_61 + estimate_male_62_64 +
      estimate_male_65_66 + estimate_male_70_74 + estimate_male_75_79 + 
      estimate_male_80_84 + estimate_male_85plus +
      estimate_female_20 + estimate_female_21 + estimate_female_22_24 + 
      estimate_female_25_29 + estimate_female_30_34 + 
      estimate_female_35_39 + estimate_female_40_44 +
      estimate_female_45_49 + estimate_female_50_54 + 
      estimate_female_55_59 + estimate_female_60_61 + 
      estimate_female_62_64 + estimate_female_65_66 +
      estimate_female_67_69 + estimate_female_70_74 + 
      estimate_female_75_79 + estimate_female_80_84 + 
      estimate_female_85plus,
    age_20plus_se = sqrt(moe_male_20 + moe_male_21 + moe_male_22_24 + 
                           moe_male_25_29 + moe_male_30_34 + moe_male_35_39 + 
                           moe_male_40_44 + moe_male_45_49 + moe_male_50_54 + 
                           moe_male_55_59 + moe_male_60_61 + moe_male_62_64 +
                           moe_male_65_66 + moe_male_70_74 + moe_male_75_79 + 
                           moe_male_80_84 + moe_male_85plus +
                           moe_female_20 + moe_female_21 + moe_female_22_24 + 
                           moe_female_25_29 + moe_female_30_34 + 
                           moe_female_35_39 + moe_female_40_44 +
                           moe_female_45_49 + moe_female_50_54 + 
                           moe_female_55_59 + moe_female_60_61 + 
                           moe_female_62_64 + moe_female_65_66 +
                           moe_female_67_69 + moe_female_70_74 + 
                           moe_female_75_79 + moe_female_80_84 + 
                           moe_female_85plus) / Z,
    # female 15 to 44
    female_age_15_44 = estimate_female_15_17 + estimate_female_18_19 + 
      estimate_female_20 + estimate_female_21 +
      estimate_female_22_24 + estimate_female_25_29 +
      estimate_female_30_34 + estimate_female_35_39 +
      estimate_female_40_44,
    female_age_15_44_se = sqrt(moe_female_15_17 + moe_female_18_19 + 
                                 moe_female_20 + moe_female_21 +
                                 moe_female_22_24 + moe_female_25_29 +
                                 moe_female_30_34 + moe_female_35_39 +
                                 moe_female_40_44) / Z
  )

#################################################################################
# Assuming equal distribution across all ages, calculate age range stats for    #
# 0 to 17,                                                                      #
# under 18,                                                                     #
# under 20                                                                      #
# females 15 to 44 years old.                                                   #
#################################################################################         

den_acs_pop2 <- den_acs_pop %>%
  mutate(# Ages 0 to 17
    pop_0_17 = under5 + age_5_9 + age_10_14 + age_15_17,
    pop_0_17_se = sqrt(under5_se^2 + age_5_9_se^2 + age_10_14_se^2 +
                         age_15_17_se^2),
    # Under 18
    pop_under18 = pop_0_17 + (0.5 * age_18_19),
    pop_under18_se = sqrt((pop_0_17_se^2) + ((age_18_19_se * 0.5)^2)),
    # Under 20
    pop_under20 = pop_0_17 + age_18_19,
    pop_under20_se = sqrt((pop_0_17_se^2) + (age_18_19_se^2)),
    # female 15 to 44
    female_pop_15_44 = female_age_15_44,
    female_pop_15_44_se = sqrt((female_age_15_44_se^2)),
    # Adults 18+
    pop_18plus = age_18plus,
    pop_18plus_se = sqrt((age_18plus_se^2)),
    # Adults 20+
    pop_20plus = age_20plus,
    pop_20plus_se = sqrt((age_20_se^2))
  )


save(file = "Outputs/den_acs_pop.RData", x = den_acs_pop2)


#################################################################################         
den_acs_pop3 <- den_acs_pop2 %>%
  dplyr::select(geometry, GEOID, total_pop, total_pop_se,
                pop_under18, pop_under18_se,
                pop_18plus, pop_18plus_se,
                pop_20plus, pop_20plus_se,
                female_pop_15_44, female_pop_15_44_se)

save(file = "Outputs/den_acs_pop3.RData",
     x = den_acs_pop3)


#################################################################################
# Get population data for GES specifically

#################################################################################
     
# Load ACS populations
load(file = "den_acs_pop3.RData")

# Rename for consistency
den_acs <- den_acs_pop3

# Set GES CT ID
load(file = "Outputs/ges_ct_nbhd.RData")

# Filter total population to GES only
ges_pop <- den_acs %>%
  dplyr::filter(GEOID %in% ges_ct_nbhd$GEOID)

# Add new column with neighborhood "ID" 
ges_pop2 <- ges_pop %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(GES = 
                  if_else(condition = (stringr::str_detect(string = GEOID,
                                                           pattern = "08031001500")),
                          true = "glbvle",
                          if_else(condition = 
                                    (stringr::str_detect(string = GEOID,
                                                         pattern = "08031003501")),
                                  true = "es_1",
                                  false = "es_2")))
# Rename columns for consistency
ges_pop3 <- ges_pop2 %>%
  select(-GEOID) %>%
  rename(total_pop = total_pop,
         under18_pop = pop_under18,
         adult18_pop = pop_18plus,
         plus20_pop = pop_20plus,
         female_15_44_pop = female_pop_15_44,
         total_se = total_pop_se,
         under18_se = pop_under18_se,
         adult18_se = pop_18plus_se,
         plus20_se = pop_20plus_se,
         female_15_44_se = female_pop_15_44_se)

# Pull Globeville only
glbvle <- as.data.frame(ges_pop3) %>%
  filter(GES == "glbvle")

# Combine Elyria-Swansea's 2 CT as 1 
es <- c("es_1", "es_2")

elsw <- as.data.frame(ges_pop3) %>%
  filter(GES %in% es)

# Mutate population to get sum of ES1 and ES2
elsw_pop <- elsw %>%
  select(ends_with("pop")) 

total_row <- colSums(elsw_pop)

total_elsw_pop <- rbind(elsw_pop, total_row)

rownames(total_elsw_pop)[nrow(total_elsw_pop)] <- "total"

# Mutate to get combined SE values for ES1 and ES2
elsw_se <- elsw %>%
  select(ends_with("_se"))

new_row <- sqrt(elsw_se[1,]^2 + elsw_se[2,]^2)

total_elsw_se <- rbind(elsw_se, new_row)

rownames(total_elsw_se)[nrow(total_elsw_se)] <- "total"

# Combine population and SE for Elyria-Swansea
elsw2 <- cbind(total_elsw_pop,
               total_elsw_se)

# Mutate to create combined ES1 and ES2
elsw3 <- elsw2 %>%
  mutate(GES = "elsw")

elsw_total_row <- elsw3["total",]

# Combine Globeville and Elyria-Swansea
ges_pop_table <- rbind(elsw_total_row,
                       glbvle)

# rename rows to match neighborhoods
rownames(ges_pop_table) <- ges_pop_table$GES

# remove column of neighborhood text
ges_pop_table2 <- ges_pop_table %>%
  select(-GES)

#################################################################################
# Calculate CIs for population data

df <- as.data.frame(ges_pop_table2)

pop_CI_results <- data.frame(matrix(ncol = 0, nrow = nrow(df))) 

# Loop through population columns to calculate 95% CI
for (pop_col in grep("_pop$", names(df), value = TRUE)) {
  
  # Debugging: Print the population column being processed
  print(paste("Processing column:", pop_col))
  
  # Use str_replace to replace "_pop" with "_se" for corresponding standard error columns
  se_col <- str_replace(pop_col, "_pop$", "_se")  # Replaces "_pop" with "_se"
  
  # Debugging: Print the corresponding standard error column
  print(paste("Corresponding SE column:", se_col))
  
  # Ensure the standard error column exists in the dataframe
  if (se_col %in% names(df)) {
    # Debugging: Print the values of the population and SE columns
    print(paste("Values of", pop_col, ":", df[[pop_col]]))
    print(paste("Values of", se_col, ":", df[[se_col]]))
    
    # Calculate the 95% CI for each row
    CI_lower <- df[[pop_col]] - (1.96 * df[[se_col]])
    CI_upper <- df[[pop_col]] + (1.96 * df[[se_col]])
    
    # Store the results in CI_results dataframe
    pop_CI_results[[paste0(pop_col, "_lower")]] <- CI_lower
    pop_CI_results[[paste0(pop_col, "_upper")]] <- CI_upper
  } else {
    warning(paste("Standard error column for", pop_col, "does not exist"))
  }
}

rownames(pop_CI_results) <- rownames(ges_pop_table)

ges_acs2021 <- cbind(ges_pop_table2,
                     pop_CI_results)

ges_acs2021_v2 <- ges_acs2021 %>%
  mutate(GES = rownames(.))

#################################################################################
# Calculate the female population at risk based on average birth rate of 36%

# From CHED data cleaning, we found that the overall birth rate in CO 
# between 2017 and 2021 is ~36% -> Multiply pop by 0.36
#################################################################################

ges_acs2021_v3 <- ges_acs2021_v2 %>%
  mutate(at_risk_fem = 0.36 * female_15_44_pop,
         at_risk_fem_lower = 0.36 * female_15_44_pop_lower,
         at_risk_fem_upper = 0.36 * female_15_44_pop_upper)


save(file = "Outputs/ges_acs2021.RData",
     x = ges_acs2021_v3)


write_csv(file = "Outputs/ges_acs2021.csv",
          x = ges_acs2021_v3)

load(file = "Outputs/ges_acs2021.RData")
