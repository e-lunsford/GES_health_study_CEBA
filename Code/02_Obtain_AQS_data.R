#-------------------------------------------------------------------------------
# title: "Obtain EPA AQS data"  
# author: E Lunsford  
# date: 2025-08-04 
#  
# This code is to obtain AQS data for PM2.5, PM10, O3, NO2, SO2, and CO from
# EPA AQS API platform
#
# Last Run: 08/04/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#
#-------------------------------------------------------------------------------

#################################################################################
# Load Libraries                                                                #
#################################################################################
library(broom)
library(ggmap)
library(ggplot2)
library(ggspatial)
library(ggthemes)
library(knitr)
library(lubridate)
library(purrr)
library(RAQSAPI)
library(raster)
library(sf)
library(skimr)
library(stringr)
library(tidycensus)
library(tidyverse)
library(tigris)
library(units)

################################################################################
# Load key-ring information to set AQS API credentials.                        #
################################################################################

library(keyring)
key_list()

server <- keyring::key_list()[2, "service"]
datamartAPI_user <- keyring::key_list()[2, "username"]

aqs_credentials(username = datamartAPI_user, 
                key = key_get(service = server,
                              username = datamartAPI_user))

#################################################################################
# Load parameters 01_setting_area file.                                         #
#################################################################################

# Years of interest
years <- c(2017:2021)

# Denver Metro counties
load(file = "Outputs/denver_metro_counties_wgs84.RData")

# City and County of Denver census tracts & neighborhoods
load(file = "Outputs/den_ct_nbhd2.RData")

# Denver metro bbox
metro_bbox <- st_bbox(metro_tract) + c(-0.005,-0.005,0.005,0.005)

# City and County of Denver bbox
den_bbox <- st_bbox(den_tract) + c(-0.005,-0.005,0.005,0.005)



#################################################################################
# title: "Obtaining EPA AQS for PM2.5"                                          #
#################################################################################

# Set air parameter of interest.
## PM2.5 - 88101


#################################################################################
# Get the Data.                                                                 #
## Obtain PM2.5 monitor data across Colorado for 2017-2021.                     #
## EPA Data Frame                                                               #
## Use the aqs daily summary function under the RAQSAPI package.                #
#################################################################################

pm25_aqs <- aqs_dailysummary_by_state(parameter="88101",#PM2.5
                                      bdate=as.Date("2017-01-01"),
                                      edate=as.Date("2021-12-31"),
                                      stateFIPS = "08")

# Save data frame as a csv file. 
write_csv(file = "Data/Colorado_PM25_2017_2021.csv", x = pm25_aqs)

pm25_aqs <- read_csv(file = "Data/Colorado_PM25_2017_2021.csv")

#################################################################################
# Data Manipulation.                                                            #  
## Manipulate data into usable format.                                          #
## 1.Filter data to Denver metropolitan counties (Adams, Arapahoe, Broomfield,  #
### Denver, Douglas, Jefferson).                                                #
## 2. Convert data frame into an sf object.                                     #
## 3. Set the CRS to WGS84.                                                     #
## 4. Mutate date into "date" format, and pull year using the lubridate package.#
## 5. Remove some of the columns for easier future use.                         #
#################################################################################

pm25_aqs2 <- pm25_aqs %>%
  filter(#pollutant_standard == "PM25 Annual 2012",
    county_code %in% denver_metro_counties_wgs84$COUNTYFP)%>% 
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
  mutate(sample_date = as.Date(date_local),
         #pollutant = "PM2.5",
         year = lubridate::year(date_local)) %>%
  rename("pm25_ugm3"="arithmetic_mean") %>%
  dplyr::select(-poc, -datum, -parameter,
                -sample_duration_code, -units_of_measure,
                -event_type, -observation_count, -observation_percent,
                -validity_indicator, 
                -first_max_value, -first_max_hour, 
                -method_code, -method)
# geometry, state_code, county_code, site_number, Date, 
# Pollutant, Concentration)

# Save as a csv file. 
write_csv(file="Data/Denver_PM25_2017_2021.csv",x = pm25_aqs2)

# Save as .RData file.
save(file = "Data/denver_pm25_2017_2021.RData", x = pm25_aqs2)

load(file = "Data/denver_pm25_2017_2021.RData")

#################################################################################
# Use the glimpse function from dplyr to view data.                             #
#################################################################################

dplyr::glimpse(pm25_aqs2)

#################################################################################
# Calculate summary statistics by county, site number.                          #
#################################################################################

pm25_sum_stats <- pm25_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, year, local_site_name, site_address) %>%
  summarize(count = n(),
            min = min(pm25_ugm3),
            max = max(pm25_ugm3),
            average = mean(pm25_ugm3),
            med = median(pm25_ugm3))

knitr::kable(pm25_sum_stats, digits = 2)

write_csv(file = "Data/pm25_locations.csv", x = pm25_sum_stats)

pm25_sum_stats2 <- pm25_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, local_site_name, site_address) %>%
  summarize(count = n(),
            min = min(pm25_ugm3),
            max = max(pm25_ugm3),
            average = mean(pm25_ugm3),
            med = median(pm25_ugm3))

knitr::kable(pm25_sum_stats2, digits = 3)

write_csv(file = "Outputs/pm25_sum_stats.csv", x = pm25_sum_stats2)

#################################################################################
## Check monitor locations for PM2.5.                                           #
#################################################################################

ggplot() +
  geom_sf(data = den_ct_nbhd2,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = pm25_aqs2)


#################################################################################
# title: "Obtaining EPA AQS for O3"                                             #
#################################################################################

# Set air parameters of interest.
## Ozone - 44201


#################################################################################
# Get the Data.                                                                 #
## Obtain Ozone monitor data across Colorado for 2017-2021.                     #
## EPA Data Frame                                                               #
## Use the aqs daily summary function under the RAQSAPI package.                #
#################################################################################

ozone_aqs <- aqs_dailysummary_by_state(parameter="44201",#O3
                                       bdate=as.Date("2017-01-01"),
                                       edate=as.Date("2021-12-31"),
                                       stateFIPS = "08")

# Save data frame as a csv file. 
write_csv(file = "Data/Colorado_O3_2017_2021.csv", x = ozone_aqs)

ozone_aqs <- read_csv(file = "Data/Colorado_O3_2017_2021.csv")

#################################################################################
# Data Manipulation.                                                            #  
## Manipulate data into usable format.                                          #
## 1.Filter data to Denver metropolitan counties (Adams, Araphaoe, Broomfield,  #
### Denver, Douglas, Jefferson).                                                #
## 2. Convert data frame into an sf object.                                     #
## 3. Set the CRS to WGS84.                                                     #
## 4. Mutate date into "date" format, and pull year using the lubridate package.#
## 5. Remove some of the columns for easier future use.                         #
#################################################################################

ozone_aqs2 <- ozone_aqs %>%
  filter(#pollutant_standard == "",
    county_code %in% denver_metro_counties_wgs84$COUNTYFP)%>% 
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
  mutate(sample_date = as.Date(date_local),
         #pollutant = "ozone",
         year = lubridate::year(date_local)) %>%
  rename("o3_mean_ppm"="arithmetic_mean",
         "o3_max_ppm" = "first_max_value") %>%
  dplyr::select(-parameter_code, -poc, -datum, -parameter,
                -sample_duration_code, -sample_duration, 
                -units_of_measure, -event_type,
                -observation_count, -observation_percent,
                -validity_indicator, -first_max_hour,
                -method_code, -method)
# geometry, state_code, county_code, site_number, 
# Date, Pollutant, Concentration)

# Save as a csv file. 
write_csv(file="Data/Denver_O3_2017_2021.csv",
          x = ozone_aqs2)

# Save as .RData file.
save(file = "Data/denver_o3_2017_2021.RData", x = ozone_aqs2)

load(file = "Data/denver_o3_2017_2021.RData")

#################################################################################
# Use the glimpse function from dplyr to view data.                             #
#################################################################################

dplyr::glimpse(ozone_aqs2)

################################################################################
# Calculate summary statistics by county, site number, and year.               #
################################################################################

ozone_sum_stats <- ozone_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, year, site_address, local_site_name) %>%
  summarize(count = n(),
            min_av = min(o3_mean_ppm),
            max_av = max(o3_mean_ppm),
            average_av = mean(o3_mean_ppm),
            med_av = median(o3_mean_ppm),
            min_mx = min(o3_max_ppm),
            max_mx = max(o3_max_ppm),
            average_mx = mean(o3_max_ppm),
            med_mx = median(o3_max_ppm))

knitr::kable(ozone_sum_stats, digits = 2)

write_csv(file = "Data/o3_locations.csv", x = ozone_sum_stats)

ozone_sum_stats2 <- ozone_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, site_address, local_site_name) %>%
  summarize(count = n(),
            min_av = min(o3_mean_ppm),
            max_av = max(o3_mean_ppm),
            average_av = mean(o3_mean_ppm),
            med_av = median(o3_mean_ppm),
            min_mx = min(o3_max_ppm),
            max_mx = max(o3_max_ppm),
            average_mx = mean(o3_max_ppm),
            med_mx = median(o3_max_ppm))

knitr::kable(ozone_sum_stats2, digits = 2)

write_csv(file = "Outputs/o3_sum_stats.csv", x = ozone_sum_stats2)

#################################################################################
## Check monitor locations.                                                     #
#################################################################################

ggplot() +
  geom_sf(data = den_ct_nbhd2,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = ozone_aqs2)



#################################################################################
# title: "Obtaining EPA AQS for PM10"                                           #
#################################################################################

# Set air parameters of interest.
## PM10 - 81102


#################################################################################
# Get the Data.                                                                 #
## Obtain PM10 monitor data across Colorado for 2017-2021.                      #
## EPA Data Frame                                                               #
## Use the aqs daily summary function under the RAQSAPI package.                #
#################################################################################

pm10_aqs <- aqs_dailysummary_by_state(parameter="81102",#PM10
                                      bdate=as.Date("2017-01-01"),
                                      edate=as.Date("2021-12-31"),
                                      stateFIPS = "08")
# Save dataframe as a csv file. 
write_csv(file = "Data/Colorado_PM10_2017_2021.csv", x = pm10_aqs)

pm10_aqs <- read_csv(file = "Data/Colorado_PM10_2017_2021.csv")

#################################################################################
# Data Manipulation.                                                            #  
## Manipulate data into usable format.                                          #
## 1.Filter data to Denver metropolitan counties (Adams, Araphaoe, Broomfield,  #
### Denver, Douglas, Jefferson).                                                #
## 2. Convert data frame into an sf object.                                     #
## 3. Set the CRS to WGS84.                                                     #
## 4. Mutate date into "date" format, and pull year using the lubridate package.#
## 5. Remove some of the columns for easier future use.                         #
#################################################################################

pm10_aqs2 <- pm10_aqs %>%
  filter(#pollutant_standard == "PM25 Annual 2012",
    county_code %in% denver_metro_counties_wgs84$COUNTYFP)%>% 
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
  mutate(sample_date = as.Date(date_local),
         #pollutant = "PM2.5",
         year = lubridate::year(date_local)) %>%
  rename("pm10_ugm3"="arithmetic_mean") %>%
  dplyr::select(-parameter_code, -poc,-datum, -parameter, 
                -sample_duration_code, -sample_duration,
                -units_of_measure, -event_type,
                -observation_count, -observation_percent,
                -validity_indicator,
                -first_max_value, -first_max_hour, 
                -method_code, -method)
# geometry, state_code, county_code, site_number, Date, Pollutant, Concentration)

# Save as a csv file. 
write_csv(file = "Data/Denver_PM10_2017_2021.csv",
          x = pm10_aqs2)

# Save as .RData file.
save(file = "Data/denver_pm10_2017_2021.RData", x = pm10_aqs2)

load(file="Data/denver_pm10_2017_2021.RData")

#################################################################################
# Use the glimpse function from dplyr to view data.                             #
#################################################################################

dplyr::glimpse(pm10_aqs2)

#################################################################################
# Calculate summary statistics by county, site number, and year.                #
#################################################################################

pm10_sum_stats <- pm10_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, year, site_address, local_site_name) %>%
  summarize(count = n(),
            min = min(pm10_ugm3),
            max = max(pm10_ugm3),
            average = mean(pm10_ugm3),
            med = median(pm10_ugm3))

knitr::kable(pm10_sum_stats, digits = 2)

write_csv(file = "Data/pm10_locations.csv", x = pm10_sum_stats)

pm10_sum_stats2 <- pm10_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, site_address, local_site_name) %>%
  summarize(count = n(),
            min = min(pm10_ugm3),
            max = max(pm10_ugm3),
            average = mean(pm10_ugm3),
            med = median(pm10_ugm3))

knitr::kable(pm10_sum_stats2, digits = 2)

write_csv(file = "Outputs/pm10_sum_stats.csv", x = pm10_sum_stats2)

#################################################################################
## Check monitor locations.                                                     #
#################################################################################

ggplot() +
  geom_sf(data = den_ct_nbhd2,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = pm10_aqs2)


#################################################################################
# title: "Obtaining EPA AQS for NO2"                                            #
#################################################################################

# Set air parameters of interest.
## NO2 - 42602

#################################################################################
# Get the Data.                                                                 #
## Obtain NO2 monitor data across Colorado for 2017-2021.                       #
## EPA Data Frame                                                               #
## Use the aqs daily summary function under the RAQSAPI package.                #
#################################################################################

no2_aqs <- aqs_dailysummary_by_state(parameter="42602",#NO2
                                     bdate=as.Date("2017-01-01"),
                                     edate=as.Date("2021-12-31"),
                                     stateFIPS = "08")
# Save dataframe as a csv file. 
write_csv(file = "Data/Colorado_NO2_2017_2021.csv", x = no2_aqs)

no2_aqs <- read_csv(file = "Data/Colorado_NO2_2017_2021.csv")

#################################################################################
# Data Manipulation.                                                            #  
## Manipulate data into usable format.                                          #
## 1.Filter data to Denver metropolitan counties (Adams, Araphaoe, Broomfield,  #
### Denver, Douglas, Jefferson).                                                #
## 2. Convert data frame into an sf object.                                     #
## 3. Set the CRS to WGS84.                                                     #
## 4. Mutate date into "date" format, and pull year using the lubridate package.#
## 5. Remove some of the columns for easier future use.                         #
#################################################################################

no2_aqs2 <- no2_aqs %>%
  filter(#pollutant_standard == "NO2 1-hour 2010",
    county_code %in% denver_metro_counties_wgs84$COUNTYFP)%>% 
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
  mutate(sample_date = as.Date(date_local),
         #pollutant = "NO2",
         year = lubridate::year(date_local)) %>%
  rename("no2_ppb"="arithmetic_mean") %>%
  dplyr::select(-parameter_code, -poc,-datum, -parameter, 
                -sample_duration_code, -sample_duration,
                -units_of_measure, -event_type,
                -observation_count, -observation_percent,
                -validity_indicator,
                -first_max_value, -first_max_hour, 
                -method_code, -method)
# geometry, state_code, county_code, site_number, 
# Date, Pollutant, Concentration)

# Save as a csv file. 
write_csv(file = "Data/Denver_NO2_2017_2021.csv",
          x = no2_aqs2)

# Save as .RData file.
save(file = "Data/denver_no2_2017_2021.RData", x = no2_aqs2)

load(file = "Data/denver_no2_2017_2021.RData")

#################################################################################
# Use the glimpse function from dplyr to view data.                             #
#################################################################################

dplyr::glimpse(no2_aqs2)

#################################################################################
# Calculate summary statistics by county, site number, and year.                #
#################################################################################

no2_sum_stats <- no2_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, year, site_address, local_site_name) %>%
  summarize(count = n(),
            min = min(no2_ppb),
            max = max(no2_ppb),
            average = mean(no2_ppb),
            med = median(no2_ppb))

knitr::kable(no2_sum_stats, digits = 2)

write_csv(file = "Data/no2_locations.csv", x = no2_sum_stats)

no2_sum_stats2 <- no2_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, site_address, local_site_name) %>%
  summarize(count = n(),
            min = min(no2_ppb),
            max = max(no2_ppb),
            average = mean(no2_ppb),
            med = median(no2_ppb))

knitr::kable(no2_sum_stats2, digits = 2)

write_csv(file = "Outputs/no2_sum_stats.csv", x = no2_sum_stats2)

#################################################################################
## Check monitor locations.                                                     #
#################################################################################

ggplot() +
  geom_sf(data = den_ct_nbhd2,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = no2_aqs2)



#################################################################################
# title: "Obtaining EPA AQS for SO2"                                            #
#################################################################################

# Set air parameters of interest.
## SO2 - 42401


#################################################################################
# Get the Data.                                                                 #
## Obtain SO2 monitor data across Colorado for 2017-2021.                       #
## EPA Data Frame                                                               #
## Use the aqs daily summary function under the RAQSAPI package.                #
#################################################################################

so2_aqs <- aqs_dailysummary_by_state(parameter="42401",#SO2
                                     bdate=as.Date("2017-01-01"),
                                     edate=as.Date("2021-12-31"),
                                     stateFIPS = "08")
# Save dataframe as a csv file. 
write_csv(file = "Data/Colorado_SO2_2017_2021.csv", x = so2_aqs)

so2_aqs <- read_csv(file = "Data/Colorado_SO2_2017_2021.csv")

#################################################################################
# Data Manipulation.                                                            #  
## Manipulate data into usable format.                                          #
## 1.Filter data to Denver metropolitan counties (Adams, Araphaoe, Broomfield,  #
### Denver, Douglas, Jefferson).                                                #
## 2. Convert data frame into an sf object.                                     #
## 3. Set the CRS to WGS84.                                                     #
## 4. Mutate date into "date" format, and pull year using the lubridate package.#
## 5. Remove some of the columns for easier future use.                         #
#################################################################################

so2_aqs2 <- so2_aqs %>%
  filter(#pollutant_standard == "SO2 1-hour 2010",
    county_code %in% denver_metro_counties_wgs84$COUNTYFP)%>% 
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
  mutate(sample_date = as.Date(date_local),
         #pollutant = "SO2",
         year = lubridate::year(date_local)) %>%
  dplyr::rename("so2_ppb"="arithmetic_mean") %>%
  dplyr::select(-parameter_code, -poc,-datum, -parameter, 
                -sample_duration_code, -sample_duration,
                -units_of_measure, -event_type,
                -observation_count, -observation_percent,
                -validity_indicator,
                -first_max_value, -first_max_hour, 
                -method_code, -method)
# geometry, state_code, county_code, site_number, 
# Date, Pollutant, Concentration)

# Save as a csv file. 
write_csv(file = "Data/Denver_SO2_2017_2021.csv",
          x = so2_aqs2)

# Save as .RData file.
save(file = "Data/denver_so2_2017_2021.RData", x = so2_aqs2)

load(file = "Data/denver_so2_2017_2021.RData")

#################################################################################
# Use the glimpse function from dplyr to view data.                             #
## Alternatively, use the skim function from the skimr package to see the       #
## distribution of data.                                                        #
#################################################################################

dplyr::glimpse(so2_aqs2)
skimr::skim(so2_aqs2)

#################################################################################
# Calculate summary statistics by county, site number, and year.                #
#################################################################################

so2_sum_stats <- so2_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, site_address, year, local_site_name) %>%
  summarize(count = n(),
            min = min(so2_ppb),
            max = max(so2_ppb),
            average = mean(so2_ppb),
            med = median(so2_ppb))

knitr::kable(so2_sum_stats, digits = 2)

write_csv(file = "Data/so2_locations.csv", x = so2_sum_stats)

so2_sum_stats2 <- so2_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, site_address, local_site_name) %>%
  summarize(count = n(),
            min = min(so2_ppb),
            max = max(so2_ppb),
            average = mean(so2_ppb),
            med = median(so2_ppb))

knitr::kable(so2_sum_stats2, digits = 2)

write_csv(file = "Outputs/so2_sum_stats.csv", x = so2_sum_stats2)

#################################################################################
## Check monitor locations.                                                     #
#################################################################################

ggplot() +
  geom_sf(data = den_ct_nbhd2,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = so2_aqs2)



#################################################################################
# title: "Obtaining EPA AQS for CO"                                             #
#################################################################################

# Set air parameters of interest.
## CO - 42101


#################################################################################
# Get the Data.                                                                 #
## Obtain CO monitor data across Colorado for 2017-2021.                        #
## EPA Data Frame                                                               #
## Use the aqs daily summary function under the RAQSAPI package.                #
#################################################################################

co_aqs <- aqs_dailysummary_by_state(parameter="42101",#CO
                                    bdate=as.Date("2017-01-01"),
                                    edate=as.Date("2021-12-31"),
                                    stateFIPS = "08")

# Save dataframe as a csv file. 
write_csv(file = "Data/Colorado_CO_2017_2021.csv", x = co_aqs)

co_aqs <- read_csv(file = "Data/Colorado_CO_2017_2021.csv")

#################################################################################
# Data Manipulation.                                                            #  
## Manipulate data into usable format.                                          #
## 1.Filter data to Denver metropolitan counties (Adams, Araphaoe, Broomfield,  #
### Denver, Douglas, Jefferson).                                                #
## 2. Convert data frame into an sf object.                                     #
## 3. Set the CRS to WGS84.                                                     #
## 4. Mutate date into "date" format, and pull year using the lubridate package.#
## 5. Remove some of the columns for easier future use.                         #
#################################################################################

co_aqs2 <- co_aqs %>%
  filter(#pollutant_standard == "CO 8-hour 1971",
    county_code %in% denver_metro_counties_wgs84$COUNTYFP)%>% 
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
  mutate(sample_date = as.Date(date_local),
         #pollutant = "CO",
         year = lubridate::year(date_local)) %>%
  rename("co_ppm"="arithmetic_mean") %>%
  dplyr::select(-parameter_code, -poc,-datum, -parameter,
                -sample_duration_code, -sample_duration,
                -units_of_measure, -event_type,
                -observation_count, -observation_percent,
                -validity_indicator,
                -first_max_value, -first_max_hour, 
                -method_code, -method)
# geometry, state_code, county_code, site_number, 
# Date, Pollutant, Concentration)

# Save as a csv file. 
write_csv(file = "Data/Denver_CO_2017_2021.csv",
          x = co_aqs2)

# Save as .RData file.
save(file = "Data/denver_co_2017_2021.RData", x = co_aqs2)

load(file = "Data/denver_co_2017_2021.RData")

#################################################################################
# Use the glimpse function from dplyr to view data.                             #
## Alternatively, use the skim function from the skimr package to see the       #
## distribution of data.                                                        #
#################################################################################

dplyr::glimpse(co_aqs2)
#skimr::skim(co_aqs2)

#################################################################################
# Calculate summary statistics by county, site number, and year.                #
#################################################################################

co_sum_stats <- co_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, site_address, year, local_site_name) %>%
  summarize(count = n(),
            min = min(co_ppm),
            max = max(co_ppm),
            average = mean(co_ppm),
            med = median(co_ppm))

knitr::kable(co_sum_stats, digits = 2)

write_csv(file = "Data/co_locations.csv", x = co_sum_stats)

co_sum_stats2 <- co_aqs2 %>%
  st_drop_geometry() %>%
  group_by(site_number, site_address, local_site_name) %>%
  summarize(count = n(),
            min = min(co_ppm),
            max = max(co_ppm),
            average = mean(co_ppm),
            med = median(co_ppm))

knitr::kable(co_sum_stats2, digits = 2)

write_csv(file = "Outputs/co_sum_stats.csv", x = co_sum_stats2)

#################################################################################
## Check monitor locations.                                                     #
#################################################################################

ggplot() +
  geom_sf(data = den_ct_nbhd2,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = co_aqs2)


