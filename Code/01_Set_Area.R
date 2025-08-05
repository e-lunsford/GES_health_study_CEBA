#-------------------------------------------------------------------------------
# title: "Setting area of interest for CEBA"  
# author: E Lunsford  
# date: 2025-07-30 
#  
# This code is to set the area of interest for further spatial coding and 
# map making.                                                              
#
# Last Run: 07/30/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#-------------------------------------------------------------------------------

#################################################################################
# Load Libraries
#################################################################################
library(broom)
library(colorspace)
library(ggspatial)
library(ggthemes)
library(ggplot2)
library(keyring)
library(knitr)
library(lubridate)
library(purrr)
library(RAQSAPI)
library(raster)
library(sf)
library(stringr)
library(tidycensus)
library(tidyverse)
library(tigris)
library(units)

#################################################################################
# Set map theme.
#################################################################################
map_theme <- theme(
  #aspect.ratio = 1,
  text  = element_text(size = 12, color = 'black'),
  #panel.spacing.y = unit(0,"cm"),
  #panel.spacing.x = unit(0.25, "lines"),
  panel.grid.minor = element_line(color = "transparent"),
  panel.grid.major = element_line(color = "transparent"),
  panel.border = element_blank(),
  panel.background=element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  # legend.position = c(0.1,0.1),
  # plot.margin = grid::unit(c(0,0,0,0), "mm"),
  legend.key = element_blank(),
  #legend.background = element_rect(fill='transparent'),
  plot.margin = unit(c(0, 0, 0, 0), "cm"),
  panel.spacing = unit(c(0, 0, 0, 0), "cm"),
  legend.position = "inside",
  legend.position.inside = c(0.84, 0.4),
  legend.title = element_text(size = 12),
  legend.background = element_rect(fill = "white"),
  legend.spacing.y = unit(0.2, 'cm'),
  legend.margin = margin(0.1,0,0,0, unit="cm")
)

den_col <- "darkgrey"
ges_col <- "#e41a1c"
oth_col <- "#4daf4a"
hwy_col <- "darkblue"

blue <- "#377eb8"
red <- "#e41a1c"
green <- "#4daf4a"
purple <- "#984ea3"
orange <- "#ff7f00"  
yellow <- "#ffff33"

#################################################################################
# Set parameter variables.
#################################################################################

## Years of interest
years <- c(2017:2021)

# Colorado Boundary
colorado_sf_nad83 <- states(year = 2021, cb = T) %>%
  filter(GEOID == "08")

save(file = "Outputs/colorado_sf_nad83.RData", 
     x = colorado_sf_nad83)

# Change CRS to WGS 84
colorado_sf_wgs84 <- colorado_sf_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "Outputs/colorado_sf_wgs84.RData", 
     x = colorado_sf_wgs84)

#################################################################################
# Get counties of interest.
#################################################################################

# Colorado Counties
colorado_counties_nad83 <- counties(state = "08",
                                    year = 2021,
                                    cb = T)

save(file = "Outputs/colorado_counties_nad83.RData", 
     x = colorado_counties_nad83)

# Transform to WGS 84
colorado_counties_wgs84 <- colorado_counties_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "Outputs/colorado_counties_wgs84.RData", 
     x = colorado_counties_wgs84)

# Denver Metro counties
## Adams (001)
## Arapaho (005)
## Broomfield (014) 
## Denver (031)
## Douglas (035)
## Jefferson (059)

counties_of_interest <- c("001", "005", "014", "031", "035", "059")

denver_metro_counties_nad83 <- colorado_counties_nad83 %>%
  filter(COUNTYFP %in% counties_of_interest)

save(file = "Outputs/denver_metro_counties_nad83.RData",
     x = denver_metro_counties_nad83)

denver_metro_counties_wgs84 <- colorado_counties_wgs84 %>%
  filter(COUNTYFP  %in%  counties_of_interest)

save(file = "Outputs/denver_metro_counties_wgs84.RData",
     x = denver_metro_counties_wgs84)


#################################################################################
# Census tract boundaries.
#################################################################################

# Colorado census tracts
colorado_ct_nad83 <- tracts(state = "08", year = 2020, cb = T)

save(file = "Outputs/colorado_ct_nad83.RData", x = colorado_ct_nad83)

colorado_ct_wgs_84 <- colorado_ct_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "Outputs/colorado_ct_wgs_84.RData", colorado_ct_wgs_84)

#############################
# Denver metropolitan census tracts
metro_tract_nad83 <- colorado_ct_nad83 %>%
  filter(COUNTYFP %in% counties_of_interest)

save(file = "Outputs/metro_tract_nad83.RData",
     x = metro_tract_nad83)

metro_tract_wgs84 <- colorado_ct_wgs_84 %>%
  filter(COUNTYFP  %in%  counties_of_interest)

save(file = "Outputs/metro_tract_wgs84.RData", 
     x = metro_tract_wgs84)

#############################
# City and County of Denver census tracts
den_tract_nad83 <- colorado_ct_nad83 %>%
  filter(COUNTYFP == "031")

save(file = "Outputs/den_tract_nad83.RData",
     x = den_tract_nad83)

den_tract_wgs84 <- colorado_ct_wgs_84 %>%
  filter(COUNTYFP == "031")

save(file = "Outputs/den_tract_wgs84.RData", 
     x = den_tract_wgs84)


#################################################################################
## Set bbox.
#################################################################################

# Colorado bbox
colorado_bbox <- st_bbox(colorado_sf_nad83)

# Denver metro bbox
metro_bbox <- st_bbox(metro_tract_wgs84) + c(-0.005,-0.005,0.005,0.005)

## City and County of Denver bbox
den_bbox <- st_bbox(den_tract_wgs84) + c(-0.005,-0.005,0.005,0.005)

#################################################################################
# Statistical Neighborhoods.
#
# Statistical Neighborhoods are typically combinations of census tracts.
# Geographic place names, such as Windsor and Mar Lee, were assigned to each
# area and reflect commonly used names of subdivisions and historical
# parts of the city.
# 
# Download Denver's statistical neighborhoods from Denver Open Data Catalog.
# https://www.denvergov.org/opendata/dataset/city-and-county-of-denver-statistical-neighborhoods
# 
#################################################################################

den_nbhd <- 
  sf::st_read("Data/statistical_neighborhoods/statistical_neighborhood.shp") %>%
  st_make_valid() %>%
  select(-TYPOLOGY, -NOTES)

plot(st_geometry(den_nbhd))

# Transform to NAD83 (CRS: 4269)
den_nbhd_nad83 <- den_nbhd %>%
  sf::st_transform(crs = 4269)

save(file = "Outputs/den_nbhd_nad83.RData",
     x = den_nbhd_nad83)

# Transform to WGS84 (CRS: 4326)
den_nbhd_wgs84 <- den_nbhd %>%
  sf::st_transform(crs = 4326)

save(file = "Outputs/den_nbhd_wgs84.RData", 
     x = den_nbhd_wgs84)

#################################################################################
## Check neighborhood overlays with census tracts.
#  
## We can see that the neighborhoods have a combination of census tracts, and 
## in some instances, census tracts cover multiple neighborhoods. For example,
## GEOID: 08031007091 is split between Lowry Field and Windsor neighborhoods.
# 
#################################################################################

# Figure
ggplot() +
  geom_sf(data = den_tract_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = den_nbhd_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "red",
          linewidth = 1.25) +
  map_theme


#################################################################################
# Join census tracts and neighborhoods.
# 
### Method 1: Don't specify join command in st_join. This results in
### multiple entries per GEOID because census tracts and nbhds don't
### perfectly align. 
# 
### Method 2: include specific join command to remove duplicated entries.
# 
#################################################################################

# Method 1:
 den_tract_nbhd <- den_tract_wgs84 %>%
  st_join(den_nbhd_wgs84)

nrow(den_tract_nbhd)
# 608
length(unique(den_tract_nbhd$GEOID))
# 178

save(file = "Outputs/den_tract_nbhd.RData", x = den_tract_nbhd)

# Method 2: 
den_ct_nbhd <- sf::st_join(x = den_tract_wgs84,
                           y = den_nbhd_wgs84,
                           join = st_intersects,
                           left = TRUE,
                           largest = TRUE)


nrow(den_ct_nbhd)
#178
length(unique(den_ct_nbhd$GEOID))
#178

save(file = "Outputs/den_ct_nbhd.RData", x = den_ct_nbhd)


#################################################################################
# Filter data to remove DIA and pull GES specifically.
#
# You can use the grep function to remove rows.
# Alternatively, you can use filter with "!"
#
#################################################################################

## DIA CT are 83.88, 83.89, 9800.01, 158
# Largest DIA area is CT 9800.01
dia_tract_id <- c("08031980001")

# Remove DIA from Denver neighborhood & tract data set.
den_ct_nbhd2 <- den_ct_nbhd %>%
  dplyr::filter(!GEOID %in% dia_tract_id)

save(file = "Outputs/den_ct_nbhd2.RData", x = den_ct_nbhd2)

# Create separate GES sf objects
ges_tract_ids <- c("08031001500", "08031003501", "08031003502")

# Pull GES from condensed data set. (3 observations)
ges_ct_nbhd <- den_ct_nbhd2 %>%
  filter(GEOID %in% ges_tract_ids)

save(file = "Outputs/ges_ct_nbhd.RData", x = ges_ct_nbhd)

# Test plot of filtered data
ggplot() +
  geom_sf(data = den_ct_nbhd2,
          color = "darkgreen",
          fill = NA,
          inherit.aes = F) +
  geom_sf(data = ges_ct_nbhd,
          color = ges_col,
          fill = NA,
          inherit.aes = F) +
  map_theme


