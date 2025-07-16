#-------------------------------------------------------------------------------
# title: "Setting area of interest for CEBA"  
# author: E Lunsford  
# date: 2025-07-16  
#  
# This code is to set the area of interest for further spatial coding and 
# map making.                                                              
#
# Last Run: 07/16/2025 and code was in working order using R 4.5.1 and RStudio 2025.05.1+513 
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
# Set map theme.                                                                #
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
# Set parameter variables.                                                      #
#################################################################################

## Years of interest
years <- c(2017:2021)

# Colorado Boundary
colorado_sf_nad83 <- states(year = 2021, 
                            cb = T) %>%
  filter(GEOID == "08")

save(file = "colorado_sf_nad83.RData", 
     x = colorado_sf_nad83)

# Change CRS to WGS 84
colorado_sf_wgs84 <- colorado_sf_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "colorado_sf_wgs84.RData", 
     x = colorado_sf_wgs84)

# Colorado bbox
colorado_bbox <- st_bbox(colorado_sf_nad83)

#################################################################################
# Get counties of interest                                                      #
#################################################################################

# Colorado Counties
colorado_counties_nad83 <- counties(state = "08",
                                    year = 2021,
                                    cb = T)

save(file = "colorado_counties_nad83.RData", x = colorado_counties_nad83)

# Transform to WGS 84
colorado_counties_wgs84 <- colorado_counties_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "colorado_counties_wgs84.RData", colorado_counties_wgs84)

# Denver Metro counties
## Adams (001)
## Arapaho (005)
## Broomfield (014) 
## Denver (031)
## Douglas (035)
## Jefferson (059)

counties_of_interest <- c("001", "005", "014", "031", "035", "059")

denver_metro_counties_wgs84 <- colorado_counties_wgs84 %>%
  filter(COUNTYFP  %in%  counties_of_interest)

save(file = "denver_metro_counties_wgs84.RData", denver_metro_counties_wgs84)




