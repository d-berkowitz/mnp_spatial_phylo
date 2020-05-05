##clean spatial data script. rectifies issues regarding data that contains UTM, 
##Lat/long, decimal degrees (DD), degree minute seconds (DMS). converts data to 
## sf object (spatial points data frame) for projection / reprojection.

#load libraries
library(dplyr)
library(tidyverse)
library(sf)
library(mapview)
library(leaflet)

# change to your project HOME directory
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/") 

# set filepath
path <- "data/semiclean/clean_species_names_with_unknowns.csv" # path to your data

my_data <- read.csv(file = path)
my_data

#drop columns in Easting, Northing due to uncertainty about their origins. 
#Will create Easting, Northing from Latitude and Longitude for all data points
#to minimize error / assumptions about data since some have Easting/Northing and some don't.
lat_long_only <- my_data 
#%>% select(-c(Easting, Northing))

#create new df with lat, long in dms only to clean and prepare for conversion to dd
dd_only <- lat_long_only %>% 
  filter(!grepl(pattern = c(' '), Latitude))

dms_only <- lat_long_only %>% 
  filter(grepl(pattern = c(' '), Latitude))

# remove ° symbol from longitudes (idiosyncracy of the data I am working with)
dd_only$Longitude <- dd_only$Longitude %>%
  gsub(pattern = '°', replacement = '')
#convert latitude, longitude columns from factor to double in preparation for creating northing, easting
dd_only$Latitude <- as.numeric(as.character(dd_only$Latitude))
dd_only$Longitude <- as.numeric(dd_only$Longitude)

#remove rows where latitude == longitude
dd_only <- dd_only %>%
  filter(!dd_only$Longitude == dd_only$Latitude)


#convert DMS to DD 
#convert latitude, longitude columns to character vectors to prepare for conversion to decimal degrees
# dms_only$Latitude <- lapply(dms_only$Latitude, as.character)
# dms_only$Longitude <- lapply(dms_only$Longitude, as.character)

# convert degree, minute, second characters in Latitude and Longitude column to 
#format required by char2dms function

#helper function to fix degree symbol
fix_deg <- function(degree_symbol) {
  gsub(pattern = "° ", replacement = "d", degree_symbol)
}

#helper function to remove spaces
fix_space <- function(space_symbol) {
  gsub(pattern = " ", replacement = "", space_symbol)
}

#helper function to insert explicit N indicating latitude
fix_north <- function(north_symbol) {
  gsub(pattern = "\"", replacement = "\"N", north_symbol)
}

#remove negatives from Longitude to prepare for inserting explicit W
remove_neg <- function(neg_symbol) {
  gsub(pattern = "-", replacement = "", neg_symbol)
}

#helper function to insert explicit W indicating longitude
fix_west <- function(west_symbol) {
  gsub(pattern = "\"", replacement = "\"W", west_symbol)
}


#function to convert latitude coordinates to decimal degrees
fix_latitude <- function(lat_coord){
  lat_char <- lapply(lat_coord, as.character)
  lat_coord1 <- fix_deg(lat_char)
  lat_coord2 <- fix_space(lat_coord1)
  lat_coord3 <- fix_north(lat_coord2)
  dd_lat <- lat_coord3 %>% sp::char2dms() %>% as.double()
  dd_lat
}

#function to convert longitude coordinates to decimal degrees
fix_longitude <- function(long_coord) {
  long_coord_pos <- remove_neg(long_coord)
  long_coord1 <- fix_deg(long_coord_pos)
  long_coord2 <- fix_space(long_coord1)
  long_coord3 <- fix_west(long_coord2)
  dd_long <- long_coord3 %>% sp::char2dms() %>% as.numeric() 
  dd_long
}

#convert DMS to DD
dms_only$Latitude_DD <- map_dbl(dms_only$Latitude, fix_latitude)
dms_only$Longitude_DD <- map_dbl(dms_only$Longitude, fix_longitude)

# drop DMS latitude, longitude columns
dms_now_dd <- dms_only %>% select(-c('Latitude', 'Longitude'))
#rename latitude_DD, longitude_DD columns to prepare for merge
dms_now_dd <- dms_now_dd %>% rename('Latitude' = 'Latitude_DD',
                    'Longitude' = 'Longitude_DD')

#remove decimal degree columns from dms_only df
dms_only <- dms_only %>% select(-c('Latitude_DD', 'Longitude_DD'))

#merge spatial dataframes back together
spatial_clean <- rbind(dms_now_dd, dd_only)

# filter out erroneous lat / long values that are outside region of interest
spatial_clean <- spatial_clean %>%
  filter(Latitude <= 35.241)

#reorder columns 
spatial_clean <- spatial_clean %>% select(
  c('Latitude', 'Longitude', 'Genus_species', 'Family', 'Cover_Class', 'Number', 'Observation'))

# spatial_no_unk <- spatial_clean %>% filter(Genus_species != 'Unknown' &
#                                              Family != 'Unknown')
# 
# unknown <- spatial_clean %>% filter(Genus_species == 'Unknown' &
#                                       Family == 'Unknown')
# nrow(unknown)
# 
# ten_seventeen <- spatial_clean %>% filter(grepl('MJV10-17.', Observation))

#plot data to figure out what field season points are from

leaflet() %>% 
  addProviderTiles(provider =  "Esri.WorldImagery") %>% 
  addCircleMarkers(data = spatial_no_unk, radius = 3, color = 'blue', opacity = 0.025, 
                   weight = .05, label = spatial_no_unk$Genus_species) %>%
  addCircleMarkers(data = unknown, radius = 3, color = 'red', opacity = 0.025,
                   weight = .05, label = unknown$Observation)

###############STOP HERE

#convert dataframe to sf spatial dataframe object
sp_clean_sf <- st_as_sf(spatial_clean, coords = c('Longitude', 'Latitude'))
sp_clean_sf #inspect min, max lat longs to see if they make sense

#project sf dataframe to UTM Zone 11N
sp_clean_utm <- sp_clean_sf %>%
  st_set_crs(4326) %>% #assume lat long to be WGS84
  st_transform(crs = 32611) #project data to WGS84 UTM zone 11N
sp_clean_utm




#do.call(st_geometry(dd_only_utm), dd_only_utm$Genus_species) %>% as_tibble() %>% setNames(c('east', 'north'))
#test_df <- as.data.frame(dd_only_utm)
#test_df$geometry
###extract northing, easting from geom column?
#test <- do.call(rbind, st_geometry(dd_only_utm)) %>% 
#  as_tibble(.name_repair = "minimal") %>% setNames(c("east","north"))
#test


#############CLEAN SPATIAL DATA END

#view data
ggplot2::ggplot(spatial_clean, aes(x = Longitude, y = Latitude)) + geom_point() + 
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank())

######summary statistics 

#find number of observations in df
df_len <- nrow(my_data)

#check number of observations that have lat/long
num_latlong <- sum(!is.na(my_data$Latitude))
#sum(is.na(my_data$Longitude)) # check to see if number that has long matches number that has lat

num_utm <- sum(!is.na(my_data$Easting))# check to see number of observations that have easting /northing
#sum(is.na(my_data$Northing)) 

prop_utm <- round((num_utm / df_len), 3)
prop_latlong <- round((num_latlong / df_len), 3)

paste0(prop_utm*100, '% of the data have northing / easting, while ', prop_latlong*100, ' % of the data have lat/long')


#######atttempt to resolve unknown moss / licken issue

example_moss <- c('gray moss', 'Moss', 'Bright Green Moss', 'dried moss', 'moss', 'gray-brown moss')
example_moss