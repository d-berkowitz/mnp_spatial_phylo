##clean spatial data script. rectifies issues regarding data that contains UTM, 
##Lat/long, decimal degrees (DD), degree minute seconds (DMS). converts data to 
## sf object (spatial points data frame) for projection / reprojection.

#load libraries
library(dplyr)
library(tidyverse)
library(sf)
library(mapview)
library(leaflet)
library(ggmap)
library(rgdal)
library(raster)
library(ggplot2)
# change to your project HOME directory
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/") 

# set filepath
path <- "data/semiclean/clean_species_names.csv" # path to your data

raw_spatial_data <- read.csv(file = path)

#create new df with lat, long in dms only to clean and prepare for conversion to dd
dd_only <- raw_spatial_data %>% 
  filter(!grepl(pattern = c(' '), Latitude))

dms_only <- raw_spatial_data %>% 
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
dms_now_dd <- dms_only %>% dplyr::select(-c('Latitude', 'Longitude'))
#rename latitude_DD, longitude_DD columns to prepare for merge
dms_now_dd <- dms_now_dd %>% dplyr::rename('Latitude' = 'Latitude_DD',
                    'Longitude' = 'Longitude_DD')

#remove decimal degree columns from dms_only df
dms_only <- dms_only %>% dplyr::select(-c('Latitude_DD', 'Longitude_DD'))

#merge spatial dataframes back together
spatial_clean <- rbind(dms_now_dd, dd_only)

# filter out erroneous lat / long values that are outside region of interest
spatial_clean <- spatial_clean %>%
  filter(Latitude <= 35.241)
spatial_clean <- spatial_clean %>%
  filter(Longitude < -115)

nrow(spatial_clean)

#remove duplicate species
spatial_clean <- spatial_clean %>% distinct(Latitude, Longitude, Genus_species, .keep_all = TRUE)
nrow(spatial_clean)
#reorder columns 
spatial_clean <- spatial_clean %>% dplyr::select(
  c('Latitude', 'Longitude', 'Genus_species', 'Family', 'Cover_Class', 'Number', 'Observation'))


spatial_clean_no_unk <- spatial_clean %>% filter(Genus_species != 'Unknown')
nrow(spatial_clean_no_unk)
unknown <- spatial_clean %>% filter(Genus_species == 'Unknown')
nrow(unknown)

#plot data to figure out what field season points are from

leaflet() %>% 
  addProviderTiles(provider =  "Esri.WorldImagery") %>% 
  addCircleMarkers(data = spatial_clean_no_unk, radius = 3, color = 'blue', opacity = 0.025, 
                   weight = .05, label = spatial_clean_no_unk$Genus_species) 
# %>%
#   addCircleMarkers(data = unknown, radius = 3, color = 'red', opacity = 0.025,
#                    weight = .05, label = unknown$Genus_species)

#extract unique taxa names for use in Genbank query
clean_spatial_taxa_list <- spatial_clean_no_unk$Genus_species %>% sort() %>% unique()
clean_spatial_taxa_list_df <- data.frame(taxa = clean_spatial_taxa_list)
#write taxa list to CSV, .txt files
write.csv(clean_spatial_taxa_list_df, file = "data/clean/clean_taxa_list.csv", row.names = FALSE) # change path relevant to your directory organization
write.table(clean_spatial_taxa_list_df, file = "data/clean/clean_taxa_list.txt", sep = " ", col.names = FALSE)


###############STOP HERE

#convert dataframe to sf spatial dataframe object
sp_clean_sf <- st_as_sf(spatial_clean_no_unk, coords = c('Longitude', 'Latitude'))
sp_clean_sf #inspect min, max lat longs to see if they make sense

#project sf dataframe to UTM Zone 11N
sp_clean_utm <- sp_clean_sf %>%
  st_set_crs(4326) %>% #assume lat long to be WGS84
  st_transform(crs = 32611) #project data to WGS84 UTM zone 11N
sp_clean_utm

#write function to split geometry column into Easting, Northing
sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

#split geometry column
sp_clean_utm_biodiverse <- sfc_as_cols(sp_clean_utm, names = c("Easting", "Northing"))
#convert to dataframe in format required for Biodiverse software
sp_clean_utm_biodiverse_df <- sp_clean_utm_biodiverse %>% 
                              as.data.frame() %>% 
                              arrange(Genus_species) %>% 
                              dplyr::select(c('Genus_species', 'Easting', 'Northing'))
#write spatial data csv for use in Biodiverse
write.csv(sp_clean_utm_biodiverse_df, file = "data/clean/spatial/clean_spatial_utm_biodiverse.csv")


##Spatial analysis

species_data <- sp_clean_utm_biodiverse %>%
  as.data.frame() %>%
  arrange(Genus_species) %>%
  dplyr::select(c('Genus_species', 'Easting', 'Northing'))



#######GET FAMILY INFORMATION

library(myTAI)
library(furrr)

get_fam <- function(sp){
            taxonomy(organism = sp, 
            db = "itis", output = "classification") %>%
            filter(rank == 'family') %>%
            dplyr::select('name') %>%
            pull()
}

sp_df <- species_data %>%
  dplyr::select('Genus_species') %>%
  distinct() %>% 
  filter(Genus_species != 'Castilleja chromosa') %>%
  filter(Genus_species != 'Cleomella arborea') %>%
  filter(Genus_species != 'Johnstonella angustifolia') %>%
  filter(Genus_species != 'Lotus strigosus') %>%
  filter(Genus_species != 'Munroa pulchella') %>%
  filter(Genus_species != 'Scutellaria mexicana') %>%
  filter(Genus_species != 'Stipa hymenoides')


plan(multiprocess)
sp_df$Family <- future_map_chr(sp_df$Genus_species, get_fam)
              
syn_error_sp <- c('Castilleja chromosa', 'Cleomella arborea', 'Johnstonella angustifolia',
               'Lotus strigosus', 'Munroa pulchella', 'Scutellaria mexicana',
               'Stipa hymenoides')

syn_error_fam <- c('Orobanchaceae', 'Cleomaceae', 'Boraginaceae',
                   'Fabaceae', 'Poaceae', 'Lamiaceae',
                   'Poaceae')

syn_error_df <- data.frame('Genus_species' = syn_error_sp,
                           'Family' = syn_error_fam)

sp_df_bound <- rbind(sp_df, syn_error_df) %>%
  arrange(Genus_species)

species_data <- left_join(species_data, sp_df_bound, by = 'Genus_species')

abundance <- species_data %>% 
  group_by(Genus_species) %>%
  summarise(Abundance = n()) %>%
  arrange(desc(Abundance))

abundance$Rank <- as.numeric(as.factor(abundance$Abundance))

abundance %>% merge(abundance %>% group_by(Rank) %>% summarise(dens = n())) -> abundance

top20 <- abundance %>% filter(Rank > 33)
top20$Abundance %>% sum() / sum(abundance$Abundance)
############### MAKE FIGURES

library(wesanderson)
library(RColorBrewer)
library(extrafont)
library(viridis)


join_test <- join_test %>% arrange(Abundance) 
#run together
set.seed(8)

ggplot(join_test, aes(x = reorder(Genus_species, -Abundance))) + 
geom_bar(stat = "count")



#run together
set.seed(8)
ggplot(abundance, aes(x = (Rank), y = Abundance)) + 
  geom_jitter(height = 0.2, width = 0.5,
              alpha = 0.2) + 
  labs(x = "Species rank (from least to most abundant)", 
       y = "Log abundance (# of individuals)",
       size = "# of species that share rank") +
  guides(size = guide_legend(override.aes = list(size = c(1, 2, 3, 4))))+
  ggtitle("Unevenness of taxa abundance across study area") +
  geom_vline(xintercept = 32.375, color = "red", size = .5) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(1, 10, 50, 100, 600)) + 
  scale_size(range = c(.75, 3)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 12, family = "Times New Roman")) 


#plot families of 10 most abundant taxa by rank
abund_10 <- abundance %>% head(10)

family_abundance <- species_data %>% 
  group_by(Family) %>%
  summarise(abundance = n()) %>%
  arrange(desc(abundance))

family_num_taxa <- species_data %>% 
  dplyr::select('Genus_species', 'Family') %>%
  distinct() %>%
  group_by(Family) %>%
  summarise(num_taxa = n()) %>%
  arrange(desc(num_taxa))

family_stats <- left_join(family_abundance, family_num_taxa, by = 'Family') %>%
  dplyr::rename(num_sp = num_taxa) %>%
  dplyr::rename(fam_abundance = abundance)


sp_fam <- species_data %>% 
  dplyr::select('Genus_species', 'Family') %>%
  distinct()



font_import()
fonts()

abundance$Family < factor(abundance$Family, levels = abundance$Rank)
abundance$Abundance %>% class()


#Plot species evenness
abundance %>% arrange(desc(Abundance)) %>%
  ggplot(aes(x= reorder(Genus_species, -Abundance), y= Abundance)) +
  geom_col() +
  theme_minimal(base_family = "Times New Roman") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x = 'Species', y = 'Abundance (# of individuals',
       title = 'Uneven Species Abundance') 

#Most abundant species in study area
abundance %>% arrange(desc(Abundance)) %>% head(6) %>%
  ggplot(aes(x = reorder(Genus_species, -Abundance), y = Abundance)) +
  geom_col(aes(fill = reorder(Family, Abundance)),
           colour = "black",
           size = 0.2) +
  scale_fill_brewer(palette = "YlOrRd") + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_minimal(base_family = "Times New Roman") + 
  labs(y = "# of individuals", x = "Genus species",
       fill = "Family") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
  ggtitle("Most Abundant Species")

bluecols <- brewer.pal(6, 'Blues')

#Diversity of families associated with most abundance
family_stats %>% arrange(desc(fam_abundance)) %>% head(6) %>%
  ggplot(aes(x = reorder(Family, -fam_abundance), y = fam_abundance)) + 
  geom_col(aes(fill = reorder(Family, fam_abundance)),
           colour = "black",
           size = 0.2) +
  scale_fill_brewer(palette = "Blues") + 
  theme_minimal(base_family = "Times New Roman") + 
  labs(y = "total # of individuals", x = "Family") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) + 
  ggtitle("Most Abundant Families") +
  theme(legend.position = "none")


  

merged_for_facet <- left_join(abundance, family_stats, by = 'Family')

#plot species + families
ggplot(merged_for_facet, aes(x = Abundance, y = Genus_))

#plot most abundant families
family_stats %>% arrange(desc(num_sp)) %>% head(6) %>%
  ggplot(aes(y = Family, x = abundance)) + 
  geom_col(fill = family_stats$num_sp)



# log normal distribution of species abundance preston
hist(abundance$Abundance, breaks = c(1, 2, 5, 50, 100, 700))

min(abundance$Abundance)
# log_plot <- ggplot(abundance, aes(x = Rank, y = Abundance)) + 
#   geom_point() + 
#   geom_line()

st_read("shapefile/35m_richness.shp")

richness <- raster(x = '35m_richness.tif')
values(richness)
plot(richness)
richness_utm <- richness %>%
                st_set_crs(4326) %>% 
                st_transform(crs = 32611)
#do.call(st_geometry(dd_only_utm), dd_only_utm$Genus_species) %>% as_tibble() %>% setNames(c('east', 'north'))
#test_df <- as.data.frame(dd_only_utm)
#test_df$geometry
###extract northing, easting from geom column?
#test <- do.call(rbind, st_geometry(dd_only_utm)) %>% 
#  as_tibble(.name_repair = "minimal") %>% setNames(c("east","north"))
#test


#############CLEAN SPATIAL DATA END

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