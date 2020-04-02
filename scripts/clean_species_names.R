#clean_species_names script, based on Viktoria Wagner's blog post on how to use the taxize package

#load libraries
library(dplyr)
library(taxize)
library(magrittr)
library(tidyverse)
library(DataCombine)
library(comprehenr)
library(sf)

# change to your project HOME directory
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/") 

# set filepath
path <- "data/semiclean" # path to your data

#list files
files <- list.files(path = path, full.names = T)

#read in data
my_data <- read.csv(files[1])
head(my_data)

#replace '#N/A' with NA across dataframe
my_data <- as.data.frame(sapply(my_data, function(x) {
  gsub('#N/A', NA, x)
}))

#split data into 2 dataframes in order to move observations from the columns genus, species to the column Genus_Species

na_gensp <- my_data %>%
  filter(is.na(Genus_Species))

not_na_gensp <- my_data %>%
  filter(!is.na(Genus_Species))

#move all values from separate genus, species columns to the Genus_Species column to prepare for taxa name spellchecking/cleaning

na_gensp$Genus_Species <- paste(na_gensp$genus, na_gensp$species)
na_gensp

#replace 'NA NA' with NA in new genus_species_column
na_gensp$Genus_Species <- sapply(na_gensp$Genus_Species, function(x) {
  gsub('NA NA', NA, x)
})


## bind the two dataframes back together
bound <- rbind(not_na_gensp, na_gensp)


#Subset df to include only useful columns
my_data_subset <- bound %>% 
                  select(Easting, Northing, Latitude, Longitude, 
                         Genus_Species, Family, Cover_Class, Number)

#test Remove unknown from gen_sp to try and retain information on type of unknown
my_data_subset$Genus_Species <- sapply(my_data_subset$Genus_Species, function(unknown_taxon) {
  gsub('Unknown', '', unknown_taxon)
})


#specify reference databases to use
sources <- c("EOL", "The International Plant Names Index", "ITIS")

#show specified sources
subset(gnr_datasources(), title %in% sources)

#submit request to global names resolver
my_data_subset$Genus_Species <- as.character(my_data_subset$Genus_Species) # convert column from factor to character vector

resolved_long <- my_data_subset$Genus_Species %>%
                 gnr_resolve(data_source_ids = c(3, 167), 
                             with_canonical_ranks=T)

#subset results to remove duplicates and to only include user supplied names, submitted names, matched names, scores 
resolved_distinct <- resolved_long %>%
                  select(user_supplied_name, submitted_name, matched_name2) %>%
                  distinct()

#join the resolved taxa df back to original
merged <- my_data_subset %>%
          left_join(resolved_distinct, by = c("Genus_Species" = "user_supplied_name")) 

#rename matched_name2 to clean taxa to explicitly indicate column with cleaned taxa names
merged <- merged %>% 
          rename(clean_taxa = matched_name2)

#remove original (erroneous) Genus_Species column, extraneous columns
pruned <- subset(merged, select = -c(submitted_name, Genus_Species))

#rename clean_taxa to Genus_species
data.table::setnames(pruned, 'clean_taxa', 'Genus_species')

#Remove var, ssp, subsp
# pruned$clean_taxa <- str_replace_all(pruned$clean_taxa, pattern = c('var. |ssp. |subsp. '), replacement = '')

#remove duplicate rows
dupl_remov <- pruned %>% 
              distinct()

#replace taxa, family NAs with unknowns to maintain data
dupl_remov$Family <- as.character(dupl_remov$Family) # convert family column to character
dupl_remov <- dupl_remov %>%
  tidyr::replace_na(list(Genus_species = 'Unknown', Family = 'Unknown'))

#remove observations where both taxa and family are unknown
cleaned_data <- dupl_remov[!(dupl_remov$Genus_species == 'Unknown' & dupl_remov$Family =='Unknown'),]

#IDIOSYNCRATIC STEP, UNECCESSARY TO REPLICATE. replace all occurrences of 'Cylindropuntia bigelovii' with 
#'Cylindropuntia echinocarpa' and also 'Cylindropuntia californica' to 'Cylindropuntia acanthocarpa' 
#'to rectify plant misidentification from data collection process
cleaned_data$Genus_species <- sapply(cleaned_data$Genus_species, function(taxon) {
  gsub("Cylindropuntia bigelovii", "Cylindropuntia echinocarpa", taxon)
})

cleaned_data$Genus_species <- sapply(cleaned_data$Genus_species, function(taxon) {
  gsub("Cylindropuntia californica", "Cylindropuntia acanthocarpa", taxon)
})

###########STILL NEED TO REMOVE ALL SINGLE WORD TAXA FROM DATA (EXCEPT PECTOCARYA & UNKNOWN)


#extract unique taxa names for use in Genbank query
taxa_list <- cleaned_data$Genus_species %>% sort() %>% unique()


#write taxa list to CSV, .txt files
write.csv(taxa_list, file = "data/clean/taxa_list.csv") # change path relevant to your directory organization
write.table(taxa_list, file = "data/clean/taxa_list.txt", sep = " ", col.names = FALSE)


###########CLEAN SPATIAL DATA START
#drop columns in Easting, Northing due to uncertainty about their origins. Will create Easting, Northing from Latitude and Longitude for all data points
#to minimize error / assumptions about data since some have Easting/Northing and some don't.
lat_long_only <- cleaned_data %>% select(-c(Easting, Northing))

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
dms_only$Latitude <- lapply(dms_only$Latitude, as.character)
dms_only$Longitude <- lapply(dms_only$Longitude, as.character)

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

#helper function to insert explicit W indicating longitude
fix_west <- function(west_symbol) {
  gsub(pattern = "\"", replacement = "\"W", west_symbol)
}

#function to convert latitude coordinates to decimal degrees
fix_latitude <- function(lat_coord){
  lat_coord1 <- fix_deg(lat_coord)
  lat_coord2 <- fix_space(lat_coord1)
  lat_coord3 <- fix_north(lat_coord2)
  dd_lat <- lat_coord3 %>% sp::char2dms() %>% as.numeric()
  dd_lat
}

#function to convert longitude coordinates to decimal degrees
fix_longitude <- function(long_coord) {
  long_coord1 <- fix_deg(long_coord)
  long_coord2 <- fix_space(long_coord1)
  long_coord3 <- fix_west(long_coord2)
  dd_long <- long_coord3 %>% sp::char2dms() %>% as.numeric()
  dd_long
}

#create test longitude vector to check if function works
test_long <- dms_only$Longitude %>% head(1)
#create test latitude vector to check if function works
test_lat <- dms_only$Latitude %>% head(1)

#see if function correctly converts longitude, it DOESNT WORK!
long_works_test <- lapply(test_long, fix_longitude)
print(test_long)
print(long_works_test)

#see if function correcty converts latitude, it WORKS!
lat_works_test <- lapply(test_lat, fix_latitude)
print(lat_works_test)
print(test_lat)

##attempt to build own DMS to DD function to see if it resolves incorrect conversion issue when using 
##sp::char2dms()
dms2dd <- function(dms_coord) {
  
}

## try to understand how to split strings, need to do this to create the deg, min, sec variables for
## dms2dd function
foo <- "35d14'4.322\"N"
deg <- str_split(foo, 'd')[[1]][1]
str_split(foo, ".*d")

l 



############LEFT OFF HERE: NEXT TIME - APPLY GSUB PATTERNS TO ENTIRE LATITUDE COLUMN, 
##LONGITUDE (BUT W NOT N)
#MAP COLUMNS IN DMS DATAFRAME TO FUNCTION, REPLACE THEM WITH DD, MERGE DATAFRAME WITH DD DATAFRAME


as.numeric(test2)
#convert dataframe to sf spatial dataframe object
dd_only_sf <- st_as_sf(dd_only, coords = c('Longitude', 'Latitude'))
dd_only_sf #inspect min, max lat longs to see if they make sense

dd_only_utm <- dd_only_sf %>%
  st_set_crs(4326) %>% #assume lat long to be WGS84
  st_transform(crs = 32611) #project data to WGS84 UTM zone 11N
dd_only_utm

do.call(st_geometry(dd_only_utm), dd_only_utm$Genus_species) %>% as_tibble() %>% setNames(c('east', 'north'))
test_df <- as.data.frame(dd_only_utm)
test_df$geometry
###extract northing, easting from geom column?
test <- do.call(rbind, st_geometry(dd_only_utm)) %>% 
  as_tibble(.name_repair = "minimal") %>% setNames(c("east","north"))
test


#############CLEAN SPATIAL DATA END

#view data
ggplot2::ggplot(spatial_bound, aes(x = Longitude, y = Latitude)) + geom_point() + 
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


