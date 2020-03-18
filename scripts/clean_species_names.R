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

###########STILL NEED TO REMOVE ALL SINGLE WORD TAXA FROM DATA (EXCEPT PECTOCARYA & UNKNOWN)


#extract unique taxa names for use in Genbank query
taxa_list <- cleaned_data$Genus_species %>% sort() %>% unique()

#write taxa list to CSV, .txt files
write.csv(taxa_list, file = "data/clean/taxa_list.csv") # change path relevant to your directory organization
write.table(taxa_list, file = "data/clean/taxa_list.txt", sep = " ", col.names = FALSE)


###########CLEAN SPATIAL DATA START
#create new df with lat, long in dms only to clean and prepare for convesion to dd
dd_only <- cleaned_data %>% 
  filter(!grepl(pattern = c(' '), Latitude))

# remove ° symbol from longitudes (idiosyncracy of the data I am working with)
dd_only$Longitude <- dd_only$Longitude %>%
  gsub(pattern = '°', replacement = '')
#convert latitude, longitude columns from factor to double in preparation for creating northing, easting
dd_only$Latitude <- as.numeric(as.character(dd_only$Latitude))
dd_only$Longitude <- as.numeric(dd_only$Longitude)

#remove rows where latitude == longitude
dd_only <- dd_only %>%
  filter(!dd_only$Longitude == dd_only$Latitude)
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

dd_only$Latitude <- as.numeric(dd_only$Latitude)
sapply(dd_only$Latitude, as.numeric)
dms_only <- cleaned_data %>% 
  filter(grepl(pattern = c(' '), Latitude))

#replace all special characters with spaces

dms_only$Latitude <- gsub('°', '', dms_only$Latitude)
dms_only$Latitude <- gsub('\'', '', dms_only$Latitude)
dms_only$Latitude <- gsub('"', '', dms_only$Latitude)
dms_only$Longitude <- gsub('°', '', dms_only$Longitude)
dms_only$Longitude <- gsub('\'', '', dms_only$Longitude)
dms_only$Longitude <- gsub('"', '', dms_only$Longitude)

view(dms_only)
#convert from DMS to DD
dms_only$Latitude <- measurements::conv_unit(dms_only$Latitude, from = 'deg_min_sec', to = 'dec_deg')
dms_only$Longitude <- measurements::conv_unit(dms_only$Longitude, from = 'deg_min_sec', to = 'dec_deg')

#join to original df by binding rows of DD and newly converted DMS to DD

spatial_bound <- rbind(dms_only, dd_only)
view(spatial_bound)
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

                        
