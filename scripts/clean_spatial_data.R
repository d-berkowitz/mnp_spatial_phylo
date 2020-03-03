#merge different CSVs into a single master dataframe, subset that dataframe and remove duplicate rows.

library(tidyverse)
library(dplyr)
library(readr)

setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo") # change to your project HOME directory directory

# load all csv files and merge into one large data table
path <- "data/raw/spatial" # path to your raw data, if different organization than this

files <- list.files(path = path, full.names = T)

# make a list of the files that are loaded
tbl <- sapply(files, read_csv, simplify=FALSE) 

# canonicalize column names so they properly join, remove erroneous/misspelled column names
for (i in seq_along(tbl)){
  colnames(tbl[[i]])[grep("Cover Class", colnames(tbl[[i]]))] <-"Cover_Class"
  colnames(tbl[[i]])[grep("LAT", colnames(tbl[[i]]))] <-"Latitude"
  colnames(tbl[[i]])[grep("LONG", colnames(tbl[[i]]))] <-"Longitude"
  colnames(tbl[[i]])[grep("Genius_Species", colnames(tbl[[i]]))] <-"Genus_Species"
  
}

# bind all data tables by column names
bound <- bind_rows(lapply(tbl, function(dtt){mutate_all(dtt, as.character)}))

# remove any fully duplicated rows
data <- distinct(bound)

# view column names to find if any are duplicates that should be renamed before bind_rows() as above
colnames(data) 

# pull out columns we're interested in:
# "" "Easting", "Northing", Longitude", "Latitude", Genus_Species", "genus", "species", , "Cover_Class, "Number" 

spatial_data <- select(data, Easting, Northing, Latitude, Longitude, Genus_Species, genus, species, Family, Cover_Class, Number)

#add subspecies column for tidy data
spatial_data <- spatial_data %>%
                  mutate(subspecies = NA)


##export to csv in directory for semiclean data
write_csv(spatial_data, path = "data/semiclean/spatial_data", na = "NA", 
          append = FALSE, col_names = TRUE, 
          quote_escape = "double")

# next tasks: 
# fix taxonomy 
# convert easting and northing to lat long, or vice versa
# convert genus_species to genus and species (separate columns), additional column for subspecies/variety
# combine this script with dean's

