library(tidyverse)
library(dplyr)
library(readr)

setwd("/Users/jennaekwealor/Documents/dean_project/mnp_spatial_phylo") # change to your project HOME directory directory

# load all csv files and merge into one large data table
path <- "data/raw/spatial" # path to your raw data, if different organization than this

files <- list.files(path = path, full.names = T)

# make a list of the files that are loaded
tbl <- sapply(files, read_csv, simplify=FALSE) 

# replace column name "Cover Class" with "Cover_Class" so they properly join
for (i in seq_along(tbl)){
  colnames(tbl[[i]])[grep("Cover Class", colnames(tbl[[i]]))] <-"Cover_Class"
}

# bind all data tables by column names
bound <- bind_rows(lapply(tbl, function(dtt){mutate_all(dtt, as.character)}))

# remove any fully duplicated rows
data <- distinct(bound)

# view column names to find if any are duplicates that should be renamed before bind_rows() as above
colnames(data) 


# pull out columns we're interested in:
# "Genus_Species", "genus", "species", "Longitude", "Latitude", "Easting", "Northing", "Cover_Class" 

spatial_data <- select(data, Easting, Northing, Latitude, Longitude, Genus_Species, genus, species, Family, Cover_Class)

# remove any rows that have NA in the family column, as as short cut to removing unknown genus and species 
spatial_data <- spatial_data[!is.na(spatial_data$Family),]


# next tasks: 

# convert easting and northing to lat long, or vice versa
# convert genus_species to genus and species (separate columns)
