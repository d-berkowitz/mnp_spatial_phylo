#aggregate spatial data script.
#merge different CSVs into a single master dataframe, subset that dataframe and remove duplicate rows.

library(tidyverse)
library(dplyr)
library(readr)

setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo") # change to your project HOME directory directory

# load all csv files and merge into one large data table
path <- "data/raw/spatial" # path to your raw data, if different organization than this

files <- list.files(path = path, full.names = T)

# make a list of the files that are loaded
my_tbl <- sapply(files, read_csv, simplify=FALSE) 

# canonicalize column names so they properly join, remove erroneous/misspelled column names
for (i in seq_along(my_tbl)){
  colnames(my_tbl[[i]])[grep("Cover Class", colnames(my_tbl[[i]]))] <-"Cover_Class"
  colnames(my_tbl[[i]])[grep("LAT", colnames(my_tbl[[i]]))] <-"Latitude"
  colnames(my_tbl[[i]])[grep("LONG", colnames(my_tbl[[i]]))] <-"Longitude"
  colnames(my_tbl[[i]])[grep("Genius_Species", colnames(my_tbl[[i]]))] <-"Genus_Species"
  
}

# bind all data tables by column names
bound <- bind_rows(lapply(my_tbl, function(dtt){mutate_all(dtt, as.character)}))

# remove any fully duplicated rows
my_data <- distinct(bound)

# view column names to find if any are duplicates that should be renamed before bind_rows() as above
my_data %>% colnames() %>% sort()

# pull out columns we're interested in
spatial_data <- my_data %>% 
  select(Latitude, Longitude, 
         Genus_Species, genus, species, Family, 
         Cover_Class, Number, Observation, `Observation #`)


##export to csv in directory for semiclean data
write_csv(spatial_data, path = "data/semiclean/spatial_data.csv", na = "NA", 
          append = FALSE, col_names = TRUE, 
          quote_escape = "double")



