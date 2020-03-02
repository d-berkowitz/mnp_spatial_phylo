library(tidyverse)
library(dplyr)
library(readr)
library(sp)

setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/") # change to your project HOME directory directory

# load all csv files and merge into one large data table
path <- "data/raw/spatial" # path to your raw data, if different organization than this

files <- list.files(path = path, full.names = T)
tbl <- sapply(files, read_csv, simplify=FALSE) %>% 
 bind_rows(.id = "id")

#check files
test_df <- read_csv(files[1])
class(test_df$Genus_Species)

# remove any fully duplicated rows
my_data <- distinct(test_df)

# pull out columns we're interested in
spatial_data <- data[]

# plot to view initial shape
coordinates(my_data) <-c('Longitude', 'Latitude')
plot(my_data)


