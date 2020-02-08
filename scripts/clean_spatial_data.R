library(tidyverse)
library(dplyr)
library(readr)

setwd("/Users/jennaekwealor/Documents/dean_project/mnp_spatial_phylo") # change to your project HOME directory directory

# load all csv files and merge into one large data table
path <- "data/raw/spatial" # path to your raw data, if different organization than this

files <- list.files(path = path, full.names = T)
tbl <- sapply(files, read_csv, simplify=FALSE) %>% 
 bind_rows(.id = "id")

# remove any fully duplicated rows
data <- distinct(data)

# pull out columns we're interested in
spatial_data <- data[]

# 