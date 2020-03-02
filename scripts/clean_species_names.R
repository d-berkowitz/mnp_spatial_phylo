#clean_species_names script, based on Viktoria Wagner's blog post on how to use the taxize package

#load libraries
library(dplyr)
library(taxize)
library(magrittr)
library(tidyverse)
library(DataCombine)

# change to your project HOME directory directory
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/") 

# load all csv files and merge into one large data table
path <- "data/raw/spatial" # path to your raw data, if different organization than this

files <- list.files(path = path, full.names = T)
#check files
my_data <- read_csv(files[1])
head(my_data$Genus_Species)

#Subset df to include only useful columns
my_data_subset <- my_data %>% select(Latitude, Longitude, Genus_Species, Family, Cover_Class, Observation, Number)
head(my_data_subset)

#specify reference databases to use
sources <- c("EOL", "The International Plant Names Index", "ITIS")

#show specified sources
subset(gnr_datasources(), title %in% sources)

#submit request to global names resolver
resolved_long <- my_data_subset$Genus_Species %>%
                 gnr_resolve(data_source_ids = c(3, 167), 
                             with_canonical_ranks=T)

#subset results to remove duplicates and to only include user supplied names, submitted names, matched names, scores 
resolved_short <- resolved_long %>%
                  select(user_supplied_name, submitted_name, matched_name2, score) %>%
                  distinct()

#join the resolved taxa df back to original
merged <- my_data_subset %>%
          left_join(resolved_short, by = c("Genus_Species" = "user_supplied_name")) 

#rename matched_name2 to clean taxa to explicitly indicate column with cleaned taxa names
merged <- merged %>% 
          rename(clean_taxa = matched_name2)
#remove extraneous columns
columns <- c('Latitude', 'Longitude', 'clean_taxa', 'Family', 'Cover_Class', 'Observation', 'Number')
pruned <- merged[columns]


#Remove var, ssp, subsp
pruned$clean_taxa <- str_replace_all(pruned$clean_taxa, pattern = c('var. |ssp. |subsp. '), replacement = '')
dupl_remov <-pruned %>% distinct()

#replace taxa, family NAs with unknowns to maintain data
df <- dupl_remov %>%
  tidyr::replace_na(list(clean_taxa = 'Unknown', Family = 'Unknown'))

#remove observations where both taxa and family are unknown
cleaned_data <- df[!(df$clean_taxa == 'Unknown' & df$Family =='Unknown'),]

#Rename columns to be more descriptive
data.table::setnames(cleaned_data, 'clean_taxa', 'Genus_species')
data.table::setnames(cleaned_data, 'Cover_Class', 'Cover_class')
view(cleaned_data)

###OUTSTANDING ISSUES: how to pull family names using taxize or other package? 
###                    family column was manually entered, incorrect
                        
                        
