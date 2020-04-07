#clean_species_names script, based on Viktoria Wagner's blog post on how to use the taxize package.
#automate cleaning of erroneous taxon names by comparing taxa to reference databases.

#load libraries
library(dplyr)
library(taxize)
library(magrittr)
library(tidyverse)
library(sf)

# change to your project HOME directory
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/") 

# set filepath
path <- "data/semiclean/spatial_data.csv" # path to your data

#read in data
my_data <- read.csv(file = "data/semiclean/spatial_data.csv")

#replace '#N/A' with NA across dataframe
my_data <- as.data.frame(sapply(my_data, function(x) {
  gsub('#N/A', NA, x)
}))

#split data into 2 dataframes in order to move observations from the columns genus, 
#species to the column Genus_Species

na_gensp <- my_data %>%
  filter(is.na(Genus_Species))

not_na_gensp <- my_data %>%
  filter(!is.na(Genus_Species))

#move all values from separate genus, species columns to the Genus_Species column 
#to prepare for taxa name spellchecking/cleaning

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
cleaned_sp_names <- dupl_remov[!(dupl_remov$Genus_species == 'Unknown' & dupl_remov$Family =='Unknown'),]

#IDIOSYNCRATIC STEP, UNECCESSARY TO REPLICATE. replace all occurrences of 'Cylindropuntia bigelovii' with 
#'Cylindropuntia echinocarpa' and also 'Cylindropuntia californica' to 'Cylindropuntia acanthocarpa' 
#'to rectify plant misidentification from data collection process
cleaned_sp_names$Genus_species <- sapply(cleaned_sp_names$Genus_species, function(taxon) {
  gsub("Cylindropuntia bigelovii", "Cylindropuntia echinocarpa", taxon)
})

cleaned_sp_names$Genus_species <- sapply(cleaned_sp_names$Genus_species, function(taxon) {
  gsub("Cylindropuntia californica", "Cylindropuntia acanthocarpa", taxon)
})

#write cleaned_species_names to .csv
write.csv(cleaned_sp_names, file = "data/semiclean/clean_species_names.csv", row.names = FALSE)

###########STILL NEED TO REMOVE ALL SINGLE WORD TAXA FROM DATA (EXCEPT PECTOCARYA & UNKNOWN)


#extract unique taxa names for use in Genbank query
taxa_list <- cleaned_data$Genus_species %>% sort() %>% unique()
taxa_list_df <- data.frame(taxa = taxa_list)
#write taxa list to CSV, .txt files
write.csv(taxa_list_df, file = "data/clean/taxa_list.csv", row.names = FALSE) # change path relevant to your directory organization
write.table(taxa_list_df, file = "data/clean/taxa_list.txt", sep = " ", col.names = FALSE)


