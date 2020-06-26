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

##repeat to move Observation # to Observation column

obs_num <- bound %>%
  filter(is.na(Observation))

obs <- bound %>%
  filter(!is.na(Observation))


#move Observation.. to Observation column
obs_num$Observation <- paste(obs_num$Observation..)

#bind two observation df together
bound2 <- rbind(obs, obs_num)

#Subset df to include only useful columns
my_data_subset <- bound2 %>% 
                  select(Latitude, Longitude, 
                         Genus_Species, Family, Cover_Class, 
                         Number, Observation)

#test Remove unknown from gen_sp to try and retain information on type of unknown
my_data_subset$Genus_Species <- sapply(my_data_subset$Genus_Species, function(unknown_taxon) {
  gsub('Unknown', '', unknown_taxon)
})

#remove underscores
my_data_subset$Genus_Species <- sapply(my_data_subset$Genus_Species, function(underscore) {
  gsub('_', ' ', underscore)
})

# keep important info (grass, lichen, moss)
#DO THIS LATER -db 05052020
  
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
merged <- dplyr::rename(merged,
                 clean_taxa = matched_name2)
#remove original (erroneous) Genus_Species column, extraneous columns
pruned <- subset(merged, select = -c(submitted_name, Genus_Species))

#rename clean_taxa to Genus_species
data.table::setnames(pruned, 'clean_taxa', 'Genus_species')

#remove duplicate rows
dupl_remov <- pruned %>% 
              distinct()

#replace taxa, family NAs with unknowns to maintain data
dupl_remov$Family <- as.character(dupl_remov$Family) # convert Family column to character
dupl_remov$Observation <- as.character(dupl_remov$Observation) # convert Observation column to character

dupl_remov <- dupl_remov %>%
  tidyr::replace_na(list(Genus_species = 'Unknown', Family = 'Unknown'))
dupl_remov <- dupl_remov %>%
  tidyr::replace_na(list(Observation = 'Unknown'))

#replace all occurrences of 'Cylindropuntia bigelovii' with 
#'Cylindropuntia echinocarpa' and also 'Cylindropuntia californica' to 'Cylindropuntia acanthocarpa' 
#'to rectify plant misidentification from data collection process

#Opuntia bigelovii has been renamed to Cylindropuntia bigelovii
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Opuntia bigelovii", "Cylindropuntia echinocarpa", taxon)
})

#Opuntia ramosissima has been renamed to Cylindropuntia ramosissima
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Opuntia ramosissima", "Cylindropuntia ramosissima", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Cylindropuntia bigelovii", "Cylindropuntia echinocarpa", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Cylindropuntia californica", "Cylindropuntia acanthocarpa", taxon)
})

#### fix typos and synonyms not caught by GNR step
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Cylindropuntia ramossissima", "Cylindropuntia ramosissima", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Lepidium aschersonii", "Lepidium andersonii", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Echinocereus triglochidiatus.+", "Echinocereus triglochidiatus", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Cleomella arborea", "Peritoma arborea", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Camissonia boothii", "Eremothera boothii", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Castilleja angustifolia", "Castilleja chromosa", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Psorodendron californicum", "Phoradendron californicum", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Chamaesyce albomarginata", "Euphorbia albomarginata", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Cryptandra angustifolia", "Cryptantha angustifolia", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Cryptanthus angustifolius", "Cryptantha angustifolia", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Dryophytes", "Bryophyta", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Larra tridentata", "Larrea tridentata", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Anemopsis$", "Anemopsis californica", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Baccharis sarothroides", "Baccharis sergiloides", taxon)
})


#Lepidium andersonii isn't a real species
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Lepidium andersonii", "Unknown", taxon)
})

#Removing this species as CalFlora only lists one occurrence which is outside of the Mojave desert
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Dalea greggii", "Unknown", taxon)
})

#Removing this species as CalFlora lists occurrences which are in Northern California
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Hieracium nudicaule", "Unknown", taxon)
})

#Swapping out for Solidago Juliae for Ericameria paniculata as it was probably a mis-identification.
#Solidago juliae isn't listed on CalFlora and therefore is assumed to not occur in California.
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Solidago juliae", "Ericameria paniculata", taxon)
})

#replace with GenBank synonym
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Dasyochloa pulchella", "Munroa pulchella", taxon)
})

#replace with GenBank synonym
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Peritoma arborea", "Cleomella arborea", taxon)
})

#Replacing Cryptantha angustifolia with Jonhstonella angustifolia, a synonym that is on GenBank
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Cryptantha angustifolia", "Johnstonella angustifolia", taxon)
})

#Replacing Acmispon maritimus with an older synonym, Lotus strigosus, which is on GenBank
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Acmispon maritimus", "Lotus strigosus", taxon)
})

##Remove all subspecies, varities #LATER REMOVE THIS STEP TO PRESERVE SSP INFO
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Bromus madritensis ssp. rubens", "Bromus madritensis", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Astragalus nuttallianus var. imperfectus", "Astragalus nuttallianus", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Senecio flaccidus var. monoensis", "Senecio flaccidus", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Eucrypta chrysanthemifolia var. bipinnatifida", "Eucrypta chrysanthemifolia", taxon)
})

####### Change this if time: Want to retain observations that have Family name in the 
#######Genus_species column by changing the Genus_species to Unknown and Moving the Family over to Family
####### Also: want to write a function that, for a Genus_species, if it has a Family that is not Unknown,
####### Apply that family to all Genus_species of that kind. To have higher level taxonomic data !!

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Didymodon", "Unknown", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Lichens", "Unknown", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Bryoria", "Unknown", taxon)
})

dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Bryophyta", "Unknown", taxon)
})
###########COLLAPSE PECTOCARYA sp. INTO PECTOCARYA?
dupl_remov$Genus_species <- sapply(dupl_remov$Genus_species, function(taxon) {
  gsub("Pectocarya$", "Unknown", taxon)
})

#step to visually check the cleaned taxa list in a dataframe to make sure there are no errors
check_taxa_list <- dupl_remov$Genus_species %>% sort() %>% unique()
check_taxa_list_df <- data.frame(taxa_check = check_taxa_list)

#rename to indicate end of clean species script
clean_species_names <- dupl_remov

#write cleaned_species_names to .csv
write.csv(clean_species_names, file = "data/semiclean/clean_species_names.csv", row.names = FALSE)

