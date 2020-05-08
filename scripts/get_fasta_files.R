# Get FASTA files
packages <- c('rentrez', 'dplyr', 'tidyverse', 'ape')
lapply(packages, library, character.only = TRUE) #load packages

#set working directory to project home folder
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/")

#set filepath to taxa list
path <- "data/clean/clean_taxa_list.csv"

#read csv
taxa_list <- read.csv(path)

#convert column to character vector
taxa_list$taxa <- map_chr(taxa_list$taxa, as.character)

#create character vector with desired genes
gene_list <- c('internal transcribed spacer', 'atpb', 'trnK', 'trnL', 'matK', 'matR', 'ndhF', 'rbcL')

#write function to create GenBank query term based on input taxon and gene
create_term <- function(taxon, gene){
  term <- paste0(taxon, "[ORGN] AND ", gene, "[WORD]")
  term
}

#initialize term_df to store terms for each gene
v <- vector(mode = "character", length = nrow(taxa_list))
term_df <- data.frame(x = v)
term_df[gene_list] <- NA # add genes as columns to df with NA as placeholders 
term_df <- term_df %>% select(-c(x)) # remove placeholder column created when initializing df

#create terms row-wise in term_df 
for (gene in colnames(term_df)) {
  term_df[[gene]] <- map2(taxa_list$taxa, gene, create_term)
}

# write function to get GenBank IDs based on desired criteria. This function only gets IDs for sequences that are
# less than 10,000 base pairs long in order for alignment to work.
get_max_ID_all_small <- function(term){
  search <- entrez_search(db = "nuccore", term = term, retmax = 25)
  if (search$count == 0){
    NA
  } else if (search$count > 0){
    max_ID = NA # variable to store ID associated with entry that has the most base pairs
    max_len_bp = 0 # number of basepairs for longest sequence
    for (ID in search$ids) {
      search_sum <- entrez_summary(db = "nuccore", 
                                   id = ID)
      if (search_sum$slen < 10000){
        if (max_len_bp < search_sum$slen){
          max_ID = ID
          max_len_bp = search_sum$slen 
        }
      }
    }
    max_ID
  }
}

#create variable df to store results of unique id
result_all_small <- taxa_list

# loop get_max_ID_all_small function columnwise across dataframe containing query terms to gather 
#unique IDs for every taxa and gene
for (gene in colnames(term_df)) {
  gene_terms <- term_df %>% select(gene)
  gene_terms_vector <- gene_terms[, gene]
  result_all_small[[gene]] <- map_chr(gene_terms_vector, get_max_ID_all_small)
}

#save IDs df result to file
write.csv(result_all_small, file = "data/genetic/ids/accession_ids_all_small.csv", row.names = FALSE)

#Retrieve FASTAs, written based on JC Santos' script found at 
#http://www.jcsantosresearch.org/Class_2014_Spring_Comparative/pdf/week_2/Jan_13_15_2015_GenBank_part_2.pdf

#read in dataframe containing IDs for sequences < 10,000 base pairs
all_small_id <- read.csv('data/genetic/ids/accession_ids_all_small.csv')

# get list of species for which no genetic data was available
taxa_no_genes <- all_small_id %>% 
  filter(is.na(internal.transcribed.spacer)) %>%
  filter(is.na(atpb)) %>%
  filter(is.na(trnK)) %>%
  filter(is.na(trnL)) %>%
  filter(is.na(matK)) %>%
  filter(is.na(matR)) %>%
  filter(is.na(ndhF)) %>%
  filter(is.na(rbcL)) %>%
  select(taxa) %>%
  dplyr::rename(taxa_no_genetic_data = taxa)

#save df as CSV to know priority species for field surveys
write.csv(taxa_no_genes, file = "data/clean/taxa_missing_genbank_data.csv")

# remove taxa column for function to work
all_small_no_taxa <- all_small_id %>% select(-c('taxa')) 

#function to write fasta files for individual genes
write_fasta_func_all_small <- function(id_df){
  for (gene in colnames(id_df)){ # loop through each column (gene) in the dataframe 
    gene_id_vector <- id_df[, gene] # isolate unqiue ID vector for a single gene 
    gene_ids_no_na <- gene_id_vector[!is.na(gene_id_vector)]
    sequences <- read.GenBank(gene_ids_no_na) # create a DNABin object for the gene, containing all sequences and their respective accession ID's
    sequences_GenBank_IDs <- paste(attr(sequences, "species"), names(sequences), sep = paste0(" | ", gene, "_")) #build a character vector with the species, GenBank accession numbers, and gene name to create more informative labels associated with each sequence in the FASTA file
    new_names <- updateLabel(sequences, old = names(sequences), new = sequences_GenBank_IDs) #replace accession IDs with the informative labels created above
    write.FASTA(x = new_names, file = paste0("data/genetic/sequences/all_small/", gene, "_all_small.fasta")) # write FASTA file to desired directory
  }
}

#write fastas
write_fasta_func_all_small(all_small_no_taxa)


################ Walkthrough of functions in 'rentrez' package. Include b/c informative or remove b/c unnecessary? 

#see list of all NCBI databases
entrez_dbs()
entrez_db_summary('nuccore') # get info about database of interest

#set up NCBI search specific to database of interest
r_search <- entrez_search(db = 'nuccore', term = 'R Language')
r_search$ids

#get list of query options for database 
entrez_db_searchable("nuccore")

sage_rbcl <- "Ambrosia dumosa[ORGN]"
sage_rbcl_search <- entrez_search(db="nuccore",
              term= sage_rbcl,
              retmax=25)
sage_rbcl_search

sage_rbcl_search$ids

#create summary for id to determine basis for selection
sage_sum <- entrez_summary(db = "nuccore", id = sage_rbcl_search$ids[1])
sage_sum$slen

#iterate through the various sequences available by sequence length (# of base pairs) to determine which one to use


#obtain fasta files from accession IDs
sage_rbcl_seqs <- entrez_fetch(db = "nuccore", id = sage_rbcl_search$ids[1], rettype = "fasta")
sage_rbcl_seqs





