# Get FASTA files
library(rentrez)
library(dplyr)
library(tidyverse)

#set working directory to project home folder
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/")

#set filepath to taxa list
path <- "data/clean/taxa_list.csv"

#read csv
taxa_list <- read.csv(path)

#subset csv to remove extra index column, rename taxa column
taxa_list <- taxa_list %>% subset(select = x) %>%
  setNames('taxa')
taxa_list$taxa <- lapply(taxa_list$taxa, as.character)

#create character vector with desired genes
gene_list <- c('internal transcribed spacer', 'atpb', 'trnK', 'trnL', 'matK', 'matR', 'ndhF', 'rbcL')

# add genes as columns to df with NA as placeholders before adding Accession IDs from GenBank
taxa_list[gene_list] <- NA

#write function to create GenBank query term based on input taxon and gene
create_term <- function(taxon, gene){
  term <- paste0(taxon, "[ORGN] AND ", gene, "[WORD]")
  term
}

#initialize term_df to store terms for each gene
v <- vector(mode = "character", length = length(taxa_list$taxa))
term_df <- data.frame(x = v)
term_df[gene_list] <- NA
term_df <- term_df[-x]

#create terms row-wise in term_df 
for (name in colnames(term_df)) {
  term_df[[name]] <- map2(taxa_list$taxa, name, create_term)
}

#define get_max_id function to query GenBank and pull the ID of the entry with the most basepairs
get_max_ID <- function(term){
  search <- entrez_search(db = "nuccore", term = term, retmax = 20)
  if (search$count == 0){
    print(NA)
  } else if (search$count == 1){
    search$ids
  } else if (search$count > 1){
    max_ID = NA # variable to store ID associated with entry that has the most base pairs
    max_len_bp = 0 # number of basepairs for longest sequence
    for (ID in search$ids) {
      search_sum <- entrez_summary(db = "nuccore", 
                                   id = ID)
      if (max_len < search_sum$slen){
        max_ID = ID
        max_len_bp = search_sum$slen 
      }
    }
    print(max_ID)
  }
}

#create test term list
test_terms_rbcL_df <- head(term_df, 20) %>% select('rbcL')
test_terms_rbcL_lst <- test_terms_rbcl_df$rbcL

accession_ids <- taxa_list %>% select(-c(taxa))
accession_ids

#######LEFT OFF HERE. NEED TO BUILD ACCESSION ID DATAFRAME BY QUERYING GENBANK FROM TERMS COLUMN-WISE 

for (gene in colnames(accession_ids)) {
  
}

##################GABE demonstrate package
search <- entrez_search(db = "nuccore", term = term, retmax = 20)
search$count
search_ids <- search$ids
search_sum <- entrez_summary(db = "nuccore", 
                             id = search_ids[1])
search_sum$slen

for (ID in search_ids){
  each_try <- entrez_summary(db = "nuccore",
                id = ID)
  id
}

################

#see list of all NCBI databases
entrez_dbs()
entrez_db_summary('nuccore') # get info about database of interest

#set up NCBI search specific to database of interest
r_search <- entrez_search(db = 'nuccore', term = 'R Language')
r_search$ids

#get list of query options for database 
entrez_db_searchable("nuccore")

sage_rbcl <- "Salvia dorrii[ORGN] AND rbcl[GENE]"
sage_rbcl_search <- entrez_search(db="nuccore",
              term= sage_rbcl,
              retmax=10)
sage_rbcl_search

sage_rbcl_search$ids[1]

#create summary for id to determine basis for selection
sage_sum <- entrez_summary(db = "nuccore", id = sage_rbcl_search$ids[1])
sage_sum$slen

#iterate through the various sequences available by sequence length (# of base pairs) to determine which one to use


#obtain fasta files from accession IDs
sage_rbcl_seqs <- entrez_fetch(db = "nuccore", id = sage_rbcl_search$ids[1], rettype = "fasta")
sage_rbcl_seqs





