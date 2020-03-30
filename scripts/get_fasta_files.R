# Get FASTA files
library(rentrez)
library(dplyr)

#set working directory to project home folder
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/")

#set filepath to taxa list
path <- "data/clean/taxa_list.csv"
#read csv
taxa_list <- read.csv(path)
#subset csv to remove extra index column, rename taxa column
taxa_list <- taxa_list %>% subset(select = x) %>%
  setNames('taxa')
taxa_list
#create character vector with desired genes
gene_list <- c('ITS', 'atpb', 'trnK', 'trnL', 'matK', 'matR', 'ndhF', 'rbcL')
#add genes as columns to df with NA as placeholders before adding Accession IDs from GenBank
taxa_list[gene_list] <- NA

taxa_six <- taxa_list %>% head()
taxa_six[,1]

gene_list[1]
taxa_six[,1][1]

agg_func <- function(df){}
  
  taxa <- taxa_six[,1]
  genes <- c('internal transcribed spacer', 'atpb', 'trnK', 'trnL', 'matK', 'matR', 'ndhF', 'rbcL')
?entrez_search()
term <- paste0(taxa[1], "[ORGN] AND ", genes[1], "[WORD]")
term

search <- entrez_search(db = "nuccore", term = term, retmax = 20)
search_ids <- search$ids
search_sum <- entrez_summary(db = "nuccore", 
                             id = search_ids[1])
id <- c()
slen <- c()

output <- vector("integer", length = nrow(taxa_six))
output
for (i in taxa_six[,1]){
  output[[i]] <- entrez_search(db)
}
taxa_six
nrow(taxa)

for (ID in search_ids){
  each_try <- entrez_summary(db = "nuccore",
                id = ID)
  id
}

lst1 <- list()

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

########33
test_taxa <- list(c("Acmispon maritimus", "Lotus strigosus"))

for (taxa in test_taxa){
  for (gene in gene_list){
    paste0(taxa, "[Organism] AND ", gene, "[Gene]")
  }
}
term




