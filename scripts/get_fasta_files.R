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

sage_rbcl_search$

#create summary for id to determine basis for selection
sage_sum <- entrez_summary(db = "nuccore", id = sage_rbcl_search$ids)
sage_sum$`952002070`$slen

#iterate through the various sequences available by sequence length (# of base pairs) to determine which one to use
for (i in sage_sum){
  print(paste("ID", i$uid, "Length", i$slen))}

#obtain fasta files from accession IDs
sage_rbcl_seqs <- entrez_fetch(db = "nuccore", id = sage_rbcl_search$ids, rettype = "fasta")
sage_rbcl_seqs
