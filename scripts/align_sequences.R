#Align multiple sequences script using DECIPHER
library('DECIPHER')

#set working directory to project home folder
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/data/genetic")

all_small_files <- list.files(path = "sequences/all_small")

#align all files
for (unaligned_gene in all_small_files){
  gene_path <- paste0("sequences/all_small/", unaligned_gene)
  dna <- readDNAStringSet(gene_path) #DNAStringSet object of unaligned sequences
  DNA <- AlignSeqs(dna) # align the sequences directly without translation
  writeXStringSet(DNA, filepath = paste0("alignments/", str_sub(gene_path, start = 11, end = -7), "_aligned_0506", ".fasta"))
}

