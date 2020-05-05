#Align multiple sequences script using DECIPHER
library('DECIPHER')

#set working directory to project home folder
setwd("/Users/deanberkowitz/Documents/mishler_lab/thesis/mnp_spatial_phylo/data/genetic")


one_big_files <- list.files(path = "sequences/one_big")
all_small_files <- list.files(path = "sequences/all_small")
one_big_files %>% class()

#test atpb alignment on big
gene_path <- paste0("sequences/one_big/", one_big_files[1])
dna <- readDNAStringSet(gene_path) #DNAStringSet object of unaligned sequences
DNA <- AlignSeqs(dna) # align the sequences directly without translation
writeXStringSet(DNA, filepath = paste0("alignments/", str_sub(gene_path, start = 11, end = -7), "_aligned_test", ".fasta"))

#test aligning on all small files
for (unaligned_gene in all_small_files){
  gene_path <- paste0("sequences/all_small/", unaligned_gene)
  dna <- readDNAStringSet(gene_path) #DNAStringSet object of unaligned sequences
  DNA <- AlignSeqs(dna) # align the sequences directly without translation
  writeXStringSet(DNA, filepath = paste0("alignments/", str_sub(gene_path, start = 11, end = -7), "_aligned_test", ".fasta"))
}



###################3
matK <- "data/genetic/fasta/matK.fasta"
dna <- readDNAStringSet(matK) #DNAStringSet object of unaligned sequences
DNA <- AlignSeqs(dna) # align the sequences directly without translation
writeXStringSet(DNA, file="data/genetic/fasta/matK_aligned.fasta")

#matk works without chloroplasts, now try ITS since its harder to align

its <- "data/genetic/fasta/internal\ transcribed\ spacer.fasta"
its_dna <- readDNAStringSet(its)
its_DNA <- AlignSeqs(its_dna) # align the sequences directly without translation
writeXStringSet(its_DNA, file="data/genetic/fasta/its_aligned.fasta")
