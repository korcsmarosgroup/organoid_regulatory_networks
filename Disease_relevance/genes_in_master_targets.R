# Script to see how many of an input gene list (UC SNP genes, CD SNP genes or UC DEGs) are regulated by the predicted master regulators.
# Hypergeometric significance testing calculates the significance of getting this many genes of interest in the target list of the master regulators
#
# Input genes_f: Tab delimited text files with the genes of interest. UC DEGs have genes as human symbols in the column 'gene' with goblet cell type
#                indicated by 'Goblet' in column 'ident'. The SNP genes have human ensembl ids in the column 'EnsemblID'.
# Input net_f: Goblet or Paneth cell network created previously with columns 'source' and 'target'. Choose which one depending on what you want to test.
# Input convert_f: Conversion file for mouse to human IDs from Inparanoid. Gene symbols in column 'Hum.gene.syb' and mouse ensembl IDs in column 'Mouse.Ensembl'.
# Input background_f: Unfiltered background network with all possibel mouse interactions - from which the goblet network was created.
# Input regs: Vector of mouse ensembl IDs for the predicted master regulators associated with goblet ro Paneth cells. Select which to use as appropriate.
# Input out_text: Text string to be added to output files to say which network/cell type/input list was analysed. 
#
# Output 1: Tab delimited text file which gives in sentences, the results from the overlap of genes of interest and the network considering the master regulators.
# Output 2: List of all genes of interest targeted by at least one of the relevant cell predicted master regualtors in the relevant network.

##### Setup #####

library(dplyr)
library(tidyr)

# Input files

# IBD genes from Smillie et al and from Jostins et al and Farh et al. Select just one to test
#genes_f <- "../Input_data/IBD_genes/goblet_inflammed_v_healthy_smillie_tableS4.txt"
genes_f <- "../Input_data/IBD_genes/All_UC_snps.txt"
#genes_f <- "../Input_data/IBD_genes/All_CD_snps.txt"

# GceE and PCeE networks - select just the one you want to test
#net_f <- "../Input_data/Pan_network_lfc1.txt"
net_f <- "../Input_data/Gob_network_lfc1.txt"

# Human to mouse ID conversion table
convert_f <- "../Input_data/InParanoid-Mus-homo-UniprotEnsembl-May2018-expanded.txt"

# Whole unfiltered network
background_f <- "../Input_data/Known_interactions_unfiltered.txt"

# Predicted master regulators - select just one you want to test
# goblet
regs <- c("ENSMUSG00000026815", "ENSMUSG00000022346", "ENSMUSG00000032035", "ENSMUSG00000024431", "ENSMUSG00000022479")
# paneth
#regs <- c("ENSMUSG00000034957", "ENSMUSG00000052684", "ENSMUSG00000020889", "ENSMUSG00000015846", "ENSMUSG00000032035", "ENSMUSG00000024431", "ENSMUSG00000022479")

# Name of input list and cell type for use in output filenames
out_text_name <- "UC_degs_goblet"

##### Load files #####

# Load files
genes <- read.csv(genes_f, sep = "\t", header=T)
net <- read.csv(net_f, sep = "\t", header = T)
convert <- read.csv(convert_f, sep = "\t", header = T)
background <- read.csv(background_f, sep = "\t", header=T)

# Empty list to hold the output text
out_text <- c()

##### Preprocessing #####

# Expand the conversion file for one human gene symbol per line
convert_1 <- separate_rows(convert, Hum.gene.syb, sep="; ")

# Get unique list of human input genes - different for degs and snps files
if ("ident" %in% colnames(genes)) {
  
  genes <- left_join(genes, convert, by = c("gene" = "Hum.gene.syb"))
  genes1 <- genes %>% filter(Mouse.Ensembl != "-") %>% 
    filter(ident == "Goblet") %>% select(Mouse.Ensembl) %>% 
    na.omit() %>% unique()
  
} else {
  
  genes1 <- genes %>% select(EnsemblID) %>% separate_rows(EnsemblID) %>% unique() %>% na.omit()
  genes1 <- left_join(genes1, convert, by = c("EnsemblID" = "Human.Ensembl"))
  genes1 <- genes1 %>% select(Mouse.Ensembl) %>% filter(Mouse.Ensembl != "") %>% na.omit() %>% unique()
  
}

# Remove genes which are not in the background network - targets specifically (so we can do hypergeomtric significance testing)
genes_fi <- genes1 %>% filter(Mouse.Ensembl %in% background$target)

# Write out the result
out_text <- append(out_text, paste0("There are ", (nrow(genes) - nrow(genes_fi)), " genes of interest which are not in the background network targets."))
out_text <- append(out_text, paste0("After removing these genes there are ", nrow(genes_fi), " genes of interest."))

##### Regulators of the DEGs #####

# First extract all the regulators of the deg list
gene_regs <- net %>% filter(target %in% genes_fi$Mouse.Ensembl) %>% select(source) %>% unique()

# Filter for targets of the master reg
master <- gene_regs %>% filter(source %in% regs)

out_text <- append(out_text, paste0("There are ", as.character(nrow(gene_regs)), " regulators in the network which target the ", as.character(nrow(genes)), " genes of interest."))
out_text <- append(out_text, paste0("There are ", nrow(master), " master regulators targetting the genes of interest."))

##### Targets of the master regulators #####

# Extract all targets of the master regulators
master_tars <- net %>% filter(source %in% regs) %>% select(target) %>% unique()
out_text <- append(out_text, paste0("There are ", nrow(master_tars), " targets of the master regs in the network."))

# How many of the input genes are in this target list?
master_tars_filt <- master_tars %>% filter(target %in% genes_fi$Mouse.Ensembl)
out_text <- append(out_text, paste0("There are ", nrow(master_tars), " genes of interest targeted by at least one of the master regulators in the network."))

##### Significance testing #####
# Significance of getting this many uc/cd/drug target associated genes in the target list of the relevant master regulators

# For the background I use all the targets in the background network
backg_targs <- background %>% select(target) %>% unique()

# Hypergeometric signi testing
sig <- 1 - phyper(nrow(master_tars_filt), nrow(master_tars), (nrow(backg_targs) - nrow(master_tars)), nrow(genes_fi))
out_text <- append(out_text, paste0(sig, " = The significance of having ", nrow(master_tars_filt), " genes of interest in the list of ", nrow(master_tars), " targets of the master regulators (in the network), given ", nrow(backg_targs), " targets in the background network."))


##### Output #####

# Out file path
out1 <- paste0(out_text_name, "_master_reg_targets_results.txt")
out2 <- paste0(out_text_name, "_master_reg_targets_list.txt")

# Output descriptive results
write.table(out_text, out1, row.names = F, quote = F, sep = "\t")

# Output list of genes of interest targetted by the predicted master regulators in the network
write.table(master_tars_filt, out2, row.names = F, quote = F, sep = "\t")

out_text_f2 <- "_targeted_by_1ormore_master_regs.txt"