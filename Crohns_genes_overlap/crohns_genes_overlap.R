# Script to get the number of Crohns associated genes in the differentially expressed gene lists, the DEG networks and the
# regulator-marker networks.
#
# Input: 1. All Crohn's associated SNPs and genes in a tab delimited text file with headers. Ensembl gene id (human) in column 'EnsemblID', mulitple IDs
#           seperated by ','.
#        2. ID conversion file in a tab delimited text format with headers. Only 1 ID per line for the input ID (long format).
#           The ensembl ID headers should be 'Human Ensembl' and 'Mouse Ensembl'.
#           The uniprot ID headers should be 'Human Uniprot' and 'Mouse Uniprot'.
#        3. Background network of all unfiltered interactions from external sources. Tab delimited text file with columns 'source'
#           and 'target' containing the source an dtarget nodes.
#        4. Differentially expressed genes tables for Paneth and goblet datasets. Tab delimited text file
#           (with headers) where each row is a differentially expressed gene (mouse). Ensembl id in column 'target_id'.
#        5. DEG networks in tab delimited text files. One for Paneth and one for goblet datasets.
#           The 'source' and 'target' column headers contain the source and target node ids respectively.
#        6. Regulator-marker networks in tab delimited text files. One for Paneth and one for goblet datasets.
#           The 'source' and 'target' column headers contain the source and target node ids respectively.
# Output: 1. Tab delimited text files saved with list of CD genes in the DEGs files and the networks files.
#         2. Tab delimited text files saved with the same format as the input files (DEGs and networks files) but filtered for the CD genes only. 
#           Due to the small number of CD genes these filtered networks are liekely to be empty as every interaction must be between 2 CD genes.
#         3. The script prints the results of the hypergeometric significance test on the Paneth and goblet networks.


######## Packages #########

library(reshape2)
library(dplyr)
library(tidyr)

####### File paths #########

# Crohn's SNP associated genes table
crohns_f <- "../Input_data/All_CD_snps.txt"

# ID conversion file
id_conversion_f <- "../Input_data/Inparanoid-Mus-homo-UniprotEnsembl-May2018-expanded.txt"

# Background network with which to filter the CD genes (to ensure same background when doing hypergeometric distribtuion test)
background_f <- "../Input_data/Known_interactions_unfiltered.txt"

# Differentially expressed gene tables
pan_degs <- "../Input_data/Differentially_expressed_genes/Pan_DEGs_lfc1.txt"
gob_degs <-  "../Input_data/Differentially_expressed_genes/Gob_DEGs_lfc1.txt"

# DEG networks
pan_deg_n <- "../Input_data/Pan_network_lfc1.txt"
gob_deg_n <-  "../Input_data/Gob_network_lfc1.txt"

# Regulator-marker networks
pan_deg_m <- "../Input_data/Regulator_marker_networks/Pan_network_pan_markers_all.txt"
gob_deg_m <-  "../Input_data/Regulator_marker_networks/Gob_network_gob_markers_all.txt"

inputs <- c(pan_degs,gob_degs,pan_deg_n,gob_deg_n,gob_deg_m,pan_deg_m)

##### Get background network genes #######

background_fun <- function(background_file){
  ## Description: Extracts a unique list of nodes/genes in the network
  ## Input: background_file - filepath to the background network file (see script header for description)
  ## Output: nodes - Dataframe with one column (header- 'nodes') containing unique 'list' of nodes/genes
  ##        in the background network file.
  
  # Open file as table
  bgd <- read.csv(background_file, sep = "\t")
  
  # Get unique list of all source and target nodes - dataframe with 1 column, header 'nodes'
  bgd_1 <- bgd %>% select(nodes = source) %>% unique()
  bgd_2 <- bgd %>% select(nodes = target) %>% unique()
  nodes <- rbind(bgd_1, bgd_2)
  nodes <- nodes %>% unique()
  
  return(nodes)
}

##### Get CD genes #######

cd_gene_fun <- function(cd_file, id_file, background_nodes){
    ## Description: Converts CD gene table to a unique list of mouse IDs and filters for those in background network
    ## Input: cd_file - filepath to CD table (see script header for description)
    ## Input: id_file - filepath to id conversion table (see script header for description)
    ## Input: background_nodes - Dataframe with one column (header- 'nodes') containing unique 'list' of nodes/genes in the background network file
    ## Output: cd_mus - dataframe with one column (header- 'Mouse.Ensembl') containing unique 'list' of CD nodes/genes in background network

    # Open CD file
    cd <- read.csv(cd_file, sep = "\t")

    # Get unique list of human CD genes
    cd <- cd %>% select(c(EnsemblID)) %>% mutate(EnsemblID = strsplit(as.character(EnsemblID), ",")) %>% 
      unnest(EnsemblID) %>%
      unique() %>% 
      filter(EnsemblID != "-")

    # Open ID file
    ids <- read.csv(id_file, sep = "\t")
    ids_filt <- ids %>% select(c(Human.Ensembl, Mouse.Ensembl))
    
    # Convert CD genes
    cd_mus <- left_join(cd, ids_filt, by = c("EnsemblID" = "Human.Ensembl")) %>% select(c(Mouse.Ensembl)) %>%
      unique() %>% na.omit()
    
    # Filter for background nodes only
    cd_mus <- cd_mus %>% filter(Mouse.Ensembl %in% background_nodes$nodes)
  
    return(cd_mus)
}


##### Overlap CD genes with input lists #######

overlap_fun <- function(cd_list, input_data){
  ## Description: Overlaps an input dataset with cd list
  ## Input: cd_list - Dataframe with one column containing unique mouse Ensembl CD gene ids
  ## Input: input_data - filepath to input data, could be network or DEGs table.
  ## Output: overlap - Dataframe with one column containinf unique mouse ensembl CD gene ids which are also in the input_data.
  ## Output: in_filt - Filtered version of the input data where all nodes/degs are CD genes. For the networks, these will be networks, 
  ##                  but they are likely to be empty due to no interactions between 2 cd genes.
  
  # Open input data
  in_data <- read.csv(input_data, sep = "\t")
  
  # Check if we have a DEGs file or a network
  # Get unique list of all genes/nodes
  if (names(in_data)[1] == "target_id"){
    genes <- in_data %>% select(target_id) %>% unique()
  } else {
    genes_1 <- in_data %>% select(target_id = source) %>% unique()
    genes_2 <- in_data %>% select(target_id = target) %>% unique()
    genes <- rbind(genes_1, genes_2)
    genes <- genes %>% unique()
  }
  
  # Get CD genes in the gene/node lists
  overlap <- cd_list %>% filter(Mouse.Ensembl %in% genes$target_id)
  
  # Filter the original input for these genes
  if (names(in_data)[1] == "target_id"){
    in_filt <- in_data %>% filter(target_id %in% overlap$Mouse.Ensembl)
  } else {
    in_filt <- in_data %>% filter((source %in% overlap$Mouse.Ensembl)&(target %in% overlap$Mouse.Ensembl))
  }
  
  # List the outputs
  out_list <- list(overlap, in_filt)
  
  return(out_list)
}

##### Hypergeometric test #######

hypergeo_fun <- function(num_pan, num_gob){
  ## Description: Calculate hypergeometric distribution statistic for enrichment of CD genes
  ##              in the Paneth and goblet networks.
  ## Input: num_pan/gob - number of cd genes in the Paneth and goblet networks calculated previously in the script.
  ## Output: list_out - list of the paneth and goblet significance test statistics
  ## Warning: Must run all previous code for this function to get it's inputs
  
  
  # Number of unique nodes in the background network
  back_net <- nrow(bgd_nodes)
  #  Number of unique nodes in the Paneth network (calculated previously)
  pan_net <- 3231
  #  Number of unique nodes in the goblet network (calculated previously)
  gob_net <- 2219
  # Number of CD genes in the background network
  num_cd <- nrow(cd_mus)
  # Number of mouse CD genes in background networks = num_pan, num_gob

  # Paneth network hypergeometric sig test
  pan_sig <- 1-phyper(num_pan,num_cd,back_net-num_cd,pan_net)
  
  # Paneth network hypergeometric sig test
  gob_sig <- 1-phyper(num_gob,num_cd,back_net-num_cd,gob_net)
  
  list_out <- list(pan_sig, gob_sig)
  
  return(list_out)
}

##### Call functions #######

# Get list of background network nodes/genes
bgd_nodes <- background_fun(background_f)

# Get list of mouse cd genes
cd_mus <- cd_gene_fun(crohns_f,id_conversion_f,bgd_nodes)

# Iterate the deg lists/networks
for (item in inputs){
  
  # Overlap cd genes and network/deg list
  out_list <- overlap_fun(cd_mus, item)
  filtered_cd <- out_list[[1]]
  filtered_in <- out_list[[2]]
  
  # Get output filename
  input_filename <- basename(item)
  output_filename1 <- paste("cd_genes_",input_filename, sep="")
  output_filename2 <- paste("list_cd_genes_",input_filename, sep="")
  
  # Get number of CD genes in Pan and gob networks for hypergeo test
  if (output_filename2 == "list_cd_genes_Pan_network_lfc1.txt") {
    num_pan = nrow(filtered_cd)
  } 
  if (output_filename2 == "list_cd_genes_Gob_network_lfc1.txt") {
    num_gob = nrow(filtered_cd)
  }
  
  # Save output
  write.table(filtered_in, output_filename1, sep = "\t", quote = F, row.names =F)
  write.table(filtered_cd, output_filename2, sep = "\t", quote = F, row.names =F)
}

# Call hypergeometric significance test
significance <- hypergeo_fun(num_pan, num_gob)

print(paste("Paneth significance = ", significance[1]))
print(paste("Goblet significance = ", significance[2]))
print(paste("Number of CD genes in background network = ",nrow(cd_mus)))
print(paste("Number of CD genes genes in Paneth network = ",num_pan))
print(paste("Number of CD genes genes in Paneth network = ",num_gob))
