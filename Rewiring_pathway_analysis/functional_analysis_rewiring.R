# Script to carry out functional analysis on the top most rewired regulators.
#
# Input: 1. List of most rewired regulators.
#        2. DEG networks in tab delimited text files. One for Paneth and one for goblet datasets.
#           The 'source' and 'target' column headers contain the source and target node ids respectively.
#        3. ID conversion file in a tab delimited text format with headers. Only 1 ID per line for the input ID (long format).
#           The ensembl ID headers should be 'Human Ensembl' and 'Mouse Ensembl'.
#           The uniprot ID headers should be 'Human Uniprot' and 'Mouse Uniprot'.
# Output: 1. One tab delimited text file with the results for the Reactome, KEGG and GO enrichment analysis.

######## Packages #########

library(dplyr)
library(ReactomePA)
library(clusterProfiler)

########### Input and filter #########

# Most rewired regulators
regulators <- c("ENSMUSG00000017724","mmu-let-7e-5p","mmu-miR-152-3p","ENSMUSG00000019982","ENSMUSG00000032238")

# Paneth and gobet DEG networks
pan_net <- read.table("Pan_network_lfc1.txt")
gob_net <- read.table("Gob_network_lfc1.txt")

# ID conversion table
human <- read.delim("InParanoid-Mus-homo-UniprotEnsembl-May2018-expanded.txt")

all2 <- data.frame(File=character(), 
                   User=character(), 
                   stringsAsFactors=FALSE) 
  
# Filter networks for only the regulators of interest and their downstream interactors
pan_net2 <- subset(pan_net, source %in% regulators)
gob_net2 <- subset(gob_net, source %in% regulators)
rm(pan_net, gob_net)

########### Get overlap lists #########

for(n in regulators){
  print(n)
  #Filter network for specific regulator
  pan_net_filt <- subset(pan_net2, source == n)
  gob_net_filt <- subset(gob_net2, source == n)
  shared <- intersect(pan_net_filt$target, gob_net_filt$target)
  pan_only <- setdiff(pan_net_filt$target, gob_net_filt$target)
  gob_only <- setdiff(gob_net_filt$target, pan_net_filt$target)
  
  
  ###### Convert lists to human #########
  
  #Keep only the first instance of a mouse ensembl id
  human_2 <- human[!duplicated(human$Mouse.Ensembl),]
  
  shared <- unique(shared)
  pan_only <- unique(pan_only)
  gob_only <- unique(gob_only)
  
  #Get human conversions
  shared_h <- subset(human_2, Mouse.Ensembl %in% shared)
  pan_only_h <- subset(human_2, Mouse.Ensembl %in% pan_only)
  gob_only_h <- subset(human_2, Mouse.Ensembl %in% gob_only)
  
  
  ############# Functional enrichment ###########
  
  # Convert to entrez ids
  shared_h2 <- bitr_kegg(shared_h$Human.Uniprot, fromType='uniprot', toType='kegg', organism='hsa')
  pan_only_h2 <- bitr_kegg(pan_only_h$Human.Uniprot, fromType='uniprot', toType='kegg', organism='hsa')
  gob_only_h2 <- bitr_kegg(gob_only_h$Human.Uniprot, fromType='uniprot', toType='kegg', organism='hsa')
  
  # Reactome pathway enrichment
  shared_reactome <- enrichPathway(gene=shared_h2$kegg,pvalueCutoff=0.1,readable=T)
  pan_only_reactome <- enrichPathway(gene=pan_only_h2$kegg,pvalueCutoff=0.1,readable=T)
  gob_only_reactome <- enrichPathway(gene=gob_only_h2$kegg,pvalueCutoff=0.1,readable=T)
  
  # Get reactome results
  shared_reactome2 <- as.data.frame(shared_reactome)
  pan_only_reactome2 <- as.data.frame(pan_only_reactome)
  gob_only_reactome2 <- as.data.frame(gob_only_reactome)
  
  # KEGG pathway enrichment
  shared_kegg <- enrichKEGG(gene=shared_h2$kegg ,organism = 'hsa',pvalueCutoff = 0.1)
  pan_only_kegg <- enrichKEGG(gene=pan_only_h2$kegg,organism = 'hsa',pvalueCutoff=0.1)
  gob_only_kegg <- enrichKEGG(gene=gob_only_h2$kegg,organism = 'hsa',pvalueCutoff=0.1)
  
  # Get KEGG results
  shared_kegg2 <- as.data.frame(shared_kegg)
  pan_only_kegg2 <- as.data.frame(pan_only_kegg)
  gob_only_kegg2 <- as.data.frame(gob_only_kegg)
  
  # GO over representation
  shared_go <- enrichGO(gene=shared_h2$kegg,OrgDb = "org.Hs.eg.db", ont= "BP", pvalueCutoff=0.1)
  pan_only_go <- enrichGO(gene=pan_only_h2$kegg,OrgDb = "org.Hs.eg.db",ont= "BP",pvalueCutoff=0.1)
  gob_only_go <- enrichGO(gene=gob_only_h2$kegg,OrgDb = "org.Hs.eg.db",ont= "BP",pvalueCutoff=0.1)
  
  # Get GO results
  shared_go2 <- as.data.frame(shared_go)
  pan_only_go2 <- as.data.frame(pan_only_go)
  gob_only_go2 <- as.data.frame(gob_only_go)
  
  
  ############# Join results ###########
  
  #Add column of regulator, column of list type and column of kegg or reactome, join together
  if (nrow(shared_reactome2) > 0){
    shared_reactome2$TargetList <- "shared"
    shared_reactome2$Pathway <- "reactome"
    shared_reactome2$Regulator <- n
    all2 <- rbind(all2, shared_reactome2)
  }
  if (nrow(pan_only_reactome2) > 0){
    pan_only_reactome2$TargetList <- "pan_only"
    pan_only_reactome2$Pathway <- "reactome"
    pan_only_reactome2$Regulator <- n
    all2 <- rbind(all2, pan_only_reactome2)
  }
  if (nrow(gob_only_reactome2) > 0){
    gob_only_reactome2$TargetList <- "gob_only"
    gob_only_reactome2$Pathway <- "reactome"
    gob_only_reactome2$Regulator <- n
    all2 <- rbind(all2, gob_only_reactome2)
  }
  if (nrow(shared_kegg2) > 0){
    shared_kegg2$Pathway <- "kegg"
    shared_kegg2$TargetList <- "shared"
    shared_kegg2$Regulator <- n
    all2 <- rbind(all2, shared_kegg2)
  }
  if (nrow(pan_only_kegg2) > 0){
    pan_only_kegg2$TargetList <- "pan_only"
    pan_only_kegg2$Pathway <- "kegg"
    pan_only_kegg2$Regulator <- n
    all2 <- rbind(all2, pan_only_kegg2)
  }
  if (nrow(gob_only_kegg2) > 0){
    gob_only_kegg2$TargetList <- "gob_only"
    gob_only_kegg2$Pathway <- "kegg"
    gob_only_kegg2$Regulator <- n
    all2 <- rbind(all2, gob_only_kegg2)
  }
  if (nrow(shared_go2) > 0){
    shared_go2$TargetList <- "shared"
    shared_go2$Pathway <- "go"
    shared_go2$Regulator <- n
    all2 <- rbind(all2, shared_go2)
  }
  if (nrow(pan_only_go2) > 0){
    pan_only_go2$TargetList <- "pan_only"
    pan_only_go2$Pathway <- "go"
    pan_only_go2$Regulator <- n
    all2 <- rbind(all2, pan_only_go2)
  }
  if (nrow(gob_only_go2) > 0){
    gob_only_go2$TargetList <- "gob_only"
    gob_only_go2$Pathway <- "go"
    gob_only_go2$Regulator <- n
    all2 <- rbind(all2, gob_only_go2)
  }
  
}

###### output #######
write.table(all2, "Functional_analysis_rewired_p0.1.txt", quote = F, sep = "\t", row.names = F)
