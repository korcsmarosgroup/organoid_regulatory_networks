# Script to get marker regulators from interactions files, split by specificity,
# convert to human and carry out functional enrichment.
#
# Input: 1. Regulator-marker networks in tab delimited text files. One for Paneth and one for goblet datasets.
#           The 'source' and 'target' column headers contain the source and target node ids respectively.
#        2. ID conversion file in a tab delimited text format with headers. Only 1 ID per line for the input ID (long format).
#           The ensembl ID headers should be 'Human Ensembl' and 'Mouse Ensembl'.
#           The uniprot ID headers should be 'Human Uniprot' and 'Mouse Uniprot'.
# Output: 1. Tab delimited text files with the significant Reactome pathways associated with specific lists of
#            regulators: paneth all, paneth specific, goblet all, goblet specific, shared.
#         2. PDF images of dot plot output from, KEGG and Reactome group analaysis of regulators split by specificity
#           (pan only, gob only, shared)
#         3. Tab delimited text files with the human and mouse IDs for each of the regulators
#           (Paneth and goblet marker regulators).

######## Packages #########

library(dplyr)
library(tidyr)
library(ReactomePA)
library(clusterProfiler)

####### Input files #########

# Regulator-marker networks
pan_f <- "Pan_network_pan_markers_all.txt"
gob_f <-  "Gon_network_gob_markers_all.txt"

# Id conversion table
human_f <- "Inparanoid-Mus-homo-UniprotEnsembl-May2018-expanded.txt"

##### Open files #######

pan <- read.delim(pan_f, header=T)
gob <- read.delim(gob_f, header=T)

convert <- read.delim(human_f)

##### Convert regulators to human #####

# Keep only the first instance of a mouse ensembl id
convert2 = convert[!duplicated(convert$Mouse.Ensembl),]

# Warning - the select doesnt work if the reactome pa or cluster profiler packages are loaded, unless 'dplyr::' specified.
gob_r <- gob %>% dplyr::select(V1) %>% unique()
pan_r <- dplyr::select(pan, V1) %>% unique()

# Join conversion table to regulator tables
pan_join <- left_join(pan_r, convert2, by = c("source" = "Mouse.Ensembl" ))
gob_join <- left_join(gob_r, convert2, by = c("source" = "Mouse.Ensembl" ))

##### Split regulators by specificity #######

# Get lists of human regulators
gob_r <- unique(gob_join$Human.Uniprot)
pan_r <- unique(pan_join$Human.Uniprot)

# Get shared and cell type specific regulators
shared_r <- intersect(gob_r, pan_r)
pan_sp_r <- setdiff(pan_r, gob_r)
gob_sp_r <- setdiff(gob_r, pan_r)

# Remove any NAs
shared_r <- shared_r[!is.na(shared_r)]
pan_r <- pan_r[!is.na(pan_r)]
gob_r <- gob_r[!is.na(gob_r)]
pan_sp_r <- pan_sp_r[!is.na(pan_sp_r)]
gob_sp_r <- gob_sp_r[!is.na(gob_sp_r)]

#### Functional Enrichment ####

# Convert to entrez ids
shared_r_e <- bitr_kegg(shared_r, fromType='uniprot', toType='kegg', organism='hsa')
pan_sp_r_e <- bitr_kegg(pan_sp_r, fromType='uniprot', toType='kegg', organism='hsa')
gob_sp_r_e <- bitr_kegg(gob_sp_r, fromType='uniprot', toType='kegg', organism='hsa')
pan_r_e <- bitr_kegg(pan_r, fromType='uniprot', toType='kegg', organism='hsa')
gob_r_e <- bitr_kegg(gob_r, fromType='uniprot', toType='kegg', organism='hsa')

# Reactome pathway enrichment
shared_reactome <- enrichPathway(gene=shared_r_e$kegg,pvalueCutoff=0.05, readable=T)
pan_sp_reactome <- enrichPathway(gene=pan_sp_r_e$kegg,pvalueCutoff=0.05, readable=T)
gob_sp_reactome <- enrichPathway(gene=gob_sp_r_e$kegg,pvalueCutoff=0.05, readable=T)
pan_reactome <- enrichPathway(gene=pan_r_e$kegg,pvalueCutoff=0.05, readable=T)
gob_reactome <- enrichPathway(gene=gob_r_e$kegg,pvalueCutoff=0.05, readable=T)

# Get results
shared_reactome2 <- as.data.frame(shared_reactome)
pan_sp_reactome2 <- as.data.frame(pan_sp_reactome)
gob_sp_reactome2 <- as.data.frame(gob_sp_reactome)
pan_reactome2 <- as.data.frame(pan_reactome)
gob_reactome2 <- as.data.frame(gob_reactome)

# KEGG pathway enrichment
#kegg_shared <- enrichKEGG(gene =shared_r_e$kegg ,organism = 'hsa',pvalueCutoff = 0.05)
#kegg_pan_sp <- enrichKEGG(gene =pan_sp_r_e$kegg ,organism = 'hsa',pvalueCutoff = 0.05)
#kegg_gob_sp <- enrichKEGG(gene =gob_sp_r_e$kegg ,organism = 'hsa',pvalueCutoff = 0.05)


#### Visualise results ####

# All 3 reactome together
list1 <- list(pan_sp_r_e$kegg, shared_r_e$kegg, gob_sp_r_e$kegg)
names(list1) <- c("Paneth specific","Shared","Goblet specific")
res <- compareCluster(list1, fun="enrichPathway")
all_reactome <- dotplot(res)

#cnetplot(pan_sp_reactome, categorySize="pvalue")

# Kegg all together
res_kegg <- compareCluster(list1, fun="enrichKEGG")
all_kegg <- dotplot(res_kegg)

#### Export results ####

write.table(shared_reactome2, "Shared_Regulators_Reactome.txt", quote = F, sep = "\t", row.names =F)
write.table(pan_sp_reactome2, "Paneth_specifc_Regulators_Reactome.txt", quote = F, sep = "\t", row.names =F)
write.table(gob_sp_reactome2, "Goblet_specific_Regulators_Reactome.txt", quote = F, sep = "\t", row.names =F)
write.table(pan_reactome2, "Paneth_all_Regulators_Reactome.txt", quote = F, sep = "\t", row.names =F)
write.table(gob_reactome2, "Goblet_all_Regulators_Reactome.txt", quote = F, sep = "\t", row.names =F)

pdf("Regulators_bySpecificity_Reactome.pdf", width = 10, height = 7)
all_reactome
dev.off()

pdf("Regulators_bySpecificity_KEGG.pdf", width = 10, height = 7)
all_kegg
dev.off()

# Output lists of regulators by specificity - mouse and human
write.table(pan_join, "Regulator_List_Paneth.txt", quote = F, sep = "\t", row.names = F)
write.table(gob_join, "Regulator_List_Goblet.txt", quote = F, sep = "\t", row.names = F)
