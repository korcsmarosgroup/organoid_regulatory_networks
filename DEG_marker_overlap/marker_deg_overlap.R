# Carry out DEG-marker overlaps
#
# Input: 1. Tab delimited text files with headers where each row represents a differentially expressed gene and the mouse
#        ensembl ID for the gene is in the column with header 'target_id', the miRNA annotation is in the column labelled 'final_annotation',
#        and the novel lncRNAs are labelled in the column 'final_annotation.y'. One file for Paneth cell data and one for goblet cell data.
#        2. Text files with 1 column lists of marker genes in mouse ensembl format - no header. All markers and high confidence markers (Haber et al.) 
#        after removing markers which are not in the Wald output of differential expression (genes with variance greater than zero among samples).
#        3. File containing the background gene set for which the hypergeometric test requires. We are using the Wald test output which has Ensembl ids in the target_id' column.
# Output: 1. Tab delimited txt file with statistics for the overlap between the degs and the markers
#         2. Tab delimited text file with heatplot style layout with the hypergeometric test results (best opened in excel).

###### Setup ##########

# Libraries
library(dplyr)

# Input paths
p_deg1_f <- "Pan_DEGs_lfc1.txt"
g_deg1_f <- "Gob_DEGs_lfc1.txt"

pan_mark_f <- "PanethMarkers_Aviv.txt"
gob_mark_f <- "GobletMarkers_Aviv.txt"
eec_mark_f <- "EnteroendoMarkers_Aviv.txt"
enteroc_mark_f <- "EnterocyteMarkers_Aviv.txt"
tuft_mark_f <-  "TuftMarkers_Aviv.txt"

pan_mark_hf <- "PanethMarkers_Aviv.txt"
gob_mark_hf <- "GobletMarkers_Aviv.txt"
eec_mark_hf <- "EnteroendoMarkers_Aviv.txt"
enteroc_mark_hf <- "EnterocyteMarkers_Aviv.txt"
tuft_mark_hf <-  "TuftMarkers_Aviv.txt"

# Output paths
overlap_tab <- "Marker-DEG-Overlap.txt"
overlap_hyper <- "Marker-DEG-Overlap-hyp.txt"

#Use wald file for background gene set (both wald goblet and paneth file share same genes, so only need to import 1)
wald_f <-  "wald_test_Paneth_vs_Control.with_gene_names.txt"

###### Load Input files ##########

p_deg1 <- read.delim(p_deg1_f)
g_deg1 <- read.delim(g_deg1_f)

pan_mark <- read.delim(pan_mark_f, header=F)
gob_mark <- read.delim(gob_mark_f, header=F)
eec_mark <- read.delim(eec_mark_f, header=F)
enteroc_mark <- read.delim(enteroc_mark_f, header=F)
tuft_mark <- read.delim(tuft_mark_f, header=F)

pan_markh <- read.delim(pan_mark_hf, header=F)
gob_markh <- read.delim(gob_mark_hf, header=F)
eec_markh <- read.delim(eec_mark_hf, header=F)
enteroc_markh <- read.delim(enteroc_mark_hf, header=F)
tuft_markh <- read.delim(tuft_mark_hf, header=F)

wald <- read.delim(wald_f)

###### File pre-processing ##########

#Remove -ve degs from files
p_deg1 <- filter(p_deg1, log2fc >=1)
g_deg1 <- filter(g_deg1, log2fc >=1)

#Remove novel lincrnas and mirnas from the DEG lists
p_deg1<- filter(p_deg1, final_annotation != "miRNA" & final_annotation.y != "potential_novel_lncRNA")
g_deg1<- filter(g_deg1, final_annotation != "miRNA" & final_annotation.y != "potential_novel_lncRNA")

#Remove MSTRG (novel lncrna) genes from wald file
wald1 <- filter(wald, !grepl("^MSTR", gene_name))
wald_count <- nrow(wald1)

###### Overlap All Markers ##########

#Paneth degs
p1_p_all <- inner_join(p_deg1, pan_mark, by = c(target_id = "V1"))
p1_g_all <- inner_join(p_deg1, gob_mark, by = c(target_id = "V1"))
p1_eec_all <- inner_join(p_deg1, eec_mark, by = c(target_id = "V1"))
p1_ec_all <- inner_join(p_deg1, enteroc_mark, by = c(target_id = "V1"))
p1_t_all <- inner_join(p_deg1, tuft_mark, by = c(target_id = "V1"))

#Paneth table
marks <- c("Paneth", "Goblet", "Entoerendocrine", "Tuft", "Enterocyte")
marks2 <- c("All","All","All","All","All")
overlap <- c(nrow(p1_p_all), nrow(p1_g_all),nrow(p1_eec_all),nrow(p1_t_all),nrow(p1_ec_all))
totalmarks <- c(nrow(pan_mark),nrow(gob_mark),nrow(eec_mark),nrow(tuft_mark),nrow(enteroc_mark))
cell <- c("Paneth_lfc1","Paneth_lfc1","Paneth_lfc1","Paneth_lfc1","Paneth_lfc1")
table <- data.frame(marks,marks2, totalmarks, cell,overlap)
colnames(table) <- c("MarkerList_Cell","MarkerList_Length","Total_Markers","Input_DEG_List","DEG_Markers")

#Add percentages
table <- mutate(table, Percentage = ((DEG_Markers/Total_Markers)*100))

#Goblet degs
g1_p_all <- inner_join(g_deg1, pan_mark, by = c(target_id = "V1"))
g1_g_all <- inner_join(g_deg1, gob_mark, by = c(target_id = "V1"))
g1_eec_all <- inner_join(g_deg1, eec_mark, by = c(target_id = "V1"))
g1_ec_all <- inner_join(g_deg1, enteroc_mark, by = c(target_id = "V1"))
g1_t_all <- inner_join(g_deg1, tuft_mark, by = c(target_id = "V1"))

#Goblet table
marks <- c("Paneth", "Goblet", "Entoerendocrine", "Tuft", "Enterocyte")
marks2 <- c("All","All","All","All","All")
overlap <- c(nrow(g1_p_all), nrow(g1_g_all),nrow(g1_eec_all),nrow(g1_t_all),nrow(g1_ec_all))
totalmarks <- c(nrow(pan_mark),nrow(gob_mark),nrow(eec_mark),nrow(tuft_mark),nrow(enteroc_mark))
cell <- c("Goblet-Ent_lfc1","Goblet-Ent_lfc1","Goblet-Ent_lfc1","Goblet-Ent_lfc1","Goblet-Ent_lfc1")
table2 <- data.frame(marks,marks2, totalmarks, cell,overlap)
colnames(table2) <- c("MarkerList_Cell","MarkerList_Length","Total_Markers","Input_DEG_List","DEG_Markers")

#Add percentages
table2 <- mutate(table2, Percentage = ((DEG_Markers/Total_Markers)*100))

#Export Paneth degs- Paneth markers overlap list
write.table(p1_p_all, "PanethDEGs_overlap_PanethMarkersAll.txt", quote = F, row.names = F, sep = "\t")

#Export Paneth degs- Paneth markers overlap list etc.
write.table(p1_p_all, "PanethDEGs_overlap_PanethMarkersAll.txt", quote = F, row.names = F, sep = "\t")
write.table(g1_g_all, "GobDEGs_overlap_GobletMarkersAll.txt", quote = F, row.names = F, sep = "\t")
write.table(g1_eec_all, "GobDEGs_overlap_EECMarkersAll.txt", quote = F, row.names = F, sep = "\t")

###### Overlap Hconf Markers ##########

#Paneth degs
p1_p_h <- inner_join(p_deg1, pan_markh, by = c(target_id = "V1"))
p1_g_h <- inner_join(p_deg1, gob_markh, by = c(target_id = "V1"))
p1_eec_h <- inner_join(p_deg1, eec_markh, by = c(target_id = "V1"))
p1_ec_h <- inner_join(p_deg1, enteroc_markh, by = c(target_id = "V1"))
p1_t_h <- inner_join(p_deg1, tuft_markh, by = c(target_id = "V1"))

#Paneth table
marks <- c("Paneth", "Goblet", "Entoerendocrine", "Tuft", "Enterocyte")
marks2 <- c("H-conf","H-conf","H-conf","H-conf","H-conf")
overlap <- c(nrow(p1_p_h), nrow(p1_g_h),nrow(p1_eec_h),nrow(p1_t_h),nrow(p1_ec_h))
totalmarks <- c(nrow(pan_markh),nrow(gob_markh),nrow(eec_markh),nrow(tuft_markh),nrow(enteroc_markh))
cell <- c("Paneth_lfc1","Paneth_lfc1","Paneth_lfc1","Paneth_lfc1","Paneth_lfc1")
table3 <- data.frame(marks,marks2, totalmarks, cell,overlap)
colnames(table3) <- c("MarkerList_Cell","MarkerList_Length","Total_Markers","Input_DEG_List","DEG_Markers")

#Add percentages
table3 <- mutate(table3, Percentage = ((DEG_Markers/Total_Markers)*100))

#Goblet degs
g1_p_h <- inner_join(g_deg1, pan_markh, by = c(target_id = "V1"))
g1_g_h <- inner_join(g_deg1, gob_markh, by = c(target_id = "V1"))
g1_eec_h <- inner_join(g_deg1, eec_markh, by = c(target_id = "V1"))
g1_ec_h <- inner_join(g_deg1, enteroc_markh, by = c(target_id = "V1"))
g1_t_h <- inner_join(g_deg1, tuft_markh, by = c(target_id = "V1"))

#Goblet table
marks <- c("Paneth", "Goblet", "Entoerendocrine", "Tuft", "Enterocyte")
marks2 <- c("H-conf","H-conf","H-conf","H-conf","H-conf")
overlap <- c(nrow(g1_p_h), nrow(g1_g_h),nrow(g1_eec_h),nrow(g1_t_h),nrow(g1_ec_h))
totalmarks <- c(nrow(pan_markh),nrow(gob_markh),nrow(eec_markh),nrow(tuft_markh),nrow(enteroc_markh))
cell <- c("Goblet-Ent_lfc1","Goblet-Ent_lfc1","Goblet-Ent_lfc1","Goblet-Ent_lfc1","Goblet-Ent_lfc1")
table4 <- data.frame(marks,marks2, totalmarks, cell,overlap)
colnames(table4) <- c("MarkerList_Cell","MarkerList_Length","Total_Markers","Input_DEG_List","DEG_Markers")

#Add percentages
table4 <- mutate(table4, Percentage = ((DEG_Markers/Total_Markers)*100))

#Export Paneth degs- Paneth markers overlap list
write.table(p1_p_h, "PanethDEGs_overlap_PanethMarkersHconf.txt", quote = F, row.names = F, sep = "\t")
write.table(g1_g_h, "GobDEGs_overlap_GobletMarkersHconf.txt", quote = F, row.names = F, sep = "\t")
write.table(g1_eec_h, "GobDEGs_overlap_EECMarkersHconf.txt", quote = F, row.names = F, sep = "\t")


###### Join tables and export ##########
all_tables <- rbind(table,table2,table3,table4)

write.table(all_tables, overlap_tab, row.names = F, quote = F, sep = "\t")

rm(table,table2,table3,table4)

###### Hypogeometric significance score ##########

# DEG list lengths (without mirnas or novel linrnas)
yp1 <- nrow(p_deg1)
yg1 <- nrow(g_deg1)

# Number of genes in wald test (doesn't nclude novel lncrnas and mirnas)
a <- wald_count

# Remove tuft and enterocyte (or if you want to keep them, need to alter the multiple testing correction of the hypergeometric score)
all_tables_2 <- filter(all_tables, MarkerList_Cell != "Tuft" & MarkerList_Cell != "Enterocyte")

# Carry out hypergeometic score calcualtion followed by multiple testing correction
all_tables_3 <- all_tables_2 %>% 
  mutate(Num_DEGs_linc_pcgs = ifelse(Input_DEG_List == "Paneth_lfc1", yp1, yg1)) %>%
  mutate(Hypergeo_pval = phyper(DEG_Markers-1, Total_Markers, a-Total_Markers, Num_DEGs_linc_pcgs, lower.tail = F, log.p = F)) %>%
  mutate(hyp_pval_corr = Hypergeo_pval*6)

# Export/save
write.table(all_tables_3, overlap_hyper, row.names = F, quote = F, sep = "\t")

