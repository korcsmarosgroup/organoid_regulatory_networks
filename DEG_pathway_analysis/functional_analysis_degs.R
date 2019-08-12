# Functional analysis on input list of human ensembl ids
#
# Input: Tab delimited text files with headers where each row represents a differentially expressed gene and the human
#        ensembl ID for the gene is in the column with header 'Human.Ensembl'. One file for Paneth cell data and one for goblet cell data.
# Output: 1. Tab delimited text file with headers where each row is one significantly enriched pathway/function. 
#         The Paneth specific, goblet specific and shared results are output in one file. 
#         The Reactome, Go and KEGG results are output in one file each.
#         2. Plot of top 10 Reactome results for each DEG list.

##### Setup #####

# Packages
library(dplyr)
library(ReactomePA)
library(clusterProfiler)
library(ggplot2)
library(forcats)


#Input files
input_g <- read.delim("Gob_DEGs_lfc1_hum.txt", header=T)
input_p <- read.delim("Pan_DEGs_lfc1_hum.txt", header=T)

##### Preprocesing files #####

# Tidy input data
input_g2 <- input_g %>% select(Human.Ensembl, log2fc, qval) %>% filter(Human.Ensembl != "" & Human.Ensembl != "-") %>% 
  mutate(Hum.Ens = gsub("\\;.*","",Human.Ensembl)) %>% select (-c(Human.Ensembl)) %>% mutate(cell = "goblet")
input_p2 <- input_p %>% select(Human.Ensembl, log2fc, qval) %>% filter(Human.Ensembl != "" & Human.Ensembl != "-") %>%
  mutate(Hum.Ens = gsub("\\;.*","",Human.Ensembl)) %>% select (-c(Human.Ensembl)) %>% mutate(cell = "paneth")

# Get overlaps
input_both = full_join(input_g2, input_p2, by=c("Hum.Ens", "Hum.Ens")) %>% mutate(cell_types = paste(cell.x, cell.y))
input_split = split(input_both, input_both$cell_types)

# Empty df's for data
all_reactome <- data.frame()
all_kegg <- data.frame()
all_go <- data.frame()

##### Functional analysis #####

# Carry out analysis on the Paneth specific, goblet specific and shared DEGs separately
for (list in input_split) {
  
  # Sort by and extract top 50 by lfc
  if (list$cell_types[1] == "goblet NA"){
    list_g <- list %>% top_n(-50,qval.x)
  } else if (list$cell_types[1] == "NA paneth") {
    list_p <- list %>% top_n(-50,qval.y)
  } else {
    list_gp <- list %>% mutate(qval.xy = qval.x + qval.y) %>% top_n(-50,qval.xy)
  }
  
  #convert to entrez ids
  input_e <- bitr(list$Hum.Ens, fromType='ENSEMBL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  #Get list name
  list_f <- list$cell_types[1]
  
  #### REACTOME ####
  
  #Reactome pathway enrichment
  genes_reactome <- enrichPathway(gene=input_e$ENTREZID,pvalueCutoff=0.05,readable=T)
  #Get reactome results
  genes_reactome2 <- as.data.frame(genes_reactome)
  genes_reactome2 <- mutate(genes_reactome2, group = list_f)
  
  #Append reactome results to previous results
  all_reactome <- rbind(all_reactome, genes_reactome2)
  
  # Seperate plots for each
  #pdf(paste(list_f,"2.pdf"))
  #print(dotplot(genes_reactome, showCategory=10))
  #dev.off() 
  
  #### KEGG ####
  
  #KEGG pathway enrichment
  genes_kegg <- enrichKEGG(gene=input_e$ENTREZID ,organism = 'hsa',pvalueCutoff = 0.05)
  #Get KEGG results
  genes_kegg2 <- as.data.frame(genes_kegg)
  genes_kegg2 <- mutate(genes_kegg2, group = list_f)
  
  #Append reactome results to previous results
  all_kegg <- rbind(all_kegg, genes_kegg2)
  
  #### GO ####
  
  #GO over representation
  genes_go <- enrichGO(gene=input_e$ENTREZID,OrgDb = "org.Hs.eg.db", ont= "BP", pvalueCutoff=0.05)
  #Get GO results
  genes_go2 <- as.data.frame(genes_go)
  genes_go2 <- mutate(genes_go2, group = list_f)
  
  #Append reactome results to previous results
  all_go <- rbind(all_go, genes_go2)
}

##### Save all Reactome plots together #####

# Change labels and get top 10 by p adjust
all_reactome_top10 <- all_reactome %>% 
  mutate(cell = ifelse(group == "goblet NA", "GCeE - specific", ifelse(group == "NA paneth", "PCeE - specific", "PCeE and GCeE"))) %>% 
  group_by(cell) %>% top_n(-10, p.adjust)
# Convert gene ratio to decimal
all_reactome_top10$GeneRatio <- sapply(all_reactome_top10$GeneRatio, function(x) eval(parse(text=x)))
# plot
p <- ggplot(all_reactome_top10, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(color = p.adjust,size = Count)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(trans = "log", high="blue", low="red", breaks= c(1e-4, 1e-8, 1e-12, 1e-17)) +
  ylab(NULL)

pdf("FunctionalEnrichment_Reactome_0.05.pdf")
p + facet_grid(cell~., scales = "free") + 
  theme(legend.key.size = unit(0.5, "cm"),legend.key.width = unit(0.5,"cm"),legend.text = element_text(size=10),
        text = element_text(size=14, color="black"),
        strip.text.y = element_text(size = 14, face = "bold"))
  
dev.off()

##### Save table outputs #####
write.table(all_reactome, "FunctionalEnrichment_Reactome_0.05.txt",quote = F, sep = "\t", row.names = F)
write.table(all_kegg, "FunctionalEnrichment_KEGG_0.05.txt",quote = F, sep = "\t", row.names = F)
write.table(all_go, "FunctionalEnrichment_GO_0.05.txt",quote = F, sep = "\t", row.names = F)
