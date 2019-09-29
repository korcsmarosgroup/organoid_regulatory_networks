# Functional analysis on the output tables from MCode module analysis in Cytoscape. 
# Converts mouse ensembl IDs to human and then will only carry out the enrichment analysis if there at least 15 gene ids in the cluster.
#
# Input: 1. input_p{g}: Mcode results tables downloaded from Cytoscape for Paneth and goblet datasets.
#        2. conversion: ID conversion file in a tab delimited text format with headers. Only 1 ID per line for the input ID (long format).
#           The ID headers should be 'Human.Ensembl' and 'Mouse.Ensembl'.
# Output: 1. mcode_clusters_functions_0.05qval.txt: Dataframe of enriched Reactome and KEGG pathways (q val <= 0.05) for each cluster in the input datasets.
#         2. Reactome{kegg}_functions_mcode_heatplot_top5.pdf: 2 PDF files of heatplots for the reactome and kegg results. 

##### Setup #####

# Packages
library(dplyr)
library(tidyr)
library(ReactomePA)
library(clusterProfiler)
library(pheatmap)
library(dendsort)
#library(RColorBrewer)
library(ggplot2)

# Input mcode files (skip the first 9 lines which had metadata in)
input_p <- read.csv("../Input_data/Mcode_Paneth.txt", header=T, skip=9, sep="\t")
input_g <- read.csv("../Input_data/Mcode_Goblet.txt", header=T, skip=9, sep ="\t")

# Input ID conversion file
conversion <- read.csv("../Input_data/InParanoid-Mus-homo-UniprotEnsembl-May2018-expanded.txt", sep = "\t")
  
##### Preprocess files #####

preprocess <- function(table, conversion_f, celltype){
  # Function to preprocess the imported mcode datatables
  #
  # Input table: Dataframe of mcode results.
  # Input conversion_f: Mouse to human ID conversion table with mouse ids in the column
  #                     'Mouse.Ensembl' and human in 'Human.Ensembl'
  # Input celltype: String representing the cell type for the imported datatable results.
  # Return clusters: List of long format dataframes. One for each cluster in the mcode results.
  #                 Contains column for human IDs and a column to represent the cell type.
  
  # Add cell type column
  table <- table %>% mutate(cell = celltype)
  
  # In the mcode data, split the gene ids and add as rows
  long <- table %>% separate_rows(Node.IDs, sep = ", ")
  
  # Add human gene IDs to the mcode data
  long <- left_join(long, conversion_f, by = c("Node.IDs" = "Mouse.Ensembl"))
  
  # Split the mcode data by clusters
  clusters <-split(long, f = long$Cluster)
  
  return(clusters)
}

# Remove cols in conversion table which are not required
conversion_f <- conversion %>% dplyr::select(Mouse.Ensembl, Human.Ensembl)

# Run preprocessing
clusters_p <- preprocess(input_p, conversion_f, "Paneth")
clusters_g <- preprocess(input_g, conversion_f, "goblet")

# List of cell type cluster tables
clusters_list <- list(clusters_p, clusters_g)

##### Functional enrichment #####

# Empty df's for data
all_functions <- data.frame()

# Function to run the Reactome enrichment analysis
reactome_analysis <- function(gene_factor){
  # Function to run Reactome enrichment analysis on an input gene list
  #
  # Input gene_factor: factor of human ensembl gene ids.
  # Return genes_reactome_df: Dataframe of enriched Reactome fucntions for the input gene ids.
  
  # Convert to entrez ids
  genes_entrez <- bitr(gene_factor, fromType='ENSEMBL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  #Reactome pathway enrichment
  genes_reactome <- enrichPathway(gene=genes_entrez$ENTREZID,qvalueCutoff=0.05,readable=T)
  
  #Get reactome results
  genes_reactome_df <- as.data.frame(genes_reactome)
  
  return(genes_reactome_df)
}

# Function to run the KEGG enrichment analysis
kegg_analysis <- function(gene_factor){
  # Function to run KEGG enrichment analysis on an input gene list
  #
  # Input gene_factor: factor of human ensembl gene ids.
  # Return genes_kegg_df: Dataframe of enriched KEGG fucntions for the input gene ids.
  
  # Convert to entrez ids
  genes_entrez <- bitr(gene_factor, fromType='ENSEMBL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  # KEGG pathway enrichment
  genes_kegg <- enrichKEGG(gene=genes_entrez$ENTREZID ,organism = 'hsa', qvalueCutoff = 0.05)
  
  # Get KEGG results
  genes_kegg_df <- as.data.frame(genes_kegg)
  
  return(genes_kegg_df)
}

# Iterate the cell types
for (cell in clusters_list){
    
  # Iterate the clusters
  for (i in cell){
    
    # Get cell type
    celltype <- i$cell[1]
    
    # Drop NAs and duplicates in the Human ensembl id column
    genes <- i %>% drop_na(Human.Ensembl) %>% distinct(Human.Ensembl)
    
    # Run Reactome enrichment if there are at least 20 genes in the cluster
    if(nrow(genes) >= 20) {
      
      # Run Reactome analysis
      reactome <- reactome_analysis(genes$Human.Ensembl)
      
      # Add name of cluster and celltype to dataframe
      reactome <- mutate(reactome, Cluster = i$Cluster[1], Pathway = "Reactome", CellType = celltype)
      
      #Append reactome results to previous results
      all_functions <- rbind(all_functions, reactome)
      
      # Run KEGG analysis
      kegg <- kegg_analysis(genes$Human.Ensembl)
      
      # Add name of cluster to dataframe
      kegg <- mutate(kegg, Cluster = i$Cluster[1], Pathway = "KEGG", CellType = celltype)
      
      #Append reactome results to previous results
      all_functions <- rbind(all_functions, kegg)
    }
  }
}

# Save the results
write.table(all_functions, "mcode_clusters_functions_0.05qval.txt", row.names = F, sep = "\t", quote = F)

##### Heatplot #####

plot_heatmap <- function(functions_df, colours) {
  # Function to create the heatplot
  #
  # Input functions_df : Dataframe with function name in the ' Description' column, 
  #                     cluster ids in the 'ClusterB' column and numerical value of 
  #                     r the pathway resource in the 'value' column.
  # Input colours : Character vector of colours to pass to the heatplot function 
  #                 (must have the same amount of categories as the input).
  #
  # Return plot : Heatplot
  
  # Create matrix of reactome terms for each cluster, with the value as the result
  mat <- acast(functions_df, list("Description", "ClusterB"), value.var="Cluster_num")
  
  # Replace nas with 0
  mat[is.na(mat)] = 0
  
  # Sort the dendograms
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse=FALSE))
  cluster_cols <- hclust(dist(t(mat)))
  cluster_cols <- sort_hclust(cluster_cols)
  cluster_rows <- sort_hclust(hclust(dist(mat)))
  cluster_rows <- sort_hclust(cluster_rows)
  
  # Create plot with all the data
  plot <- pheatmap(mat, 
                   #color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),
                   color = colours,
                   cellwidth = 11,
                   cellheight = 11,
                   fontsize_col = 12,
                   fontsize_row = 12,
                   show_colnames = TRUE,
                   show_rownames=TRUE,
                   legend = FALSE,
                   treeheight_col= 15,
                   treeheight_row = 15,
                   border_color = "lightgrey",
                   #annotation_row = ,
                   #main = ,
                   cluster_cols = cluster_cols,
                   cluster_rows = cluster_rows)
  return(plot)
}

# Edit the cluster ids so that the goblet cell clusters are 4,5,6 (distict numerically from the paneth ones, for the purpose of colouring the heatplots)
functions <- all_functions %>% mutate(Cluster_num = ifelse(CellType == "goblet", (Cluster+3), Cluster))

# Edit the cluster ids to include the celltype
functions <- functions %>% mutate(ClusterB = ifelse(CellType == "Paneth", paste0(Cluster, "_P"), paste0(Cluster, "_G")))

# Get the Reactome and KEGG results seperately
functions_r <- functions %>% filter(Pathway == "Reactome")
functions_k <- functions %>% filter(Pathway == "KEGG")

# Get the top 5 functions for each cluster
functions_r <- functions_r %>% group_by(Cluster_num) %>% top_n(-5, qvalue)
functions_k <- functions_k %>% group_by(Cluster_num) %>% top_n(-5, qvalue)

# Run the heatplot function
plot_r <- plot_heatmap(functions_r, c("white", "#FF6600", "#FF8432", "#FFBB8E", "#9900CC", "#BB44E5", "#C35DE5"))
plot_k <- plot_heatmap(functions_k, c("white", "#FF6600", "#FF8432", "#FFBB8E", "#9900CC", "#BB44E5"))

# Save the results
ggsave(plot_r, f = "Reactome_functions_mcode_heatplot_top5.pdf", device = "pdf",  width = 10, height = 11) 
ggsave(plot_k, f = "KEGG_functions_mcode_heatplot_top5.pdf", device = "pdf",  width = 8, height = 11)

