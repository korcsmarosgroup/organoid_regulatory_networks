# Script to draw 'heatplot' of mark-regulator relationship using an input file
#
# Input: 1. Regulator-marker networks in tab delimited text files. One for Paneth and one for goblet datasets.
#           The 'source' and 'target' column headers contain the source and target node ids respectively.
#        2. File with gene name annotations for each ensembl gene id. Tab delimited text file with first column
#           labelled 'gene_id' and the second 'gene_name'.
#        3. Tab delimited text files with details of the differentially expressed genes for the Paneth
#           and goblet datasets. Gene id in the 'target_id' column and log fold change in the 'log2fc' column.
#        4. Lists of high confidence markers for Paneth and goblet cell types (the input networks are for all markers).
#           Text file with one column, no header.
# Output: 1. Regulator-marker networks with gene names and specificity of the regulator. Tab delimited text files.
#         2. Heatplots generated but not saved automatically (plots split by marker type or by regulator type).

######## Packages #########

library(pheatmap)
library(reshape2)
library(dplyr)
library(dendsort)
library(grid)

####### Input files #########

# Regulator-marker networks
p_f <- "Pan_network_pan_markers_all.txt"
g_f <- "Gob_network_gob_markers_all.txt"

# Gene name annotations
annotation_f <- "gene_TPM_AllTogether.txt"

# Differentially expressed genes
p_lfc <- "Pan_DEGs_lfc1.txt"
g_lfc <- "Gob_DEGs_lfc1.txt"

# High conf markers
p_hc <- "PanethMarkers_Aviv_highc.txt"
g_hc <- "GobletMarkers_Aviv_highc.txt"

####### Output files #########

# Networks labelled with regulator specificity
pan_net_out <- "MarkerRegulators_Pan_Specificity.txt"
gob_net_out <- "MarkerRegulators_Gob_Specificity.txt"

###### Load Data #########

p <- read.delim(p_f, header = T)
g <- read.delim(g_f, header = T)
ann <- read.delim(annotation_f, header = T)
lfcp <- read.delim(p_lfc, header = T)
lfcg <- read.delim(g_lfc, header = T)
p_h<- read.delim(p_hc, header = F)
g_h<- read.delim(g_hc, header = F)

####### Preprocess Data #########

# Add discrete and cell type column, remove interaction type column
p1 <- p %>% select(source, target) %>%
  mutate(CellType = "Paneth") %>%
  mutate(Interaction = "1")
g1 <- g %>% select(source, target) %>%
  mutate(CellType = "Goblet") %>%
  mutate(Interaction = "1")

# Remove unnecessary columns of tpm file
ann <- ann %>% select(1,2)

# Get gene names instead of gene ids
p1 <- left_join(p1,ann, by = c("source" = "gene_id"))
p1 <- left_join(p1,ann, by = c("target" = "gene_id"))
p1 <- p1 %>% rename(Reg_gene = gene_name.x, Targ_gene = gene_name.y)
g1 <- left_join(g1,ann, by = c("source" = "gene_id"))
g1 <- left_join(g1,ann, by = c("target" = "gene_id"))
g1 <- g1 %>% rename(Reg_gene = gene_name.x, Targ_gene = gene_name.y)

# Remove mmu- from mirnas
p1$Reg_gene <- gsub('mmu-', '', p1$Reg_gene)
p1$Targ_gene <- gsub('mmu-', '', p1$Targ_gene)
g1$Reg_gene <- gsub('mmu-', '', g1$Reg_gene)
g1$Targ_gene <- gsub('mmu-', '', g1$Targ_gene)

####### Regulator categories #########

# Get all regulators + cell specificity
all_reg <- rbind(p1,g1)
all_reg <- all_reg %>% select(source,Reg_gene,CellType) %>%
  unique() %>%
  group_by(Reg_gene, source) %>%
  summarize(CellType = paste0(CellType, collapse = "-"))

# Append lfc value columns for each cell type
lfcp <- lfcp %>% select(target_id,log2fc_p = log2fc)
lfcg <- lfcg %>% select(target_id,log2fc_g = log2fc)
all_reg <- left_join(all_reg, lfcp, by = c("source" = "target_id"))
all_reg <- left_join(all_reg, lfcg, by = c("source" = "target_id"))

# categorise regulators  - if shared use different or same lfc
all_reg2 <- all_reg %>% mutate("Regulator specificity" = ifelse(CellType == "Paneth", "Paneth",
                                                  ifelse(CellType == "Goblet", "Goblet", 
                                                         ifelse((log2fc_p > 0 & log2fc_g >0), "Paneth & Goblet", 
                                                                ifelse((log2fc_p < 0 & log2fc_g <0), "Paneth & Goblet","Shared_dif")))))
# Remove other columns, row annotations
all_reg2 <- all_reg2 %>% ungroup() %>% select(Reg_gene,"Regulator specificity")

# Pan regulator annotations
p1_mat2 <- dcast(p1, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = length)
pan_regs <- p1_mat2 %>% select(Reg_gene)
pan_regs <- left_join(pan_regs, all_reg2)
pan_regs <- na.omit(pan_regs)
row.names(pan_regs) <- pan_regs$Reg_gene

# Gob regulator annotations
g1_mat2 <- dcast(g1, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = length)
gob_regs <- g1_mat2 %>% select(Reg_gene)
gob_regs <- left_join(gob_regs, all_reg2)
gob_regs <- na.omit(gob_regs)
row.names(gob_regs) <- gob_regs$Reg_gene

# Append regulators and edit Interaction to show specificity
p1 <- left_join(p1, pan_regs)
p1 <- p1 %>% rename(RegSpec = "Regulator specificity") %>%
  mutate(Interaction = ifelse(RegSpec == "Paneth", 1, ifelse(RegSpec == "Paneth & Goblet", 2,0)))
g1 <- left_join(g1, gob_regs)
g1 <- g1 %>% rename(RegSpec = "Regulator specificity") %>%
  mutate(Interaction = ifelse(RegSpec == "Goblet", 1, ifelse(RegSpec == "Paneth & Goblet", 2,0)))

pan_regs <- pan_regs %>% select("Regulator specificity")
gob_regs <- gob_regs %>% select("Regulator specificity")

# Save interactions
write.table(p1, pan_net_out, quote = F, row.names = F, sep = "\t")
write.table(g1, gob_net_out, quote = F, row.names = F, sep = "\t")

####### Reformat Data #########

p1_mat <- acast(p1, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = sum)
g1_mat <- acast(g1, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = sum)

####### Heatplot prep #########

#Sort the dendograms
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse=FALSE))
p1_cluster_cols <- hclust(dist(t(p1_mat)))
p1_cluster_cols <- sort_hclust(p1_cluster_cols)
p1_cluster_rows <- sort_hclust(hclust(dist(p1_mat)))
g1_cluster_cols <- hclust(dist(t(g1_mat)))
g1_cluster_cols <- sort_hclust(g1_cluster_cols)
g1_cluster_rows <- sort_hclust(hclust(dist(g1_mat)))

# Overwrite default draw_colnames in the pheatmap package.
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
#draw_colnames_45 <- function (coln, gaps, ...) {
#  coord <- pheatmap:::find_coordinates(length(coln), gaps)
#  x     <- coord$coord - 0.5 * coord$size
#  res   <- grid::textGrob(
#    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
#    vjust = 0.75, hjust = 1, rot = 90, gp = grid::gpar(...)
#  )
#  return(res)
#}
#assignInNamespace(
#  x = "draw_colnames",
#  value = "draw_colnames_45",
#  ns = asNamespace("pheatmap")
#)

####### Heatplot by Marker Type/Network #########

# Create plots
pplot<-pheatmap(p1_mat, 
            color = c("white", "darkorange", "dodgerblue"),
            cellwidth = 11,
            cellheight = 11,
            fontsize_col = 12,
            fontsize_row = 12,
            show_colnames = TRUE,
            show_rownames=TRUE,
            legend = FALSE,
            treeheight_col=0,
            treeheight_row = 0,
            border_color = "lightgrey",
            #annotation_row = pan_regs,
            #main = "Regulators of Paneth Markers - Paneth Network",
            cluster_cols=p1_cluster_cols,
            cluster_rows =p1_cluster_rows)
pplot

gplot<-pheatmap(g1_mat, 
                color = c("white",  "darkorange", "dodgerblue"),
                cellwidth = 5,
                cellheight = 5,
                fontsize_col=6,
                fontsize_row =6,
                show_colnames = TRUE,
                show_rownames=TRUE,
                legend = FALSE,
                treeheight_col=0,
                treeheight_row =0,
                border_color = "lightgrey",
                #annotation_row = gob_regs,
                #main = "Regulators of Goblet Markers - Goblet Network",
                cluster_cols=g1_cluster_cols,
                cluster_rows =g1_cluster_rows)
gplot

####### Heatplot by Regulator type #########

# Join dataframes together
all <- rbind(p1,g1)

# Add regulator specificity information
all <-  left_join(all, all_reg2)

# Split by regulator specificity
split_reg <- split(all, all$`Regulator specificity`)
gob_reg_sp <- split_reg[[1]]
pan_reg_sp <- split_reg[[2]]
sh_reg_sp <- split_reg[[3]]

# Get marker annotation
mark_ann <- sh_reg_sp %>% select(CellType,Targ_gene) %>%
  group_by(Targ_gene) %>%
  summarise(CellType = paste0(CellType, collapse = "-")) %>%
  mutate('Marker Cell Specificity' = ifelse(CellType == "Paneth-Goblet", "Paneth & Goblet", ifelse(grepl('^P', .$CellType), 'Paneth', 'Goblet'))) %>%
  select(Targ_gene,'Marker Cell Specificity')

# Remove duplicated row from shared (due to mmp7)
sh_reg_sp2 <- sh_reg_sp %>% select(Reg_gene, Targ_gene, Interaction) %>%
  unique()
# Get annotation from the exact dataset
sh_reg_sp2 <- left_join(sh_reg_sp2, mark_ann)
annotation <- sh_reg_sp2 %>% select(Targ_gene, `Marker Cell Specificity`) %>%
  unique()
row.names(annotation) <- annotation$Targ_gene
annotation <- annotation %>% select(`Marker Cell Specificity`)

# Reformat
p_mat_r <- acast(pan_reg_sp, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = length)
g_mat_r <- acast(gob_reg_sp, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = length)
sh_mat_r <- acast(sh_reg_sp2, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = length)

# Sort the dendograms
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse=FALSE))
pr_cluster_cols <- hclust(dist(t(p_mat_r)))
pr_cluster_cols <- sort_hclust(pr_cluster_cols)
pr_cluster_rows <- sort_hclust(hclust(dist(p_mat_r)))
gr_cluster_cols <- hclust(dist(t(g_mat_r)))
gr_cluster_cols <- sort_hclust(gr_cluster_cols)
gr_cluster_rows <- sort_hclust(hclust(dist(g_mat_r)))
sr_cluster_cols <- hclust(dist(t(sh_mat_r)))
sr_cluster_cols <- sort_hclust(sr_cluster_cols)
sr_cluster_rows <- sort_hclust(hclust(dist(sh_mat_r)))

# Create plots
pr_plot<-pheatmap(p_mat_r, 
                color = c("white", "slategray3"),
                cellwidth = 10,
                cellheight = 10,
                fontsize_col = 10,
                fontsize_row = 10,
                show_colnames = TRUE,
                show_rownames=TRUE,
                legend = FALSE,
                treeheight_col=0,
                treeheight_row = 0,
                main = "Regulator-marker interactions : Paneth specific markers, Paneth network",
                cluster_cols=pr_cluster_cols,
                cluster_rows =pr_cluster_rows)
gr_plot<-pheatmap(g_mat_r, 
                  color = c("white", "slategray3"),
                  cellwidth = 5,
                  cellheight = 5,
                  fontsize_col = 6,
                  fontsize_row = 6,
                  show_colnames = TRUE,
                  show_rownames=TRUE,
                  legend = FALSE,
                  treeheight_col=0,
                  treeheight_row = 0,
                  main = "Regulator-marker interactions : Goblet specific markers, Goblet network",
                  cluster_cols=gr_cluster_cols,
                  cluster_rows =gr_cluster_rows)
sr_plot<-pheatmap(sh_mat_r, 
                  color = c("white", "slategray3"),
                  cellwidth = 4.5,
                  cellheight = 4.5,
                  fontsize_col = 4.5,
                  fontsize_row = 4.5,
                  show_colnames = TRUE,
                  show_rownames=TRUE,
                  legend = FALSE,
                  treeheight_col=0,
                  treeheight_row = 0,
                  annotation_col = annotation,
                  annotation_legend = F,
                  annotation_names_col = F,
                  main = "Regulator-marker interactions : Shared markers",
                  cluster_cols = sr_cluster_cols,
                  cluster_rows =sr_cluster_rows)

####### Heatplot High Conf Markers #########

# filter for high confidnce markers
p1_hc <- inner_join(p1, p_h, by = c("target" = "V1"))
g1_hc <- inner_join(g1, g_h, by = c("target" = "V1"))

# Pan regulator annotations - annotations refer to regualtors of all markers
p1_hc_mat2 <- dcast(p1_hc, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = length)
pan_regs_hc <- p1_hc_mat2 %>% select(Reg_gene)
pan_regs_hc <- left_join(pan_regs_hc, all_reg2)
row.names(pan_regs_hc) <- pan_regs_hc$Reg_gene
pan_regs_hc <- pan_regs_hc %>% select("Regulator specificity")

# Gob regulator annotations - annotations refer to regualtors of all markers
g1_hc_mat2 <- dcast(g1_hc, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = length)
gob_regs_hc <- g1_mat2 %>% select(Reg_gene)
gob_regs_hc <- left_join(gob_regs_hc, all_reg2)
row.names(gob_regs_hc) <- gob_regs_hc$Reg_gene
gob_regs_hc <- gob_regs_hc %>% select("Regulator specificity")

# Reformat
p1_mat_hc <- acast(p1_hc, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = sum)
g1_mat_hc <- acast(g1_hc, Reg_gene~Targ_gene, value.var="Interaction", fun.aggregate = sum)

# Sort the dendograms
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse=FALSE))
p1_hc_cluster_cols <- hclust(dist(t(p1_mat_hc)))
p1_hc_cluster_cols <- sort_hclust(p1_hc_cluster_cols)
p1_hc_cluster_rows <- sort_hclust(hclust(dist(p1_mat_hc)))
g1_hc_cluster_cols <- hclust(dist(t(g1_mat_hc)))
g1_hc_cluster_cols <- sort_hclust(g1_hc_cluster_cols)
g1_hc_cluster_rows <- sort_hclust(hclust(dist(g1_mat_hc)))

# Create plots
pplot<-pheatmap(p1_mat_hc, 
                color = c("white",  "darkorange", "dodgerblue"),
                cellwidth = 14,
                cellheight = 14,
                fontsize_col = 12,
                fontsize_row = 12,
                show_colnames = TRUE,
                show_rownames=TRUE,
                legend = FALSE,
                treeheight_col=0,
                treeheight_row = 0,
                border_color = "lightgrey",
                #annotation_row = pan_regs_hc,
                #main = "Regulators of Paneth Markers - High confidence markers",
                cluster_cols=p1_hc_cluster_cols,
                cluster_rows =p1_hc_cluster_rows)
pplot

gplot<-pheatmap(g1_mat_hc, 
                color = c("white", "slategray3"),
                cellwidth = 7.5,
                cellheight = 7.5,
                fontsize_col=7,
                fontsize_row =7,
                show_colnames = TRUE,
                show_rownames=TRUE,
                legend = FALSE,
                treeheight_col=0,
                treeheight_row =0,
                annotation_row = gob_regs_hc,
                main = "Regulators of Goblet Markers - High confidence markers",
                cluster_cols=g1_hc_cluster_cols,
                cluster_rows =g1_hc_cluster_rows)
gplot
