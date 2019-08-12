# Filter the cell type marker lists to exclude genes not in the wald file list of genes (genes with variance greater than zero among samples)
# This is done so that the background is standardised to calculate hypergeometric significance on the DEG-marker overlaps.
#
# Input: 1. Wald test output which has Ensembl ids in the 'target_id' column.
#        2. Text files with 1 column lists of marker genes in mouse ensembl format - no header. All markers and high confidence markers (Haber et al.).
# Output: 1. Tab delimited txt file with statistics for the overlap between the degs and the markers
#         2. Tab delimited text file with heatplot style layout with the hypergeometric test results (best opened in excel).

###### Setup ##########

# Libraries
library(dplyr)

# Input
wald_f <- "wald_test_Paneth_vs_Control.with_gene_names.txt" # Only one needed as the Paneth and goblet files have the same genes in.

# Marker lists
pan_mark_f <- "/All Markers - Haber/PanethMarkers_Aviv.txt"
gob_mark_f <- "/All Markers - Haber/GobletMarkers_Aviv.txt"
eec_mark_f <- "All Markers - Haber/EnteroendoMarkers_Aviv.txt"
enteroc_mark_f <- "/All Markers - Haber/EnterocyteMarkers_Aviv.txt"
tuft_mark_f <-  "/All Markers - Haber/TuftMarkers_Aviv.txt"

pan_mark_hf <- "//High confidence Markers - Haber/PanethMarkers_Aviv.txt"
gob_mark_hf <- "/High confidence Markers - Haber/GobletMarkers_Aviv.txt"
eec_mark_hf <- "/High confidence Markers - Haber/EnteroendoMarkers_Aviv.txt"
enteroc_mark_hf <- "/High confidence Markers - Haber/EnterocyteMarkers_Aviv.txt"
tuft_mark_hf <-  "/High confidence Markers - Haber/TuftMarkers_Aviv.txt"

###### Input files ##########

wald <- read.delim(wald_f)
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

###### Remove genes ######

# Remove markers which are not present in wald list
pan_mark1 <- semi_join(pan_mark, wald, by = c("V1" = "target_id"))
gob_mark1 <- semi_join(gob_mark, wald, by = c("V1" = "target_id"))
eec_mark1 <- semi_join(eec_mark, wald, by = c("V1" = "target_id"))
enteroc_mark1 <- semi_join(enteroc_mark, wald, by = c("V1" = "target_id"))
tuft_mark1 <- semi_join(tuft_mark, wald, by = c("V1" = "target_id"))

pan_markh1 <- semi_join(pan_markh, wald, by = c("V1" = "target_id"))
gob_markh1 <- semi_join(gob_markh, wald, by = c("V1" = "target_id"))
eec_markh1 <- semi_join(eec_markh, wald, by = c("V1" = "target_id"))
enteroc_markh1 <- semi_join(enteroc_markh, wald, by = c("V1" = "target_id"))
tuft_markh1 <- semi_join(tuft_markh, wald, by = c("V1" = "target_id"))

###### Output ######

# Output new marker files
write.table(pan_mark1, "/All Markers - Haber/After Wald Filtering/PanethMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(gob_mark1, "/All Markers - Haber/After Wald Filtering/GobletMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(eec_mark1, "/All Markers - Haber/After Wald Filtering/EnteroendoMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(enteroc_mark1, "/All Markers - Haber/After Wald Filtering/EnterocyteMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(tuft_mark1, "/All Markers - Haber/After Wald Filtering/TuftMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)

write.table(pan_markh1, "/High confidence Markers - Haber/After Wald Filtering/PanethMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(gob_markh1, "/High confidence Markers - Haber/After Wald Filtering/GobletMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(eec_markh1, "/High confidence Markers - Haber/After Wald Filtering/EnteroendoMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(enteroc_markh1, "/High confidence Markers - Haber/After Wald Filtering/EnterocyteMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
write.table(tuft_markh1, "/High confidence Markers - Haber/After Wald Filtering/TuftMarkers_Aviv.txt", row.names = F, quote = F, sep = "\t", col.names = F)
