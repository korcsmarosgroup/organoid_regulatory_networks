# Organoid regulatory networks

This repository holds the code and necessary input files relating to the following paper:
 
 * Treveil A, Sudhakar P, Matthews Z J, Wrzesinski T, Jones E J, Brooks J, Olbei M, Hautefort I, Hall L J, Carding S R, Mayer U, Powell P P, Wileman T, Di Palma F, Haerty W, Korcsmáros T. Regulatory network analysis of Paneth cell and goblet cell enriched gut organoids using transcriptomics approaches. 

Specifically the contained code is used to carry out scripted analyses of differentially expressed genes using network approaches. These gene lists were obtained from transcriptomics data of Paneth cell enriched, goblet cell enriched and conventionally differentiated enteroids.

This repository is organised into analysis stages, which can be run separately from each other. For each stage, there is a directory containing the relevant scripts - details of which are outlined below.  All the necessary input files are given in the *Input_data/* directory. Please note that some manual/non-scripted analysis was carried out alongside these scripts, for example the analysis of most rewired regulators was carried out using a Cytoscape app called DyNet.

## Pathway analysis of differentially expressed genes

This stage takes the two lists of differentially expressed genes (DEGs - Paneth enriched vs conventionally differentiated expressed enteroids, goblet enriched vs conventionally differentiated enteroids), converts the gene IDs to human (from mouse) and carries out pathway enrichment using Reactome, KEGG and Gene Ontology. 

#### Output
The scripts output a table of results and a plot of the top 10 Reactome results for each DEG list.

#### Code
The code is found in the directory *DEG_pathway_analysis/*. The script *mouse_human_id_conversion.py* should be run first.


## Differentially expressed marker genes

This stage identifies the number of cell type specific markers (from Haber et al., Nature, 2017) in each of the DEG lists. First the marker lists (previously converted to Ensembl IDs) are filtered for only those which are present in the list of genes with variance greater than zero among samples (Wald file output by analysis upstream of DEGs). This is done to ensure that the background for each of the datasets is the same, enabling hypergeometric significance testing. Then the markers are identified in the input DEG lists and the hypergeometric significant test is applied to determine if the markers are significantly overrepresented in the DEG lists. 

#### Output
The scripts output a table with the results of all the overlaps, as well as another text file with the results of the hypergeometric significance tests. NB. miRNAs and potential novel lncRNAs are removed from this analysis, as they are not present in the marker gene lists.

#### Code
The code is found in the directory *DEG_marker_overlap/*. The script *marker_filtering_wald.R* should be run first.


## Generating regulatory networks

This stage identifies possible regulatory connections between the DEGs in the Paneth and goblet conditions using an input set of interactions. These 'background' interactions were collated previously from open source databases as described in the manuscript (Treveil A, Regulatory network analysis of Paneth cell and goblet cell enriched gut organoids using transcriptomics approaches). 

#### Output
The output of the stage is two regulatory networks (one for each DEG dataset) in a tab delimited text format where each row represents an interaction.

#### Code
The code is found in the directory *Creating_regulatory_networks/*.

## Cluster analysis

Here the nodes of network clusters are analysed for enriched pathway/functional associations. Prior to running this code, the Paneth and goblet networks are imported into Cytoscape and the MCode app is used to identify clusters of highly interconnected nodes using default settings. The results of this analysis are saved as text files named *Mcode_Paneth.txt* and *Mcode_Goblet.txt*.

#### Output

The script outputs a dataframe of enriched Reactome and KEGG pathways (q val <= 0.05) for each cluster in the input datasets, as well as heatplots saved as PDF files for the Reactome and KEGG results.

#### Code
The code is found in the directory *Cluster_analysis/*.

## Filtering networks for markers and their regulators

This stage takes the networks output by the *Generating regulatory networks* stage and filters them for cell type specific markers and their regulators. NB. Only the markers also present in the Wald file are included here to standardise the analysis with the **Differentially expressed marker genes** stage.

#### Output
The output of the stage is two regulatory networks (one for each DEG dataset) in a tab delimited text format where each row represents an interaction. All target nodes are markers of the relevant cell type based on the markers from Haber el al. All nodes are DEGs in the relevant dataset.

Additionally, the code creates heatplots (and associated data tables) indicating the regulator-marker interactions and the overlap of regulators between the Paneth cell and goblet cell datasets.

#### Code
The code is found in the directory *Regulator_marker_network/*. The script *filter_network_target.py* should be run first.

## Pathway analysis of regulators

This stage carries out pathway enrichment analysis on the regulators of the cell type specific markers using the output networks of the **Filtering networks for markers and their regulators** stage. The regulators are categorised depending on their presence in only one or both of the networks and then converted to human genes IDs to allow for functional analysis.

#### Output

This script outputs tab delimited text files containing the significant Reactome pathways and PDF images of dot plot output from KEGG and Reactome analysis of regulators split by specificity. It will also output a text file with the human and mouse IDs for each regulator.

#### Code
The code is found in the directory *Regulator_pathway_analysis/*. 

## Pathway analysis of rewired regulators

This stage does the following: 

  a) takes a list of most rewired regulators (determined using Cytoscape app DyNet)
  
  b) identifies their regulatory targets in the networks output from the **Generating regulatory networks** stage
  
  c) carries out pathway enrichment analysis on the targets (grouped by which networks the targets are found in) - using Reactome, KEGG and Gene Ontology.

#### Output
The script outputs a tab delimited text file with the results for the Reactome, KEGG and GO enrichment analysis.

#### Code
The code is found in the directory *Rewiring_pathway_analysis/*.

## Analysis of Crohn's disease and Ulcerative Colitis associated genes

This stage gets the number of Crohn's / UC associated genes present in the DEG lists, the regulatory networks and the regulator-marker networks. It requires the output from the stages **Generating regulatory networks** and **Filtering networks for markers and their regulators**. I then calculates the hypergeometric significance statistic to determine the significance of the enrichment of Crohn's genes in the regulatory networks.

The Crohn's / UC associated genes are taken from the papers Jostins et al., 2012 and Farh et al., 2015.

#### Output
The script outputs tab delimited text files listing the Crohn's / UC genes identified in the DEG lists/networks. It prints the results of the hypergeometric significance test.

#### Code
The code is found in the directory *Crohns_uc_genes_overlap/*.

## Dependencies

R (>= 3.4.0) with the following packages:
* clusterProfiler (>= 3.4.0)
* dendsort (>= 0.3.3)
* dplyr (>= 3.2.0)
* forcats (>= 0.4.0)
* ggplot2 (>= 3.2.0)
* grid (NA)
* pheatmap (>= 1.0.12)
* reactomePA (>= 1.28.0)
* reshape2 (>= 1.4.3)
* tidyr (>= 0.8.3)

Python (>= 3) with the following packages:
* Pandas (>= 0.24.2)

## Authors

All code in this repository was written and applied by Agatha Treveil in 2018/9 - [agathat](https://https://github.com/agathat)

## Acknowledgments

* The other authors of the paper and members of the Korcsmaros lab for their support and contributions to the project. Particularly to Padhmanand Sudhakar and Tomasz Wrezesinski for the collation of background interactions and analysis of raw sequencing data used in this repo.
* Cytoscape - Shannon P et al., Genome Research, 2003
* Mcode - Bader GD and Hogue CW, BMC Bioinformatics, 2003
* DyNet - Goenawan IH, Bryan K and Lynn DJ, Bioinformatics, 2016
* Marker data - Haber AL et al., Nature, 2017
* Crohn's / UC associated gene lists - Jostins L et al., Nature. 2012 and Farh et al., Nature. 2015.
* Mouse-human ID conversion - InParanoid resource: O’Brien KP, Remm M and Sonnhammer ELL., Nucleic Acids Res. 2005.