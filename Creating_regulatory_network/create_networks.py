"""
Add regulatory connections between the differentially expressed genes of the Paneth and goblet datasets using 
input known interaction network.

Inputs:
    Known interactions network - tab delimited text file with headers. Columns 'source' and 'target'
    contain source and target nodes (only 1 id per node).
    
    Table differentially expressed genes for the Paneth and goblet datasets. Tab delimited text files where
    the id for the gene is in the column with header 'target_id'.
     
Output:
    Saves 2 Tab delimited text files - one for Paneth and one for goblet - with the filtered interactions.
    Same columns and headers as the input network.

"""

import pandas as pd

# Known interactions
network_path = "Known_interactions_unfiltered.txt"
# Tables of DEGs
pan_degs_path = "Pan_DEGs_lfc1.txt"
gob_degs_path = "Gob_DEGs_lfc1.txt"
# Ouput file paths
pan_out_path = "Pan_network_lfc1.txt"
gob_out_path = "Gob_network_lfc1.txt"

def filter(deg_path, network_path):
    """ Function to filter input network for interactions where both nodes are in first column of input dataframe.

    :param deg_path: Path (string) to the file containing the differentially expressed genes (DEGs). The file is tab
                     delimited with the gene id in the column with header 'target_id'.
    :param network_path: Path (string) to file containing the interactions to filter. The source nodes are in column
                    with header 'source' and the target nodes in column with header 'target'.
    :return net_filt: Pandas dataframe containing the filtered interactions with same columns and headers as
                    the input network.
    """

    degs = pd.read_csv(deg_path, delimiter = "\t", skip_blank_lines=True)
    net = pd.read_csv(network_path, delimiter = "\t", skip_blank_lines=True)

    # Create filters for the soruce and target nodes
    filt1 = net.source.isin(degs.target_id)
    filt2 = net.target.isin(degs.target_id)

    # Apply both filters to the input network
    net_filt = net[filt1 & filt2]

    return net_filt

# Create network using Paneth DEGs and save to file
pan_network = filter(pan_degs_path, network_path)
pan_network.to_csv(pan_out_path, index=False, sep = "\t")

# Create network using goblet DEGs and save to file
gob_network = filter(gob_degs_path, network_path)
gob_network.to_csv(gob_out_path, index=False, sep = "\t")
