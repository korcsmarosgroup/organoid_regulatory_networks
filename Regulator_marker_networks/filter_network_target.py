"""
Filter a network for interactions where the target node is in an input list of gene ids eg. marker list.

Inputs:
    Input network - tab delimited text file with headers. Columns 'source' and 'target'
    contain source and target nodes (only 1 id per node).

    Text file containing list of gene ids to filter for. No header.

Output:
    Saved a tab delimited text file with the filtered interactions. Same columns and headers as the input network.

"""

import pandas as pd

# Input network
network_path = "Gob_network_lfc1.txt"
# List to filter with
list_path = "GobletMarkers_Aviv.txt"
# Ouput file path
out_path = "Gob_network_gob_markers.txt"

def filter(list_path, network_path):
    """ Function to filter input network for interactions where the target nodes is in the first column of input data.

    :param list_path: Path (string) to txt file containing list of gene ids to filter for. No header.
    :param network_path: Path (string) to file containing the interactions to filter. The source nodes are in column
                    with header 'source' and the target nodes in column with header 'target'.
    :return net_filt: Pandas dataframe containing the filtered interactions with same columns and headers as
                    the input network.
    """

    list = pd.read_csv(list_path, skip_blank_lines=True, header=None)
    net = pd.read_csv(network_path, delimiter = "\t", skip_blank_lines=True)

    print(list)
    # Create filters for the target nodes
    filt = net.target.isin(list[0])

    # Apply filter to the input network
    net_filt = net[filt]

    return net_filt

# Filter network using input list and save to file
filt_network = filter(list_path, network_path)
filt_network.to_csv(out_path, index=False, sep = "\t")
