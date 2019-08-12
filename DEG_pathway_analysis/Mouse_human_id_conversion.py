"""
Converts mouse ensembl/uniprot IDs to human ensembl/uniprot IDs

Input:
Tab delimited text file (with headers) where each row is a differentially expressed gene (mouse). One column
 must contain a mouse Ensembl ID or a Uniprot ID. The header name of this column is specified as an input variable.

ID conversion file in a tab delimited text format with headers. Only 1 ID per line for the input ID (long format).
The ensembl ID headers should be 'Human Ensembl' and 'Mouse Ensembl'.
The uniprot ID headers should be 'Human Uniprot' and 'Mouse Uniprot'.
The id type (ensembl or uniprot) for the input and output is specified as a variable.

Output:
Tab delimited text file. Same as input data but with 2 additional columns
 - the to (mouse) and from (human) converted IDs.
 Where multiple human IDs exist for 1 mouse ID, multiple rows will be added to the output file.

"""

# Packages
import pandas as pd
import sys

# Variables
input_f = "../Input_data/Pan_DEGs_lfc1.txt"
inparanoid = "Inparanoid-Mus-hom-UniprotEnsembl-May2018-expanded.txt"
output_f ="Pan_DEGs_lfc1_hum.txt"
id_type = "ensembl" # "ensembl' or 'uniprot' : Used to get the right columns of the conversion file
id_type_input = "ensembl" # "ensembl' or 'uniprot' : Used to get the right columns of the input file file
col_input = "target_id" # In the input file, what is the column header where the id to convert is (base 1)

##### Functions #####

# Function to convert ID
def conversion(conversion,input, mus):
    """
    Joins a dataframe with an ID conversion table, adding the ID conversion columns to the dataframe.

    :param conversion: pandas df with 2 column, first column is mouse id, second is human
    :param input: dataframe of input mouse ids to convert, where the ids are in the column labelled 'col_input'
    :param mus: name of mouse id column in convert file
    :return: Pandas df of mouse and human ids
    """
    # Left join the 2 dataframes by mouse id, put nan where no conversion
    output = input.merge(conversion, left_on=col_input, right_on=mus, how="left")
    return output


##### Import/Prep conversion file #####

# Get columns names for id conversion
if id_type == "ensembl":
    hum = 'Human Ensembl'
elif id_type == "uniprot":
    hum = "Human Uniprot"
else:
    sys.exit("Cannot determine the output id type.")
if id_type_input == "ensembl":
    mus = 'Mouse Ensembl'
elif id_type == "uniprot":
    mus = "Mouse Uniprot"
else:
    sys.exit("Cannot determine the input id type.")

# Open the conversion file and extract just the columns required, based on Id type
convert = pd.DataFrame.from_csv(inparanoid, sep="\t", index_col=None)
convert_f = convert[[mus,hum]]

##### Import/Prep id input file #####

# Open input file
input = pd.DataFrame.from_csv(input_f, sep="\t", index_col=None)

# Run function to get id conversion
output_id = conversion(convert_f, input, mus)

# Save the output file : the input file with added column for human id
output_id.to_csv(output_f, sep = "\t", quoting=None, index=False)

