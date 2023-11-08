import Bio
import pandas as pd
import blosum as bl
import math
from Bio import SubsMat
from Bio import Alphabet
from math import log
import sys

output_file = "entropies_correlations_output.txt"
sys.stdout = open(output_file, "w")

# Import the custom matrix, in this case as an excel file
torppa90_excel = "finalissimo/torppa90_final.xlsx"

# Read the Excel file into a DataFrame
df = pd.read_excel(torppa90_excel, header=0, index_col=0)

# Initialize an empty dictionary to store the data
torppa90_dict = {}

# Iterate through the rows and columns of the DataFrame
for row_label, row in df.iterrows():
    for col_label, value in row.items():
        # Create a 2-tuple key (row_label, col_label) and store the value
        torppa90_dict[(row_label, col_label)] = value

# Now, the custom matrix is in a format usable in Bio.SubsMat in biopython 1.75

# Import the other matrices and transform into dictionaries

blosum62 = bl.BLOSUM(62)
blosum90 = bl.BLOSUM(90)

def matrix_trim(bl_matrix):
    # Remove the last 5 rows and columns as they correspond to ambiguous amino acids B, J, Z, X, * which is not accounted for in torppa90
    # Residues to remove
    residues_to_remove = {'B', 'J', 'Z', 'X', '*'}
    # Create a new dictionary without the specified residues
    matrix_trimmed = {key: value for key, value in bl_matrix.items() if key not in residues_to_remove}
    # Remove the specified residues from the inner dictionaries
    for key, value in matrix_trimmed.items():
        matrix_trimmed[key] = {inner_key: inner_value for inner_key, inner_value in value.items() if inner_key not in residues_to_remove}

    return matrix_trimmed

def dict_generator(matrix):
    # Return an int dictionary that can be easily used to compare each blosum matrix with torppa90
    matrix_dict = {}
    for aa1 in matrix:
        for aa2 in matrix[aa1]:
            matrix[aa1][aa2] = float(matrix[aa1][aa2])
            matrix_dict[(aa1, aa2)] = matrix[aa1][aa2]
    return matrix_dict

blosum62_dict = dict_generator(matrix_trim(blosum62))
blosum90_dict = dict_generator(matrix_trim(blosum90))

# Now we have the blosum matrices in the format usable by biopython 1.75 Bio.SubsMat

'''print("\ntorppa90 matrix dictionary\n")
print(torppa90_dict)

print("\nblosum62 matrix dictionary\n")
print(blosum62_dict)

print("\nblosum90 matrix dictionary\n")
print(blosum90_dict)'''

# Now we can calculate the relative entropies

rel_ent_62 = SubsMat.two_mat_relative_entropy(torppa90_dict, blosum62_dict)
rel_ent_90 = SubsMat.two_mat_relative_entropy(torppa90_dict, blosum90_dict)
rel_ent_bl90vsbl62 = SubsMat.two_mat_relative_entropy(blosum90_dict, blosum62_dict)


# And the linear correlation coefficient 

def two_mat_correlation_edit(mat_1, mat_2):
    """Return linear correlation coefficient between two matrices."""
    try:
        import numpy
    except ImportError:
        raise ImportError("Please install Numerical Python (numpy) if you want to use this function")
    values = []
    # Edit below to put the following code in quotes since our new dictionaries have no "ab_list" attribute
    """assert mat_1.ab_list == mat_2.ab_list"""
    for ab_pair in mat_1:
        try:
            values.append((mat_1[ab_pair], mat_2[ab_pair]))
        except KeyError:
            raise ValueError("%s is not a common key" % ab_pair)
    correlation_matrix = numpy.corrcoef(values, rowvar=0)
    correlation = correlation_matrix[0, 1]
    return correlation

corr_62 = two_mat_correlation_edit(torppa90_dict, blosum62_dict)
corr_90 = two_mat_correlation_edit(torppa90_dict, blosum90_dict)
corr_bl90vsbl62 = two_mat_correlation_edit(blosum90_dict, blosum62_dict)


print("----------------------------------------------------\nRelative entropies:")

print("\nTORPPA90 vs BLOSUM62: ", rel_ent_62)
print("\nTORPPA90 vs BLOSUM90: ", rel_ent_90)
print("\nBLOSUM90 vs BLOSUM62: ", rel_ent_bl90vsbl62)

print("----------------------------------------------------\nLinear correlations:")

print("\nTORPPA90 vs BLOSUM62: ", corr_62)
print("\nTORPPA90 vs BLOSUM90: ", corr_90)
print("\nBLOSUM90 vs BLOSUM62: ", corr_bl90vsbl62)