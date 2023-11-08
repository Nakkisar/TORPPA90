import pandas as pd
import numpy as np
import sys

output_file = "D:\\BioinfoHons\\PROJECT\\sub_mats\\torppa90_finalissimo\\matrix_maker_output.txt"
sys.stdout = open(output_file, "w")

# Read the long-format CSV file with the unnormalized true and predicted residue pair values
df = pd.read_csv("D:\\BioinfoHons\\PROJECT\SeqPredNN\\models\\oversamples_20_09_2023\\200_20\\unnormalised_conf_matrix.csv")

# Now we need to pivot the long format df into a wide format
raw_conf = df.pivot(index=' True residue', columns='Predicted residue', values=' Value')

# The resulting raw_conf is a 20-by-20 square DataFrame with True residues as rows and Predicted residues as
# columns
print("Raw confusion matrix:\n", raw_conf)

# Next we can generate a matrix showing the expected frequency for the amino acid true-predicted pairs
# Make a blank copy of true_norm_conf, this will become the expected frequency matrix.
exp_freq = raw_conf.copy()
exp_freq[:] = 0

# Loop through every value in true_norm_conf and calculate the expected frequency value for each of the pairs
# Iterate through rows and columns to calculate expected frequencies
for true_residue in raw_conf.index:
    for predicted_residue in raw_conf.columns:
        true_count = raw_conf.loc[true_residue, predicted_residue]
        row_sum = raw_conf.loc[true_residue, :].sum()
        col_sum = raw_conf.loc[:, predicted_residue].sum()
        total_sum = raw_conf.values.sum()
        expected_frequency = (true_count/row_sum)/(col_sum/total_sum)
        exp_freq.loc[true_residue, predicted_residue] = expected_frequency

# The "exp_freq" DataFrame now contains the expected frequencies for each pair of amino acids
print("Expected frequencies:\n", exp_freq)


# Now we can make the dataframe for the log-odds matrix
log_odds_matrix = raw_conf.copy()
log_odds_matrix[:] = 0

# We'll use a small positive constant "epsilon" to avoid the possibility of dividing by zero
epsilon = 1e-10

# Calculate the log-odds matrix, the same expected frequency calculation is used here as from earlier
for true_residue in raw_conf.index:
    for predicted_residue in raw_conf.columns:
        # The number of times a certain residue was predicted (columns) instead of the true residue present (rows)
        true_count = raw_conf.loc[true_residue, predicted_residue]
        # Sum of the row associated with that cell
        row_sum = raw_conf.loc[true_residue, :].sum()
        # Sum of the column associated with that cell
        col_sum = raw_conf.loc[:, predicted_residue].sum()
        # Sum of all values in the raw conf matrix
        total_sum = raw_conf.values.sum()

        # Calculate the expected frequency
        expected_frequency = (true_count / row_sum) / (col_sum / total_sum)

        # Calculate the log-odds score, add scaling factor of 2
        log_odds_score = 2 * np.log2(expected_frequency + epsilon)
        log_odds_matrix.loc[true_residue, predicted_residue] = log_odds_score


# Round the log-odds matrix to 2 decimal places
torppa90 = np.round(log_odds_matrix, 2)

# The rounded_log_odds_matrix is now log-odds scoring matrix with the original BLOSUM62 scaling factor
print("Initial unrounded log_odds_matrix:\n", log_odds_matrix)

# We would like the single letter abbreviations for the amino acid labels
# Define a mapping from full names to single-letter abbreviations
name_to_abbr = {
    'CYS': 'C', 'SER': 'S', 'THR': 'T', 'ALA': 'A', 'GLY': 'G',
    'PRO': 'P', 'ASP': 'D', 'GLU': 'E', 'GLN': 'Q', 'ASN': 'N',
    'HIS': 'H', 'ARG': 'R', 'LYS': 'K', 'MET': 'M', 'ILE': 'I',
    'LEU': 'L', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y', 'PHE': 'F'
}

# Rename the row and column labels using the mapping
torppa90.rename(index=name_to_abbr, columns=name_to_abbr, inplace=True)
torppa90.index.name = None
torppa90.columns.name = None

# The log-odds matrix now has single-letter abbreviations as labels, and is rounded to 2 decimal places to yield TORPPA90
print("TORPPA90:\n", torppa90)

# Print the stats for torppa90 (interesting, and useful for normalizing it next)
torppa90_range_value = torppa90.max().max() - torppa90.min().min()
torppa90_max_value = torppa90.max().max()
torppa90_min_value = torppa90.min().min()
print("torppa90 \nRange: ", torppa90_range_value, "\nMax: ", torppa90_max_value, " Min: ", torppa90_min_value)
torppa90_median_value = torppa90.stack().median()
torppa90_mean_value = torppa90.mean().mean()
print("Median: ", torppa90_median_value)
print("Mean: ", torppa90_mean_value)

# Write out the csv and txt files for torppa90
torppa90.to_csv("D:\\BioinfoHons\\PROJECT\\sub_mats\\torppa90_finalissimo\\torppa90_final.csv")
torppa90.to_csv("D:\\BioinfoHons\\PROJECT\\sub_mats\\torppa90_finalissimo\\torppa90_final.txt", sep=',', index=False)

# Generate an int version of the matrix, write them out too
torppa90_int = np.round(torppa90, 0).astype(int)
torppa90_int.to_csv("D:\\BioinfoHons\\PROJECT\\sub_mats\\torppa90_finalissimo\\torppa90_final_int.csv")
torppa90_int.to_csv("D:\\BioinfoHons\\PROJECT\\sub_mats\\torppa90_finalissimo\\torppa90_final_int.txt", sep=',', index=False)