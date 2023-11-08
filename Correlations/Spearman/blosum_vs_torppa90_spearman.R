library(gt)
library(gtExtras)
library(dplyr)
library(scales)
library(tibble)
library(webshot2)

setwd("D:/BioinfoHons/PROJECT/sub_mats")

# read in our new torppa90 substitution matrix

torppa90 <- read.csv2("torppa90_finalissimo/torppa90_final.csv", sep = ",", row.names = 1 )
torppa90[] <- lapply(torppa90, as.numeric)

# for each of the existing matrices which we load in for the comparison, the last residues in the respective matrices (which don't represent any specific amino acid, but rather ambiduities that may appear in the aligned sequences) are trimmed off to get BLOSUM matrices similar in dimension to our new matrix.
# these include: "B", "J", "Z", "X", "*", yielding a 20-by-20 matrix with values for each amino acid substitution pair. The already trimmed BLOSUM matrices

blosum62 <- read.csv2("blosum62_20x20.csv", sep = ",", row.names = 1)
blosum62[] <- lapply(blosum62, as.numeric)

blosum90 <- read.csv2("blosum90_20x20.csv", sep = ",", row.names = 1)
blosum90[] <- lapply(blosum90, as.numeric)

#View(sub_matrix)
#View(blosum62)
#View(blosum90)

# Generate the correlation matrices comparing blosum and torppa90 matrices
correlation_matrix_62 <- cor(torppa90, blosum62, method = "spearman")
correlation_matrix_90 <- cor(torppa90, blosum90, method = "spearman")

print(correlation_matrix_90)

#View(correlation_matrix_62)
#View(correlation_matrix_62_norm)

correlation_figure <- function(correlation_matrix, mat_name, file_name) {
  
  # Visualise correlation as a gt table, coloured by value
  # Convert the correlation matrix to a data frame
  correlation_df <- as.data.frame(correlation_matrix)
  
  
  # Visualize the correlation matrix as a gt table
  correlation_values_table <- correlation_df %>%
    rownames_to_column(var = "Amino Acid") %>%
    gt() %>%
    data_color(
      columns = 2:21,
      colors = scales::col_numeric(
        palette = c("red", "white", "green"),
        domain = c(-1, 1)
      )
    ) %>%
    fmt_number(columns = 2:21, decimals = 1) %>%
    cols_align(align = 'center', columns = everything()) %>%
    tab_header(
      title = sprintf('%s vs TORPPA90', mat_name),
      subtitle = 'Spearman Correlation'
    )
  
  # Print the table
  correlation_values_table
  
  # save table for report
  correlation_values_table %>% 
    gtsave_extra(filename=file_name, path="D:/BioinfoHons/PROJECT/sub_mats/torppa90_finalissimo/spearman_correlations")
}

# Generate and save the figures
correlation_figure(correlation_matrix_62, "BLOSUM62", "blosum62_vs_torppa90_spearman.png")
correlation_figure(correlation_matrix_90, "BLOSUM90", "blosum90_vs_torppa90_spearman.png")
