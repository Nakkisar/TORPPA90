import pandas as pd
import numpy as np
import sys

output_file = "D:\\BioinfoHons\\PROJECT\\sub_mats\\torppa90_finalissimo\\6KYE_D_nat_vs_pred.txt"
sys.stdout = open(output_file, "w")

# 6KYE_D

native = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
seqpred_predict = "LKLSPEDRGTFKAFLKKLDWEHLGGTAAGAFAAAFPWCLKFHPWWGDIETPEDFMSCPELAAEGREAVKEILYFLAHPDNPEEAFKEECEKEVEKCCVHPELALYLALIMLMAEANKYGEEASPEEIAALMLIYGTIAKGMMKCQC"

print("Native and prediction equal length: ", len(native) == len(seqpred_predict))

match_count = 0
for aa1, aa2 in zip(native, seqpred_predict):
    if aa1 == aa2:
        match_count += 1

similarity_percentage_nat_seqpred = (match_count / len(native)) * 100



torppa90 = pd.read_csv("D:\\BioinfoHons\\PROJECT\sub_mats\\torppa90_finalissimo\\torppa90_final.csv", index_col=0,
                            header=0)
blosum90 = pd.read_csv("D:\\BioinfoHons\\PROJECT\\sub_mats\\blosum90_20x20.csv", index_col=0, header=0)

blosum62 = pd.read_csv("D:\\BioinfoHons\\PROJECT\\sub_mats\\blosum62_20x20.csv", index_col=0, header=0)

print("\n6KYE chain D:\n", "\nNative sequence: ", native, "\nSeqPredNN prediction: ", seqpred_predict)
print(f"Percentage Similarity, native sequence vs SeqPredNN prediction: {similarity_percentage_nat_seqpred:.2f}%")
# print(f"Percentage Similarity, manual prediction vs SeqPredNN prediction: {similarity_percentage_man_seqpred:.2f}%")

'''print("\nTORPPA90: \n", torppa90)
print("\nBLOSUM90: \n", blosum90)
print("\nBLOSUM62: \n", blosum62)'''

torppa90_alignment_score_seqpred = 0
torppa90_alignment_score_native = 0

blosum90_alignment_score_seqpred = 0
blosum90_alignment_score_native = 0

blosum62_alignment_score_seqpred = 0
blosum62_alignment_score_native = 0


for i in range(len(native)):
    amino_acid1 = native[i]
    amino_acid3 = seqpred_predict[i]

    torppa90_score_seqpred = torppa90.loc[amino_acid1, amino_acid3]
    torppa90_score_native = torppa90.loc[amino_acid1, amino_acid1]
    torppa90_alignment_score_seqpred += torppa90_score_seqpred
    torppa90_alignment_score_native += torppa90_score_native

    blosum90_score_seqpred = blosum90.loc[amino_acid1, amino_acid3]
    blosum90_score_native = blosum90.loc[amino_acid1, amino_acid1]
    blosum90_alignment_score_seqpred += blosum90_score_seqpred
    blosum90_alignment_score_native += blosum90_score_native

    blosum62_score_seqpred = blosum62.loc[amino_acid1, amino_acid3]
    blosum62_score_native = blosum62.loc[amino_acid1, amino_acid1]
    blosum62_alignment_score_seqpred += blosum62_score_seqpred
    blosum62_alignment_score_native += blosum62_score_native


print("\nTORPPA90 alignment score, native vs seqpred prediction: ", torppa90_alignment_score_seqpred)
print("TORPPA90 alignment score, native vs native: ", torppa90_alignment_score_native)

print("\nBLOSUM90 alignment score, native vs seqpred prediction: ", blosum90_alignment_score_seqpred)
print("BLOSUM90 alignment score, native vs native: ", blosum90_alignment_score_native)

print("\nBLOSUM62 alignment score, native vs seqpred prediction: ", blosum62_alignment_score_seqpred)
print("BLOSUM62 alignment score, native vs native: ", blosum62_alignment_score_native)