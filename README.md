# TORPPA90
Substitution matrix based on SeqPredNN (https://github.com/falategan/SeqPredNN)

TORPPA90 is a substitution matrix applicable for use in BLASTP searches. Values in TORPPA90, should be transposed into the table files of substitution matrices within BLAST+ source code in order to be used. An already transformed "blosum90" table is provided in the BLAST+ directory

This substitution matrix is asymmetrical and reflects inequality of forward and reverse substitutions with respect to the maintenance of local 3D structures within protein sequences.

BLAST+ pdbaa searches were performed using the preformatted database avalaible here: https://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz
