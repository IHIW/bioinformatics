# filename: Run_Validation.sh
# date: 2020-11-17
# version: 1.0
# author: Ben Matern <B.M.Matern@umcutrecht.nl>

# description:
# 	This Script runs the Reference Sequence Generator. This chooses valid reference sequences to be used in HML documents for the 18th Workshop.
#   The --validate flag is set in order to chose the best matching reference for each full-length sequence.

#	Inputs: 
#		an IPD-IMGT/HLA release version.
#	Outputs:
#		X.XX.X_Reference_Alleles.txt (A list of reference alleles)
#		X.XX.X_Reference_Sequences.fasta (The actual reference sequences)
# Options:
#   --validate
#     Validate the reference sequences by aligning full length sequences to the chosen references using BLAST.
#     Output the validation results to ReferenceFinderValidationResults.csv
#     This CSV shows A mapping of this release's full-length sequences to their closest chosen references
#   --supplemental
#     Keep all supplemental files. Including:
#       IMGT/HLA zip and XML files.
#       A fasta containing the full-length allele sequences from this release
#       X.XX.X_Missing_Reference_Alleles.txt (Alleles from the previous reference sequence that were not found in this release)
#   --blast
#     Keep all BLAST alignment XMLS
#     Warning, this generates several gigabytes of output data):
#   --verbose
#     more output text

# 	Prerequisites (other versions likely work):
#		Python 3.6
#		Biopython 1.78
#		NCBI Blast commandline 2.6.0+
#   requests
#   zipfile

# My biopython etc. packages are installed in a virtual environment. Change or delete this if your venv is different.
source ../../venv/bin/activate

# Iterate a list of IPD-IMGT/HLA releases, separated by spaces. Or just chose a single release.
#releases="3.42.0 3.41.0 3.40.0 3.39.0 3.38.0 3.37.0 3.36.0 3.35.0 3.34.0 3.33.0 3.32.0 3.31.0 3.30.0 3.29.0 3.28.0 3.27.0 3.26.0"
releases="3.43.0"
version="1.3"
for release in $releases; do
    echo $release
    python GenerateReferences.py --release=$release --output="../../reference_alleles/"$release"_catalog" --version=$version --supplementary --validate --threads=5
done

# Deactivating my virtual environment
deactivate