# filename: Run_GenerateReferences.sh
# date: 2020-11-17
# version: 1.0
# author: Ben Matern <B.M.Matern@umcutrecht.nl>
# description:
# 	This Script runs the Reference Sequence Generator. This chooses valid reference sequences to be used in HML documents for the 18th Workshop.
#	Inputs: 
#		hla.xml (an IPD-IMGT/HLA release)
#		17IHIWS_Reference_Alleles.txt (Previous Reference Sequences , for backwards compatibility)
# 		I manually made input folders containing each release of unzipped IPD-IMGT/HLA input xml files (https://github.com/ANHIG/IMGTHLA/blob/Latest/xml/hla.xml.zip)
# 		For example: /reference_alleles/generate_references/input/3.42.0/hla.xml
#		TODO: Could probably script/automate the downloading and unzipping of hla.xml files from github. What a good idea.
#	Outputs:
#		X.XX.X_Reference_Alleles.txt (A list of reference alleles)
#		X.XX.X_Reference_Sequences.fasta (The actual reference sequences)
#		X.XX.X_Missing_Reference_Alleles.txt (A list of sequences from the previous set that could not be found in the current release)
#		ReferenceFinderValidationResults.csv (A mapping of this release's full-length sequences to their closest chosen references)
# 	Prerequisites (other versions likely work):
#		Python 3.6
#		Biopython 1.78
#		NCBI Blast commandline 2.6.0+

PrevRefSeqs="../../reference_alleles/3.25.0_catalog/17IHIWS_Reference_Alleles.txt"

# My biopython etc. packages are installed in a virtual environment. Change or delete this if your venv is different.
source ../../venv/bin/activate

python GenerateReferences.py --validate --allelelist=$PrevRefSeqs --xml=./input/3.42.0/hla.xml --output=../../reference_alleles/3.42.0_catalog --threads=4
python GenerateReferences.py --validate --allelelist=$PrevRefSeqs --xml=./input/3.41.0/hla.xml --output=../../reference_alleles/3.41.0_catalog --threads=4
python GenerateReferences.py --validate --allelelist=$PrevRefSeqs --xml=./input/3.40.0/hla.xml --output=../../reference_alleles/3.40.0_catalog --threads=4
python GenerateReferences.py --validate --allelelist=$PrevRefSeqs --xml=./input/3.39.0/hla.xml --output=../../reference_alleles/3.39.0_catalog --threads=4
python GenerateReferences.py --validate --allelelist=$PrevRefSeqs --xml=./input/3.38.0/hla.xml --output=../../reference_alleles/3.38.0_catalog --threads=4
python GenerateReferences.py --validate --allelelist=$PrevRefSeqs --xml=./input/3.37.0/hla.xml --output=../../reference_alleles/3.37.0_catalog --threads=4

# Deactivating my virtual environment
deactivate


