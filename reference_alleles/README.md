# Reference alleles for the 18th IHIWS
This repository contains selected reference sequences to be used in creating HLA Genotyping data in NGS HLA typing reports for the 18th International HLA and Immunogenetics Workshop. 

All 18th International HLA and Immunogenetics Workshop database functions used a subset of the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) sequence information.

Robinson, J, Barker, DJ, Georgiou, X, Cooper, MA, Flicek, P, Marsh, SGE\
The IPD-IMGT/HLA Database\
Nucleic Acids Research (2020) 43:D948-D955

These reference sequences were selected using [a set of python scripts](https://github.com/IHIW/bioinformatics/tree/master/reference_alleles/generate_references), which choose the lowest-numbered full-length sequence available in a given IPD-IMGT/HLA release. 

For these purposes, "full length" refers to the 'genomic' sequences curated by the IPD-IMGT/HLA Database, which include sequence for all exons and introns, and the 5' and 3' untranslated regions (UTRs). This reference allele catalog was assembled for use in generating and aligning consensus sequences, and describing novel sequence variation. 

A list of designated reference sequences for each major release, starting with 3.26.0.  Reference sequences will be provided for for each IPD-IMGT/HLA release up to, and soon after the 18th IHIWS.

These allele lists are hosted at the [Anthony Nolan HLA Informatics Group Github repository.](https://github.com/ANHIG/IMGTHLA/tree/Latest/ihiw/hml)

## Loci
Reference alleles are provided for each HLA locus (HLA-A, -B, -C, -E, -F, -G, -DPA1, -DPB1, -DQA1, -DQB1, -DRA, -DRB1, -DRB3, -DRB4 and -DRB5) When available, sequences are also provided for the MICA and MICB loci, as well as HLA-DMA, -DMB, -DOA, -DOB. For each release starting with 3.26.0, this catalog is available in this repository as [A list of alleles, ordered by locus and allele name. A fasta file containing the nucleotide sequences is also available for download.](https://github.com/ANHIG/IMGTHLA/tree/Latest/ihiw/hml) Reference sequence lists for historical releases of the IPD-IMGT/HLA Database can be found in other branches in that repository.

The designated reference sequences are hosted by the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/), and can be downloaded from [Anthony Nolan HLA Informatics Group Github repository.](https://github.com/ANHIG/IMGTHLA/tree/Latest/ihiw/hml) 

## Reference Sequences and HML documents for the 18th IHIWS
More guidelines, including database versions and sequences reference sequences in HML documents can be found in the file in [18IHIWS_Vendor_Genotyping_Requirements](https://github.com/IHIW/bioinformatics/blob/master/reference_alleles/18IHIWS_Vendor_Genotyping_Requirements.md) 

## Usage

###Inputs:
IPD-IMGT/HLA Database release version.

###Outputs:
X.XX.X_Reference_Alleles.txt (A list of reference alleles)\
X.XX.X_Reference_Sequences.fasta (The sequence of each designated reference allele)

###Options/Flags:
--validate\
Validate the reference sequences by aligning full length sequences to the chosen references using BLAST.
Output the validation results to ReferenceFinderValidationResults.csv
This CSV shows A mapping of this release's full-length sequences to their closest chosen references

--supplemental\
Keep all supplemental files. Including:
IMGT/HLA zip and XML files.
A fasta containing the full-length allele sequences from this release
X.XX.X_Missing_Reference_Alleles.txt (Alleles from the previous reference sequence that were not found in this release)

--blast\
Keep all BLAST alignment XMLS
Warning, this generates several gigabytes of output data):
--verbose
more output text

###Prerequisites:
(other versions likely work)\
Python 3.6\
pip\
Biopython 1.78 (biopython)\
NCBI Blast commandline 2.6.0+\
requests\
zipfile