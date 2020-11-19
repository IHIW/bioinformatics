# filename: README.md
# date: 2020-11-17
# version: 1.0
# author: Ben Matern <B.M.Matern@umcutrecht.nl>
# text: 

https://github.com/IHIW/bioinformatics/tree/master/reference_alleles

This repository contains selected reference sequences for the 18th International HLA and Immunogenetics Workshop. All 18th International HLA and Immunogenetics Workshop database functions used [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) sequence information.

Starting with version 3.37.0, A set of reference HLA alleles for which "full length" consensus sequence was available in IPD-IMGT/HLA Database were selected to facilitate the use of this specific database version in NGS HLA Typing Reports.

Fore these purposes, "full length" refers to the 'genomic' sequences curated by the IPD-IMGT/HLA Database, which include 
sequence for all exons and introns, and the 5' and 3' untranslated regions (UTRs). This reference allele catalog was assembled
for use in generating and aligning consensus sequences, and describing novel sequence variation. 


##
Reference alleles are provided for each HLA locus (HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1, -DRB3, -DRB4 and -DRB5)
, when a full length allele is available and appropriate. Some alleles serve as both locus and allele family references, and some serve as references for
multiple allele families. No allele family reference alleles are provided for the HLA-DPB1, HLA-DQA1 and HLA-DQB1 loci.

This catalog is available in this repository in three files: 
- [A list of alleles, ordered by locus and allele name](https://github.com/IHIW/bioinformatics/blob/master/reference_alleles/3.25.0_catalog/17IHIWS_Reference_Alleles.txt)
- [A list of alleles, ordered by category](https://github.com/IHIW/bioinformatics/blob/master/reference_alleles/3.25.0_catalog/17IHIWS_Reference_Categories.txt) -- locus references are shown first, followed by allelic family references
- [A ReadMe file](https://github.com/IHIW/bioinformatics/blob/master/reference_alleles/3.25.0_catalog/17IHIWS_Reference_README.txt) that provides more detail about the first two files
  



This file is meant to accompany and explain two other files:
17IHIWS_Reference_Alleles.txt is a tab-delimited text file that identifies all 75 17th IHIWS reference alleles by their locus, IMGT/HLA Database version 3.25.0 accession number and allele name, and reference category.
17IHIWS_Reference_Categories.txt is a tab-delimited text file that identifies the reference alleles in each of 99 17th IHIWS reference categories by their locus, IMGT/HLA Database version 3.25.0 accession number and allele name.

Each reference allele has a full-length genomic sequence in IMGT/HLA Database release version 3.25.0. Even though an allele name may include only two or three fields (e.g., HLA-B*73:01 or HLA-B*07:02:01), nucleotide sequence is available for all introns, exons and UTRs for that allele.

Each reference allele can belong to multiple reference categories. Some reference alleles are used by the IMGT/HLA Database as a locus-level reference, and can also serve as an allele-family reference. For example, HLA-A*01:01:01:01 is the IMGT/HLA Database’s reference allele for the HLA-A locus, and can serve as the reference allele for all HLA-A*01 alleles. 

Not all of the locus-level reference alleles used by the IMGT/HLA Database have full-length genomic sequences; in the cases of the HLA-DPA1 and HLA-DPB1 loci, IHIWS Locus Reference alleles with full-length genomic sequences have been identified.

In some cases (e.g., for alleles in the DRB1*10 allele family), no allele-family reference has been identified. In these cases, the locus-level reference should be used for these alleles.

No allele-family reference alleles are identified for the HLA-DPB1, HLA-DQA1 and HLA-DQB1 loci. The locus-level reference should be used for alleles at these loci.

Please use these reference alleles and categories for 17th IHIWS consensus sequence reporting; no other reference alleles or categories should be used to report consensus sequence for the 17th IHIWS.

Version Notes: 

