# Reference alleles for the 18th IHIWS
This repository contains selected reference sequences for the 18th International HLA and Immunogenetics Workshop. All 18th International HLA and Immunogenetics Workshop database functions used [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) sequence information to facilitate the use of this specific database version in NGS HLA Typing Reports..

These reference sequences were selected using [a set of python scripts](https://github.com/IHIW/bioinformatics/tree/master/reference_alleles/generate_references), which choose the lowest-numbered full-length sequence available in a given IPD-IMGT/HLA release. The [sequences selected for the 17th IHIWS](https://github.com/IHIW/bioinformatics/tree/master/reference_alleles/3.25.0_catalog) were included, for backwards compatibility.

For these purposes, "full length" refers to the 'genomic' sequences curated by the IPD-IMGT/HLA Database, which include sequence for all exons and introns, and the 5' and 3' untranslated regions (UTRs). This reference allele catalog was assembled for use in generating and aligning consensus sequences, and describing novel sequence variation. 

A list of Sequences, and a .fasta were provided, for each major release, starting with 3.26.0.  Reference sequences will be provided for for each IPD-IMGT/HLA release up to, and soon after the 18th IHIWS.

## Loci
Reference alleles are provided for each HLA locus (HLA-A, -B, -C, -E, -F, -G, -DPA1, -DPB1, -DQA1, -DQB1, -DRA, -DRB1, -DRB3, -DRB4 and -DRB5) Sequences are also provided for the MICA and MICB loci, as well as HLA-DMA, -DMB, -DOA, -DOB.

For each release tarting with 3.26.0, this catalog is available in this repository in two files: 
- [A list of alleles, ordered by locus and allele name](https://github.com/IHIW/bioinformatics/blob/master/reference_alleles/3.42.0_catalog/3.42.0_Reference_Alleles.txt)
- [A fasta containing the reference sequences](https://github.com/IHIW/bioinformatics/blob/master/reference_alleles/3.42.0_catalog/3.42.0_Reference_Alleles.txt)

## Reference Sequences and HML documents for the 18th IHIWS
More guidelines, including database versions and sequences reference sequences in HML documentscan be found in the file in [18IHIWS_Vendor_Genotyping_Requirements](https://github.com/IHIW/bioinformatics/blob/18WorkshopUpdates/reference_alleles/18IHIWS_Vendor_Genotyping_Requirements.md)

