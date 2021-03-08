# Vendor requirements for generation of Genotyping Result documents.
When reporting HLA genotyping for the [18th IHIWS](https://www.ihiw18.org/), it is necessary to be consistent and predictable, to facilitate ease of interpretation and reusability. There are many potential difficulties and ambiguities that may be encountered when creating a genotyping document. The result should be easy to interpret, both by humans and by computers. It is valuable to address these difficulties and find the best approach based on community consensus. This document contains a set of guidelines for vendors on how we recommend to provide genotyping data.

These guidelines are for the community; Your questions, comments, complaints, and suggestions are welcomed and encouraged. Please send them to Ben Matern (B.M.Matern@umcutrecht.nl) or Eric Spierings (E.Spierings@umcutrecht.nl)

## HML
HLA genotyping data should be provided in the [HML 1.0.1 format](schemas.nmdp.org) ([Milius et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4674307/)). This format is described by an [XML Schema](http://schemas.nmdp.org/spec/hml/1.0.1/hml-1.0.1.xsd), and facilitates the transmission of the most important data related to the protocols, sequence, and interpretation involved in creating the genotyping result. The schema is strict in it's data format, but allows a wide variety of data to be included.

HML Documents that are uploaded to the [workshop database](https://data.ihiws.org/) are validated according to some sets of standards.

## Standard Validation Requirements
As a general rule, HML documents should follow these guidelines for validation:

* Documents must conform to the [HML Schema](http://schemas.nmdp.org/spec/hml/1.0.1/hml-1.0.1.xsd).
* Documents should follow the [MIRING](miring.b12x.org) guidelines ([Mack et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4674382/)).
* Genotypes should be reported as valid GLStrings ([Milius et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23849068/)).
* HML documents are validation against the NMDP gateway, which provides further guidelines for high quality data.

## Reference Sequences          
Sequence data should be presented as an alignment against a known sequence. Sequences must be from a subset of full-length alleles chosen from a given release of IPD-IMGT/HLA Database. Sequences chosen for the 18th IHIWS can be found in the [Anthony Nolan HLA Informatics Group Github repository.](https://github.com/ANHIG/IMGTHLA/tree/Latest/ihiw/hml)

Robinson, J, Barker, DJ, Georgiou, X, Cooper, MA, Flicek, P, Marsh, SGE\
The IPD-IMGT/HLA Database\
Nucleic Acids Research (2020) 43:D948-D955

* Submission of known sequences that match known alleles, over the full length of the submitted sequence can be reported as matching the known allele sequence from IPD-IMGT/HLA
* The reference sequence should be chosen based on homology to the submitted sequence. Choosing a sequence based on either exon homology or full-length sequence homology is acceptable.
* Sequences with Novel variants should be reported against a reference from the provided list. Variants should be reported using indices relative to the chosen reference sequence.
* If the sequence is typed as a known allele, but the allele is not on the provided list of reference sequences, allele variants should be reported relative to one of the provided full-length IPD-IMGT/HLA sequences. The interpreted genotype can be provided as a GLString.
* If the submittd sequence is reported outside the region of the reference sequence, such as in extended UTRs, the extra sequence can be reported as an insertion at the beginning or end of the reference.

## IPD-IMGT/HLA Database Releases
HML Documents will specify the release version of the IPD-IMGT/HLA database. The chosen release should be recent, to hopefully provide some consistency in data interpretation.

* If possible, use release version 3.39.0 or later for the submission of HML documents for the 18th IHIWS
	* We will support future IPD-IMGT/HLA releases up to and including the versions released during the Workshop in May 2022.
* Reference sequences are also provided for each IPD-IMGT/HLA release since 3.26.0, mostly for historical and analysis reasons. Versions 3.26.0-3.38.0 are discouraged for the 18th IHIWS.
* Designated reference sequence lists for historical releases of the IPD-IMGT/HLA Database can be found in other branches in the [Anthony Nolan HLA Informatics Group Github repository.](https://github.com/ANHIG/IMGTHLA/tree/Latest/ihiw/hml).

## Potential Issues:
There are many strategies on the correct way to present data and choose reference sequences, and many strategies have merit. There are also potential difficulties and ambiguities.

* How should I report a novel allele in a clinically relevant way? 
	* It is valuable to use a reference sequence that is clinically similar to the reported sequence.
	* The submitter may use a reference sequence, from the provided list of full-length references, that is likely to be clinically similar to the submitted sequence.
		* Likely, this would match in exon sequences.

* Genotyping results are "simpler" if they're reported against the closest matching sequence, regardless of clinical outcomes.
	* It may be more relevant to use a full-length common allele as a reference, and perhaps easier for reporting
	* Ideally,  "variants" are reported against confirmed genotypes in glstring.
	* We would like to know where the variants are in the closest alleles. 
	* A tool is in development to provide the closest-matching sequence to a given consensus sequence. 
		* Vendors are encouraged to use their own strategy for choosing a reference sequence, from the provided list of full-length references.

* It is critical that results can be translated back to IPD-IMGT/HLA.
	* Variants are defined against the reference sequence.
	* This is why it is required to use a subset of selected alleles from a subset of database releases.
	* The use of internal vendor databases is strongly discouraged, as it makes interpretation of sequence results very difficult.


