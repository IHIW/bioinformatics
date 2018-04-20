# Workshop XML NGS HLA Typing Report Format -- *IHIW XML*

As described by [Chang et al. 2017](https://doi.org/10.1016/j.humimm.2017.12.004), an NGS HLA Typing Report contains both NGS genotyping data and associated meta-data. 
_HLA Genotyping data_ are defined as the HLA genotype and associated consensus sequences, and _meta-data_ are the information that allow the genotyping data to be put in the appropriate reference and experimental context. 

Meta-data include laboratory, report and specimen identifiers; documentation of the instrumentation and software used to generate the genotyping data; reference, alignment and phasing information; the availability and location of primary data (FASTQ files); quality metrics pertinent to the consensus sequences; and annotations providing context for novel sequence variation.
[Mack et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4674382/) provide additional discussion and definition of meta-data and NGS typing report requirements.

##
IHIW XML was developed for the collection of NGS HLA Typing Reports generated as part of the 17th International HLA and 
Immunogenetics Workshop effort. 
The [IHIW XML Typing Report Schema (v1.4) ](https://github.com/IHIW/bioinformatics/blob/master/typing_report_formats/IHIW_XML/IHIW_XML_Schema_v1.4.txt) and associated [IHIW XML XSD](https://github.com/IHIW/bioinformatics/blob/master/typing_report_formats/IHIW_XML/current_ws_xml.xsd) are available in this repository.


