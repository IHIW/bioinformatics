# IHIWS 17th Database
The IHIWS 17th database contains two main parts
1. the web application in https://ihiws17.stanford.edu/ 
2. the sftp server hidpl.stanford.edu

The functions are described in the paper and the scripts in this folder "db_related" could be used to re-build the database with minor changes.

## Web application
The web application was built with oracle SQL database (12c Standard Edition) and APEX 5.0. Before running the scripts, installation a compatible version of oracle SQL database and APEX is required.

After the installation and the successful login as a privileged user into Oracle APEX,

1. Edit the file oracle_sql/ihiws17_schema20171120.mssql so that the strings "/home/oracle/scripts/" are replaced with the folder containing scripts in the folder "scripts".
2. The web app allows users to execute pre-defined shell scripts with functions like merging multiple HML files, converting typing reports to IHIWS XML format with DBMS_SCHEDULER. Please set corresponding privileges for the schema holding the data.
3. Go to SQL Workshop -> SQL scripts -> Upload -> Choose File to upload oracle_sql/ihiws17_schema20171120.mssql and execute it to create the tables, functions, packages, views, procedures and etc used for the web app. 
4. Go to Application Builder -> Import -> Database Application, Page or Component Export to import the file oracle_sql/ihiws_apex_app.mssql and the application should be created.

## data

## scripts

## hlaPoly

Install package:

Change folder to db_related/

In R
```
>install.packages("hlaPoly",repos = NULL, type="source")
```

Usage:
```
>library("hlaPoly")
>hladata<-hlaDataInit(hlaVersion="3.34.0")
>seq="TACCGATAACTAACTGAGTAGTTAATATGGTCAGGCGCTATTCTGAGGATTTACATTTATTAACTCACTTTATTCTCACACATAGTCTTTGAGGTAGGTACTATTATTTTCACTATTTCACATGAGAGATACTTACATCTTTTTACATACACAGAGACTTTAAGCACTTTGATCAAGTTCCCACAGCTATGAAGTAGTAGGGCTAGCTTCCAATCCAGAAAGTCTGGATCCAAGACTGTTTATCCACTGTCCTATTCACCCTATTTTGTGAAGGAAAAGACCAAGTTCAAATTCTCCAGAGTCCATTGCCAAATAATGGAGTCAGATCTATATTTCTATACATAATTACAACACAGTGTGGTGGGTGCCTGTAACTACTTACTGTCTCTACTTGGACTCATTCCATGGCAATGTTCACACAAAAAATGCC"
>hlaGetNovel(hladata,"HLA-DPA1*01:03:01:01",seq)
```

To use a specific version of IMGT/HLA

1. Create a folder with the version name in hlaPoly/inst/hlalib/, e.g. 
```
mkdir hlaPoly/inst/hlalib/3.34.0
```
2. Copy two files: hla_gen.fasta and hla_nuc.fasta from https://github.com/ANHIG/IMGTHLA to the new folder

3. Download hla.xml.zip from https://github.com/ANHIG/IMGTHLA/tree/Latest/xml or from a specific version and unzip it.

```
unzip hla.xml.zip
```
4. Execute scripts/parseHlaxml.pl in the folder with hla.xml and move the output file to the version folder
```
./scripts/parseHlaxml.pl > hlaxml.tsv
mv hlaxml.tsv hlaPoly/inst/hlalib/3.34.0/
```

5. Copy the file hlaRefAllele.tsv from previous version folder to the latest version
```
cp hlaPoly/inst/hlalib/3.30.0/hlaRefAllele.tsv hlaPoly/inst/hlalib/3.34.0/
```

6. Go to the folder of ther version. Execute checkHlaRef.pl in scripts
```
perl ../../../../scripts/checkHlaRef.pl
```
The script will check the hlaRefAllele.tsv and hla_gen.fasta to see if the reference alleles have sequences in hla_gen.fasta.
If a reference doesn't have a sequence in hla_gen.fasta, please edit hlaRefAllele.tsv manually to change the reference allele.

7. In the end, there will be four files in the version folder

hla_nuc.fasta: fasta file with nucleotide sequences for HLA alleles

hla_gen.fasta: fasta file with genomic sequences

hlaxml.tsv: The locations and annotations of HLA gene features

hlaRefAllele.tsv: Indicate which reference allele to use for alleles with only nucleotide sequences 




