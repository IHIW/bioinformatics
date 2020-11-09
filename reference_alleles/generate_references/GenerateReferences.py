import argparse
from os.path import join, isdir
from os import makedirs
from datetime import datetime

from ParseXml import parseXmlFile, clusterSequences, parsePreviousReferences

def printSequences(alleleSequences=None, outputFilename=None, verbose=False):
    if(verbose):
        print('Writing ' + str(len(alleleSequences)) + ' allele sequences to:' + str(outputFilename))

    outputFile = open(outputFilename,'w')
    for alleleSequence in alleleSequences:
        outputFile.write('>' + str(alleleSequence.alleleName) + '\n')
        outputFile.write(str(alleleSequence.getSequence()) + '\n')
    outputFile.close()

def printSequenceList(alleleSequences=None, databaseVersion=None, outputDirectory=None, fileVersion='1.0', verbose=False):
    outputFileNameShort = databaseVersion + '_Reference_Alleles.txt'
    outputFileNameFull = join(outputDirectory, outputFileNameShort)
    if(verbose):
        print('Writing a list of ' + str(len(alleleSequences)) + ' allele sequences to:' + str(outputFileNameFull))
    outputFile=open(outputFileNameFull,'w')

    outputFile.write('# filename: ' + str(outputFileNameShort) + '\n')
    outputFile.write('# date: ' + datetime.today().strftime('%Y-%m-%d') + '\n')
    outputFile.write('# version: ' + str(fileVersion) + '\n')
    outputFile.write('# author: ' + str('Ben Matern <B.M.Matern@umcutrecht.nl>') + '\n')
    outputFile.write('IMGT/HLA Database ' + str(databaseVersion) + ' Accession Number\tLocus\tIMGT/HLA Database ' + str(databaseVersion) + ' Allele Name\tDescription\n')

    #Cluster to sort them
    alleleClusters=clusterSequences(alleleSequences)

    for locus in sorted(alleleClusters.keys()):
        #print('Finding reference for locus ' + str(locus))
        for alleleGroup in sorted(alleleClusters[locus].keys()):
            for allele in alleleClusters[locus][alleleGroup]:
                currentLocus, nomenclatureFields = allele.alleleName.split('*')
                nomenclatureTokens = nomenclatureFields.split(':')
                currentGroup = str(nomenclatureTokens[0])
                outputFile.write(str(allele.accessionNumber) + '\t' + currentLocus
                    + '\t' + str(allele.alleleName) + '\t' + currentLocus + '*' + currentGroup + ' Reference\n')

    outputFile.close()



def createReferenceSequences(clusteredFullLenAlleleSequences=None, previousReferenceSequences=None, verbose=False):
    clusteredPreviousReferences = clusterSequences(alleleSequences=previousReferenceSequences)

    print('Creating Reference Sequences for ' + str(len(clusteredFullLenAlleleSequences.keys()))
          + ' loci, previous reference sequences contained ' + str(len(clusteredPreviousReferences.keys()))
          + ' loci.')

    referenceSequences = []

    # Start by looping previous reference sequences
    for locus in sorted(clusteredPreviousReferences.keys()):
        #print('Finding reference for locus ' + str(locus))
        for alleleGroup in sorted(clusteredPreviousReferences[locus].keys()):
            #alleleGroupFull = locus + '*' + alleleGroup
            #print('Finding reference for group ' + alleleGroupFull)
            for previousReferenceSequence in clusteredPreviousReferences[locus][alleleGroup]:
                #print('Finding reference for allele ' + previousReferenceSequence.alleleName)
                # Is there a full-length seq available? Look in the newest list.
                try:
                    matchFound = False
                    for fullLenSequence in clusteredFullLenAlleleSequences[locus][alleleGroup]:
                        if(previousReferenceSequence.alleleName in fullLenSequence.alleleName):
                            #print('Match:' + previousReferenceSequence.alleleName + '/' + fullLenSequence.alleleName)
                            referenceSequences.append(fullLenSequence)
                            matchFound = True
                            break
                except Exception as e:
                    matchFound=False

                if(not matchFound and verbose):
                    print('Warning, Could not find a matching full-length sequence for allele ' + str(previousReferenceSequence.alleleName))


    # Loop full-length sequences from current IPD-IMGT/HLA XML.
    # For each Loci
    for locus in sorted(clusteredFullLenAlleleSequences.keys()):
        #print('Finding reference for locus ' + str(locus))
        for alleleGroup in sorted(clusteredFullLenAlleleSequences[locus].keys()):
            alleleGroupFull = locus + '*' + alleleGroup
            #print('Finding reference for group ' + alleleGroupFull)

            # Do we already have a reference sequence for this group? From Previous Reference Sequences.
            alreadyHaveReference=False
            for newReferenceSequence in referenceSequences:
                if (alleleGroupFull in newReferenceSequence.alleleName):
                    alreadyHaveReference = True
                    break

            if not alreadyHaveReference:
                # Take the first sequence available. Skip DPB1 and add other exceptions if I need it.
                # TODO: Am I handling MICA and MICB right? I don't know actually
                if ( (locus=='HLA-DPB1' and alleleGroup != '02')
                    or (locus=='HLA-MICA' and alleleGroup != '001')
                    or (locus=='HLA-MICB' and alleleGroup != '002') ):
                    if(verbose):
                        print('Not adding a reference sequence for allele group ' + alleleGroupFull)
                else:
                    referenceSequences.append(clusteredFullLenAlleleSequences[locus][alleleGroup][0])
                    if(verbose):
                        print('adding sequence ' + clusteredFullLenAlleleSequences[locus][alleleGroup][0].alleleName
                            + ' as a reference for group ' + str(alleleGroupFull))

    return referenceSequences

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-x", "--xml", required=True, help="path to IPD-IMGT/HLA xml file", type=str)
    parser.add_argument("-o", "--output", required=True, help="Output Directory", type=str)
    parser.add_argument("-a", "--allelelist", required=False, help="Previous Allele List", type=str)

    args = parser.parse_args()
    verbose = args.verbose

    outputDirectory = args.output
    if not isdir(outputDirectory):
        makedirs(outputDirectory)

    print('Generating Reference Sequences')

    if verbose:
        print("Running in verbose mode.")

    # TODO: Get the database version from the xml file
    alleleSequences=parseXmlFile(xmlFile=args.xml,fullLengthOnly=True, verbose=verbose)
    databaseVersion='3.42.0'
    printSequences(alleleSequences=alleleSequences, outputFilename=join(outputDirectory,'FullLengthSequences.fasta'), verbose=verbose)
    alleleSequenceClusters=clusterSequences(alleleSequences=alleleSequences)
    previousReferenceSequences = parsePreviousReferences(referenceSequenceFileName=args.allelelist, verbose=verbose)
    newReferenceSequences = createReferenceSequences(clusteredFullLenAlleleSequences=alleleSequenceClusters, previousReferenceSequences=previousReferenceSequences, verbose=verbose)

    # Print Reference Sequences to fasta.
    printSequences(alleleSequences=newReferenceSequences, outputFilename=join(outputDirectory, 'ReferenceSequences.fasta'), verbose=verbose)

    # Print Reference Sequence List to txt.
    printSequenceList(alleleSequences=newReferenceSequences, databaseVersion=databaseVersion, outputDirectory=outputDirectory, verbose=verbose)

    # Align Reference Sequences
        # For each reference sequence
            # Find corresponding full-length sequences
            # Align them.

    print('Done. Reference Sequences were written to ' + str(outputDirectory))

