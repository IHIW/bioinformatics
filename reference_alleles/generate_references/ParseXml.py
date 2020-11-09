import xml.etree.ElementTree as ElementTree
from os.path import isfile

from AlleleSequence import AlleleSequence


def parseXmlFile(xmlFile=None, fullLengthOnly=True, verbose=False):
    print('Parsing this IPD-IMGT/HLA input file:' + str(xmlFile))
    if(xmlFile is None or not isfile(xmlFile)):
        raise Exception ('Specify the input IPD-IMGT/HLA input file:' + str(xmlFile))

    alleleSequences=[]
    documentRoot = ElementTree.parse(xmlFile).getroot()
    if(verbose):
        print(str(len(documentRoot)) + ' allele nodes found.')
    for index, alleleNode in enumerate(documentRoot):
        currentAllele = AlleleSequence()
        currentAllele.alleleName = alleleNode.get('name')
        currentAllele.accessionNumber =alleleNode.get('id')
        if(currentAllele.accessionNumber is None or len(currentAllele.accessionNumber)<1):
            raise Exception('No id/accession number found for allele ' + str(currentAllele.alleleName))

        sequenceList = alleleNode.findall('{http://hla.alleles.org/xml}sequence')

        # Should be exactly one sequence node for each allele.
        if (len(sequenceList) == 0):
            if(verbose):
                print('Warning, Zero sequences for the allele ' + str(currentAllele.alleleName) + ', this allele was probably deleted or renamed.')
        elif (len(sequenceList) > 1):
            raise Exception('Oops! More than one sequence node found for an allele.' + str(currentAllele.alleleName))

        for sequenceNode in sequenceList:
            nucSequenceList = sequenceNode.findall('{http://hla.alleles.org/xml}nucsequence')

            # One nucleotide sequence per sequence node.
            if (len(nucSequenceList) != 1):
                raise Exception('Error: More than one nucsequence node on ' + str(currentAllele.alleleName) + ' found: ' + str(len(nucSequenceList)))

            fullSequence = nucSequenceList[0].text

            # Store each feature Sequence
            featureList = sequenceNode.findall('{http://hla.alleles.org/xml}feature')
            for featureNode in featureList:
                featureName = featureNode.get('name').replace('\'','').replace(' ','')
                featureType = featureNode.get('featuretype')

                if (featureType in ['Intron', 'Exon', 'UTR']):
                    coordinates = featureNode.findall('{http://hla.alleles.org/xml}SequenceCoordinates')
                    if (len(coordinates) != 1):
                        raise Exception('Error: More than one coordinates node on ' + str(currentAllele.alleleName) + ' found: ' + str(len(coordinates)))

                    coordinate=coordinates[0]
                    # IPD-IMGT/HLA XML uses 1-based indexing
                    featureStart = int(coordinate.get('start')) - 1
                    featureEnd = int(coordinate.get('end'))
                    featureSequence = fullSequence[featureStart:featureEnd]
                    currentAllele.featureSequences[featureName] = featureSequence
                elif(featureType=='Protein'):
                    # TODO: can parse protein sequence if I want to...
                    pass
                else:
                    raise Exception('I do not understand this feature type for allele ' + str(currentAllele.alleleName) + ' :' + str(featureType))

            # Sanity check
            if(fullSequence != currentAllele.getSequence()):
                if(verbose):
                    print('Warning, Full sequence from XML file does not match sequence from annotated features, probably has un-annotated UTR sequence:'
                          + str(currentAllele.alleleName) + ', FeatureKeys=' + str(list(currentAllele.featureSequences.keys())))
                    #print('Full Length from XML :' + str(fullSequence))
                    #print('Assembled Full Length:' + str(currentAllele.getSequence()))

        if(not fullLengthOnly or currentAllele.isFullLength()):
            alleleSequences.append(currentAllele)

    return alleleSequences

def clusterSequences(alleleSequences=None):
    print('Clustering ' + str(len(alleleSequences)) + ' Allele Sequences by Locus and Group.')
    # Create Dictionary clusteredSequences[locus][group]
    # containing list of allele sequences
    clusteredSequences = {}
    for alleleSequence in alleleSequences:
        alleleName = str(alleleSequence.alleleName)
        currentLocus, nomenclatureFields = alleleName.split('*')

        nomenclatureTokens = nomenclatureFields.split(':')
        currentGroup = str(nomenclatureTokens[0])

        if(currentLocus not in clusteredSequences.keys()):
            clusteredSequences[currentLocus] = {}
        if (currentGroup not in clusteredSequences[currentLocus].keys()):
            clusteredSequences[currentLocus][currentGroup] = []

        clusteredSequences[currentLocus][currentGroup].append(alleleSequence)

    return clusteredSequences

def parsePreviousReferences(referenceSequenceFileName=None, verbose=False):
    if(referenceSequenceFileName is None or not isfile(referenceSequenceFileName)):
        print('No reference sequence filename was provided.')
        return []
    else:
        print('Parsing this file for previous reference sequence names:' + str(referenceSequenceFileName))
        referenceSequenceNames = []

        referenceSequenceFile = open(referenceSequenceFileName, 'r')
        textLines = referenceSequenceFile.readlines()
        for textLine in textLines:
            if(textLine.startswith('#')):
                pass
            elif(textLine.startswith('IMGT/HLA Database')):
                pass
            elif(textLine.strip().replace('\n','')==''):
                pass
            else:
                accession, locus, alleleName, description = textLine.strip().split('\t')
                referenceAllele = AlleleSequence()
                referenceAllele.accession=accession
                referenceAllele.alleleName=alleleName

                referenceSequenceNames.append(referenceAllele)

        if(verbose):
            print('I found ' + str(len(referenceSequenceNames)) + ' reference sequences in file ' + str(referenceSequenceFileName))

        return referenceSequenceNames