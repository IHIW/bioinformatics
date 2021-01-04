import argparse
from os.path import join, isdir, isfile
from os import makedirs, remove, rmdir
from datetime import datetime
from requests import get
from zipfile import ZipFile
from random import shuffle

from ParseXml import parseXmlFile, clusterSequences, parsePreviousReferences
from FindBestReference import validateFullLengthSequences

def printSequences(alleleSequences=None, outputFilename=None, verbose=False):
    if(verbose):
        print('Writing ' + str(len(alleleSequences)) + ' allele sequences to:' + str(outputFilename))

    outputFile = open(outputFilename,'w')
    alleleClusters = clusterSequences(alleleSequences=alleleSequences, verbose=verbose)

    for locus in sorted(alleleClusters.keys()):
        # print('Finding reference for locus ' + str(locus))
        for alleleGroup in sorted(alleleClusters[locus].keys()):
            for alleleSequence in alleleClusters[locus][alleleGroup]:
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

    locusReferences = getLocusReferences()

    #Cluster to sort them
    alleleClusters=clusterSequences(alleleSequences=alleleSequences, verbose=verbose)

    for locus in sorted(alleleClusters.keys()):
        #print('Finding reference for locus ' + str(locus))
        for alleleGroup in sorted(alleleClusters[locus].keys()):
            for allele in alleleClusters[locus][alleleGroup]:
                currentLocus, nomenclatureFields = allele.alleleName.split('*')
                nomenclatureTokens = nomenclatureFields.split(':')
                currentGroup = str(nomenclatureTokens[0])

                allele.description = ''
                # Is it a locus reference?
                for locusReference in list(set(locusReferences)):
                    if(allele.alleleName in locusReference or locusReference in allele.alleleName):
                        allele.description = currentLocus + ' Locus Reference;'

                #if(allele.description is None or len(allele.description) < 1):
                #    allele.description = currentLocus.replace('HLA-','') + '*' + currentGroup + ' Reference'
                allele.description += currentLocus.replace('HLA-', '') + '*' + currentGroup + ' Reference'



                outputFile.write(str(allele.accessionNumber) + '\t' + currentLocus
                    + '\t' + str(allele.alleleName) + '\t' + str(allele.description) + '\n')

    outputFile.close()

def printMissingSequences(missingSequences=None, databaseVersion=None, outputDirectory=None, fileVersion='1.0', verbose=False):
    outputFileNameShort = databaseVersion + '_Missing_Reference_Alleles.txt'
    outputFileNameFull = join(outputDirectory, outputFileNameShort)
    if(verbose):
        print('Writing a list of ' + str(len(missingSequences)) + ' missing allele sequences to:' + str(outputFileNameFull))
    outputFile=open(outputFileNameFull,'w')

    outputFile.write('# filename: ' + str(outputFileNameShort) + '\n')
    outputFile.write('# date: ' + datetime.today().strftime('%Y-%m-%d') + '\n')
    outputFile.write('# version: ' + str(fileVersion) + '\n')
    outputFile.write('# author: ' + str('Ben Matern <B.M.Matern@umcutrecht.nl>') + '\n')

    alleleClusters=clusterSequences(alleleSequences=missingSequences,verbose=verbose)

    for locus in sorted(alleleClusters.keys()):
        for alleleGroup in sorted(alleleClusters[locus].keys()):
            for allele in alleleClusters[locus][alleleGroup]:
                outputFile.write(str(allele.alleleName) + '\n')

    outputFile.close()

def createReferenceSequences(clusteredFullLenAlleleSequences=None, verbose=False, imgtReleaseVersion=None):
    print('Creating Reference Sequences for ' + str(len(clusteredFullLenAlleleSequences.keys()))
          + ' loci.')

    referenceSequences = []
    missingSequences = []
    locusReferences=getLocusReferences()

    # Add hard-coded alleles
    includeSequences, alleleDescriptionLookup = getIncludeSequenceList(imgtReleaseVersion=imgtReleaseVersion)
    for locus in sorted(clusteredFullLenAlleleSequences.keys()):
        for alleleGroup in sorted(clusteredFullLenAlleleSequences[locus].keys()):
            # print('Finding reference for group ' + alleleGroupFull)
            for fullLengthSequence in clusteredFullLenAlleleSequences[locus][alleleGroup]:
                # print('Finding reference for allele ' + previousReferenceSequence.alleleName)
                # Skip the sequence if it's in the excludeSequenceList
                if (fullLengthSequence.alleleName in (includeSequences)):
                    referenceSequences.append(fullLengthSequence)
                    if verbose:
                        print('Forced adding hard-coded reference:' + str(fullLengthSequence.alleleName))

    # Add locus references if they don't already exist
    for locus in sorted(clusteredFullLenAlleleSequences.keys()):
        for alleleGroup in sorted(clusteredFullLenAlleleSequences[locus].keys()):
            for allele in clusteredFullLenAlleleSequences[locus][alleleGroup]:
                for locusReference in list(set(locusReferences)):
                    if(allele.alleleName in locusReference or locusReference in allele.alleleName):
                        #print('adding locus reference:' + str(locusReference))
                        # This reference be already in the list...
                        addReference = True
                        for referenceSequence in referenceSequences:
                            if(referenceSequence.alleleName in locusReference or locusReference in referenceSequence.alleleName):
                                #print('Skipping adding locus reference:' + str(locusReference) + ' because it is already in the list.')
                                addReference=False
                                break
                        if(addReference):
                            #print('indeed, adding locus reference:' + str(locusReference))
                            referenceSequences.append(clusteredFullLenAlleleSequences[locus][alleleGroup][0])

    # Add groupwise seq references
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
                if ( (locus=='HLA-DPB1' and alleleGroup != '01')
                    or (locus in ['HLA-MICA','MICA'] and alleleGroup != '001')
                    or (locus in ['HLA-MICB','MICB'] and alleleGroup != '004') ):
                    if(verbose):
                        print('Not adding a reference sequence for allele group ' + alleleGroupFull)
                else:
                    referenceSequences.append(clusteredFullLenAlleleSequences[locus][alleleGroup][0])
                    if(verbose):
                        print('adding sequence ' + clusteredFullLenAlleleSequences[locus][alleleGroup][0].alleleName
                            + ' as a reference for group ' + str(alleleGroupFull))

    return referenceSequences, missingSequences

def downloadImgtXml(outputDirectory=None, release=None, verbose=False):
    print('Downloading IPD-IMGT/HLA xml file for release ' + str(release))

    zipLocalFileName = join(outputDirectory, 'hla.xml.zip')
    xmlLocalFileName = join(outputDirectory, 'hla.xml')
    zipRemoteFileName = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/' + release.replace('.','') + '/xml/hla.xml.zip'

    # This is a hack. Because the 3.35.0 release was uploaded to github using lfs, the normal url doesn't work.
    # Todo: check for this, somehow, and respond in a smart way.
    if(release=='3.35.0'):
        zipRemoteFileName = 'https://github.com/ANHIG/IMGTHLA/raw/3350/xml/hla.xml.zip'

    if(verbose):
        print('Local xml filename:' + str(xmlLocalFileName))
        print('Remote zip filename:' + str(zipRemoteFileName))

    requestData = get(zipRemoteFileName, allow_redirects=True)
    open(zipLocalFileName, 'wb').write(requestData.content)

    with ZipFile(zipLocalFileName, 'r') as zipObj:
        zipObj.extract('hla.xml', path=outputDirectory)

    if(isfile(xmlLocalFileName)):
        return xmlLocalFileName
    else:
        raise Exception('Problem when downloading and unzipping IMGT XML:' + str(zipRemoteFileName))

def cleanupSupplementalFiles(keepSuppFiles=False, supplementalFileDirectory=None):
    if(not keepSuppFiles):
        # delete the zip file and the xml file.
        zipFileName = join(supplementalFileDirectory, 'hla.xml.zip')
        xmlFileName = join(supplementalFileDirectory, 'hla.xml')
        remove(zipFileName)
        remove(xmlFileName)
        rmdir(supplementalFileDirectory)

def getLocusReferences():
    # A list of locus-level references. These hopefully do not change over IMGT HLA releases.
    # They also hopefully have longest UTR sequences, not in all cases.
    locusReferences = [
        'HLA-A*01:01:01:01'
        ,'HLA-B*07:02:01:01'
        ,'HLA-C*01:02:01:01'
        ,'HLA-DPA1*01:03:01:01'
        ,'HLA-DPB1*01:01:01:01'
        ,'HLA-DQA1*01:01:01:01'
        ,'HLA-DQB1*05:01:01:01'
        ,'HLA-DRB1*01:01:01:01'
        ,'HLA-DRB3*01:01:02:01'
        ,'HLA-DRB4*01:01:01:01'
        ,'HLA-DRB5*01:01:01:01'
        ,'HLA-E*01:01:01:01'
        ,'HLA-F*01:01:01:01'
        ,'HLA-G*01:01:01:01'
        ,'HLA-DMA*01:01:01:01'
        ,'HLA-DMB*01:01:01:01'
        ,'HLA-DOA*01:01:01'
        ,'HLA-DOB*01:01:01:01'
        ,'HLA-DRA*01:01:01:01'
        ,'MICA*001'
        ,'MICB*004:01:01'
    ]
    return locusReferences

def getIncludeSequenceList(imgtReleaseVersion=None):
    # These are hard-coded alleles that are included in the reference lists.
    # Try to only use full allele names here.
    includeSequences = []
    alleleDescriptionLookup = {}

    # DRB1*15:01:01:02 has longer UTRs. Better to keep for continuity.
    includeSequences.append('HLA-DRB1*15:01:01:02')
    includeSequences.append('HLA-B*56:01:01:02')

    # Kazu suggests to use  B*15:10:01:01 as a B71 reference, it differs greatly in exon 2 and this is useful.
    # SHould still include B*15:01 if possible.
    includeSequences.append('HLA-B*15:10:01:01')
    includeSequences.append('HLA-B*15:01:01:01')
    alleleDescriptionLookup['HLA-B*15:10:01:01'] = 'B71 Serotype Reference'

    # Kazu's suggestion, indeed this sequence has longest UTRs.
    includeSequences.append('HLA-B*27:05:02:01')

    # 08:01 and 08:03(from 17th IHIW) have shorter UTRs. Use 08:02.
    includeSequences.append('HLA-DRB1*08:02:01:01')

    # Kazu's suggestion, Historical Reference. Full-length appears in version 3.29.0. Otherwise use the 17th ref.
    if (imgtReleaseVersion in ['3.25.0', '3.26.0','3.27.0', '3.28.0']):
        includeSequences.append('HLA-DRB1*14:05:01')
    else:
        includeSequences.append('HLA-DRB1*14:54:01:01')

    # Kazu's suggestion, 03:03 and 03:04 have longer UTRs than 03:02(17th IHIW) or 03:01(not full length)
    includeSequences.append('HLA-C*03:03:01:01')

    # In 3.25 and 3.26 we didn't have a good full-length A*74:01 reference. Using this one for continuity and clarity.
    if(imgtReleaseVersion in ['3.25.0','3.26.0']):
        includeSequences.append('HLA-A*74:02:01:02')

    # # Starting in release 3.27.0 there is a full length DPB1*01 reference. Before that we use DPB1*02:01:02 for continuity
    if(imgtReleaseVersion in ['3.25.0','3.26.0']):
        includeSequences.append('HLA-DPB1*02:01:02')

    return sorted(list(set(includeSequences))), alleleDescriptionLookup

# TODO: Since we're no longer reading the 17th ws list, this isn't necessary. This is only alleles excluded from 17th list.
'''
def getExcludeSequenceList(imgtReleaseVersion=None):
    excludeSequenceList = []

    # These are hard-coded alleles that are excluded from the reference lists.

    # Starting in release 3.26.0 there is a full length HLA-B*41:01:01 reference,
    # So, B*40:305 is not needed.
    if(imgtReleaseVersion not in ['3.25.0']):
        excludeSequenceList.append('HLA-B*40:305')

    # Starting in release 3.27.0 there is a full length DPB1*01 reference
    # So, DPB1*02:01:01 are not needed.
    if(imgtReleaseVersion not in ['3.25.0','3.26.0']):
        excludeSequenceList.append('HLA-DPB1*02:01:02')

    # Starting in release 3.29.0, there is a full-length DPA1:01:03:01:01, used as standard locus reference.
    # Stop using DPA1:01:03:01:02
    if (imgtReleaseVersion not in ['3.25.0', '3.26.0', '3.27.0', '3.28.0']):
        excludeSequenceList.append('HLA-DPA1*01:03:01:02')


    return list(set(excludeSequenceList))
'''

def printSequenceDetails(alleleSequences=None, outputFilename=None, verbose=False, delimiter='\t', imgtReleaseVersion=None):
    if(verbose):
        print('Creating Sequence Details file:' + str(outputFilename))

    outputFile = open(outputFilename, 'w')

    outputFile.write(imgtReleaseVersion + ' Allele Name' + delimiter
        + imgtReleaseVersion + ' Sequence Length' + delimiter
        + imgtReleaseVersion + ' 5\'UTR Length' + delimiter
        + imgtReleaseVersion + ' 3\'UTR Length' + delimiter
        + imgtReleaseVersion + ' CWD Status' + '\n')

    # Loop sequences:
    alleleClusters = clusterSequences(alleleSequences=alleleSequences, verbose=verbose)

    for locus in sorted(alleleClusters.keys()):
        # print('Finding reference for locus ' + str(locus))
        for alleleGroup in sorted(alleleClusters[locus].keys()):
            for alleleSequence in alleleClusters[locus][alleleGroup]:
                # TODO: This will break if the feature is missing. Only running on full-len for now.
                outputFile.write(alleleSequence.alleleName + delimiter
                    + str(len(alleleSequence.getSequence())) + delimiter
                    + str(len(alleleSequence.featureSequences['5UTR'])) + delimiter
                    + str(len(alleleSequence.featureSequences['3UTR'])) + delimiter
                    + str(alleleSequence.cwdStatus) + '\n')

    outputFile.close()

def printSequenceCountsPerLocus(alleleSequences=None, outputFilename=None, verbose=None, imgtReleaseVersion=None, delimiter='\t'):
    if(verbose):
        print('Creating sequence counts per locus:' + str(outputFilename))

    outputFile = open(outputFilename, 'w')

    outputFile.write('Locus' + delimiter
        + imgtReleaseVersion + ' Reference Count' + '\n')

    # Loop sequences:
    alleleClusters = clusterSequences(alleleSequences=alleleSequences, verbose=verbose)

    for locus in sorted(alleleClusters.keys()):
        # print('Finding reference for locus ' + str(locus))
        seqCount = 0
        for alleleGroup in sorted(alleleClusters[locus].keys()):
            seqCount += len(alleleClusters[locus][alleleGroup])

        outputFile.write(locus + delimiter
            + str(seqCount) + '\n')

    # Write TOtal
    outputFile.write('Total' + delimiter
                     + str(len(alleleSequences)) + '\n')
    outputFile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-s", "--supplementary", help="keep supplementary files", action="store_true", default=False)
    parser.add_argument("-b", "--blast", help="keep blast output files (Very big)", action="store_true", default=False)
    parser.add_argument("-x", "--version", help="export allele list version number, default 1.0", default='1.0', type=str)
    parser.add_argument("-V", "--validate", help="validate full-length sequences by aligning against References", action="store_true")
    parser.add_argument("-r", "--release", required=True, help="IPD-IMGT/HLA release version", type=str)
    parser.add_argument("-o", "--output", required=True, help="Output Directory", type=str)
    parser.add_argument("-t", "--threads", required=False, help="Processor Threads", type=int, default=1)

    args = parser.parse_args()
    verbose = args.verbose

    outputDirectory = args.output
    supplementalFileDirectory = join(outputDirectory,'supplemental_files')

    if not isdir(outputDirectory):
        makedirs(outputDirectory)
    if not isdir(supplementalFileDirectory):
        makedirs(supplementalFileDirectory)

    print('Generating Reference Sequences')

    if verbose:
        print("Running in verbose mode.")

    xmlFileLocation = downloadImgtXml(outputDirectory=supplementalFileDirectory, release=args.release, verbose=verbose)
    alleleSequences, databaseVersion = parseXmlFile(xmlFile=xmlFileLocation,fullLengthOnly=True, verbose=verbose)
    # TODO: The database version from the XML file may be slightly different than the provided release number, due to minor versioning.
    #  I am naming files based on the "0" version but there might be 3.42.1 for example.
    #  Not completely accurate but its more consistent this way. May cause confusion.
    if(databaseVersion != args.release):
        print('Warning! the latest IMGT/HLA xml file shows a different (newer?) release date ('
              + str(databaseVersion) + ') than the provided release version (' + str(args.release) + ')')
    if(args.supplementary):
        printSequences(alleleSequences=alleleSequences, outputFilename=join(supplementalFileDirectory,str(args.release) + '_FullLengthSequences.fasta'), verbose=verbose)
    alleleSequenceClusters=clusterSequences(alleleSequences=alleleSequences,verbose=verbose)
    newReferenceSequences, missingSequences = createReferenceSequences(clusteredFullLenAlleleSequences=alleleSequenceClusters, verbose=verbose, imgtReleaseVersion=args.release)
    printSequences(alleleSequences=newReferenceSequences, outputFilename=join(outputDirectory, str(args.release) + '_ReferenceSequences.fasta'), verbose=verbose)
    printSequenceList(alleleSequences=newReferenceSequences, databaseVersion=args.release, outputDirectory=outputDirectory, verbose=verbose, fileVersion=args.version)
    if (args.supplementary):
        printSequenceDetails(alleleSequences=alleleSequences, outputFilename=join(supplementalFileDirectory, 'FullLengthSequenceDetails.csv'), verbose=verbose, imgtReleaseVersion=args.release)
        printSequenceDetails(alleleSequences=newReferenceSequences, outputFilename=join(supplementalFileDirectory, 'ReferenceSequenceDetails.csv'), verbose=verbose, imgtReleaseVersion=args.release)
        printSequenceCountsPerLocus(alleleSequences=alleleSequences, outputFilename=join(supplementalFileDirectory, 'FullLengthSequenceCountsPerLocus.csv'), verbose=verbose, imgtReleaseVersion=args.release)
        printSequenceCountsPerLocus(alleleSequences=newReferenceSequences, outputFilename=join(supplementalFileDirectory, 'ReferenceSequenceCountsPerLocus.csv'), verbose=verbose, imgtReleaseVersion=args.release)
        printMissingSequences(missingSequences=missingSequences, databaseVersion=args.release, outputDirectory=supplementalFileDirectory, verbose=verbose)
    if(args.validate):
        validationSet = alleleSequences
        validateFullLengthSequences(referenceSequences=newReferenceSequences, fullLengthSequences=validationSet
            , outputDirectory=outputDirectory, verbose=verbose, threadCount=args.threads, keepBlastFiles=args.blast)
    cleanupSupplementalFiles(keepSuppFiles=args.supplementary, supplementalFileDirectory=supplementalFileDirectory)
    print('Done. Reference Sequences were written to ' + str(outputDirectory))

