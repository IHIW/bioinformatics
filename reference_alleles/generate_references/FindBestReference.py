from os.path import join, isdir
from os import makedirs, remove, rmdir
from Bio import pairwise2
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio.Blast.NCBIXML import parse as parseNCBIXML
from Bio.pairwise2 import format_alignment
import multiprocessing
from subprocess import call
import time

def currentMillis():
    return int(round(time.time()))

def printSequences(alleleSequences=None, outputFilename=None, verbose=False):
    if(verbose):
        print('Writing ' + str(len(alleleSequences)) + ' allele sequences to:' + str(outputFilename))

    outputFile = open(outputFilename,'w')
    for alleleSequence in alleleSequences:
        outputFile.write('>' + str(alleleSequence.alleleName) + '\n')
        outputFile.write(str(alleleSequence.getSequence()) + '\n')
    outputFile.close()

def findBestReferenceSequence(referenceFileName=None, batchFileName=None, verbose=False):
    if(verbose):
        print('Finding best reference for query sequence (Blast):' + str(batchFileName))
    alignmentScores = {}

    cline = NcbiblastnCommandline(query=str(batchFileName)
        , db=str(referenceFileName)
        , out=str(batchFileName) + '.xml'
        , outfmt=5
        , num_threads=1)

    stdout, stderr = cline()
    if(verbose):
        print ('Blast Commandline:\n' + str(cline))
        print ('Output:' + str(stdout))
        print ('Errors?:' + str(stderr))

    blastRecords = parseNCBIXML(open(str(batchFileName) + '.xml'))
    #print('blastRecords:\n' + str(blastRecords))
    for blastRecordIndex, blastRecord in enumerate(blastRecords):
        queryName = blastRecord.query
        if(len(blastRecord.alignments)<1):
            # No alignments for this query!
            hitName = 'None'
            hitScore = 0
        else:
            # take the first alignment. It has the highest score.
            hitName = blastRecord.alignments[0].hit_def
            hitScore = blastRecord.descriptions[0].score
            if(verbose):
                print('Query ' + queryName + ' has hit ' + hitName + ' with score ' + str(hitScore))

        alignmentScores[queryName] = (hitName, hitScore)

    return alignmentScores


def validateFullLengthSequences(referenceSequences=None, fullLengthSequences=None, outputDirectory=None, threadCount=1
    , verbose=False, delimiter='\t', keepSuppFiles=False):
    validateFullLengthSequencesUsingBlast(referenceSequences=referenceSequences
        , fullLengthSequences=fullLengthSequences
        , outputDirectory=outputDirectory
        , threadCount=threadCount
        , verbose=verbose
        , delimiter=delimiter
        ,keepSuppFiles=keepSuppFiles)


def cleanupBlastOutputFiles(blastDirectory=None, referenceFileName=None, queryBatches=None):
    # Obviously, Be careful with deleting files. Only files that I made, don't search the directory or anything like that.
    print('Cleaning up BLAST output files.')

    # Delete Blast Reference files
    remove(referenceFileName)
    remove(referenceFileName + '.nhr')
    remove(referenceFileName + '.nin')
    remove(referenceFileName + '.nsq')

    for batchIndex, batch in enumerate(queryBatches):
        batchFastaFileName = join(blastDirectory, 'Batch' + str(batchIndex) + 'Sequences.fasta')
        batchXmlFileName = join(blastDirectory, 'Batch' + str(batchIndex) + 'Sequences.fasta.xml')
        remove(batchFastaFileName)
        remove(batchXmlFileName)

    # Delete Directory
    rmdir(blastDirectory)

    pass


def validateFullLengthSequencesUsingBlast(referenceSequences=None, fullLengthSequences=None, outputDirectory=None
        , threadCount=1, batchSize=50, verbose=False, delimiter='\t', keepSuppFiles=False):
    # TODO: Blast is hopefully faster than pairwise alignments.
    #  But it's only doing local alignments.
    print('Validating ' + str(len(fullLengthSequences)) + ' sequences against ' + str(len(referenceSequences))
          + ' Reference Sequences using Blast Alignments (threads=' + str(threadCount) + ')' )
    # Start a thread pool
    queryBatches=[]
    batchResults=[]

    pool = multiprocessing.Pool(threadCount)
    before = currentMillis()

    # Create Blast Reference
    blastDirectory = join(outputDirectory,'blast_results')
    if(not isdir(blastDirectory)):
        makedirs(blastDirectory)
    referenceFileName = join(blastDirectory,'BlastReference.fasta')
    printSequences(alleleSequences=referenceSequences, outputFilename=referenceFileName, verbose=verbose)
    cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=referenceFileName)
    stdout, stderr = cline()
    if(verbose):
        print ('MakeDB Commandline:\n' + str(cline))
        print ('Output:' + str(stdout))
        print ('Errors?:' + str(stderr))

    # Split Query Sequences into Batches
    if(verbose):
        print ('Splitting ' + str(len(fullLengthSequences)) + ' sequences into batches of size ' + str(batchSize))
    newBatch=[]
    for sequenceIndex, sequence in enumerate(fullLengthSequences):
        newBatch.append(sequence)
        if(len(newBatch) >= batchSize or sequenceIndex==len(fullLengthSequences)-1):
            # Done with this batch. Write it to file
            batchFileName = join(blastDirectory, 'Batch' + str(len(queryBatches)) + 'Sequences.fasta')
            printSequences(alleleSequences=newBatch, outputFilename=batchFileName, verbose=verbose)

            queryBatches.append(newBatch)
            newBatch=[]
    if(verbose):
        print ('Found ' + str(len(queryBatches)) + ' batches of size <= ' + str(batchSize))

    # For each Batch
    for batchIndex, batch in enumerate(queryBatches):
        batchFileName = join(blastDirectory, 'Batch' + str(batchIndex) + 'Sequences.fasta')
        # Start thread to run blast against references
        if (threadCount > 1):
            batchResults.append(pool.starmap_async(findBestReferenceSequence, [[referenceFileName, batchFileName, verbose]]))
        else:
            batchResults.append(findBestReferenceSequence(referenceFileName=referenceFileName, batchFileName=batchFileName, verbose=verbose))

    pool.close()
    pool.join()

    # Delete blast output files
    if(not keepSuppFiles):
        cleanupBlastOutputFiles(blastDirectory=blastDirectory, referenceFileName=referenceFileName, queryBatches=queryBatches)

    # Create output file
    sequenceValidationResultsFile = open(join(outputDirectory,'ReferenceFinderValidationResults.csv'), 'w')
    sequenceValidationResultsFile.write('Query_Name' + delimiter + 'Best_Reference' + delimiter + 'Alignment_Score\n')

    # Each batch result should be a dictionary. Take those results and write them
    for batchResult in batchResults:
        if (threadCount > 1):
            # If it's multi threaded we need to "get" the value
            currentBatchResults = batchResult.get()[0]
        else:
            currentBatchResults = batchResult

        for queryAlleleName in currentBatchResults.keys():
            bestReferenceName, alignmentScore = currentBatchResults[queryAlleleName]
            sequenceValidationResultsFile.write(str(queryAlleleName) + delimiter + str(bestReferenceName) + delimiter + str(alignmentScore) + '\n')

    sequenceValidationResultsFile.close()

    if(verbose):
        after = currentMillis()
        print('Finding References ' + str(len(fullLengthSequences)) + ' sequences took ' + str((after-before)) + ' seconds.')

# TODO: Pairwise might work the best, it finds optimum alignments.
#  But it does not work well for longer sequences with lots of gaps,
#  Aligning DPB1 alleles overflows my memory, might work better on a machine with lots of memory.
#  For now this code is deprecated, because BLAST works better.
'''
def validateFullLengthSequencesUsingPairwise(referenceSequences=None, fullLengthSequences=None, outputDirectory=None, threadCount=1, verbose=False, delimiter='\t'):

    print('Validating ' + str(len(fullLengthSequences)) + ' sequences against ' + str(len(referenceSequences)) + ' Reference Sequences using Pairwise Alignments')
    # Start a thread pool
    resultSets = {}
    pool = multiprocessing.Pool(threadCount)
    before = currentMillis()

    for allele in fullLengthSequences:
        locus=allele.getLocus()
        if (threadCount > 1):
            resultSets[allele.alleleName]=pool.starmap_async(findBestReferenceSequencePairwise, [[allele, referenceSequences, locus, verbose]])
        else:
            resultSets[allele.alleleName]=findBestReferenceSequencePairwise(querySequence=allele, referenceSequences=referenceSequences, locusFilter=locus, verbose=verbose)

    pool.close()
    pool.join()

    # Create output file
    sequenceValidationResultsFile = open(join(outputDirectory,'ReferenceFinderResultsPairwise.csv'), 'w')
    sequenceValidationResultsFile.write('Query_Name' + delimiter + 'Best_Reference' + delimiter + 'Alignment_Score\n')

    for alleleName in resultSets.keys():
        if (threadCount > 1):
            # If it's multi threaded we need to "get" the value
            bestReferenceName, alignmentScore = resultSets[alleleName].get()[0]
        else:
            bestReferenceName, alignmentScore = resultSets[alleleName]

        # write Query Name, Best Matching Reference, Score
        sequenceValidationResultsFile.write(str(alleleName) + delimiter + str(bestReferenceName) + delimiter + str(alignmentScore) + '\n')


    sequenceValidationResultsFile.close()

    if(verbose):
        after = currentMillis()
        print('Finding references ' + str(len(fullLengthSequences)) + ' sequences took ' + str((after-before)) + ' seconds.')

def findBestReferenceSequencePairwise(querySequence=None, referenceSequences=None, locusFilter=None, verbose=False):
    if(verbose):
        print('Finding best reference for query sequence:' + str(querySequence.alleleName))
        print('Locus Filter:' + str(locusFilter))

    alignmentScores = {}

    # for each reference sequence
    for referenceSequence in referenceSequences:
        if (locusFilter is None or locusFilter == referenceSequence.getLocus()):
            # alignmentScore
            #print('Aligning against ' + str(referenceSequence.alleleName))
            # Match=2, Mismatch=-1, GapOpen=-0.5, GapExtend=-0.1
            alignments = pairwise2.align.globalms(querySequence.getSequence(), referenceSequence.getSequence(),  2, -1, -0.5, -0.1)
            alignmentScore = int(alignments[0][2])
            #print('Found a score of ' + str(alignmentScore))
            alignmentScores[referenceSequence.alleleName]=alignmentScore

    bestReferenceName = max(alignmentScores, key=alignmentScores.get)
    if(verbose):
        print('Best reference for query sequence:' + str(querySequence.alleleName) + ' is sequence ' + bestReferenceName + ' with a score of ' + str(alignmentScores[bestReferenceName]))

    return (bestReferenceName, alignmentScores[bestReferenceName])
'''



# TODO: Make a main method.
#  Pass Fasta sequences (Or other?) for query and reference.
#  Output a csv file of results.
