import argparse

delimiter = '\t'
newline = '\n'

def readAlleleNames(inputFileName='CWDInputFiles/AlleleList.csv'):
    print('Reading Allele Names:' + str(inputFileName))
    alleleNames = []
    cwd2Lookup = {}
    with open(inputFileName, 'r') as inputFile:
        for row in inputFile:
            #print('Allele Names input row:' + str(row))
            row = row.replace('\n', '')
            # Columns are [AlleleName, CWD2.0 status]
            cols = row.split(delimiter)
            alleleNames.append(cols[0])
            cwd2Lookup[cols[0]]=cols[1]

    return alleleNames,cwd2Lookup

def readCwd1InputFile(inputFileName='CWDInputFiles/CWD1.0Lookup.csv'):
    print('Reading CWD1.0 input:' + str(inputFileName))
    # A list of alleles. One column
    cwd1Lookup = []
    with open(inputFileName, 'r') as inputFile:
        for row in inputFile:
            #print('CWD1.0 input:' + str(row))
            row = row.replace('\n', '')
            cols = row.split(delimiter)
            cwd1Lookup.append(cols[0])
    return cwd1Lookup

def readCwd3InputFile(inputFileName='CWDInputFiles/CWD3.0Lookup.csv'):
    print('Reading CWD3.0 input:' + str(inputFileName))
    cwd3Lookup = {}
    with open(inputFileName, 'r') as inputFile:
        for row in inputFile:
            #print('CWD3.0 input:' + str(row))
            row = row.replace('\n', '')
            # Columns are [AlleleName, CWD2.0 status]
            cols = row.split(delimiter)
            cwd3Lookup[cols[0]]=cols[1]

    return cwd3Lookup

def readCwdEfiInputFile(inputFileName='CWDInputFiles/CWD.EFI.csv'):
    print('Reading CWDEfiLookup input:' + str(inputFileName))
    cwdEfiLookup = {}
    with open(inputFileName, 'r') as inputFile:
        for row in inputFile:
            #print('cwdEfiLookup input:' + str(row))
            row = row.replace('\n', '')
            # Columns are [AlleleName, CWDEfi status]
            cols = row.split(delimiter)
            cwdEfiLookup[cols[0]]=cols[1]

    return cwdEfiLookup

def readCwdDeInputFile(inputFileName='CWDInputFiles/CWD.DE.csv'):
    print('Reading CWDDeLookup input:' + str(inputFileName))
    cwdDeLookup = {}
    with open(inputFileName, 'r') as inputFile:
        for row in inputFile:
            #print('cwdDeLookup input:' + str(row))
            row = row.replace('\n', '')
            # Columns are [AlleleName, CWDEfi status]
            cols = row.split(delimiter)
            cwdDeLookup[cols[0]]=cols[1]

    return cwdDeLookup

def readCwdChinaInputFile(inputFileName='CWDInputFiles/CWD.China.csv'):
    print('Reading CWDChinaLookup input:' + str(inputFileName))
    cwdChinaLookup = {}
    with open(inputFileName, 'r') as inputFile:
        for row in inputFile:
            #print('CWDChinaLookup input:' + str(row))
            row = row.replace('\n', '')
            # Columns are [locus,AlleleName, count, CWDChina status]
            cols = row.split(delimiter)
            cwdChinaLookup[cols[1]]=cols[3]

    return cwdChinaLookup



def getCwd1Status(alleleName=None, cwd1Lookup=None):
    #print('cwd1Lookup:' + str(cwd1Lookup))
    # CWD1.0 is a list of alleles. There are no colons for proper nomenclature.
    # Easy case, check full allele name
    if(alleleName.replace(':','') in cwd1Lookup):
        return 'CWD'
    # Check 3 field and 2 field
    while(alleleName.count(':')>1):
        alleleName = alleleName[0:alleleName.rindex(':')]
        #print('Checking this shortened allele:' + str(alleleName))
        if (alleleName.replace(':', '') in cwd1Lookup):
            return 'CWD'
    return 'Not CWD defined'


def getCwd3Status(alleleName=None, cwd3Lookup=None):
    #print('cwd3Lookup:' + str(cwd1Lookup))
    # CWD3.0 is a dictionary lookup. They're just missing 'HLA-' at the beginning
    alleleName=alleleName.replace('HLA-','')

    cwdStatus=''
    # Easy case, check full allele name
    #print('Checking allele name ' + str(alleleName) + ' against keys ' + str(cwd3Lookup.keys()))
    if(alleleName in cwd3Lookup.keys()):
        if(cwd3Lookup[alleleName].strip()==''):
            # This one is blank, no status found "yet"
            pass
        else:
            cwdStatus = cwd3Lookup[alleleName]
    # Check 3 field and 2 field
    while(alleleName.count(':')>1):
        alleleName = alleleName[0:alleleName.rindex(':')]
        #print('Checking this shortened allele:' + str(alleleName))
        if (alleleName in cwd3Lookup.keys()):
            if (cwd3Lookup[alleleName].strip() == ''):
                # This one is blank, no status found "yet"
                pass
            else:
                cwdStatus = cwd3Lookup[alleleName]
        # Check 3 field and 2 field

    if (cwdStatus.strip() == ''):
        return 'Not CWD defined'
    else:
        return cwdStatus


def getCwdEfiStatus(alleleName=None, cwdEfiLookup=None):
    # This lookup only uses 2 field allele names, fairly straight forward
    # remove "('"HLA-"
    # Trim down to two fields:
    #print('Checking EFI allele:' + str(alleleName))
    alleleName = alleleName.replace('HLA-','')
    while (alleleName.count(':') > 1):
        alleleName = alleleName[0:alleleName.rindex(':')]
    #print('Checking EFI short name:' + str(alleleName))

    if (alleleName in cwdEfiLookup.keys()):
        return cwdEfiLookup[alleleName]
    else:
        return 'Not CWD defined'


def getCwdDEStatus(alleleName=None, cwdDeLookup=None):
    # This lookup only uses 2 field allele names, fairly straight forward
    # remove "('"HLA-"
    # Trim down to two fields:
    #print('Checking DE allele:' + str(alleleName))
    alleleName = alleleName.replace('HLA-', '')
    while (alleleName.count(':') > 1):
        alleleName = alleleName[0:alleleName.rindex(':')]
    #print('Checking DE short name:' + str(alleleName))

    if (alleleName in cwdDeLookup.keys()):
        return cwdDeLookup[alleleName]
    else:
        return 'Not CWD defined'

def getCwdChinaStatus(alleleName=None, cwdChinaLookup=None):
    # This lookup only uses 2 field allele names, fairly straight forward
    # remove "('"HLA-"
    # Trim down to two fields:
    # print('Checking China allele:' + str(alleleName))
    alleleName = alleleName.replace('HLA-', '')
    while (alleleName.count(':') > 1):
        alleleName = alleleName[0:alleleName.rindex(':')]
    # print('Checking China short name:' + str(alleleName))

    if (alleleName in cwdChinaLookup.keys()):
        return cwdChinaLookup[alleleName]
    else:
        return 'Not CWD defined'


def writeOutputFile(outputFileName=None, alleleNames=None, cwd1Lookup=None
                    , cwd2Lookup=None, cwd3Lookup=None, cwdEfiLookup=None, cwdDeLookup=None
                    , cwdChinaLookup=None):
    with open(outputFileName,'w') as outputFile:
        # Header
        # TODO: Add more info in the header here. DOI #, Loci represented. Description of what the columns mean.
        outputFile.write('Allele Name'
            + delimiter + '"CWD 1.0\nhttps://doi.org/10.1016/j.humimm.2007.01.014\nLoci:HLA-A,-B,-C,-DRB1,-DRB3,-DRB4,-DRB5,-DQA1,-DQB1,-DPB1"'
            + delimiter + '"CWD 2.0\nhttps://doi.org/10.1111/tan.12093\nLoci:HLA-A,-B,-C,-DRB1,-DRB3,-DRB4,-DRB5,-DQA1,-DQB1,-DPA1,-DPB1"'
            + delimiter + '"CWD 3.0\nhttps://doi.org/10.1111/tan.13811\nLoci:HLA-A,-B,-C,-DRB1,-DRB3,-DRB4,-DRB5,-DQB1,-DPB1"'
            + delimiter + '"CWD European (EFI)\nhttps://doi.org/10.1111/tan.12956\nLoci:HLA-A,-B,-C,-DRB1,-DQA1,-DQB1,-DPB1"'
            + delimiter + '"CWD German Stem Cell Donors\nhttps://doi.org/10.1111/tan.13378\nLoci:HLA-A,-B,-C,-DRB1,-DQB1,-DPB1"'
            + delimiter + '"CWD China\nhttps://doi.org/10.1111/tan.13358\nLoci:HLA-A,-B,-C,-DRB1,-DQB1"'
            + newline)

        for alleleName in alleleNames:
            outputFile.write(alleleName
                + delimiter + getCwd1Status(alleleName=alleleName, cwd1Lookup=cwd1Lookup)
                + delimiter + cwd2Lookup[alleleName]
                + delimiter + getCwd3Status(alleleName=alleleName, cwd3Lookup=cwd3Lookup)
                + delimiter + getCwdEfiStatus(alleleName=alleleName, cwdEfiLookup=cwdEfiLookup)
                + delimiter + getCwdDEStatus(alleleName=alleleName, cwdDeLookup=cwdDeLookup)
                + delimiter + getCwdChinaStatus(alleleName=alleleName, cwdChinaLookup=cwdChinaLookup)
                + newline)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    # Read Allele List File (Columns Allele Name, CWD2.0 status)
    alleleNames, cwd2Lookup = readAlleleNames()
    cwd1Lookup = readCwd1InputFile()
    cwd3Lookup = readCwd3InputFile()
    cwdEfiLookup = readCwdEfiInputFile()
    cwdDeLookup = readCwdDeInputFile()
    cwdChinaLookup = readCwdChinaInputFile()

    # Write everything.
    outputFileName='CwdStatus.csv'
    writeOutputFile(outputFileName=outputFileName, alleleNames=alleleNames
        , cwd1Lookup=cwd1Lookup, cwd2Lookup=cwd2Lookup, cwd3Lookup=cwd3Lookup, cwdEfiLookup=cwdEfiLookup
        , cwdDeLookup=cwdDeLookup, cwdChinaLookup=cwdChinaLookup)



    print('Done. Reference Sequences were written to ' + str(outputFileName))

