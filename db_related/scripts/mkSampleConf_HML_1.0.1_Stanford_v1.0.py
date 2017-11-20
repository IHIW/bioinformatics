#!/usr/bin/python

import os
import sys
import re
import argparse
import platform

def parse_args():
        parser = argparse.ArgumentParser(description='Creates sample_metadata.conf file based on hla typing result files.')
        parser.add_argument('-dir','--fastq_dir', help='Folder where the results are located.', required=False, action="store", dest="fdir", default=os.curdir)
        parser.add_argument('-o', '--output', help='Output file name (optional, default value is sample_metadata.conf).', required=False, action="store", dest="output_file", default="sample_metadata.conf")
        parser.add_argument('-c', '--centerCode', help='Center code',  required=False, action="store", dest="ccode", default="")
        parser.add_argument('-id', '--testid', help='Test id (optional).',  required=False, action="store", dest="test_id", default="")
        parser.add_argument('-hmlIdExt', '--hmlIdextension', help='HML id extension code (optional). Unique document identifier.',  required=False, action="store", dest="hmlid", default="")
        parser.add_argument('-uri', '--fileUrl', help='Url for R1 input file.',  required=False, action="store", dest="uri", default="")
        parser.add_argument('-s', '--source', help='Test id source (optional).',  required=False, action="store", dest="test_id_source", default="")
        parser.add_argument('-form', '--format', help='File format (default fastq).',  required=False, action="store", dest="format", default="FASTQ")
        parser.add_argument('-p', '--paired', help='Read file pairing (paired=true or paired=false), default value is true.',  required=False, action="store", dest="paired", default="true")
        parser.add_argument('-pool', '--pooled', help='Read file pooling (pooled=true or pooled=false), default value is true.',  required=False, action="store", dest="pooled", default="true")
        parser.add_argument('-avail', '--availabilty', help='Availabilty of read files (default value is private).',  required=False, action="store", dest="avail", default="private")
        parser.add_argument('-adTrim', '--adapterTrimmed', help='Sample is adapter trimmed (default value is true).',  required=False, action="store", dest="adapt", default="true")
        parser.add_argument('-qualTrim', '--qualityTrimmed', help='Sample is primer trimmed (default value is true).',  required=False, action="store", dest="qual", default="true")
        args = parser.parse_args()
        return args

def check_platform():
#Checks os to get an idea about path formats that will likely be used
    print 'Checking OS...'
    os=platform.system()
    print "\tDetected OS: " + os
    if os is 'Windows':
        delim= "\\\\"
    elif os is 'Linux':
        delim = '/'
    else:
        delim = "/"
    print "\tThe following delimiter will be used between the read file url and the filename: " + delim
    return delim


#Gets list of htr files using filenames in current or specified folder
def get_sample_list(fdir):
#Get list of htr files
    file_list = os.listdir(fdir)

    read_list = []
    for f in file_list:
        if re.search('.*htr', f):
            read_list.append(f)
    if not read_list:
        print "Please specify a folder containing result files using the \"-dir\" flag! If a folder was specified, check permissions!"
        print "Expected filename format: ID.S([0-9].*_L001_R1_001.*htr"
        print "Example filename mathching the expected format: BC-5_S32_L001_R1_001_2015-11-19_12-32-57.htr"
        sys.exit()
    else:     
        return read_list

def get_id(sample):
    print "HTR filename is: " + sample
    s = re.search('.S([0-9].*_L001_R1_001.*htr)', sample)
    print "Extracting sample ID..."
    if not s:
        print "Filename did not match expected pattern!"
        sys.exit()
    else:
        sample_id = sample[0:s.start()]
        print "The following ID was extracted: " + sample_id          
        return sample_id

def write_json(read_list,output_file,test_id,ccode,test_id_source,hmlid,uri,format,paired,pooled,avail,adapt,qual,delim):
#Get current time

# Creates JSON formatted file using the previously filtered lines.
        conf = open(output_file,'w')
        conf.write("hla {\n\texport {\n\t\thml {\n\t\t\tsamples = [")
        for sample in read_list:
            sample_id=get_id(sample)
            filename=sample
            conf.write("\n\t\t\t\t {")
#The .htr extension must be removed from the filenames, otherwise Twin won't fill out the id and centerCode fields.
            filename_reformatted = re.sub('.htr', '', filename)
            fastq_base = re.sub(r'_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', '', filename_reformatted)
            conf.write('\n\t\t\t\t\tfilename = "%s"'% (filename_reformatted))            
            conf.write('\n\t\t\t\t\tid = "%s"' % (sample_id))
            conf.write('\n\t\t\t\t\tcenterCode = "%s"' % (ccode))
            if test_id != "":
                conf.write('\n\t\t\t\t\ttestId = "%s"' % (test_id))
            else:
                conf.write('\n\t\t\t\t\ttestId= ""')
            if hmlid != "":
                conf.write('\n\t\t\t\t\thmlIdExtension = "%s-%s"' % ("hml-1.0.1",hmlid))
            else:
                conf.write('\n\t\t\t\t\thmlIdExtension = "hml-1.0.1"')
            if test_id_source != "":
                conf.write('\n\t\t\t\t\ttestIdSource = "%s"' % (test_id_source))
            else:
                conf.write('\n\t\t\t\t\ttestIdSource = ""')
            if uri != "":
                conf.write('\n\t\t\t\t\turi = "%s%s%s.fastq.gz"' % (uri,delim,fastq_base))
            else:
                conf.write('\n\t\t\t\t\turi = ""')
            if format != "FASTQ":
                conf.write('\n\t\t\t\t\tformat = "%s"' % (format))
            else:
                conf.write('\n\t\t\t\t\tformat = "FASTQ"')
            if paired != "true":
                conf.write('\n\t\t\t\t\tpaired = "%s"' % (paired))
            else:
                conf.write('\n\t\t\t\t\tpaired = "true"')
            if paired != "true":
                conf.write('\n\t\t\t\t\tpooled = "%s"' % (pooled))
            else:
                conf.write('\n\t\t\t\t\tpooled = "true"')
            if avail != "private":
                conf.write('\n\t\t\t\t\tavailability = "%s"' % (avail))
            else:
                conf.write('\n\t\t\t\t\tavailability = "private"')
            if adapt != "true":
                conf.write('\n\t\t\t\t\tadapterTrimmed = "%s"' % (adapt))
            else:
                conf.write('\n\t\t\t\t\tadapterTrimmed = "true"')
            if qual != "true":
                conf.write('\n\t\t\t\t\tqualityTrimmed = "%s"' % (qual))
            else:
                conf.write('\n\t\t\t\t\tqualityTrimmed = "true"')
            conf.write("\n\t\t\t\t },")
        conf.write("\n\t\t\t\t]\n\t\t\t}\n\t\t}\n}")
        return conf


def main():
        args = parse_args()
        delim = check_platform()
        filelist=get_sample_list(args.fdir)
        files = write_json(filelist,args.output_file,args.test_id,args.ccode,args.test_id_source,args.hmlid,args.uri,args.format,args.paired,args.pooled,args.avail,args.adapt,args.qual,delim)
        if os.path.isfile(args.output_file):
            print "Config file successfully created."

if __name__ == "__main__":
    sys.exit(main())
