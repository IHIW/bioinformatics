#!/bin/bash
/home/oracle/scripts/mkSampleConf_HML_1.0.1_Stanford_v1.0.py -dir "/w_u/$1/upload/$2" -c "$1" -id "$3" -o "/w_u/$1/upload/$3.sample_metadata.conf" -uri "sftp://hidpl.stanford.edu/$1/upload/$2"
/home/oracle/scripts/myCreateBatch.pl /home/oracle/scripts/generic_metadata_original.conf $1 "$3" "$4" > "/w_u/$1/upload/$3.generic_metadata.conf"
