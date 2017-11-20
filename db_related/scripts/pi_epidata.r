#!/usr/bin/env Rscript
source("/home/oracle/scripts/parse_epitope2.r")
source("/home/oracle/scripts/insert_epidata2.r")
args <- commandArgs(TRUE)
epidata=read.fusion_xls(args[2])
insert_epidata(args[1],epidata)
