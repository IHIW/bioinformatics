hlaDataInit <-
function(hlaVersion="3.25.0"){
	genfa<-fasta.index(system.file("hlalib",hlaVersion,"hla_gen.fasta",package="hlaPoly"),seqtype="DNA")
	genfa$name<-paste("HLA-",sapply(strsplit(genfa$desc," ") ,"[[",2),sep="")

	nucfa<-fasta.index(system.file("hlalib",hlaVersion,"hla_nuc.fasta",package="hlaPoly"),seqtype="DNA")
	nucfa$name<-paste("HLA-",sapply(strsplit(nucfa$desc," ") ,"[[",2),sep="")

	features<-read.table(system.file("hlalib",hlaVersion,"hlaxml.tsv",package="hlaPoly"),sep="\t",quote="\"",header=F,colClasses = c("factor", "factor","factor","factor","character","integer","integer"))
	colnames(features)<-c("name","accnum","version","feature","fid","start","end")
	

	reftable<-read.table(system.file("hlalib",hlaVersion,"hlaRefAllele.tsv",package="hlaPoly"),sep="\t",header=F)
	refAllele<-reftable$V2
	names(refAllele) = reftable$V1
	refAllele<-refAllele
	result<-list(genfa=genfa, nucfa=nucfa, features=features, refAllele=refAllele)
	class(result) = "hlaData.init"
	return(result)

}
