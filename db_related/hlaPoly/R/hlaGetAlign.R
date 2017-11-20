hlaGetAlign <-
function(hlaData,genAllele,nucAllele,conseq)#,sposition=0,substitutionMatrix="Gonnet")
{
	typeIsGen=T #indicate there is genomic sequence in the database
	genIndex=which(hlaData$genfa$name==genAllele)
	nucIndex=which(hlaData$genfa$name==nucAllele)
	if(length(nucIndex)==0){
		nucIndex=which(hlaData$nucfa$name==nucAllele)
		typeIsGen=F
	}
	
	
	genFeature=hlaData$features[hlaData$features$name==genAllele,]
	genFeature<-genFeature[order(genFeature$start),]
	genFeature$ind=as.integer(sub(".*\\.","",genFeature$fid))
	nucFeature=hlaData$features[hlaData$features$name==nucAllele,]
	nucFeature<-nucFeature[order(nucFeature$start),]
	nucFeature$ind=as.integer(sub(".*\\.","",nucFeature$fid))
	
	#add pipe | to pair for the msq
	 genseq <- readDNAStringSet(hlaData$genfa[genIndex, ])
	 if(typeIsGen){
	 	nucseq<-readDNAStringSet(hlaData$genfa[nucIndex, ])
	 }
	 else{
	 	nucseq<-readDNAStringSet(hlaData$nucfa[nucIndex, ])
	 }
	conseq=gsub(" ","",gsub("\n","",gsub("\r","",conseq)))
#	if(typeIsGen){
#		nuclen=width(nucseq)
#		if(sposition+nchar(conseq)>nuclen && sposition<nuclen){
#			substr(conseq,nuclen-sposition+1,nchar(conseq))=paste(rep("+",nchar(conseq)-nuclen+sposition ),collapse="")
#		}


#	}
	conseq=gsub("N+$","",gsub("^N+","",conseq))
#	 pipeGenseq = DNAStringSet(".")
#	 for(i in 1:nrow(genFeature)){
#	 	pipeGenseq=xscat(pipeGenseq,subseq(genseq,start=genFeature$start[i],end=genFeature$end[i]),".")
#	 }
#	 pipeNucseq = DNAStringSet(".")
#	 for(i in 1:nrow(nucFeature)){
#	 	pipeNucseq=xscat(pipeNucseq,subseq(nucseq,start=nucFeature$start[i],end=nucFeature$end[i]),".")
#	 }
	
	 
	if(!typeIsGen && nrow(nucFeature)>1){
		nucNseq=subseq(nucseq,start=nucFeature$start[1],end=nucFeature$end[1])
	 for(i in 2:nrow(nucFeature)){
			 if((nucFeature$ind[i]-nucFeature$ind[i-1])==2 && length(which(genFeature$ind==(nucFeature$ind[i]-1)))==1){
				indinGen=which(genFeature$ind==(nucFeature$ind[i]-1))
				introLen=genFeature$end[indinGen]-genFeature$start[indinGen]+1
				nucNseq=xscat(nucNseq,paste(rep("N",introLen),collapse = ""))

			  }
				nucNseq=xscat(nucNseq,subseq(nucseq,start=nucFeature$start[i],end=nucFeature$end[i]))
	 }
		nucseq=nucNseq
	} 
	 conseqDNA<-DNAStringSet(conseq)
	 trio<-append(genseq,c(nucseq,conseqDNA))
#	 trio<-append(pipeGenseq,c(pipeNucseq,conseqDNA))
	# names(trio)=c("gen","nuc","input")
	 names(trio)=c(genAllele,nucAllele,"input")
	 
	 #
	 triomsa<-AlignSeqs(trio,verbose=F,processors = NULL)
#	 triomsa<-AlignSeqs(trio,verbose=F)
#	 triomsa<-msa(trio,method="ClustalOmega",order ="input")

#	 triomsa<-msaMuscle(trio,order ="input",cluster=substitutionMatrix,gapExtension=1)
	result = hlaData
	result$triomsa=triomsa
	result$genFeature = genFeature
	result$nucFeature = nucFeature
	result$conseq = conseq
	class(result) = "hlaData.msa"
	return(result)
}
