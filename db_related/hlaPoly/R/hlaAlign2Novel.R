hlaAlign2Novel <-
function(hlaData){
	 msamatrix<-t(Biostrings::as.matrix(hlaData$triomsa))
	 colnames(msamatrix)=c("gen","nuc","input")
	 trioIndex<-c(0,0,0)
	 names(trioIndex)=colnames(msamatrix)
	 
	mismatches <- data.frame(reftype=factor(levels=c("gen","nuc","begin","end","ext")),
	                 refseq=character(),
	                 conseq=character(),
	                 start=integer(),
	                 end=integer(),
	                 alnstart=integer(),
	                 alnend=integer(),
					stringsAsFactors=FALSE ) 
	 
	 for(i in 1:nrow(msamatrix)){
		if(msamatrix[i,2]=="N") msamatrix[i,2]="-"
	 	for(j in 1:ncol(msamatrix)){
	 		if(msamatrix[i,j]!='-'){
	 			trioIndex[j]=trioIndex[j]+1
	 		}
	 	}
	 	if(trioIndex["input"]>=1 && trioIndex["input"]<=nchar(hlaData$conseq)){ #only the region in input is considered
	 		if(msamatrix[i,"input"]!=msamatrix[i,"nuc"]){  #mismatch with nuc
	 			if(msamatrix[i,"nuc"]=="-"){ #insertion ?
	 				if(trioIndex["nuc"]>0 && !(trioIndex["nuc"] %in% hlaData$nucFeature$end)){ # real insertion 
	  				 mismatches[nrow(mismatches)+1,]<-data.frame(reftype="nuc",refseq="-",conseq=msamatrix[i,"input"],start=trioIndex["nuc"],end=trioIndex["nuc"],alnstart=i,alnend=i,stringsAsFactors=FALSE)
	 				}
	 				else{	#ouside or intron of nuc
	 					if(msamatrix[i,"input"]!=msamatrix[i,"gen"]){
							if(msamatrix[i,"gen"]=="-"){  #insertion or outside of ref gene
								if(trioIndex["gen"]>0 && !(trioIndex["gen"] %in% hlaData$genFeature$end)){ # with gene
	 								mismatches[nrow(mismatches)+1,]<-data.frame(reftype="gen",refseq=msamatrix[i,"gen"],conseq=msamatrix[i,"input"],start=trioIndex["gen"],end=trioIndex["gen"],alnstart=i,alnend=i,stringsAsFactors=FALSE)
								}
								else{
									if(trioIndex["gen"]==0){ #begin
	 								mismatches[nrow(mismatches)+1,]<-data.frame(reftype="begin",refseq="-",conseq=msamatrix[i,"input"],start=-1,end=-1,alnstart=i,alnend=i,stringsAsFactors=FALSE)
									}
									else{ #end
	 								mismatches[nrow(mismatches)+1,]<-data.frame(reftype="end",refseq="-",conseq=msamatrix[i,"input"],start=-1,end=-1,alnstart=i,alnend=i,stringsAsFactors=FALSE)
									}

								}
	 						}
							else{	#mismatch or deletion
	 							mismatches[nrow(mismatches)+1,]<-data.frame(reftype="gen",refseq=msamatrix[i,"gen"],conseq=msamatrix[i,"input"],start=trioIndex["gen"],end=trioIndex["gen"],alnstart=i,alnend=i,stringsAsFactors=FALSE)


							}
	 					}
						else{ # for ext
	 								mismatches[nrow(mismatches)+1,]<-data.frame(reftype="gen",refseq="",conseq="",start=trioIndex["gen"],end=trioIndex["gen"],alnstart=i,alnend=i,stringsAsFactors=FALSE)
						}
	 				}
	 			}
	 			else{
					#mismatch or deletion
	 				mismatches[nrow(mismatches)+1,]<-data.frame(reftype="nuc",refseq=msamatrix[i,"nuc"],conseq=msamatrix[i,"input"],start=trioIndex["nuc"],end=trioIndex["nuc"],alnstart=i,alnend=i,stringsAsFactors=FALSE)
	 				
	 			}
	 			
	 		}
			else if(trioIndex["nuc"]<=0 ||  (trioIndex["nuc"]>=max(hlaData$nucFeature$end))){
					if(msamatrix[i,"nuc"]=="-"){ # for ext
	 							mismatches[nrow(mismatches)+1,]<-data.frame(reftype="gen",refseq=msamatrix[i,"gen"],conseq=msamatrix[i,"input"],start=trioIndex["gen"],end=trioIndex["gen"],alnstart=i,alnend=i,stringsAsFactors=FALSE) 
					}


			}
	 		if(trioIndex["input"]==nchar(hlaData$conseq)){
	 			break
	 		}
	 	}
	 }
	# after this step, the locaitons of mismatches are OK
	hlaData$nucFeature$reftype="nuc"
	hlaData$genFeature$reftype="gen"
	refFeature=rbind(hlaData$nucFeature,hlaData$genFeature)
	refFeature$feature=as.character(refFeature$feature)
	if(nrow(mismatches)==0){
		return(NA)
	}
	mismatches$fid=""
	mismatches$feature_start=0
	mismatches$feature=""
	# get features
	attach(refFeature)
	for(i in 1:nrow(mismatches)){
		feature_pos<-refFeature[reftype==mismatches$reftype[i] & mismatches$start[i]>=start &  mismatches$end[i] <=end,c("fid","start","feature")]
		if(nrow(feature_pos)>0){
			mismatches[i,c("fid","feature_start","feature")] = feature_pos 
		}
	}
	detach(refFeature)
	
	reformMismatches <- data.frame(reftype=factor(levels=c("gen","nuc","begin","end","ext")),
	                 refseq=character(),
	                 conseq=character(),
	                 start=integer(),
	                 end=integer(),
	                 alnstart=integer(),
	                 alnend=integer(),
	                 fid=character(),
	                 feature_start = integer(),
					feature=character(),
	                 stringsAsFactors=FALSE )
	 

	if(mismatches$refseq[1]=="-"){
            mismatches$start[1] = mismatches$start[1] +1
	}
	 reformMismatches[nrow(reformMismatches)+1,]<- mismatches[1,]

	#merge
	if(nrow(mismatches)>1){               
	for(i in 2:nrow(mismatches)){
		if(mismatches$refseq[i]=="-"){ #insertion
			mismatches$start[i] = mismatches$start[i] +1
			if(reformMismatches$refseq[nrow(reformMismatches)]=="-" & reformMismatches$start[nrow(reformMismatches)]==mismatches$start[i] & reformMismatches$fid[nrow(reformMismatches)]==mismatches$fid[i]){ #merge continuous refseq "-"
				reformMismatches$conseq[nrow(reformMismatches)]=paste(reformMismatches$conseq[nrow(reformMismatches)],mismatches$conseq[i],sep="")
				reformMismatches$alnend[nrow(reformMismatches)]=mismatches$alnend[i]
				next
			} 
		}
		else if(mismatches$conseq[i]=="-"){ #deletion
			if(reformMismatches$conseq[nrow(reformMismatches)]=="-" & reformMismatches$end[nrow(reformMismatches)]==mismatches$start[i]-1 & reformMismatches$fid[nrow(reformMismatches)]==mismatches$fid[i]){
				reformMismatches$refseq[nrow(reformMismatches)]=paste(reformMismatches$refseq[nrow(reformMismatches)],mismatches$refseq[i],sep="")
				reformMismatches$end[nrow(reformMismatches)] = mismatches$end[i]
				reformMismatches$alnend[nrow(reformMismatches)]=mismatches$alnend[i]
				next
			}	
	
		
		}
		else if(mismatches$conseq[i]=="N"){	 # unknown
			if(reformMismatches$conseq[nrow(reformMismatches)]=="N" & reformMismatches$end[nrow(reformMismatches)]==mismatches$start[i]-1 & reformMismatches$fid[nrow(reformMismatches)]==mismatches$fid[i]){
				reformMismatches$refseq[nrow(reformMismatches)]=paste(reformMismatches$refseq[nrow(reformMismatches)],mismatches$refseq[i],sep="")
				reformMismatches$end[nrow(reformMismatches)] = mismatches$end[i]
				reformMismatches$alnend[nrow(reformMismatches)]=mismatches$alnend[i]
				next
			}	
	
		
		}
		else if(mismatches$conseq[i]==""){ #extension
			if(reformMismatches$conseq[nrow(reformMismatches)]=="" & reformMismatches$end[nrow(reformMismatches)]==mismatches$start[i]-1 & reformMismatches$fid[nrow(reformMismatches)]==mismatches$fid[i]){
				reformMismatches$refseq[nrow(reformMismatches)]=""
				reformMismatches$end[nrow(reformMismatches)] = mismatches$end[i]
				reformMismatches$alnend[nrow(reformMismatches)]=mismatches$alnend[i]
				next
			}	
	
		
		}
		reformMismatches[nrow(reformMismatches)+1,]<- mismatches[i,]  
	}
	}
	                 
	 reformMismatches$reftype[reformMismatches$conseq==""]="ext"


	############################ deal with ext
	tmpmis=reformMismatches[reformMismatches$reftype!="nuc",]
	if(nrow(tmpmis)>0){
	tmpmis2<- data.frame(reftype=factor(levels=c("gen","nuc","begin","end","ext")),
	                 refseq=character(),
	                 conseq=character(),
	                 start=integer(),
	                 end=integer(),
	                 alnstart=integer(),
	                 alnend=integer(),
	                 fid=character(),
	                 feature_start = integer(),
					feature=character(),
	                 stringsAsFactors=FALSE )
	 

	 tmpmis2[nrow(tmpmis2)+1,]<- tmpmis[1,]
	if(nrow(tmpmis)>1){               
		for(i in 2:nrow(tmpmis)){
#		print(tmpmis2$end[nrow(tmpmis2)])
#		print(tmpmis$start[i])
#		print("gogo")
			if(tmpmis2$end[nrow(tmpmis2)]+1==tmpmis$start[i]){
				tmpmis2$end[nrow(tmpmis2)]=tmpmis$end[i]
				tmpmis2$alnend[nrow(tmpmis2)]=tmpmis$alnend[i]
			}
			else{
	 			tmpmis2[nrow(tmpmis2)+1,]<- tmpmis[i,]

			}
		}
	}
	
	tmpmis2$reftype[]="ext"
	tmpmis2$conseq=""
	tmpmis2$refseq=""
	tmpmis=reformMismatches[reformMismatches$reftype!="ext",]
	reformMismatches=rbind(tmpmis,tmpmis2)

	}

	#############################

	reformMismatches$start=reformMismatches$start-reformMismatches$feature_start
	reformMismatches$end=reformMismatches$end-reformMismatches$feature_start+1
	
	#replace sequence with "Entire feature" if it is so
	for(i in 1:nrow(reformMismatches)){
		feature=hlaData$features[hlaData$features$fid==reformMismatches$fid[i],]
		if(reformMismatches$start[i]==0 && reformMismatches$end[i]==(feature$end-feature$start+1) && (reformMismatches$conseq[i]=="-" || reformMismatches$conseq[i]=="N")){
			reformMismatches$refseq[i]="Entire feature"
		}
	

	}
	levels(reformMismatches$reftype)=c("W","C","U","D","E")
	# "17ws reference","closest allele","upstream unreferenced","downstream unreferenced","seq extension"
	reformMismatches
}
