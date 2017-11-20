hlaGetNovel <-
function(hlaData,nucAllele,conseq,sposition=0)
{
	locus = sub("\\*.*","",nucAllele)
	allele2d = sub(":.*","",nucAllele)
	if(is.na(hlaData$refAllele[allele2d])){
       genAllele = as.character(hlaData$refAllele[locus])
	}
	else{
       genAllele = as.character(hlaData$refAllele[allele2d])

	}
	hlaData<-hlaGetAlign(hlaData,genAllele,nucAllele,conseq)
	poly = hlaAlign2Novel(hlaData)
	return(poly)
	#get sequence according to input
}
