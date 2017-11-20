a<-read.table("tmp.txt",sep="\t",header=T,stringsAsFactors=F)
a3250=a$X3250[!is.na(a$X3250)]


ruleas=which(is.na(a$X3250) & !is.na(a$X3280))

for(i in ruleas){
s=a$X3280[i]
news=NA
while(grepl(":",s)){
	s=sub("[^:]+$","",s,perl=T)
	s1=sub("\\*","\\\\*",s)
	r=grep(s1,a3250)
	if(length(r)>0){
		news=a3250[r[1]]
		break
	}
	s=sub(":$","",s)
}
print(news)
a$X3250[i]=news
}

rulebs=which(is.na(a$X3250) & grep("DPB1",a$X3280) )
a$X3250[rulebs]="HLA-DPB1*02:01:02"
