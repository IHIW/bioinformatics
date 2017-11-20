require("plyr")
parse_epitope<-function(epifolder){
	#prepare files
	report_name=regmatches(epifolder,regexpr("[\\w\\s]+[/]?$", epifolder, perl=TRUE))
	main_report_file<-file.path(epifolder,"Output.csv")
	fusion_files<-file.path(epifolder,list.files(epifolder,pattern="*.[Xx][Ll][Ss][xX]?$"))
	report_lines<-readLines(main_report_file,warn=F)
	
	#find locations
	datatype_pos<-NULL
	sep=","
	if(grepl(";",report_lines[1])){
		sep=";"
	}
	stage=0
	breakPoint1=0
	cal_pos=0
	con_pos=0
	ass_pos=0
	for(i in 1:length(report_lines)){
		if(stage==0 && grepl(":[;,\"]",report_lines[i])){
			stage=stage+1
			breakPoint1=i;

		}
		if(grepl("CALInfo:",report_lines[i])){
			cal_pos=i
		}
		else if(grepl("CONInfo:",report_lines[i])){
			con_pos=i
		}
		else if(grepl("AssayLotInfo:",report_lines[i])){
			ass_pos=i
		}
		else if(grepl("DataType:",report_lines[i])){
			datatype_pos<-c(datatype_pos,i)
		}
		
	}
	datatype_pos<-c(datatype_pos,length(report_lines))
	
	#main_report
	tmp<-read.csv(text=report_lines[1:(breakPoint1-1)],header=F,stringsAsFactors=F,sep=sep)
	tmp=tmp[tmp$V1!="",]
	rownames(tmp)=tmp$V1
	rownames(tmp)[tmp$V1=="Date"]="REPORT_DATE"
	rownames(tmp)[tmp$V1=="Session"]="REPORT_SESSION"
	# tmp$V2[tmp$V1=="Date"]=paste(tmp$V2[tmp$V1=="Date"],tmp$V3[tmp$V1=="Date"],sep=" ")
	tmp$V2[tmp$V1=="Date"]=tmp$V2[tmp$V1=="Date"]

	main_report<-as.data.frame(t(tmp[,1:2]))[2,]
	colnames(main_report)=sub(" ","_",colnames(main_report))
	
	#infos
	if(cal_pos>0 && con_pos>0 && ass_pos>0){
		calinfo<-read.csv(text=report_lines[(cal_pos+1):(con_pos-1)],header=T,stringsAsFactors=F,sep=sep)
		coninfo<-read.csv(text=report_lines[(con_pos+1):(ass_pos-1)],header=T,stringsAsFactors=F,sep=sep)
	}
	else{
		calinfo=NULL
		coninfo=NULL
	}
	
	#raw data
	rawdata<-NULL
	runs<-NULL
	for(i in 1:(length(datatype_pos)-1)){
		datatype=gsub("\"","",sub(paste(".*",sep,sep=""),"", report_lines[datatype_pos[i]]))
		tmp=read.csv(text=report_lines[(datatype_pos[i]+1):(datatype_pos[i+1]-1)],header=T,stringsAsFactors=F,sep=sep)
		if(nrow(tmp)>0 && colnames(tmp)[1]=="Location" && colnames(tmp)[2]=="Sample" && length(grep("^X",colnames(tmp)[3]))>0){
			tmp$datatype=datatype
			tmp$runid=1:nrow(tmp)
			tmp<-tmp[tmp$Sample!="",]
			rawdata<-rbind.fill(rawdata,tmp)
		}
		if(i==1){
			runs<-tmp[,c("runid","Sample")]
			rownames(runs)=runs$Sample
		}
	
	}
	rawdata$labsample<-rawdata$Sample #regmatches(rawdata$Sample,regexpr("\\d+", rawdata$Sample, perl=TRUE))
	for(col in grep("^X",colnames(rawdata))){
		colnames(rawdata)[col]=paste("X",sprintf("%03d",as.integer(sub("X","",colnames(rawdata)[col]))),sep="")
		
 		rawdata[,col]=sub(",",".",rawdata[,col])
 	}
	if(!is.null(calinfo)){
	for(col in 1:ncol(calinfo)){
 		calinfo[,col]=sub(",",".",calinfo[,col])
 	}
	}
	if(!is.null(coninfo)){
	for(col in 1:ncol(coninfo)){
 		coninfo[,col]=sub(",",".",coninfo[,col])
 	}
	}

	#fusion data
	fusion_data<-NULL
	fusion_profile<-NULL
	for(i in 1:length(fusion_files)){
		tmp = read.fusion_xls(fusion_files[i])
		runid=runs[tmp$fusion_profile$sampleid,"runid"]
		if(is.null(runid)){
			next
		}
		tmp$fusion_profile$runid=runid
		tmp$fusion_data$runid=runid
		fusion_profile<-rbind(fusion_profile,tmp$fusion_profile)
		fusion_data<-rbind(fusion_data,tmp$fusion_data)
	}
	if(!is.null(coninfo)){
	coninfo<-coninfo[,colSums(is.na(coninfo))<nrow(coninfo)]
	}
	if(!is.null(calinfo)){
	calinfo<-calinfo[,colSums(is.na(calinfo))<nrow(calinfo)]
	}
	rawdata<-rawdata[,colSums(is.na(rawdata))<nrow(rawdata)]
	
	epitopeDataframe<-list(runs=runs,report_name=report_name,main_report=main_report,calinfo=calinfo,coninfo=coninfo,rawdata=rawdata,fusion_data=fusion_data,fusion_profile=fusion_profile)
	epitopeDataframe
}

read.fusion_xls<-function(fusion_file){
	require("xlsx")
	tmp<-read.xlsx(fusion_file,1,header=F,stringsAsFactors=F)
	sampleid=tmp$X2[1]
	epitope_session<-tmp$X7[3]
	catelog<-tmp$X7[4]
	nc_bead<-tmp$X7[5]
	test_date<-tmp$X39[3]
	formula<-tmp$X39[4]
	threshold<-tmp$X39[5]
	pos<-tmp$X32[3]
	nc_sample<-tmp$X18[5]
	eversion<-sub("™","",tmp$X8[nrow(tmp)-1])
	fusion_profile<-data.frame(sampleid=sampleid,epitope_session=epitope_session,catelog=catelog,nc_bead=nc_bead,
	test_date=test_date,formula=formula,threshold=threshold,pos=pos,nc_sample=nc_sample,eversion=eversion,
	stringsAsFactors=F)
	fusion_data<-tmp[c(8:(nrow(tmp)-3)),c(1,4,9,12,13,16,21,23,27,30,36)]
	colnames(fusion_data)=c("beadid","raw_value","sample_nc","ns_raw","nsnc","normal","ratio","rxn","count","specificity","allele_specificity")
	fusion_data$rxn[fusion_data$rxn=="NC"]=0
	fusion_data$rxn[fusion_data$rxn=="PC"]=16
	fusion_data$ratio=sub(",",".",fusion_data$ratio)
	fusion_data$raw_value=sub(",",".",fusion_data$raw_value)
	fusion_data$sample_nc=sub(",",".",fusion_data$sample_nc)
	fusion_data$ns_raw=sub(",",".",fusion_data$ns_raw)
	fusion_data$nsnc=sub(",",".",fusion_data$nsnc)
	fusion_data$normal=sub(",",".",fusion_data$normal)
	fusion_data$beadid=paste("X",fusion_data$beadid,sep="")
	list(fusion_data=fusion_data,fusion_profile=fusion_profile)
}
