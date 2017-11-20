#!/usr/bin/env Rscript 
require("ORE")
require("plyr")
watchfolders=c("hml","gendx","epitopes")
root="/w_u"
labcodes<-dir(root)
ore.connect(user = "hladb", sid = "hd", host = "localhost", password = "oracle12345")
if(!ore.is.connected()){
	print("cannot connect");
	q()
}
ore.attach()
ore.sync(table = c("SFTP_DATALOG"))

while(TRUE){
labcodes<-dir(root)
datalog<-ore.pull(SFTP_DATALOG)
for(labcode in labcodes){
    for(watchfolder in watchfolders){
	path=file.path(root,labcode,"upload",watchfolder)
	files=list.files(path,pattern="*.(xml|hml|xls)$",full.names=T,ignore.case =T,recursive=T)
#	files=c(files,list.files(path,pattern="*.xml",full.names=T))
	
	for(file in files){
		fext=tolower(sub(".*\\.","",file))
		finfo<-file.info(file)
		mtime=as.character(finfo$mtime)
		if(finfo$isdir){
			next
		}
				
		if(as.numeric(Sys.time()- as.POSIXct(finfo$mtime), units = "mins")<5){
			next
		}	
			
		if(length(which(datalog$LABCODE==labcode &  datalog$FOLDER==file & datalog$MTIME==mtime))>0){
			next
		}
		if(watchfolder=="hml" & (fext=="xml" | fext=="hml")){
			prog=paste("/home/oracle/scripts/parse_hml.pl",labcode,paste("\"",file,"\"",sep=""),"2>&1",sep=" ")
			dtype="H"
		}
		else if(watchfolder=="gendx" & fext=="xml" ){
			prog=paste("/home/oracle/scripts/parse_gendx.pl",labcode,paste("\"",file,"\"",sep=""),"2>&1",sep=" ")
			dtype="G"
		}
		else if(watchfolder=="epitopes" & fext=="xls" & substring(basename(file),1,1)=="H"){
			prog=paste("/home/oracle/scripts/pi_epidata.r",labcode,paste("\"",file,"\"",sep=""),"2>&1",sep=" ")
			dtype="E"
		}
		else{
			next
		}
			
		print(prog)
		result=try(system(prog,intern=T))
		if(is.null(attr(result,"status"))){
			result="OK"
		}
		else{
			result=paste(result,collapse=" ")
			result=substr(result, nchar(result)-200+1, nchar(result))
			result=gsub("'","\"",result)

		}
		
		cmd=paste("insert into SFTP_DATALOG (datatype,labcode,folder,mtime,result) values ('",dtype,"','",labcode,"','",file,"','",mtime,"','",result,"')",sep="")
		print(cmd)
		ore.exec(cmd)
		ore.exec("commit")
	}
    }
}
Sys.sleep(300)
}   
