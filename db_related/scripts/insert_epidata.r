insert_epidata<-function(labcode,epidata){
require("ORE")
require("plyr")
	colnames(epidata$main_report)=toupper(colnames(epidata$main_report))
	colnames(epidata$rawdata)=toupper(colnames(epidata$rawdata))
	if(!is.null(epidata$coninfo) ){
	colnames(epidata$coninfo)=toupper(colnames(epidata$coninfo))
	}
	if(!is.null(epidata$calinfo) ){
	colnames(epidata$calinfo)=toupper(colnames(epidata$calinfo))
	}
	colnames(epidata$fusion_data)=toupper(colnames(epidata$fusion_data))
	colnames(epidata$fusion_profile)=toupper(colnames(epidata$fusion_profile))
	
	
	ore.connect(user = "hladb", sid = "hd", host = "localhost", password = "#######")
	if(!ore.is.connected()){
		print("cannot connect");
		q()
	}
	ore.attach()
	ore.sync( 
	      query=c("EKEYS" = sprintf("select * from EPITOPE_KEY WHERE labcode like '%s' and reportname like '%s'",labcode,epidata$report_name)
	      )
	)
	
	if(nrow(EKEYS)>0){
		mykey<-ore.pull(EKEYS) 
		ore.exec(sprintf("delete from EPITOPE_REPORT where EREPORTID=%d",mykey$EREPORTID[1]))
		ore.exec(sprintf("delete from EPITOPE_RAW where EREPORTID=%d",mykey$EREPORTID[1]))
		ore.exec(sprintf("delete from EPITOPE_RAW where EREPORTID=%d",mykey$EREPORTID[1]))
		ore.exec(sprintf("delete from EPITOPE_CALINFO where EPITOPE_CALID=%d",mykey$CALINFO1[1]))
		ore.exec(sprintf("delete from EPITOPE_CALINFO where EPITOPE_CALID=%d",mykey$CALINFO2[1]))
		ore.exec(sprintf("delete from EPITOPE_CONINFO where EPITOPE_CONID=%d",mykey$CONINFO1[1]))
		ore.exec(sprintf("delete from EPITOPE_CONINFO where EPITOPE_CONID=%d",mykey$CONINFO2[1]))
		ore.exec(sprintf("delete from FUSION_DATA where EREPORTID=%d",mykey$EREPORTID[1]))
		ore.exec(sprintf("delete from FUSION_PROFILE where EREPORTID=%d",mykey$EREPORTID[1]))
	
	}
	else{
		ore.exec(sprintf("insert into EPITOPE_KEY (labcode,reportname) values ('%s','%s')",labcode,epidata$report_name))
	
	}
	ore.sync( 
	      query=c("ILAB" = sprintf("select * from ihiw_lab where labcode like '%s'",labcode),
	              "EREPORT" = "select * from EPITOPE_REPORT where rownum=0",              
	              "ERAW" = "select * from EPITOPE_RAW where rownum=0",
	              "ECONINFO" = "select * from EPITOPE_CONINFO where rownum=0",
	              "ECALINFO" = "select * from EPITOPE_CALINFO where rownum=0",
	               "EFUSIOND" = "SELECT * FROM FUSION_DATA WHERE rownum=0",
	               "EFUSIONP" = "SELECT * FROM FUSION_PROFILE WHERE rownum=0",
	               "CTABLE" = "select * from codetable where typeid=12")
	               
	            
	              
	      			
	      
	      )
	ore.drop("TMP_EREPORT")
	ore.drop("TMP_RAW")
	ore.drop("TMP_ECONIN")
	ore.drop("TMP_ECALIN")
	ore.drop("TMP_EFUSIOND")
	ore.drop("TMP_EFUSIONP")
	 mylab<-ore.pull(ILAB)  
	 mykey<-ore.pull(EKEYS)   
	tmp<-ore.pull(EREPORT)
	class(tmp$REPORT_DATE)="character"
	tmp2<-rbind.fill(tmp,epidata$main_report)
	tmp<-tmp2[,colnames(tmp2) %in% colnames(tmp)]
	tmp$LABID=mylab$LABID
	tmp$REPORTNAME=epidata$report_name
	tmp$EREPORTID=mykey$EREPORTID
	tmp$CALINFO1=mykey$CALINFO1
	tmp$CALINFO2=mykey$CALINFO2
	tmp$CONINFO1=mykey$CONINFO1
	tmp$CONINFO2=mykey$CONINFO2
	ore.create(tmp,"TMP_EREPORT")
	ore.exec("insert into EPITOPE_REPORT select * from TMP_EREPORT")
	
	colnames(epidata$rawdata)[colnames(epidata$rawdata)=="TOTAL.EVENTS"]="TOTAL_EVENTS"
	tmp<-ore.pull(ERAW)
	tmp<-rbind.fill(tmp,epidata$rawdata)
	tmp$EREPORTID=mykey$EREPORTID
	myctable<-ore.pull(CTABLE)
	rownames(myctable)=myctable$NAME
	tmp$DATATYPE<-myctable[tmp$DATATYPE,"MTHID"]
	ore.create(tmp,"TMP_RAW")
	ore.exec("insert into EPITOPE_RAW select * from TMP_RAW")
	
	if(!is.null(epidata$coninfo) ){
	tmp<-ore.pull(ECONINFO)
	tmp<-rbind.fill(tmp,epidata$coninfo)
	tmp$EPITOPE_CONID[1]=mykey$CONINFO1
	tmp$EPITOPE_CONID[2]=mykey$CONINFO2
	ore.create(tmp,"TMP_ECONIN")
	ore.exec("insert into EPITOPE_CONINFO select * from TMP_ECONIN")
	}
	if(!is.null(epidata$calinfo) ){
	tmp<-ore.pull(ECALINFO)
	tmp<-rbind.fill(tmp,epidata$calinfo)
	tmp$EPITOPE_CALID[1]=mykey$CALINFO1
	tmp$EPITOPE_CALID[2]=mykey$CALINFO2
	ore.create(tmp,"TMP_ECALIN")
	ore.exec("insert into EPITOPE_CALINFO select * from TMP_ECALIN")
	}
	
	tmp<-ore.pull(EFUSIOND)
	tmp<-rbind.fill(tmp,epidata$fusion_data)
	tmp$EREPORTID=mykey$EREPORTID
	ore.create(tmp,"TMP_EFUSIOND")
	ore.exec("insert into FUSION_DATA select * from TMP_EFUSIOND")

	tmp<-ore.pull(EFUSIONP)
	class(tmp$TEST_DATE)="character"
	tmp<-rbind.fill(tmp,epidata$fusion_profile)
	tmp$EREPORTID=mykey$EREPORTID
	ore.create(tmp,"TMP_EFUSIONP")
	ore.exec("insert into FUSION_PROFILE select * from TMP_EFUSIONP")
	ore.detach()
	ore.drop("TMP_EREPORT")
	ore.drop("TMP_RAW")
	ore.drop("TMP_ECONIN")
	ore.drop("TMP_ECALIN")
	ore.drop("TMP_EFUSIOND")
	ore.drop("TMP_EFUSIONP")
	ore.exec("commit")
	ore.disconnect()
	return(1)
}
