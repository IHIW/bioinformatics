insert_epidata<-function(labcode,epidata){
require("ORE")
require("plyr")
	colnames(epidata$fusion_data)=toupper(colnames(epidata$fusion_data))
	colnames(epidata$fusion_profile)=toupper(colnames(epidata$fusion_profile))
	
	
	ore.connect(user = "hladb", sid = "hd", host = "localhost", password = "##########")
	if(!ore.is.connected()){
		print("cannot connect");
		q()
	}
	ore.attach()
	ore.exec(sprintf("delete from FUSION_DATA2 where sampleid='%s'",epidata$fusion_profile$SAMPLEID))
	ore.exec(sprintf("delete from FUSION_PROFILE2 where sampleid='%s'",epidata$fusion_profile$SAMPLEID))
	
	ore.sync( 
	      query=c("ILAB" = sprintf("select * from ihiw_lab where labcode like '%s'",labcode),
	               "EFUSIOND" = "SELECT * FROM FUSION_DATA2 WHERE rownum=0",
	               "EFUSIONP" = "SELECT * FROM FUSION_PROFILE2 WHERE rownum=0")
	      )
	ore.drop("TMP_EFUSIOND")
	ore.drop("TMP_EFUSIONP")
	 mylab<-ore.pull(ILAB)  
	
	
	
	tmp<-ore.pull(EFUSIOND)
	tmp<-rbind.fill(tmp,epidata$fusion_data)
	tmp$LABID=mylab$LABID
	ore.create(tmp,"TMP_EFUSIOND")
	ore.exec("insert into FUSION_DATA2 select * from TMP_EFUSIOND")

	tmp<-ore.pull(EFUSIONP)
	class(tmp$TEST_DATE)="character"
	tmp<-rbind.fill(tmp,epidata$fusion_profile)
	tmp$LABID=mylab$LABID
	ore.create(tmp,"TMP_EFUSIONP")
	ore.exec("insert into FUSION_PROFILE2 select * from TMP_EFUSIONP")

	ore.detach()
	ore.drop("TMP_EFUSIOND")
	ore.drop("TMP_EFUSIONP")
	ore.exec("commit")
	ore.disconnect()
	return(1)
}
