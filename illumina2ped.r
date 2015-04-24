illu_ped.raw<-function(fnm, sep='\t'  ,snpc=snp_name, n=NULL, p=NULL, subc=labid,getKey=TRUE){
    ### Read data	
	snpdata<-read.table(file=fnm, header=TRUE, sep='\t', stringsAsFactors=FALSE )
	### identify columns with snps and subjects
    snpc<- which(colnames(snpdata)==snpc)
    subc<- which(colnames(snpdata)==subc)
    ###start final matrix and define n and p
    snpdata[,snpc]<-as.character(snpdata[,snpc])
    if(is.null(n)) n<-length(ids<-  unique(sort(as.numeric(snpdata[,subc])  )) )
    if(is.null(p)) p<-length(snpid<-unique(snpdata[,snpc]))
    ped.raw  <- matrix(NA, nrow=n, ncol=p)
    rownames(ped.raw)<-ids
    colnames(ped.raw)<-snpid
    if(getKey){ idKey <-NULL
       mafA<-mafB<-character(p)
	}    
	sid<-raceCode<-character(n)
    for(i in 1:p){ 
	    if(i%%5==0) print(paste(round(i/p*100,0), '%',sep=''))   #i=1
        ###get data for one snp
        snp.j   <-snpid[i]
        aims.one<- snpdata[snpdata[,snpc]==snp.j,]
        aims.one<- aims.one[order(aims.one[,subc]),]
        tmp <- aims.one[,grep(patter='allele', x=colnames(aims.one))]
        rownames(tmp) <- aims.one[,subc]
        tmp<-tmp[order(as.numeric(rownames(tmp))),]
        tmp[,1]<-ifelse(tmp[,1]=='-', NA,tmp[,1]) 
        tmp[,2]<-ifelse(tmp[,2]=='-', NA,tmp[,2]) 
        colnames(tmp)<-c(paste(snp.j, '-', 1:2, sep=''))
        if(is.null(idKey) & nrow(aims.one)==n) idKey<-aims.one[,c(1:4,16:18)]
        idSnp.one<-rownames(tmp)
        (count<-table(c(tmp[,1], tmp[,2])))
        if(length(count)>1){
            fstmax <- count[1]>count[2]
            if(!fstmax) count <- count[2:1] 
            mafA[i]<-names(count[2]); mafB[i]<-names(count[1])
            tmp$x.maf<- as.numeric(tmp[,1]==mafA[i])+as.numeric(tmp[,2]==mafA[i])
            for(j in 1:nrow(tmp)){         #j=1
           		idtmp<-rownames(tmp)[j] 
           		cid<-which(rownames(ped.raw)==idtmp)
           		if(length(cid)>0) ped.raw[cid,i]<-tmp$x.maf[j]
        }} 
        if(length(count)==1) { 
     	    mafA[i]	<-names(count[1]); mafB[i]<- NA        
     	    rtmp <- which(rownames(ped.raw) %in% rownames(tmp))
     	    ped.raw[rtmp,i]<- 0;	
}}}
