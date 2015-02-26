#this script reads in a series of fasta files and spits out a series of PDF's with AA alignemnts with highlighted convgergent substitutions
library("plyr")
library("ggplot2")
#get file names
fastaFiles<-dir("C:/Projects/FlyEvolution/TopCleanAlignmentsV3/",full.names = TRUE)
#limit to just aa alignments
fastaFiles<-fastaFiles[grepl("V.aa_",fastaFiles)]
#i<-81
for(i in 1:length(fastaFiles)){
  #read in fasta file
  fasta<-readLines(con=fastaFiles[i])
  #loop through file and segregate aa's and headers
  h<-aaStart<-aaStop<-vector()
  for(j in 1:length(fasta)){
    if(grepl(">",fasta[j])){
      h<-c(h,sub(".+,(\\w{4})","\\1",fasta[j]))
      aaStart<-c(aaStart,j+1)
      if(j>1){
        aaStop<-c(aaStop,j-1)
      }
    }
    
  }
  #Add final stop postion for aa alignments
  aaStop<-c(aaStop,length(fasta))
  #bind the aa's together
  AAs<-character()
  for(k in 1:length(h)){
    AAs[k]<-paste(fasta[aaStart[k]:aaStop[k]],collapse="")
  }
  #split AA's to make columns of AA's
  columns<-list()
  for(l in 1:length(h)){
    columns[[l]]<-laply(1:nchar(AAs[l]), function(p) substr(AAs[l], p, p))
  }
  #create data frame from columns
  df<-data.frame(do.call(cbind,columns))
  names(df)<-h
  df$convergence<-with(df,ifelse(as.character(Dsuz)==as.character(Sfla) & as.character(Dsuz)!=as.character(Dpse) & as.character(Dsuz)!=as.character(Dbia) & as.character(Dsuz)!=as.character(Dmel) & as.character(Dsuz)!=as.character(Dyak) & as.character(Dsuz)!=as.character(Dmoj) & as.character(Dsuz)!=as.character(Dgri),"Yes","No"))
  df$position<-factor(1:dim(df)[1])
  longdf<-reshape(df,varying=c("Dgri","Dmoj","Dmel","Dbia","Dyak","Dpse","Dsuz","Sfla"),times=c("Dgri","Dmoj","Dmel","Dbia","Dyak","Dpse","Dsuz","Sfla"),direction="long",v.names="AA")
  longdf$time<-factor(longdf$time,levels=c("Sfla","Dsuz","Dgri","Dmoj","Dmel","Dbia","Dyak","Dpse"))
  plot1<-ggplot(data=longdf,aes(x=position,y=time))+geom_tile(aes(fill=convergence),alpha=.5)+ylab("Species")+xlab("AA Postion")+ggtitle(sub(".+(FBgn\\d+)V.+","\\1",fastaFiles[i],perl=TRUE))+scale_fill_manual(values=c("grey","red"))+geom_text(aes(label=AA))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="bottom")
  
  ggsave(filename=paste0("C:/Projects/FlyEvolution/NewPlots/",sub(".+(FBgn\\d+)V.+","\\1",fastaFiles[i],perl=TRUE),".pdf"),plot=plot1,width=as.integer(.3*length(levels(df$position))),height=4,units="in",limitsize=FALSE)
}
