genes<-read.table("C:/Projects/FlyEvolution/finalGeneSSLSresults.txt",sep="\t")
genelist<-genes[,1]
write.csv(genelist,"C:/Projects/FlyEvolution/finalGeneList.csv",row.names=FALSE,col.names=FALSE)
