rcl<-read.table("C:\\Users\\DominicLaptop\\Documents\\Work\\FlyEvolution\\FastaFilesForAlignmentV1\\CombinedRecipClust.key")
library(MASS)
truehist(rcl$V1,xlab="Number of genes/file after adding clusters")
sum(rcl$V1>60)
spnk<-read.csv("C:\\Users\\DominicLaptop\\Documents\\Work\\FlyEvolution\\FastaFilesForAlignmentV1\\SppNumberKey.csv",header=FALSE)
funq<-read.table("C:\\Users\\DominicLaptop\\Documents\\Work\\FlyEvolution\\FastaFilesForAlignmentV1\\FastaUniqueSPPV4",header=FALSE,sep="\n")
for (i in 1:dim(funq)[1]){
  funq$num[i]<-strsplit(x=as.character(funq$V1)[i],split=";")[[1]][1]
}
funq$num<-as.numeric(as.character(funq$num))
sum(funq$num==8)
truehist(funq$num)
#create indicators for each spp
funq$dmel<-grepl("Dmel",funq$V1)
funq$dbia<-grepl("Dbia",funq$V1)
funq$dgri<-grepl("Dgri",funq$V1)
funq$dpse<-grepl("Dpse",funq$V1)
funq$dyak<-grepl("Dyak",funq$V1)
funq$dsuz<-grepl("Dsuz",funq$V1)
funq$dmoj<-grepl("Dmoj",funq$V1)
funq$slfa<-grepl("Sfla",funq$V1)
summary(funq)#missing Sfla in the majority of cases

ffkey<-read.csv("C:\\Users\\DominicLaptop\\Documents\\Work\\FlyEvolution\\FastaFilesForAlignmentV1\\FastaFilesKeyV4NoSnot.txt",header=FALSE)
ffkey$V2<-as.numeric(as.character(ffkey$V2))
truehist(ffkey$V2,xlab="Number of genes/fasta file",main="Recip combined with cluster")
sum(ffkey$V2>=8)

ffkey3<-read.csv("C:\\Users\\DominicLaptop\\Documents\\Work\\FlyEvolution\\FastaFilesForAlignmentV1\\FastaFilesKeyV3NoSnot.txt",header=FALSE)
ffkey3$V2<-as.numeric(as.character(ffkey3$V2))
truehist(ffkey3$V2,xlab="Number of genes/fasta file",main="Recip w/out Cluster",prob=T)
sum(ffkey3$V2==8)
