#This script takes a initial HMMER results between pairs of spp. which must be concatenated together into a single file using the 799 script
library(MASS)
library(linkcomm)
dat<-read.delim("C:/Projects/FlyEvolution/CombinedHMMER.csv",sep=";",header=TRUE)
d1<-dim(dat)[1]
dat<-dat[as.character(dat$id1)!=as.character(dat$id2),]
#subset data for testing
dat$weight<-1-dat$evalue
dat1<-dat[sample(dim(dat)[1],10000),]
lc<-getLinkCommunities(dat1[,c(1,2,4)])
rm(dat)
save.image(file="C:/Projects/FlyEvolution/LinkClusterstest.Rdata")

linkmat<-matrix("",length(lc$clusters),1)
for(i in 1:length(lc$clusters)){
  linkmat[i,1]<-paste(getNodesIn(lc,clusterids=i),collapse=";")
}
write.csv(linkmat,"ClusterMatrix.csv",)
save.image(file="/xdisk/rlapoint/LinkMatrix.Rdata")

plot(lc, type = "graph", layout = layout.fruchterman.reingold,vlabel=F)
