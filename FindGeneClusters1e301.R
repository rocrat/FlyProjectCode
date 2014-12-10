#This script takes a initial HMMER results between pairs of spp. which must be concatenated together into a single file using the 799 script
library(MASS)
library(linkcomm)
dat<-read.delim("C:/Projects/FlyEvolution/CombinedHMMER.csv",sep=";",header=TRUE)
d1<-dim(dat)[1]
dat<-dat[as.character(dat$id1)!=as.character(dat$id2),]
#subset data for efficiency
dat$weight<-1-dat$evalue
#selected the cut-off of 1e-301 based on a plot of the ecdf of e-values
dat1<-dat[dat$evalue<=1e-301,]
system.time(lc<-getLinkCommunities(dat1[,c(1,2,4)]))
rm(dat)
save.image(file="C:/Projects/FlyEvolution/LinkClusters1e301.Rdata")

linkmat<-matrix("",length(lc$clusters),1)
for(i in 1:length(lc$clusters)){
  linkmat[i,1]<-paste(getNodesIn(lc,clusterids=i),collapse=";")
}
write.csv(linkmat,"ClusterMatrix.csv",)
save.image(file="/xdisk/rlapoint/LinkMatrix1e301.Rdata")
