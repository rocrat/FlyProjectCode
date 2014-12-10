#This script takes a initial HMMER results between pairs of spp. which must be concatenated together into a single file using the 799 script
library(MASS)
library(linkcomm)
dat<-read.delim("/xdisk/rlapoint/CombinedHMMER.csv",sep=";",header=TRUE)
dat<-dat[dat$id1!=dat$id2]
#subset data for testing
dat$weight<-1-dat$evalue
lc<-getLinkCommunities(dat[,c(1,2,4)])
rm(dat)
save.image(file="/xdisk/rlapoint/LinkClusters.Rdata")

linkmat<-matrix("",length(lc$clusters),1)
for(i in 1:length(lc$clusters)){
  linkmat[i,1]<-paste(getNodesIn(lc,clusterids=i),collapse=";")
}
write.csv(linkmat,"ClusterMatrix.csv")
save.image(file="/xdisk/rlapoint/LinkMatrix.Rdata")