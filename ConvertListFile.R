#this script takes a .list file and adds 4 character species names before each geneID
dat<-read.delim("C:\\Users\\DominicLaptop\\Documents\\Work\\FlyEvolution\\FastaFilesForAlignmentV1\\RecipBesthitsV3NoSnot.list",header=FALSE,sep=";")
out.dat<-data.frame(paste0("Dmel.",dat$V1),
                    paste0("Dbia.",dat$V2),
                    paste0("Dgri.",dat$V3),
                    paste0("Dmoj.",dat$V4),
                    paste0("Dpse.",dat$V5),
                    paste0("Dsuz.",dat$V6),
                    paste0("Dyak.",dat$V7),
                    paste0("Sfla.",dat$V8))
library(MASS)
write.table(out.dat,file="C:\\Users\\DominicLaptop\\Documents\\Work\\FlyEvolution\\FastaFilesForAlignmentV1\\ConvertedRecipListV4.list",row.names=FALSE,col.names=FALSE,quote=FALSE,sep=";")
