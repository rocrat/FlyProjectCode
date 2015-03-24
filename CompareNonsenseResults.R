#compare convergent genes between Ha and nonsense Ha
#read in both gene results files
non<-read.csv("C:/Projects/FlyProjectCode/treeSimPvalueNonsenseV1.csv")
ha<-read.table("C:/Projects/FlyEvolution/driftTreeNoSnothaV3.txt",header=T,sep="\t")
#reduce to only those genes that pass results
inc_non<-with(non,which(pvalue<0.05 & TreePvalue<0.05))
inc_ha<-with(ha,which(pvalue<0.05 & TreePvalue<0.05 & DriftPvalue<0.05))
