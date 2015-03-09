
### Parameters
# sitesDir <- "~/projects/flyCE/garliNoSnothaV3/"
#sitesDir <- "~/projects/flyCE/garliNoSnothaAAV3/"
sitesDir <- "C:/Projects/FlyEvolution/garliNonsenseV1/"
###load packages
library(snowfall)
library(parallel)

#read in files
sitesFiles <- dir(sitesDir,full.names = T)
sitesFiles <- sitesFiles[which(grepl("sitelikes",sitesFiles))]
##### Compute P value for each gene
## function testing parameters
#tmpFile <- sitesFiles[1]

compute.gene.pvalue <- function(tmpFile){
    tmpData <- read.table(tmpFile,header = T,sep = "\t",fill=T)
    ## compute deltaSSLS for given gene
    numBP <- max(tmpData$Site,na.rm = T)
    ## first check equal number of sites
    if (numBP == length((numBP+2):(nrow(tmpData)-1))) {
        ## produce p value for wilcoxon signed rank test
        ## may not be valid as pairs are NOT independent, bp are likely correlated
        pvalue <- wilcox.test(x=tmpData$X.lnL.1[1:numBP],y=tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)], alternative = "greater",paired = T)$p.value
        myMed <- median(tmpData$X.lnL.1[1:numBP]-tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)])
        myMean <- mean(tmpData$X.lnL.1[1:numBP]-tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)])
        mySD <- sd(tmpData$X.lnL.1[1:numBP]-tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)])
        return(c(pvalue=pvalue,medDSSLS=myMed,meanDSSLS=myMean,sdDSSLS=mySD))
    
    } else return(c(pvalue=NA,medDSSLS=NA,meanDSSLS=NA,sdDSSLS=NA))
}

## find p value for each gene
## test on smaller set
#set.seed(44)
#n <- 100
#tmpFiles <- sitesFiles[sample(1:length(sitesFiles),n)] 
# tmpFiles<-tmpFiles[which(grepl("sitelikes",tmpFiles))]
start.time <- proc.time()

## include starting parallel
sfInit(parallel=TRUE, cpus=detectCores())
sfExport("compute.gene.pvalue")

pvalues <- sfSapply(sitesFiles, fun=compute.gene.pvalue)    
# pvalues <- sfSapply(tmpFiles, fun=compute.gene.pvalue)    

sfStop()
(end.time <- proc.time() - start.time)

###rename and check out distribution
# newNames <- unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(pvalues),"//"),function(tmp) tmp[2])),"\\."),function(tmp) tmp[1]))
newNames <- gsub(".+(FBgn\\d+)V.+","\\1",colnames(pvalues))
# newNames <- gsub("V$","",newNames)
colnames(pvalues) <- newNames

library(ggplot2)
###check out negative pvalues
p.data <- data.frame(t(pvalues))

write.csv(p.data,file="allGenesDSSLSstats.csv")

qplot(p.data$pvalue)

## unadjusted threshold 289 genes, 258 for AA
nrow(p.data[p.data$pvalue < 0.05,])

## adjusted pvalues (bonferroni)
adj.pvalues <- p.data
adj.pvalues$bon.pvalue <- adj.pvalues$pvalue*nrow(p.data)
adj.pvalues$bon.pvalue[adj.pvalues$bon.pvalue > 1] <- 1 
adj.pvalues$fdr.pvalue <- p.adjust(adj.pvalues$pvalue,method="BH")

for (i in c(0.05, 0.01, 0.001)) {
  print(paste("There are",nrow(adj.pvalues[adj.pvalues$bon.pvalue < i,]),"genes at bon",i))
  print(paste("There are",nrow(adj.pvalues[adj.pvalues$fdr.pvalue < i,]),"genes at fdr",i))
}

adj.pvalues <- adj.pvalues[order(adj.pvalues$pvalue),]


## reshape
adj.data <- adj.pvalues
adj.data <- cbind(row.names(adj.data),adj.data)
names(adj.data) <- c("GeneID",names(adj.data)[-1])

### write output
#write.table(adj.pvalues,file="wilcoxPvaluesNoSnothaV3.txt",sep="\t",quote=F,col.names = F)

#write.table(adj.data,file="wilcoxPvaluesNoSnothaV3.txt",sep="\t",quote=F,row.names=F)
#write.table(adj.data,file="wilcoxPvaluesNoSnothaAAV3.txt",sep="\t",quote=F,row.names=F)
##write significant
#write.table(adj.data[1:length(adj.pvalues[1,][adj.pvalues[1,] < 0.05]),],file="wilcoxSignifPvaluesNoSnothaV3.txt",sep="\t",quote=F,row.names=F)

#write.table(adj.data[1:length(adj.pvalues[1,][adj.pvalues[1,] < 0.05]),],file="wilcoxSignifPvaluesNoSnothaAAV3.txt",sep="\t",quote=F,row.names=F)

##write significant and filter
#sig.data <- adj.data[1:length(adj.pvalues[1,][adj.pvalues[1,] < 0.05]),]
#sig.data <- sig.data[sig.data$medDSSLS >= 0.1,]

#write.table(sig.data,file="wilcoxFilteredGenesNoSnothaV3.txt",sep="\t",quote=F,row.names=F)
write.csv(adj.data,file="nonsenseSSLSresults.csv",quote=F,row.names=F)
write.csv(adj.data[which(adj.data$pvalue<0.05),],file="FinalnonsenseSSLSresults.csv",quote=F,row.names=F)
write.csv(adj.data[which(adj.data$pvalue<0.05),1],file="FinalnonsenseGenelist.csv",quote=F,row.names=F)

##take a look at a couple best hits
##### Compute P value for each gene
## function testing parameters
spp<-row.names(adj.pvalues)[1]
tmpFile <- sitesFiles[which(grepl(spp,sitesFiles))]

tmpData <- read.table(tmpFile,header = T,sep = "\t",fill=T)
## compute deltaSSLS for given gene
numBP <- max(tmpData$Site,na.rm = T)
newData <- data.frame(H0=tmpData$X.lnL.1[1:numBP],H1=tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)])



p <- ggplot(data=newData,aes(x=-H0,y=-H1))
p + geom_abline(slope=1,color="red") + geom_point() + labs(x="Ln Likelihood Under the Species Tree",
                                                           y="Ln Likelihood Under the CE Tree",title=paste("Each Site of",spp ,"\nAbove the Reference Line is in Favor of CE"))

qplot(newData$H0-newData$H1)

summary(newData$H0-newData$H1)
## look at deltaSSLS


qplot(x=as.numeric(row.names(newData)),newData$H0-newData$H1) + geom_hline(yintercept=0,color="red")+ labs(x="Site (NT position)",
                                                                                                           y="Difference in Site Specific Likelihood Score (SSLS)",title=paste("For Gene ",spp,"\nAbove the Reference Line is in Favor of CE"))

## problem with very little variation genes
## produce medians

sapply(siteFiles, function(tmpFile) {
  tmpFile <- sitesFiles[grep(names(adj.pvalues)[1],sitesFiles)]
  
  tmpData <- read.table(tmpFile,header = T,sep = "\t",fill=T)
  ## compute deltaSSLS for given gene
  numBP <- max(tmpData$Site,na.rm = T)
  myMed <- median(tmpData$X.lnL.1[1:numBP]-tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)])c
  mySD <- sd(tmpData$X.lnL.1[1:numBP]-tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)])
  str(tmpData)
})



### combine with drift, tree perm results

full.data <- read.table("~/projects/flyCE/driftTreeNoSnothaV3.txt",header = T,sep="\t")
numTrees <- 150
numDrift <- 100
full.data$TreePvalue[full.data$TreePvalue == 0] <- 1/(numTrees+1)
full.data$DriftPvalue[full.data$DriftPvalue == 0] <- 1/(numDrift+1)

### combine with full gene results
str(adj.data)
p.data <- merge(adj.data,full.data[,c(1,5,6)])
write.table(p.data,file="finalGeneSSLSresults.txt",sep="\t",quote=F,row.names=F)
