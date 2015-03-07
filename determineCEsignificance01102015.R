
### Parameters
sitesDir <- "~/projects/flyCE/garliNoSnothaV3/"
#sitesDir <- "~/projects/flyCE/garliNoSnothaAAV3/"

###load packages
require(snowfall)
require(parallel)

#read in files
sitesFiles <- dir(sitesDir,full.names = T)

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
start.time <- proc.time()

## include starting parallel
sfInit(parallel=TRUE, cpus=detectCores())
sfExport("compute.gene.pvalue")

pvalues <- sfSapply(sitesFiles, fun=compute.gene.pvalue)    

sfStop()
(end.time <- proc.time() - start.time)

###rename and check out distribution
newNames <- unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(pvalues),"//"),function(tmp) tmp[2])),"\\."),function(tmp) tmp[1]))
newNames <- gsub("V$","",newNames)
colnames(pvalues) <- newNames

qplot(pvalues[1,])

## unadjusted 289 genes, 258 for AA
length(pvalues[pvalues[1,] < 0.05])

## adjusted pvalues (bonferroni)
adj.pvalues <- pvalues
adj.pvalues[1,] <- pvalues[1,]*length(pvalues[1,])
length(adj.pvalues[1,][adj.pvalues[1,] < 0.05])
length(adj.pvalues[1,][adj.pvalues[1,] < 0.01])
length(adj.pvalues[1,][adj.pvalues[1,] < 0.001])

## 93 at bonferroni 5%, 80 at 1%, 68 at .1%
## for AA, 90, 82, 69
adj.pvalues <- adj.pvalues[,order(adj.pvalues[1,])]
## cap at 1
adj.pvalues[1,][adj.pvalues[1,] > 1] <- 1
## reshape
adj.data <- data.frame(t(adj.pvalues))
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
sig.data <- adj.data[1:length(adj.pvalues[1,][adj.pvalues[1,] < 0.05]),]
sig.data <- sig.data[sig.data$medDSSLS >= 0.1,]

#write.table(sig.data,file="wilcoxFilteredGenesNoSnothaV3.txt",sep="\t",quote=F,row.names=F)
write.table(sig.data,file="wilcoxFilteredGenesNoSnothaAAV3.txt",sep="\t",quote=F,row.names=F)

##take a look at a couple best hits
##### Compute P value for each gene
## function testing parameters
tmpFile <- sitesFiles[grep(names(adj.pvalues)[1],sitesFiles)]

tmpData <- read.table(tmpFile,header = T,sep = "\t",fill=T)
## compute deltaSSLS for given gene
numBP <- max(tmpData$Site,na.rm = T)
newData <- data.frame(H0=tmpData$X.lnL.1[1:numBP],H1=tmpData$X.lnL.1[(numBP+2):(nrow(tmpData)-1)])



p <- ggplot(data=newData,aes(x=-H0,y=-H1))
p + geom_abline(slope=1,color="red") + geom_point() + labs(x="Ln Likelihood Under the Species Tree",
                                   y="Ln Likelihood Under the CE Tree",title="Each Site of FBgn0040290\nAbove the Reference Line is in Favor of CE")

qplot(newData$H0-newData$H1)

summary(newData$H0-newData$H1)
## look at deltaSSLS


qplot(x=as.numeric(row.names(newData)),newData$H0-newData$H1) + geom_hline(yintercept=0,color="red")+ labs(x="Site (NT position)",
                                   y="Difference in Site Specific Likelihood Score (SSLS)",title="For Gene FBgn0040290\nAbove the Reference Line is in Favor of CE")

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
    
    
    

