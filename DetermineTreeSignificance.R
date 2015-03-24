
sitesDir<-"C:/Projects/FlyEvolution/simTreeNonsenseV1/"
sitesFiles <- dir(sitesDir,full.names = T)
sitesFiles <- sitesFiles[which(grepl("sitelikes",sitesFiles))]

######Compute deltas for sim tree
likeFiles <- sitesFiles
#need to think about where to store simTreeDelta
outDir <- paste("C:/Projects/FlyEvolution/simTreeNonsenseV1/",sep="")
numTreesTotal <- 151 #include Ho tree
#list already computed genes
doneGene <- dir(outDir,"csv")
doneGene <- unlist(lapply(strsplit(doneGene,split="SimTree"),function(temp) temp[1]))
for (ii in likeFiles){
    gene <- strsplit(x=ii,split="/")
    gene <- strsplit(gene[[1]][length(gene[[1]])],split="simTree")[[1]][1]
 #first check if already done
    if (gene %in% doneGene){
        print(paste(gene,"is already done"))
    } else {
  if (readLines(ii)[1] == "Tree\t-lnL\tSite\t-lnL"){
      
    simTree <- read.csv(file=ii,sep="\t",header=T)
    #rename and remove variables
    simTree <- simTree[,c(1,3,4)]
    names(simTree)[3] <- "Like"
    #remove garbage rows (no likelihoods)
    simTree <-simTree[!is.na(simTree$Like),]
    simTree <-simTree[!is.na(simTree$Site),]
    #find max site
    maxSite <- max(simTree$Site)
    #determine number of approximate trees
    numTrees <- nrow(simTree) / maxSite
    if (numTrees == numTreesTotal) {
      
      ##create tree labels
      treeName <- NULL
      for (i in 0:(numTrees-1)){
        treeName <- c(treeName,paste("H",rep(as.character(i),maxSite),sep=""))
      }
      simTree$Tree <- treeName
      #null tree values
      h0 <- simTree[simTree$Tree =="H0","Like"]
      #store mean deltas and p values
      meanDelta=data.frame()
      for (j in 1:(numTrees-1)){
        tempH <- simTree[simTree$Tree == paste("H",j,sep=""),"Like"]
        tempMean <- mean(h0-tempH)
        tempPvalue <- wilcox.test(h0,tempH,alternative = "greater",paired = T)$p.value
        tempData <- data.frame(Tree=paste("H",j,sep=""),pvalue=tempPvalue,meanDelta=tempMean)
        meanDelta <- rbind(meanDelta,tempData)
      }
      write.table(x=meanDelta,file=paste(outDir,gene,"SimTreeAnalysisNonsenseV1.csv",sep=""),sep=",",row.names=F,quote=F)
    } else {
      print(paste("Not right num of trees for",gene))
      write.table(x=ii,file=paste(outDir,"errorLog.txt",sep = ""),append=T,quote=F,row.names=F)
    } 
  } else {
    write.table(x=ii,file=paste(outDir,"errorLog.txt",sep = ""),append=T,quote=F,row.names=F)
  }
}
}


###############################################################################
##Compare to observed p values
## get significant p values list
nt.data <- read.csv("C:/Projects/FlyProjectCode/FinalnonsenseSSLSresults.csv",header=T,stringsAsFactors=F)

# Find simulated means under tree sims
simDeltaDir <- dir("C:/Projects/FlyEvolution/simTreeNonsenseV1/",pattern="SimTreeAnalysisNonsenseV1",full.names=T)
for (ii in nt.data$GeneID){
  simFile <-  simDeltaDir[grepl(pattern=as.character(ii),simDeltaDir)]
  if (length(simFile)>0){
    simDelta <- read.csv(simFile)
    observed <- nt.data[nt.data$GeneID == ii,"pvalue"]
    sims <- simDelta[,"pvalue"]
    pvalue <- length(sims[sims <= observed])/length(sims)
    nt.data[nt.data$GeneID == ii,"TreePvalue"] <- pvalue
  }
}
write.csv(x=nt.data[order(nt.data$TreePvalue),],file="C:/Projects/FlyEvolution/simTreeNonsenseV1/treeSimPvalueNonsenseV1.csv",col.names=T,row.names=F,quote=F)

### also compute AA
# aa.tree.dir <- paste(homePath,"simTreeOutputNoSnothaAAV3/",sep="")
# treeDir <- aa.tree.dir
######Compute deltas for sim tree
# likeFiles <- dir(treeDir,pattern="site",full.names = T)
#need to think about where to store simTreeDelta
# outDir <- paste(mainDir,"simTreeAnalysisNoSnothaAAV3/",sep="")
# numTreesTotal <- 151 #include Ho tree
#list already computed genes
# doneGene <- dir(outDir,"csv")
# doneGene <- unlist(lapply(strsplit(doneGene,split="SimTree"),function(temp) temp[1]))
# for (ii in likeFiles){
#     gene <- strsplit(x=ii,split="/")
#     gene <- strsplit(gene[[1]][length(gene[[1]])],split="simTree")[[1]][1]
#  #first check if already done
#     if (gene %in% doneGene){
#         print(paste(gene,"is already done"))
#     } else {
#   if (readLines(ii)[1] == "Tree\t-lnL\tSite\t-lnL"){
#       
#     simTree <- read.csv(file=ii,sep="\t",header=T)
#     #rename and remove variables
#     simTree <- simTree[,c(1,3,4)]
#     names(simTree)[3] <- "Like"
#     #remove garbage rows (no likelihoods)
#     simTree <-simTree[!is.na(simTree$Like),]
#     simTree <-simTree[!is.na(simTree$Site),]
#     #find max site
#     maxSite <- max(simTree$Site)
#     #determine number of approximate trees
#     numTrees <- nrow(simTree) / maxSite
#     if (numTrees == numTreesTotal) {
#       
#       ##create tree labels
#       treeName <- NULL
#       for (i in 0:(numTrees-1)){
#         treeName <- c(treeName,paste("H",rep(as.character(i),maxSite),sep=""))
#       }
#       simTree$Tree <- treeName
#       #null tree values
#       h0 <- simTree[simTree$Tree =="H0","Like"]
#       #store mean deltas and p values
#       meanDelta=data.frame()
#       for (j in 1:(numTrees-1)){
#         tempH <- simTree[simTree$Tree == paste("H",j,sep=""),"Like"]
#         tempMean <- mean(h0-tempH)
#         tempPvalue <- wilcox.test(h0,tempH,alternative = "greater",paired = T)$p.value
#         tempData <- data.frame(Tree=paste("H",j,sep=""),pvalue=tempPvalue,meanDelta=tempMean)
#         meanDelta <- rbind(meanDelta,tempData)
#       }
#       write.table(x=meanDelta,file=paste(outDir,gene,"SimTreeAnalysisNoSnothaAAV3.csv",sep=""),sep=",",row.names=F,quote=F)
#     } else {
#       print(paste("Not right num of trees for",gene))
#       write.table(x=ii,file=paste(outDir,"errorLog.txt",sep = ""),append=T,quote=F,row.names=F)
#     } 
#   } else {
#     write.table(x=ii,file=paste(outDir,"errorLog.txt",sep = ""),append=T,quote=F,row.names=F)
#   }
# }
# }

## now for AA
# aa.data <- read.table("/xdisk/rlapoint/rScripts/wilcoxSignifPvaluesNoSnothaAAV3.txt",header=T,sep="\t",stringsAsFactors=F)

#Find simulated means under tree sims
# simDeltaDir <- "/xdisk/rlapoint/simTreeAnalysisNoSnothaAAV3/"
# for (ii in aa.data$GeneID){
#   simFile <-  dir(simDeltaDir,pattern=as.character(ii),full.names = T)
#   if (length(simFile)>0){
#     simDelta <- read.csv(simFile)
#     observed <- aa.data[aa.data$GeneID == ii,"pvalue"]
#     sims <- simDelta[,"pvalue"]
#     pvalue <- length(sims[sims <= observed])/length(sims)
#     aa.data[aa.data$GeneID == ii,"TreePvalue"] <- pvalue
#   }
# }
# write.table(x=aa.data[order(aa.data$TreePvalue),],file=paste(simDeltaDir,"treeSimPvalueNoSnothaAAV3.txt",sep = ""),col.names=T,row.names=F,sep="\t",quote=F)

###pick up here: 11 Nov 2014

###############################################################################################
##2. Drift Simulation 
#Pull only candidate genes that passed tree sims
candDelta <- read.table(paste(sourceDir,"treeSimPvalueNoSnothaV3.txt",sep=""), header=T,sep="\t")
candDelta <- candDelta[order(candDelta$pvalue),]
cand4Drift <- candDelta[candDelta$TreePvalue < 0.05,]

#create fasta files from ali files for each cand gene
driftDir <- "/xdisk/rlapoint/SimsV3/SimAlignments/"
fasDir <- paste(mainDir,"simDriftV3/",sep = "")
candGenes <- as.character(cand4Drift[,"GeneID"])

###Function covert ali to fasta
#######Convert ALI to FASTA for input to GARLI
aliToFasta <- function(gene,aliFiles,outDir,idType = c("bird","fly")){
  for (ii in aliFiles){
      if (idType == "bird") {
          sampleNum <- strsplit(strsplit(ii,split="_")[[1]][4],"\\.")[[1]][1]
      } else if (idType == "fly"){
          sampleNum <- strsplit(strsplit(ii,split="_")[[1]][4],"\\.")[[1]][1]
      }
    tempText <- readLines(ii)
    #remove first line
    tempText <- tempText[2:length(tempText)]
    #reduce/replace whitespace
    tempText <- gsub(pattern="     ","\n",tempText)
    #split on newline
    tempText <- strsplit(tempText,split="\n")
      # retain the short name
    #tempText <- lapply(tempText,function(temp){
    #  name <- as.character(nameData$common[which(temp[1] == substring(text=nameData$common,first=1,last=4))])
    #  temp[1] <- name
    #  temp
    #})
    #collapse list
    tempText <- unlist(tempText)
    #species are odd indeces
    oddIndex <- seq(from=1,to=length(tempText),by=2)
    #browser()
    tempText[oddIndex] <- paste(">",tempText[oddIndex],sep="")
    #browser()
    writeLines(text=tempText,con=paste(outDir,gene,"_sim",sampleNum,".fasta",sep=""),sep="\n")
  }
}


for (zz in candGenes){
    #create directory for alignments
    outDir <- paste(fasDir,zz,"/",sep = "")
    suppressWarnings(dir.create(paste(fasDir,zz,sep = "")))
    aliFiles <- dir(driftDir,zz,full.names = T)
    aliFiles <- aliFiles[grep("ali",aliFiles)]
    #sample 100 of the 1000 for now - revisit if needed
    set.seed(44)
    if (length(aliFiles) >= 100) {
        print(paste("Converting ali to fasta for", zz))
        aliFiles <- sample(aliFiles,100)
        aliToFasta(zz,aliFiles,outDir,idType="fly")
    } else {
        print(paste("Not enough alignments for",zz))
        write.table(x=zz,file=paste(fasDir,"SimDriftErrorLog.txt",sep = ""),append=T,quote=F,row.names=F,col.names = F)
    }
}

#Sim Drift Error log
if (file.exists(paste(fasDir,"SimDriftErrorLog.txt",sep=""))) {
    redo <- readLines(con=paste(fasDir,"SimDriftErrorLog.txt",sep=""))
} else redo <- NULL

##########Write Garli Conf for Sim Drift

template <- suppressWarnings(readLines(con="/xdisk/rlapoint/garliConfFiles/garliTemplateNoSnothaV3.conf"))
garliDir <- "garliOutputNoSnothaV3/"

for (ii in candGenes){
    print(paste("working on sim conf",ii))
    inDir <- paste(fasDir,ii,"/",sep = "")
    tmp.garliDir <- paste(fasDir,ii,"/",garliDir,sep = "")
    suppressWarnings(dir.create(tmp.garliDir))
    fasFiles <- dir(inDir)
    for (jj in fasFiles){
        conf <- template
        conf[2] <- paste(conf[2],inDir,jj,sep="")
        conf[6] <- paste(conf[6],tmp.garliDir,jj,sep="")
writeLines(text=conf,con=paste(mainDir,"simDriftV3/garliConfNoSnothaV3/",jj,"Garli.conf",sep=""),sep="\n")
    }
}

## create pbs scripts and submit
dir.create(paste(mainDir,"simWorkingDirNoSnothaV3/",sep=""))
dir.create(paste(mainDir,"simWorkingDirNoSnothaV3/pbsScripts",sep=""))
setwd(paste(mainDir,"simWorkingDirNoSnothaV3/pbsScripts",sep=""))
template <- readLines(paste(mainDir,"pbsScripts/simTemplate.pbs",sep=""))

confFiles <- dir(paste(mainDir,"simDriftV3/garliConfNoSnothaV3/",sep = ""),full.names = T)
count <- 1
for (ff in confFiles){
    shortName <- strsplit(strsplit(ff,"/")[[1]][7],"\\.")[[1]][1]
    conf <- template
    conf[2] <- paste(conf[2],count,sep="")
    conf[16] <- paste("~/bin/Garli-2.01 ", ff,sep="")
    writeLines(conf,paste("simDrift",shortName,".pbs",sep=""))
    count <- count + 1
}

#submit
# note that the name was cut, make shorter next time
setwd(paste(mainDir,"simWorkingDirNoSnothaV3/",sep=""))
pbsFiles <- dir("./pbsScripts",full.names=T)
for (pp in pbsFiles) {
    system(paste("qsub",pp))
}

###############################################################
###############################################################
###############################################################
### OLD CODE

############compute Delta SSLS from Garli Output

computeGarliDelta <- function(sitesPath,outDir,multTree=F){
  #browser()
  #geneID
  gene <- strsplit(x=sitesPath,split="/")
  gene <- strsplit(gene[[1]][length(gene[[1]])],split="\\.")[[1]][1]
  sites <- read.csv(sitesPath,header=T,sep="\t")
  num <- nrow(sites)
  h0 <- sites[1:(num/2-1),4]
  h1 <- sites[(num/2+1):(num-1),4]
  pvalue <- wilcox.test(h0,h1,"greater",paired=T)$p.value
  write(x=paste(gene,pvalue,sep="\t"),append=T,file=paste(outDir,"pvaluesWilcoxon.txt",sep=""))
}

#test one gene
#sitesPath <- "/xdisk/rlapoint/simDriftV3/FBgn0034691/garliOutputNoSnothaV3/FBgn0034691_sim339.fasta.sitelikes.log"

#outDir <- "/xdisk/rlapoint/simDriftV3/FBgn0036764/pvalueV3/"
#dir.create(outDir)
#computeGarliDelta(sitesPath,outDir)

computeGarliDir <- function(inDir,outDir){
  #browser()
  gFiles <- dir(paste(inDir,"/",sep=""),"*.sitelikes.log$")
  for (gg in gFiles){
      #first check if already done
      if (!(unlist(lapply(strsplit(gg,"\\."), function(temp) temp[1])) %in% unlist(lapply(strsplit(dir(outDir),"delta"), function(temp) temp[1])))) {
          computeGarliDelta(paste(inDir,"/",gg,sep=""),outDir)
          print(paste("Computing pvalues for SSLS for",gg))
      } else print(paste("Pvalues already computed for",gg))
  }
}

#############Run this after Garli finishes
driftDir <- "/xdisk/rlapoint/simDriftV3/"
for (ii in candGenes) {
    inDir <- paste(driftDir,ii,"/garliOutputNoSnothaV3/",sep = "")
    outDir <- paste(driftDir,ii,"/pvalueV3/",sep = "")
    suppressWarnings(dir.create(outDir))
    if (file.exists(paste(outDir,"pvaluesWilcoxon.txt",sep=""))) file.remove(paste(outDir,"pvaluesWilcoxon.txt",sep=""))
    computeGarliDir(inDir,outDir)
}

###########################################################################
##Compare to observed pvalue
candDelta <- read.table(paste(sourceDir,"treeSimPvalueNoSnothaV3.txt",sep=""), header=T,sep="\t")
candDelta <- candDelta[order(candDelta$pvalue),]
#cand4Drift <- candDelta[candDelta$TreePvalue < 0.05,]

#read candidate list
#candDelta2 <- cand4Drift
#candDelta <- cand4Drift
#Find simulated means under drift
candDelta$DriftPvalue <- NA
#not which need to be redo
#candDelta$Notes <- NA
#candDelta$Notes[candDelta$geneID %in% redo] <- "NeedToRedoDriftSims"

for (ii in candGenes){
    file <- paste(driftDir,ii,"/pvalueV3/pvaluesWilcoxon.txt",sep="")
    tmp <- read.table(file)
    observed <- candDelta[candDelta$GeneID == ii,2]
    sims <- tmp[,2]
    pvalue <- length(sims[sims < observed])/length(sims)
    candDelta[candDelta$GeneID == ii,"DriftPvalue"] <- pvalue

}
sortedDelta <- candDelta[with(candDelta, order(DriftPvalue,TreePvalue,pvalue,-medDSSLS)),]
write.table(x=sortedDelta,file="/xdisk/rlapoint/rScripts/driftTreeNoSnothaV3.txt",col.names=T,row.names=F,sep="\t",quote=F)

