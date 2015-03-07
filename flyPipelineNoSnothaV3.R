############Convergent Evolution Analysis Pipeline Wrapper

fastaFile <- "~/mount/rlapoint/FastaFilesForAlignmentNoSnothaV3/FBgn0000008V3.fasta"
readLines(fastaFile,n=16)
########################################################
##Parameters will come from control file in the future
nameData <- data.frame(common=c("Dmelanogaster","Dyakuba","Dbiarmipes","Dsuzukii","Dpseudoobscura","Dgrimshawi","Sflava","Dmojavensis"),short=c("Dmel", "Dyak", "Dbia", "Dsuz", "Dpse", "Dgri", "Sfla","Dmoj"),stringsAsFactors = F)
#reorder to match FASTA file
nameData <- nameData[c(1,3,6,8,5,4,2,7),]

######Specify source code, input, output directory
mainDir <- "/xdisk/rlapoint/"
#mainDir <- "~/mount/rlapoint/"
outputDir <- paste(mainDir,"pipelineOutputNoSnothaV3/",sep = "")
fastaDir <- paste(mainDir,"FastaFilesForAlignmentNoSnothaV3/",sep = "")
sourceDir <- paste(mainDir,"rScripts/",sep = "")

###Find all gene ids
geneKey <-unlist(lapply(strsplit(dir(fastaDir,"fasta"),"\\."), function(temp) temp[1]))
### there are 11,285 genes

formatHMMERforMSA <- function(dir,outDir) {
 #Input is HMMER reciprocal hits FASTA dir
 #Output is properly formated FASTA as input to ClustalW2
    files <- dir(dir,"*.fasta")
    for (file in files) {
        geneID <- strsplit(file,"\\.")[[1]][1]
        tempText <- readLines(con=paste(dir,file,sep = "/"))
        #insert newID
        for (i in 1:nrow(nameData)){
            tempLine <- grep(nameData$short[i],tempText)
            tempText[tempLine] <- paste(">",nameData$common[i],sep = "") 
        }
        writeLines(text=tempText,con=paste(outDir,geneID,"codingDNA.fas",sep=""),sep="\n")
        print(paste(geneID,"has been formatted for MSA"))
    }
}

outDir <- paste(mainDir,"FastaFilesPreMuscleNoSnothaV3/",sep = "")
formatHMMERforMSA(dir = fastaDir,outDir=outDir)

#####################################################################################################
###make commandFile MSA
### Need to split since too many genes
halfIndex <- floor(length(geneKey) / 2)
commandText <- NULL
for (ii in geneKey[1:halfIndex]){
  tempText <- paste("muscle -in /xdisk/rlapoint/FastaFilesPreMuscleNoSnothaV3/",ii,"codingDNA.fas -out /xdisk/rlapoint/FastaFilesPreMuscleNoSnothaV3/",ii,"msa.afa",sep="")
  commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileMuscleNoSnothaV3part1.txt",sep = ""),sep="\n")

commandText <- NULL
for (ii in geneKey[(halfIndex+1):length(geneKey)]){
  tempText <- paste("muscle -in /xdisk/rlapoint/FastaFilesPreMuscleNoSnothaV3/",ii,"codingDNA.fas -out /xdisk/rlapoint/FastaFilesPreMuscleNoSnothaV3/",ii,"msa.afa",sep="")
  commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileMuscleNoSnothaV3part2.txt",sep = ""),sep="\n")

#Write and submit pbs command files from pipeline output dir
#qsub ../pbsScripts/msaMuscleNoSnothaV3part1.pbs
#qsub ../pbsScripts/msaMuscleNoSnothaV3part2.pbs

#musPart1 446644[].service2
#musPart2 446645[].service2


#done msaMuscleV1.pbs

######################
###make commandFile MSA refinement
### Need to split since too many genes
###Find all gene ids
geneKey <-unlist(lapply(strsplit(dir(paste(mainDir,"FastaFilesPreMuscleNoSnothaV3/",sep = ""),"afa"),"\\."), function(temp) temp[1]))
geneKey <- gsub("msa","",geneKey)
### there are 11214 genes (lost 71)

##get new gene count!!
halfIndex <- floor(length(geneKey) / 2)
commandText <- NULL
for (ii in geneKey[1:halfIndex]){
    tempText <- paste("muscle -in /xdisk/rlapoint/FastaFilesPreMuscleNoSnothaV3/",ii,"msa.afa -out /xdisk/rlapoint/FastaFilesPostMuscleNoSnothaV3/",ii,"refinedMSA.afa -refine",sep="")
    commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileMuscleRefineNoSnothaV3part1.txt",sep = ""),sep="\n")

commandText <- NULL
for (ii in geneKey[(halfIndex+1):length(geneKey)]){
    tempText <- paste("muscle -in /xdisk/rlapoint/FastaFilesPreMuscleNoSnothaV3/",ii,"msa.afa -out /xdisk/rlapoint/FastaFilesPostMuscleNoSnothaV3/",ii,"refinedMSA.afa -refine",sep="")
    commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileMuscleRefineNoSnothaV3part2.txt",sep = ""),sep="\n")

#Write and submit pbs command files
#qsub ../pbsScripts/msaMuscleRefineNoSnothaV3part1.pbs
#qsub ../pbsScripts/msaMuscleRefineNoSnothaV3part2.pbs

#qstat -t 446707[].service2 #msaMuscleRefineNoSnothaV3part1.pbs started @ 4:07pm 9/18/14
#qstat -t 446708[].service2 #msaMuscleRefineNoSnothaV3part1.pbs started @ 4:07pm 9/18/14


###########################################################
###########################################################
## Convert pir to fasta
pirDirToFASTA <- function(inDir,outDir){
  ###input is clustal codons in PIR format
  #browser()
  pirFiles <- dir(inDir,"pir")
  for (ii in pirFiles){
    tempText <- readLines(con=file.path(inDir,ii))
    nameIndex <- grep(">",tempText)
    tempText[nameIndex] <- unlist(lapply(strsplit(tempText[nameIndex],"_"), function(s) s[length(s)]))
    tempText <- gsub("DL;","",tempText)
    tempText <- gsub("\\*","",tempText)
    tempText <- tempText[!nchar(tempText) == 0]
    outfile <- paste(strsplit(ii,split="\\.")[[1]][1],".fas",sep="")
    print(paste("Writing file",outfile))
    writeLines(tempText,file.path(outDir,outfile),sep="\n")
  }
}

inDir <- paste(mainDir,"FastaFilesPreClustalV1/",sep = "")
outDir <- paste(mainDir,"msaFASTAfilesV1/",sep = "")
pirDirToFASTA(inDir,outDir)



#############################################
##make commandFile Gblocks
#setwd("C:\\Users\\Grant\\Dropbox\\EcolProject\\EcolProjectData\\programsForHPC")
commandText <- NULL
for (ii in geneKey){
  tempText <- paste("Gblocks /xdisk/rlapoint/msaFASTAfilesV1/",ii,"msa.fas"," -t=c -p=n",sep="")
  commandText <- c(commandText,tempText)
}
writeLines(text=commandText,con="/xdisk/rlapoint/cmdFilesLazyArray/cmdFileGBlocksV1.txt",sep="\n")

#submit pbs script, run under subdirectory to capture output.
#raise walltime
#completed 3/11/14 6018 genes left


###### 22 Oct 2014
###### pick here since TranslatorX  also ran Gblocks

#### skip this step to keep everything in frame!

###########Remove STOP codons
#load from Linux lib and load NameData
require(stringr)
inDir <- paste(mainDir,"msaFASTAfilesV1/",sep = "")
outDir <- paste(mainDir,"removeStopFASTAfilesV1/",sep = "")

removeStop <- function(file,outfile,paml=T){
  #####################################################
  ##Function to Remove Stop
  ##Input is FASTA Gblocks cleaned DNA codons
  #No Output but PAML format is written to file for input to SLR
  #browser()
  foundStop <- F
  stopCodon <- c("TAA","TAG","TGA")
  tempText <- readChar(con=file,file.info(file)$size)
  #remove \r tags
  tempText <- gsub("\\\r","",tempText)
  #pull species name
  textSpecies <- unlist(strsplit(x=tempText,split="\n>",))
  textSpecies <- gsub(pattern="[>]",replacement="",x=textSpecies)
  textSpecies <- gsub(pattern=" ",replacement="",x=textSpecies)
  listSpecies <- strsplit(textSpecies,"\n")
  listStop <- lapply(listSpecies,function(temp){
    #paste together coding DNA
    tempSeq <- paste(temp[2:length(temp)],collapse="")
    #check each codon for a STOP, Replace with XXX
    for (i in seq(from=1,to=nchar(tempSeq),by=3)){
      if (substring(text=tempSeq,first=i,last=i+2) %in% stopCodon){
        #browser()
        tempSeq <- paste(substring(text=tempSeq,first=1,last=i-1),"XXX",substring(text=tempSeq,first=i+3,last=nchar(tempSeq)),sep="")
      }
    }
    names(tempSeq) <- temp[1]
    tempSeq
  })
  #Actually remove STOP
  charStop <- unlist(listStop)
  #all the indeces
  stopIndex <- unique(unlist(lapply(str_locate_all(string=charStop,pattern="XXX"),function(temp){
    temp[,1]
  })))
  if (length(stopIndex)>0) {foundStop <- T}
  finalSeq <- lapply(listStop,function(temp){
    if (length(stopIndex)>0){
      for (i in 1:length(stopIndex)){
        tempIndex <- stopIndex[i]
        temp <- paste(substring(text=temp,first=1,last=tempIndex-1),"XXX",substring(text=temp,first=tempIndex+3,last=nchar(temp)),sep="")
      }
    }
    temp
  })
  names(finalSeq) <- names(charStop)
  #remove all XXX tags
  finalSeq <- gsub(pattern="XXX","",x=finalSeq)
  #PAML format for SLR
  if (paml) {
    num <- nchar(finalSeq[[1]][1])
    numSpecies <- length(finalSeq)
    firstLine <- paste(numSpecies," ",num,sep="")
    toWrite <- paste(names(charStop),unlist(finalSeq),sep="\n")
    toWrite <- c(firstLine,toWrite)
  } else {
    toWrite <- paste(">",names(charStop),"\n",unlist(finalSeq),sep="")  
  }
  #browser()
  strsplit(file,split="/")[[1]] -> newFile
  #path <- strsplit(file,split="/")[[1]]
  newFile <- paste(strsplit(newFile[length(newFile)],split="\\.")[[1]][1],"_Cleaned.fas",sep="")
  #path <- strsplit(file,split="/")[[1]]
  #path <- paste(path[2:(length(path)-1)],collapse="/")
  #browser()
  if (foundStop) print(paste("Stops Removed for",newFile)) else print(paste("No Stops Found for",newFile))
  writeLines(toWrite,file.path(outDir,newFile),sep="\n")
}

removeDirStop <- function(inDir,outDir,inPattern="*.fas-gb",paml=F){
  for (jj in dir(path=inDir,pattern=inPattern)){
    tempText <- readChar(con=file.path(inDir,jj),file.info(file.path(inDir,jj))$size)
    tempText <- gsub("\\\r|\\\n|>","",tempText)
    tempText <- gsub(paste(as.character(nameData$common),collapse="|"),"",tempText)
    if (tempText != ""){
      print(removeStop(file.path(inDir,jj),file.path(outDir,jj),paml=F))
    }
  }
}


removeDirStop(inDir,outDir,paml = F)

#########################################################


###get only genes with 8 species
#file.info(dir("/xdisk/rlapoint/FastaFilesForAlignmentNoSnothaV3", full.names=T)[1])
speciesKey <- readLines(dir("/xdisk/rlapoint/FastaFilesForAlignmentNoSnothaV3", full.names=T)[1])
#any(unlist(strsplit(speciesKey[11274],","))=="")
logicKey <- sapply(speciesKey,function(tmp) {
    char <- unlist(strsplit(tmp,","))
    !any(char=="")
})
geneKey <- unlist(lapply(strsplit(names(logicKey)[logicKey],","), function(tmp) tmp[1]))

#########################################################################
### GARLI on NT alignments
####find MSAs with full 8 species
inDir <- "/xdisk/rlapoint/CleanAlignmentsV3/"

####################Write Garli Conf for 2 tree 2 search 6 rate NT model
###write Garli conf files
#read in original Garli conf without geneID
garliDir <- "garliNoSnothaV3/"
inDir <- "CleanAlignmentsV3/"
mainDir <- "/xdisk/rlapoint/"
outDir <- paste(mainDir,"garliConfNoSnothaV3/",sep="")
#doneKey <- unlist(lapply(strsplit((dir(outDir)),"garli"),function(tmp) tmp[1]))
#toDoKey <- geneKey[!(geneKey %in% doneKey)]

fasFiles <- dir(paste(mainDir,inDir,sep=""))
for (ii in fasFiles) {
    print(ii)
    tmpText <- readLines(paste(mainDir,inDir,ii,sep=""))
    tmpText <- gsub(">[A-z0-9]{1,50},",">",tmpText)
    writeLines(tmpText,paste(mainDir,"CleanAlignmentsV3shortName/",ii,sep=""))
}


####first get standard names in fasta files to make
inDir <- "CleanAlignmentsV3shortName/"
for (ii in geneKey){
  print(ii)
  #remove 3 from V3 for some reason
  ii <- gsub("V3","V",ii)
  conf <- readLines(con="/xdisk/rlapoint/garliConfFiles/garliTemplateNoSnothaV3.conf")
  conf[2] <- paste(conf[2],mainDir,inDir,ii,".nt_cleanali.fasta",sep="")
  conf[6] <- paste(conf[6],mainDir,garliDir,ii,sep="")
  writeLines(text=conf,con=paste(mainDir,"garliConfNoSnothaV3/",ii,"garli.conf",sep=""),sep="\n")
}



#################################Write Command File for Garli
#stopFiles <- dir(paste(mainDir,"removeStopFASTAfilesV1/",sep = ""))
#geneTotal <- unlist(lapply(strsplit(stopFiles,split="msa_Cleaned\\.fas"),function(temp){temp[1]}))

halfIndex <- floor(length(geneKey) / 2)
commandText <- NULL
for (ii in geneKey[1:halfIndex]){
    ii <- gsub("V3","V",ii)
    tempText <- paste("Garli-2.01 /xdisk/rlapoint/garliConfNoSnothaV3/",ii,"garli.conf",sep="")
    print(ii)
    commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileGarliNoSnothaV3part1.txt",sep = ""),sep="\n")

commandText <- NULL
for (ii in geneKey[(halfIndex+1):length(geneKey)]){
    ii <- gsub("V3","V",ii)
    tempText <- paste("Garli-2.01 /xdisk/rlapoint/garliConfNoSnothaV3/",ii,"garli.conf",sep="")
    print(ii)
    commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileGarliNoSnothaV3part2.txt",sep = ""),sep="\n")


#Write and submit pbs command files
#qsub ../pbsScripts/garliNoSnothaV3part1.pbs
#qsub ../pbsScripts/garliNoSnothaV3part2.pbs

#qstat -t 183650[].service1
#garli1 started @ 4:45pm 10/22/14
#qstat -t 183651[].service1
#garli2 @ 4:45pm 10/22/14

#########################################################################
### GARLI on AA alignments
####find MSAs with full 8 species
#inDir <- "/xdisk/rlapoint/TranslatorCleanAA/"

####################Write Garli Conf for 2 tree 2 search 6 rate NT model
###write Garli conf files
#read in original Garli conf without geneID
garliDir <- "garliNoSnothaAAV3/"
inDir <- "TranslatorCleanAA/"
mainDir <- "/xdisk/rlapoint/"
outDir <- paste(mainDir,"garliConfNoSnothaAAV3/",sep="")
#doneKey <- unlist(lapply(strsplit((dir(outDir)),"garli"),function(tmp) tmp[1]))
#toDoKey <- geneKey[!(geneKey %in% doneKey)]

fasFiles <- dir(paste(mainDir,inDir,sep=""))
for (ii in fasFiles) {
    print(ii)
    tmpText <- readLines(paste(mainDir,inDir,ii,sep=""))
    tmpText <- gsub(">[A-z0-9]{1,50},",">",tmpText)
    writeLines(tmpText,paste(mainDir,"CleanAlignmentsV3shortName/",ii,sep=""))
}


####first get standard names in fasta files to make
inDir <- "CleanAlignmentsV3shortName/"
### forgot to change the name to CleanAlignmentsAAV3shortName/ # 24 Oct 2014 1:27pm

conf <- readLines(con="/xdisk/rlapoint/garliConfFiles/garliTemplateNoSnothaAAV3.conf")
for (ii in geneKey){
  print(ii)
  #remove 3 from V3 for some reason
  ii <- gsub("V3","V",ii)
  tempConf <- conf
  tempConf[2] <- paste(tempConf[2],mainDir,inDir,ii,".aa_cleanali.fasta",sep="")
  tempConf[6] <- paste(tempConf[6],mainDir,garliDir,ii,sep="")
  writeLines(text=tempConf,con=paste(mainDir,"garliConfNoSnothaAAV3/",ii,"garli.conf",sep=""),sep="\n")
}



#################################Write Command File for Garli AA
#stopFiles <- dir(paste(mainDir,"removeStopFASTAfilesV1/",sep = ""))
#geneTotal <- unlist(lapply(strsplit(stopFiles,split="msa_Cleaned\\.fas"),function(temp){temp[1]}))

halfIndex <- floor(length(geneKey) / 2)
commandText <- NULL
for (ii in geneKey[1:halfIndex]){
    ii <- gsub("V3","V",ii)
    tempText <- paste("Garli-2.01 /xdisk/rlapoint/garliConfNoSnothaAAV3/",ii,"garli.conf",sep="")
    print(ii)
    commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileGarliNoSnothaAAV3part1.txt",sep = ""),sep="\n")

commandText <- NULL
for (ii in geneKey[(halfIndex+1):length(geneKey)]){
    ii <- gsub("V3","V",ii)
    tempText <- paste("Garli-2.01 /xdisk/rlapoint/garliConfNoSnothaAAV3/",ii,"garli.conf",sep="")
    print(ii)
    commandText <- c(commandText,tempText)
}

outDir <- paste(mainDir,"cmdFilesLazyArray/",sep = "")
writeLines(text=commandText,con=paste(outDir,"cmdFileGarliNoSnothaAAV3part2.txt",sep = ""),sep="\n")


#Write and submit pbs command files
qsub ../pbsScripts/garliNoSnothaAAV3part1.pbs
qsub ../pbsScripts/garliNoSnothaAAV3part2.pbs

#qstat -t 490449[].service2
#garliAA1 started @ 4:45pm 10/22/14
#qstat -t 490450[].service2
#garliAA2 @ 4:45pm 10/22/14


###transfer to local computer and run screenGenes.R

###Access Alternative hypotheses
#################################################################################
################################################################



###Tree Sims
homePath <- "/xdisk/rlapoint/"
###find and move significant cleaned alignments
inDir <- paste(homePath,"CleanAlignmentsV3shortName/",sep="")
## get significant p values list
nt.data <- read.table("/xdisk/rlapoint/rScripts/wilcoxSignifPvaluesNoSnothaV3.txt",header=T,sep="\t",stringsAsFactors=F)
nt.files <- paste(homePath,"CleanAlignmentsV3shortName/",nt.data$GeneID,"V.nt_cleanali.fasta",sep="")
## now for AA
aa.data <- read.table("/xdisk/rlapoint/rScripts/wilcoxSignifPvaluesNoSnothaAAV3.txt",header=T,sep="\t",stringsAsFactors=F)
aa.files <- paste(homePath,"CleanAlignmentsV3shortName/",aa.data$GeneID,"V.aa_cleanali.fasta",sep="")

nt.tree.dir <- paste(homePath,"simTreeOutputNoSnothaV3/",sep="")
## grab template
conf <- readLines(con=paste(homePath,"garliConfFiles/garliTree.conf",sep=""))
###Write NT Conf file
for (ii in nt.files){
    tmpGeneID <- unlist(strsplit(unlist(lapply(strsplit(ii,"V\\."), function(temp) temp[1])),"/"))[5]
    tmp.conf <- conf
    tmp.conf[2] <- paste(tmp.conf[2],ii,sep="")
    tmp.conf[6] <- paste(tmp.conf[6],nt.tree.dir,tmpGeneID,"simTree",sep="")
    writeLines(text=tmp.conf,con=paste(homePath,"simTreeConfNoSnothaV3/",tmpGeneID,"simTree.conf",sep=""),sep="\n")
}

## grab template
conf <- readLines(con=paste(homePath,"garliConfFiles/garliTreeAA.conf",sep=""))
aa.tree.dir <- paste(homePath,"simTreeOutputNoSnothaAAV3/",sep="")
###Write aa Conf file
for (ii in aa.files){
    tmpGeneID <- unlist(strsplit(unlist(lapply(strsplit(ii,"V\\."), function(temp) temp[1])),"/"))[5]
    tmp.conf <- conf
    tmp.conf[2] <- paste(tmp.conf[2],ii,sep="")
    tmp.conf[6] <- paste(tmp.conf[6],aa.tree.dir,tmpGeneID,"simTree",sep="")
    writeLines(text=tmp.conf,con=paste(homePath,"simTreeConfNoSnothaAAV3/",tmpGeneID,"simTree.conf",sep=""),sep="\n")
}


##Write Command File for NT Sim Trees
treeConfDir <- paste(homePath,"simTreeConfNoSnothaV3/",sep="")
confFiles <- dir(treeConfDir,pattern=".conf")
toWrite <- NULL
for (jj in confFiles){
  tempCommand <- paste("Garli-2.01 ",treeConfDir,jj,sep="")
  toWrite <- c(toWrite,tempCommand)
}
writeLines(text=toWrite,con=paste(homePath,"cmdFilesLazyArray/cmdFileSimTreesNoSnothaV3.txt",sep = ""),sep="\n")

## Next copy
##Write Command File for AA Sim Trees
treeConfDir <- paste(homePath,"simTreeConfNoSnothaAAV3/",sep="")
confFiles <- dir(treeConfDir,pattern=".conf")
toWrite <- NULL
for (jj in confFiles){
  tempCommand <- paste("Garli-2.01 ",treeConfDir,jj,sep="")
  toWrite <- c(toWrite,tempCommand)
}
writeLines(text=toWrite,con=paste(homePath,"cmdFilesLazyArray/cmdFileSimTreesNoSnothaAAV3.txt",sep = ""),sep="\n")

#write and submit pbs
#/xdisk/rlapoint/pbsScripts/simTreeV3nt.pbs
#/xdisk/rlapoint/pbsScripts/simTreeV3aa.pbs
# qstat -t 203921[].service1
# qstat -t 204046[].service1

mainDir <- homePath
nt.tree.dir <- paste(homePath,"simTreeOutputNoSnothaV3/",sep="")
treeDir <- nt.tree.dir
######Compute deltas for sim tree
likeFiles <- dir(treeDir,pattern="site",full.names = T)
#need to think about where to store simTreeDelta
outDir <- paste(mainDir,"simTreeAnalysisNoSnothaV3/",sep="")
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
      write.table(x=meanDelta,file=paste(outDir,gene,"SimTreeAnalysisNoSnothaV3.csv",sep=""),sep=",",row.names=F,quote=F)
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
nt.data <- read.table("/xdisk/rlapoint/rScripts/wilcoxSignifPvaluesNoSnothaV3.txt",header=T,sep="\t",stringsAsFactors=F)

#Find simulated means under tree sims
simDeltaDir <- "/xdisk/rlapoint/simTreeAnalysisNoSnothaV3/"
for (ii in nt.data$GeneID){
  simFile <-  dir(simDeltaDir,pattern=as.character(ii),full.names = T)
  if (length(simFile)>0){
    simDelta <- read.csv(simFile)
    observed <- nt.data[nt.data$GeneID == ii,"pvalue"]
    sims <- simDelta[,"pvalue"]
    pvalue <- length(sims[sims <= observed])/length(sims)
    nt.data[nt.data$GeneID == ii,"TreePvalue"] <- pvalue
  }
}
write.table(x=nt.data[order(nt.data$TreePvalue),],file=paste(simDeltaDir,"treeSimPvalueNoSnothaV3.txt",sep = ""),col.names=T,row.names=F,sep="\t",quote=F)

### also compute AA
aa.tree.dir <- paste(homePath,"simTreeOutputNoSnothaAAV3/",sep="")
treeDir <- aa.tree.dir
######Compute deltas for sim tree
likeFiles <- dir(treeDir,pattern="site",full.names = T)
#need to think about where to store simTreeDelta
outDir <- paste(mainDir,"simTreeAnalysisNoSnothaAAV3/",sep="")
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
      write.table(x=meanDelta,file=paste(outDir,gene,"SimTreeAnalysisNoSnothaAAV3.csv",sep=""),sep=",",row.names=F,quote=F)
    } else {
      print(paste("Not right num of trees for",gene))
      write.table(x=ii,file=paste(outDir,"errorLog.txt",sep = ""),append=T,quote=F,row.names=F)
    } 
  } else {
    write.table(x=ii,file=paste(outDir,"errorLog.txt",sep = ""),append=T,quote=F,row.names=F)
  }
}
}

## now for AA
aa.data <- read.table("/xdisk/rlapoint/rScripts/wilcoxSignifPvaluesNoSnothaAAV3.txt",header=T,sep="\t",stringsAsFactors=F)

#Find simulated means under tree sims
simDeltaDir <- "/xdisk/rlapoint/simTreeAnalysisNoSnothaAAV3/"
for (ii in aa.data$GeneID){
  simFile <-  dir(simDeltaDir,pattern=as.character(ii),full.names = T)
  if (length(simFile)>0){
    simDelta <- read.csv(simFile)
    observed <- aa.data[aa.data$GeneID == ii,"pvalue"]
    sims <- simDelta[,"pvalue"]
    pvalue <- length(sims[sims <= observed])/length(sims)
    aa.data[aa.data$GeneID == ii,"TreePvalue"] <- pvalue
  }
}
write.table(x=aa.data[order(aa.data$TreePvalue),],file=paste(simDeltaDir,"treeSimPvalueNoSnothaAAV3.txt",sep = ""),col.names=T,row.names=F,sep="\t",quote=F)

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

