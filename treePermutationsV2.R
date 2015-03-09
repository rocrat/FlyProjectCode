#get random trees for not all possible combo route

notSpecial<-c("Dmelanogaster","Dyakuba","Dbiarmipes","Dpseudoobscura","Sflava","Dmojavensis")
special<-c("Dgrimshawi","Dsuzukii")
#individual functions to produce permuted trees for each of 6 templates
#swapping the the two species would create an additional 6 templates but 
#I am not sure this is necessary
tree1<-function(spplist){#pass a single row from the permutation list
  tempTree <- paste("tree 'h",countTree,"' = ",sep="")
  tempTree <- paste(tempTree,"((((",spplist[1],",",spplist[2],"),(",spplist[3],",","Dsuzukii",")),",spplist[4],"),((",spplist[5],",","Dgrimshawi","),",spplist[6],"));",sep="") 
  return(tempTree)
}
tree2<-function(spplist){#pass a single row from the permutation list
  tempTree <- paste("tree 'h",countTree,"' = ",sep="")
  tempTree <- paste(tempTree,"((((",spplist[1],",",spplist[2],"),(",spplist[3],",",spplist[4],")),","Dsuzukii","),((",spplist[5],",","Dgrimshawi","),",spplist[6],"));",sep="") 
  return(tempTree)
}
tree3<-function(spplist){#pass a single row from the permutation list
  tempTree <- paste("tree 'h",countTree,"' = ",sep="")
  tempTree <- paste(tempTree,"((((",spplist[1],",","Dsuzukii","),(",spplist[2],",",spplist[3],")),",spplist[4],"),((",spplist[5],",","Dgrimshawi","),",spplist[6],"));",sep="") 
  return(tempTree)
}
tree4<-function(spplist){#pass a single row from the permutation list
  tempTree <- paste("tree 'h",countTree,"' = ",sep="")
  tempTree <- paste(tempTree,"((((",spplist[1],",",spplist[2],"),(",spplist[3],",","Dsuzukii",")),",spplist[4],"),((",spplist[5],",",spplist[6],"),","Dgrimshawi","));",sep="") 
  return(tempTree)
}
tree5<-function(spplist){#pass a single row from the permutation list
  tempTree <- paste("tree 'h",countTree,"' = ",sep="")
  tempTree <- paste(tempTree,"((((",spplist[1],",",spplist[2],"),(",spplist[3],",",spplist[4],")),","Dsuzukii","),((",spplist[5],",",spplist[6],"),","Dgrimshawi","));",sep="") 
  return(tempTree)
}
tree6<-function(spplist){#pass a single row from the permutation list
  tempTree <- paste("tree 'h",countTree,"' = ",sep="")
  tempTree <- paste(tempTree,"((((",spplist[1],",","Dsuzukii","),(",spplist[2],",",spplist[3],")),",spplist[4],"),((",spplist[5],",",spplist[6],"),","Dgrimshawi","));",sep="") 
  return(tempTree)
}

funs<-list(tree1,tree2,tree3,tree4,tree5,tree6)#list of functions which move the special spp in non-convergent ways

#enumerate all possible constrained permutations
require(combinat)
sppPermuts<-permn(notSpecial)#list of all possible mig combs
#function to paste all possible combinations of the spp

#Write in Nexus format
countTree <- 1
toWrite <- c("#NEXUS","begin trees;")
#Write Ho tree first
toWrite <- c(toWrite, "tree 'h0' = ((((Dmelanogaster,Dyakuba),(Dbiarmipes,Dsuzukii)),Dpseudoobscura),((Dgrimshawi,Sflava),Dmojavensis));")
for (i in 1:length(sppPermuts)){
  splist<-sppPermuts[[i]]
  #browser()
  for (j in 1:length(funs)){
      tempTree<-funs[[j]](splist)
      toWrite <- c(toWrite,tempTree)
      countTree <- countTree + 1
    
    #browser()
  }
}
toWrite <- c(toWrite,"end;")

##Sample 150 trees for now
set.seed(44)
ind <- sample(4:(length(toWrite)-1),150)
ind <- c(1:3,ind,length(toWrite))
# writeLines(text=toWrite[ind],con="/xdisk/rlapoint/garliConfFiles/simTreesV1.nexus",sep="\n")
writeLines(text=toWrite[ind],con="simTreesNonsenseV1.nexus",sep="\n")

###############testing
tree4(notSpecial)
#Begin trees;
#tree 1 = ((((Dmelanogaster,Dyakuba),(Dbiarmipes,Dsuzukii)),Dpseudoobscura),(((Dgrimshawi,Sflava)),Dmojavensis));
#END;
