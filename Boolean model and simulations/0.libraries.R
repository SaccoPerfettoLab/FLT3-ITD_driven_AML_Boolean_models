library(CellNOptR)
library(tidyverse)
library(CNORfeeder)
library(pheatmap)
library("grid") 
library("gridExtra")
library(ggrepel)
library("ggdendro")
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library('ggfortify')

# --- FUNCTIONS --- #

# this script contains all the custom functions necessary for the analysis

# Name: equation of Hill function
# Hill function
# input: 
hill <- function(x, n, K){return(x^n / (K+x^n))}

v_hill <- function(v, n, K){return(sapply(v, function(x){hill(x,n,K)}))}

# Name: Normalize raw Luminex Data
# Input: path to midas file
# Output: dataframe with normalized values
normalization <- function(midas_path){
  
  CNOlist <- CNOlist(midas_path)
  
  # creating raw data df
  df <- from_midas_to_df(CNOlist)
  
  proteins <- colnames(df) #n° of analytes
  conditions<-row.names(df)
  # create normalized df container
  df_normalized <- as.data.frame(matrix(nrow=nrow(df), ncol = length(proteins)))
  names(df_normalized) <- proteins
  row.names(df_normalized) <- conditions
  
  #normalization 
  for (protein in proteins){
    #setting parameters
    maxT = 0.999 #teoric maximum
    minT = 0.001 #teoric minimum
    maxS = max(df[,protein])  #experimental maximum
    minS = min(df[,protein])  #experimental minimum
    
    log_arg = ((1-minT)*maxT)/((1-maxT)*minT)
    b = maxS/minS
    
    # calculation of n and K Hill parameters
    n = round(log(log_arg, base = b))
    K = ((1 - maxT)/maxT)*maxS^n
    
    #apply Hill function
    df_normalized[,protein] <- unlist(lapply(df[,protein], 'hill', n, K))
  }
  
  return(df_normalized)
}

# Name: From MIDAS to dataframe
# Input: CNOlist
# Output: dataframe with values 
from_midas_to_df <- function(CNOlist){
  
  #CNOlist <- rawdata
  
  # get condition from MIDAS
  cond <- CNOlist@cues
  row.names(cond) <- LETTERS[1:nrow(cond)]
  colnames(cond)[ colnames(cond) %in% colnames(CNOlist@inhibitors)] <- paste0(  colnames(cond)[ colnames(cond) %in% colnames(CNOlist@inhibitors)], 'i')
  
  for(i in c(1:nrow(cond))){
    rownames(cond)[i] <- paste0(colnames(cond)[cond[i,] == 1], collapse = '+')
  }
  conditions <- rownames(cond)
  
  
  # creating raw data df
  df <- as.data.frame(getSignals(CNOlist)$'90') #taking T90 time
  df[nrow(df)+1,] <- getSignals(CNOlist)$'0'[1,] #adding a new row for T0 (CTRL)
  row.names(df) <- c(conditions, 'CTRL')
  
  return(df)
}

# Name: Plot normalization
# Input: path of raw midas, normalized midas, path of output folder, cellline
# Output: write png files of normalized data

plot_normalization <- function(midas_path_raw, midas_path_norm,
                               graphs_path = './graphs/', 
                               cellline){
  
  # midas_path_raw <- paste0("./input/RAW_", cellline, "_MD.csv")
  # midas_path_norm <- paste0("./input/MD_", cellline, '_NORM.csv')
  # midas_path_raw<-"./input/TKD_MD.csv"
  # midas_path_norm<-"./results/MD_TKD_NORM.csv"
  # cellline = "TKD" 
  # graphs_path = './graphs/'
  
  rawdata<-CNOlist(midas_path_raw)
  normalizeddata<-CNOlist(midas_path_norm)
  raw_df <- from_midas_to_df(rawdata)
  norm_df <- from_midas_to_df(normalizeddata)
  
  rownames(raw_df) <- rownames(norm_df)
  colnames(raw_df) <- colnames(norm_df)
  
  proteins <- colnames(raw_df) 
  conditions <- rownames(norm_df)
  
  for (protein in proteins){
    print(protein)
    #protein<-"NFKB"
    x <- seq(1,length(conditions),1) #11 condizioni + 1 ctrl
    y <- sort(raw_df[,protein]) #sorted raw data of a desired protein
    z <- sort(norm_df[,protein]) #sorted normalized data of a desired protein
    
    (xlabs_sorted <- conditions[order(raw_df[,protein])]) #sorting condition labels
    
    graph_path <- paste(graphs_path, 
                        cellline, '_', protein, '_raw_vs_normalized_united', '.pdf', sep = '')
    pdf(graph_path,
        width = 6 ,height= 6, #units = "cm",
        pointsize = 12)#creating png file of the protein
    
    par(mar=c(7,7,5,5) + .3) #define margins
    plot(x, y, pch = 17, xaxt = 'n', yaxt ='n', col='lightsalmon2',
         xlab="", ylim = c(0,max(y)), ylab="",
         main=paste(protein,'Data'),
         cex.axis=0.1)
    
    axis(side = 1, at = 1:length(xlabs_sorted), labels = xlabs_sorted, cex.axis = 0.8, las=2) # X AXIS
    
    if(protein=="NFKB"&cellline=="TKD"){digits=5}else{digits=3}
    
    axis(side = 2, seq(0,max(y)+round((max(y)-0)/10, digits = digits),
                       by=round((max(y)-0)/10,digits=digits)), las=2)
    
    par(new = TRUE)
    plot(x,z, axes = FALSE, pch=15, xlab = "", ylab="", col='steelblue3', ylim = c(0,1))
    axis(side = 4, seq(0,1,by=0.1), col= 'black', las = 2) #Z AXIS
    
    #text(par("usr")[2]+1.5,mean(par("usr")[3:4]+.2), "Normalized Data", srt = -90, xpd = TRUE, pos = 4)
    #text(par("usr")[1]-2.3, mean(par("usr")[3:4]+.05), "Raw Data", srt = 90, pos = 1, xpd=TRUE)
    
    legend("topleft",legend=c("raw","normalized"),
           col=c("lightsalmon2","steelblue3"),
           pch = c(17,15),cex = 0.8)
    
    dev.off()
  }
  return('Files created')
}

RemoveMultiplicity <- function(df, idx_sw, idx_multiplicity){
  
  # df <- df_phospho_sign
  # idx_sw = 4
  # idx_multiplicity = 7
  seqs <- as.vector(unlist(unique(df[,idx_sw])))
  
  for (seq in seqs){
    
    
    
    pos <- which(df[,idx_sw] == seq)
    if (length(pos)==1){
      next
    }else{
      subset <- df[pos,]
      pos_to_remove <- pos[ which(unlist(subset[, idx_multiplicity]) != min(unlist(subset[, idx_multiplicity]))) ]
      
      if (length(pos_to_remove) >= 1){
        df <- df[-c(pos_to_remove),]
      }else{next}
    }
  }
  return(df)
}



writeScaffold<-function(
    modelComprExpanded,
    optimResT1,
    optimResT2,
    modelOriginal,
    CNOlist,
    tag=NULL){
  
  
  #get the stuff that I need for the sif file
  sif<-getSifInfo(modelComprExpanded=modelComprExpanded,
                  optimResT1=optimResT1,
                  optimResT2=optimResT2,
                  modelOriginal=modelOriginal)
  
  #get the stuff that I need for the dot file
  dot<-getDotInfo(
    modelComprExpanded=modelComprExpanded,
    modelOriginal=modelOriginal,
    CNOlist=CNOlist,
    sifFile=sif$sifFile)
  
  #this bit writes the dot, the sif and the sif edges attributes
  writeScaffoldW(
    dN=dot$dN,
    dM=dot$dM,
    modelComprExpanded=modelComprExpanded,
    sifFile=sif$sifFile,
    EApresent=sif$EApresent,
    EAweights=sif$EAweights,
    tag=tag)
  
}


###########################################################################
#######these are the functions used above
###########################################################################
#this function writes the sif file, the edge attributes, and the dot file
writeScaffoldW<-function(
    dN,
    dM,
    modelComprExpanded,
    sifFile,
    EApresent,
    EAweights,
    tag=NULL){
  
  
  create_filename<-function(x, tag=NULL){
    if (is.null(tag)){
      return(x)
    }
    else{
      return(paste(c(tag, "_", x), collapse=""))
    }
  }
  
  # writeDot(
  #     dotNodes=dN,
  #     dotMatrix=dM,
  #     model=modelComprExpanded,
  #     filename=create_filename("Scaffold.dot", tag=tag))
  
  write.table(
    sifFile[,1:3],
    file=create_filename("Scaffold.sif", tag=tag),
    row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  
  write.table(
    EApresent,
    file=create_filename("TimesScaffold.EA", tag=tag),
    row.names=FALSE,col.names="Times",quote=FALSE,sep="\t")
  
  write.table(
    EAweights,
    file=create_filename("weightsScaffold.EA", tag=tag),
    row.names=FALSE,col.names="Weights",quote=FALSE,sep="\t")
  
}

######
#this function computes the stuff that is needed for the dot file

getDotInfo<-function(modelComprExpanded,modelOriginal,CNOlist,sifFile){
  
  dM<-sifFile
  nodesCompr<-modelComprExpanded$speciesCompressed
  indexes<-indexFinder(CNOlist,modelOriginal)
  nodesSig<-modelOriginal$namesSpecies[indexes$signals]
  nodesInh<-modelOriginal$namesSpecies[indexes$inhibited]
  nodesStim<-modelOriginal$namesSpecies[indexes$stimulated]
  nodesNCNO<-findNONC(modelOriginal,indexes)
  nodesAttrNames<-nodesSig
  nodesAttr<-rep("signal",length(nodesSig))
  
  if(length(nodesInh) != 0){
    nodesAttrNames<-c(nodesAttrNames,nodesInh)
    nodesAttr<-c(nodesAttr,rep("inhibited",length(nodesInh)))
  }
  
  if(length(nodesStim) != 0){
    nodesAttrNames<-c(nodesAttrNames,nodesStim)
    nodesAttr<-c(nodesAttr,rep("stimulated",length(nodesStim)))
  }
  
  if(length(nodesNCNO) != 0){
    nodesAttrNames<-c(nodesAttrNames,nodesNCNO)
    nodesAttr<-c(nodesAttr,rep("ncno",length(nodesNCNO)))
  }
  
  if(length(nodesCompr) != 0){
    nodesAttrNames<-c(nodesAttrNames,nodesCompr)
    nodesAttr<-c(nodesAttr,rep("compressed",length(nodesCompr)))
  }
  
  dN<-cbind(nodesAttrNames,nodesAttr)
  return(list(dN=dN,dM=dM))
}

#####
#this function computes the stuff that is needed for the sif file
getSifInfo<-function(modelComprExpanded,
                     optimResT1,
                     optimResT2,
                     modelOriginal){
  
  bString1<-optimResT1$bString
  
  if(is.na(optimResT2[1])){
    bString2<-optimResT1$bString[which(optimResT1$bString == 0)]
  }else{
    bString2<-optimResT2$bString
  }
  
  #BStimes is a string containing 0,1 and 2 depending on whether the interaction is
  #absent, present at t1 or present at t2
  BStimes<-bString1
  BStimes[which(BStimes == 0)]<-bString2*2
  
  #create: weightsE is a string that holds the weights of the interactions
  
  if(!is.null(dim(optimResT1$stringsTol))){
    bW1<-apply(optimResT1$stringsTol,2,mean)
  }else{
    bW1<-bString1
  }
  
  if(!is.na(optimResT2[1])){
    
    if(!is.null(dim(optimResT2$stringsTol))){
      bW2<-apply(optimResT2$stringsTol,2,mean)
    }else{
      bW2<-bString2
    }
    
    weightsE<-bW1
    weightsE[which(optimResT1$bString == 0)]<-weightsE[which(optimResT1$bString == 0)]+bW2
    
  }else{
    weightsE<-bW1
  }
  
  #These mini functions are used to find the inputs and output of reactions
  findOutput<-function(x){
    sp<-which(x == 1)
    sp<-modelComprExpanded$namesSpecies[sp]
  }
  
  reacOutput<-apply(modelComprExpanded$interMat,2,findOutput)
  
  findInput<-function(x){
    sp<-which(x == -1)
    sp<-modelComprExpanded$namesSpecies[sp]
  }
  
  reacInput<-apply(modelComprExpanded$interMat,2,findInput)
  
  #This mini function is used to create a reaction label as used in a cystoscape edge attribute file
  createReac<-function(x){
    r<-paste(x[1]," (",x[2],") ",x[3],sep="")
    return(r)
  }
  
  #if the class of reacInput is not a list, then there are no AND gates
  if(class(reacInput) != "list"){
    isNeg<-function(x){
      isNegI<-any(x == 1)
      return(isNegI)
    }
    inpSign<-apply(modelComprExpanded$notMat,2,isNeg)
    inpSign<-!inpSign
    inpSign[inpSign]<-1
    inpSign[!inpSign]<--1
    sifFile<-cbind(reacInput,inpSign,reacOutput)
    EApresent<-apply(sifFile,1,createReac)
    EApresent<-cbind(EApresent,BStimes)
    EAweights<-cbind(EApresent,weightsE)
    
    # add a fourth and fifth column as expected in writeDot (bug report 39)
    sifFile<-cbind(sifFile,BStimes)
    sifFile<-cbind(sifFile,weightsE)
  }
  else{
    #in this case there are AND gates and so we need to create dummy "and#" nodes
    sifFile<-matrix(0,nrow=4*length(reacOutput),ncol=5)
    nR<-1
    nANDs<-1
    
    for(i in 1:length(reacOutput)){
      if(length(reacInput[[i]]) == 1){
        sifFile[nR,1]<-reacInput[[i]]
        sifFile[nR,3]<-reacOutput[i]
        sifFile[nR,2]<-ifelse(any(modelComprExpanded$notMat[,i] == 1),-1,1)
        sifFile[nR,4]<-BStimes[i]
        sifFile[nR,5]<-weightsE[i]
        nR<-nR+1
      }
      else{
        for(inp in 1:length(reacInput[[i]])){
          sifFile[nR,1]<-reacInput[[i]][inp]
          sifFile[nR,3]<-paste("and",nANDs,sep="")
          temp_indices = which(reacInput[[i]][inp]==rownames(modelComprExpanded$notMat))
          sifFile[nR,2]<-ifelse(modelComprExpanded$notMat[temp_indices, i]==1,-1,1)
          sifFile[nR,4]<-BStimes[i]
          sifFile[nR,5]<-weightsE[i]
          nR<-nR+1
        }
        sifFile[nR,1]<-paste("and",nANDs,sep="")
        sifFile[nR,3]<-reacOutput[i]
        sifFile[nR,2]<-1
        sifFile[nR,4]<-BStimes[i]
        sifFile[nR,5]<-weightsE[i]
        nANDs<-nANDs+1
        nR<-nR+1
      }
    }
    sifFile<-sifFile[1:(nR-1),]
    EApresent<-apply(sifFile[,1:3],1,createReac)
    EAweights<-cbind(EApresent,sifFile[,5])
    EApresent<-cbind(EApresent,sifFile[,4])
  }
  #this mini function makes edge attributes in the format e1 (sign) e2 = attr
  
  makeEA<-function(x){
    ea<-paste(x[1],"=",x[2])
    return(ea)
  }
  EApresent<-apply(EApresent,1,makeEA)
  EAweights<-apply(EAweights,1,makeEA)
  
  #this is the edge attribute matrix that contains, for each edge, whether it is
  #absent from the model (0), present at t1(1) or present at t2(2)
  return(list(EApresent=EApresent,EAweights=EAweights,sifFile=sifFile))
}




