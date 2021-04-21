## Robin Rabier, Adèle Erlichman, Loïc Lesobre, Alexandre Robert##
################## script created by Adèle Erlichman ####################
#########################################################################

## script used in order to simulate source and captive populations


library(dplyr)
library(stringr)
library(optiSel)
library(ggplot2)
library(data.table)
library(proto)
library(gsubfn)
library(RSQLite)
library(lubridate)
library(FactoMineR)
library(modeest)
library(miceadds)
  
setwd("directory")

# Functions ----
  # Pedigree and gene drop simulation functions based on the function simul.pedigree (synbreed package ; Wimmer et al. 2012)
    mysimul.pedigree <- function(generations=10,ids=10,genedrop=FALSE,foundergen,comparegenedrop=FALSE){
      # 1. Pedigree simulation
        # If only one value is given, set same number for all generations
        if(length(ids)==1) ids <- rep(ids,times=generations)
        
        # Initialization
          # Number of generations
          gener <- rep(1:generations,times=ids)
          # Individuals
          ID <- 1:sum(ids)
          # Parents
          Par1 <- Par2 <- rep(0,length(ID))
        
        # Define sex for the first generation (0 = female, 1 = male)
        sex <- rep(0,length(ID))
        sex[gener==1] <- sample(rep(c(0,1),length=sum(gener==1)),sum(gener==1),replace=FALSE)
        
        for (i in 1:generations) {
          if(i==1){ 
            # First generation (founders)
            founders <- data.frame(ID=ID[gener==1],Par1=Par1[gener==1],Par2=Par2[gener==1],generations=i,sex=sex[gener==1])
            
          } else if(i>1) {
            # Define which individuals from generation i-1 are dams/sires based on sex
            dams <- ID[which(gener==(i-1) & sex==0)] 
            sires <- ID[which(gener==(i-1) & sex==1)] 
            
            # Randomly attribute a dam and a sire to each individual of generation i (parents can have multiple offspring)
            Par1 <- sample(sires, ids[i], replace=TRUE) 
            Par2 <- sample(dams, ids[i], replace=TRUE) 
            
            # Define the sex of each offspring
            sex[gener==i] <- sample(rep(c(0,1),length=sum(gener==i)),sum(gener==i),replace=FALSE)
            
            # Compute pedigree of generation i
            ped <- data.frame(ID=ID[gener==i],
                              Par1=Par1,
                              Par2=Par2,
                              generations=i,
                              sex=sex[gener==i])
            
            if(i==2) { # Bind founder pedigree to second generation
              pedig <- rbind(founders, ped)
            } else if (i>2) { # Add all next generations
              pedig <- rbind(pedig, ped)
            }
          }
        } 
        
        # Prepare pedigree for optiSel
        pedig$Par1[pedig$Par1==0]<-NA
        pedig$Par2[pedig$Par2==0]<-NA
        pedig$sex[pedig$sex==1]<-2
        pedig$sex[pedig$sex==0]<-1
        colnames(pedig) <-c ("Indiv", "Sire", "Dam", "Born", "Sex")
        pedigree <- optiSel::prePed(pedig)
        
      # 2. Gene Drop simulation (MacCluer et al. 1986)
      if(genedrop==TRUE) {
        # 2.1. Define founder genome
          # Randomly attribute alleles for founders (if no genome has been provided in foundergen)
          # Founders are considers unrelated and not inbred (ie. unique alleles are provided for each individual)
          if(missing(foundergen)) {
            # Define two lists of alleles filled with zeros for the whole pedigree
            all1 <- all2 <- rep(0,length(pedigree$Indiv))
            
            # Define list of unique alleles (total number of individuals times 2 for allele 1 and allele 2)
            fall <- 1:(2*length(pedigree$Indiv[pedigree$Born==1]))
            
            # Define allele 1 and 2 for each founder
            all1[pedigree$Born==1] <- fall[1:length(pedigree$Indiv[pedigree$Born==1])]
            all2[pedigree$Born==1] <- fall[(length(pedigree$Indiv[pedigree$Born==1])+1):max(fall)]
            
            # Bind founder genome to the pedigree
            pedigree <- cbind(pedigree,all1)
            pedigree <- cbind(pedigree,all2)
            
          } else {
            # If a pre-defined genome has been provided for founders in foundergen
            # Check that the length of foundergen matches the number of founders
            if(length(pedigree$Indiv[pedigree$Born==1])==length(foundergen[,1])) {
              # Define two lists of alleles filled with zeros for the whole pedigree
              all1 <- all2 <- rep(0,length(pedigree$Indiv))
              
              # Define allele 1 and 2 for each founder
              all1[pedigree$Born==1] <- foundergen[,1]
              all2[pedigree$Born==1] <- foundergen[,2]
              
              # Bind founder genome to the pedigree
              pedigree <- cbind(pedigree,all1)
              pedigree <- cbind(pedigree,all2)
              
              # Comparegenedrop is used if a pre-defined genome has been provided for founders in foundergen but one wants to still simulate a founder genome where founders are considers unrelated and not inbred (wf : Without Founder genome)
              if(comparegenedrop==TRUE) {
                # Define two lists of alleles filled with zeros for the whole pedigree
                all1wf <- all2wf <- rep(0,length(pedigree$Indiv))
                
                # Define list of unique alleles
                fallwf <- 1:(2*length(pedigree$Indiv[pedigree$Born==1]))
                
                # Define allele 1 and 2 for each founder
                all1wf[pedigree$Born==1] <- fallwf[1:length(pedigree$Indiv[pedigree$Born==1])]
                all2wf[pedigree$Born==1] <- fallwf[(length(pedigree$Indiv[pedigree$Born==1])+1):max(fallwf)]
                
                # Bind unrelated founder genome to the pedigree
                pedigree <- cbind(pedigree,all1wf)
                pedigree <- cbind(pedigree,all2wf)
              }
            } else {
              cat("Warning: founder genome length is different than the number of founders. Please correct it.\n")
            }
          } 
          
        }
      
      # 2.1. Gene drop simulation
        pedigree$Sire <- as.numeric(pedigree$Sire)
        pedigree$Dam <- as.numeric(pedigree$Dam)
        
        # For each individual
        for(j in (length(pedigree$Born[pedigree$Born==1])+1):length(pedigree$Indiv)) {
          # Parents of each individual
          sireind <- pedigree$Sire[pedigree$Indiv==j]
          damind <- pedigree$Dam[pedigree$Indiv==j]
          
          # Both sire alleles
          all1sireind <- pedigree$all1[sireind]
          all2sireind <- pedigree$all2[sireind]
          allsire <- c(all1sireind,all2sireind)
          
          # Both dam alleles
          all1damind <- pedigree$all1[damind]
          all2damind <- pedigree$all2[damind]
          alldam <- c(all1damind,all2damind)
          
          # Sample 1 allele for each parent to define offspring genome
          pedigree$all1[pedigree$Indiv==j] <- sample(allsire,1)
          pedigree$all2[pedigree$Indiv==j] <- sample(alldam,1)
        }
        
        # This is the case where one founder genome has been given in foundergen and one has been added in the comparegenedrop part
        # Same process as above
        if(comparegenedrop==TRUE) {
          for(j in (length(pedigree$Born[pedigree$Born==1])+1):length(pedigree$Indiv)) {
            sireind <- pedigree$Sire[pedigree$Indiv==j]
            damind <- pedigree$Dam[pedigree$Indiv==j]
            
            all1sireind <- pedigree$all1wf[sireind]
            all2sireind <- pedigree$all2wf[sireind]
            allsire <- c(all1sireind,all2sireind)
            
            all1damind <- pedigree$all1wf[damind]
            all2damind <- pedigree$all2wf[damind]
            alldam <- c(all1damind,all2damind)
            
            pedigree$all1wf[pedigree$Indiv==j] <- sample(allsire,1)
            pedigree$all2wf[pedigree$Indiv==j] <- sample(alldam,1)
          }
        }
      return(pedigree)
    }
    
  # This function transforms a kinship matrix into a table
      transforkinmmatrix <- function(kin) {
        ut <- upper.tri(kin, diag = TRUE)
        row = rownames(kin)[row(kin)[ut]]
        column = rownames(kin)[col(kin)[ut]]
        kin  = (kin)[ut]
        data.frame(row,column,kin)
      }
  
# Initialization ----
  # Parameters for the source population
    g = 25 # Number of generations
    n = 100 # Number of individuals per generations
    r = 500 # Number of iterations
    a = (g*n-n)+1 # First individual of the last generation
    b = g*n # Last individual of the last generation
    nfounders = 5 # Number of founders drawn from the source population
    founders <- 1:nfounders # Founders ID
    
  # Parameters for the captive population
    g = 25 # Number of generations
    m = 100 # Number of individuals per generation
  
  # Objects
    compteur1 = 0
    compteur2 = 0
    ped1 = list() # Pedigree of the source population
    kin1 = list() # Kinship matrix of the last generation of the source population
    IndInbreeding1 = list() # Values of inbreeding for all individuals of the source population
    inds = list()
    fkin = list()
    fg = list()
    tablefkin = list()
    ped2 = list()
    kin2 = list()
    meankin2 = list()
    kin_corrected = list()
    meankin_corrected = list()
    IndInbreeding2 = list()

# Simulation 1 : source population ----
  T1 <- Sys.time() # Set timer

  # This simulates the pedigree and genome, calculates kinship and meankinship
  for(i in 1:r) {
    # Pedigree and genome simulation
    ped1[[i]] <- mysimul.pedigree(generations=g, ids=n, genedrop=TRUE)
    names(ped1)[i] = i
    
    # Calculates pedigree based inbreeding
    IndInbreeding1[[i]] <- pedInbreeding(ped1[[i]])
    names(IndInbreeding1)[i] = i
    
    # Calculates pedigree based kinship matrix
    kin1[[i]]<-pedIBD(ped1[[i]])
    kin1[[i]]<-as(kin1[[i]], "matrix")
    names(kin1)[i] = i
    kin1[[i]] <- round(kin1[[i]], 3)

    # Keeping only the last generation
    kin1[[i]] <- kin1[[i]][c(a:b),c(a:b)]
    ped1[[i]] <- ped1[[i]][ped1[[i]]$Born==g,]
    
    print(i)
  }
  
  T2 <- Sys.time()
  print(T2-T1)

# Saving files
  # save(kin1,file="path/kin1_simulID.Rdata")
  # save(ped1,file="path/ped1_simulID.Rdata")
  # save(IndInbreeding1,file="path/IndInbreeding1_simulID.Rdata")
  
# Founders : Random draw of individuals from the last generation of the source population ----
  # Loading files (optionnal)
    # load.Rdata(filename = "ped1_simulID.Rdata", "ped1")
    # load.Rdata(filename = "kin1_simulID.Rdata", "kin1")
    # load.Rdata(filename = "IndInbreeding1_simulID.Rdata", "IndInbreeding1")

  for(i in 1:r) {
    # Random draw of individuals
    inds[[i]] <- rownames(kin1[[i]])
    inds[[i]] <- sample(inds[[i]], nfounders)
    names(inds)[i] = i
    
    # Founder kinship
    fkin[[i]] <- data.frame(kin1[[i]][rownames(kin1[[i]]) %in% inds[[i]] == "TRUE",
      colnames(kin1[[i]]) %in% inds[[i]] == "TRUE"])
    names(fkin)[i] = i
    
    # Inbreeding of founders
    IndInbreeding1[[i]] <- data.frame(Indiv=IndInbreeding1[[i]]$Indiv[IndInbreeding1[[i]]$Indiv %in% inds[[i]] == "TRUE"], Inbr=IndInbreeding1[[i]]$Inbr[IndInbreeding1[[i]]$Indiv %in% inds[[i]] == "TRUE"])
    
    # Founder genome
    fg[[i]] <- ped1[[i]][,c(1,8,9)]
    fg[[i]] <- fg[[i]][fg[[i]]$Indiv %in% inds[[i]] == "TRUE",]
    fg[[i]] <- fg[[i]][,2:3]
    names(fg)[i] = i
    
    # Founders' kinship matrix as a table
    tablefkin[[i]] <- transforkinmmatrix(fkin[[i]])
    names(tablefkin)[i] = i
    tablefkin[[i]]$repet <- i
  }
  
  # Summary table of all individuals drawn to be founders
  tablefkin <- do.call(rbind, tablefkin)
  
  # Summary table of founders' alleles drawn
  tablefg <- do.call(rbind, fg)
  ids <- rep(nfounders,times=r) 
  repet <- rep(1:r,times=ids)
  tablefg <- cbind(tablefg,repet)
  
  # save(fkin,file="path/fkin_simulID.Rdata")
  # save(IndInbreeding1,file="path/IndInbreeding1_simulID.Rdata") 

# Simulation 2 : captive population ----
  # load.Rdata(filename = "fkin_simulID.Rdata", "fkin")
  # load.Rdata(filename = "ped2_simulID.Rdata", "ped2")
  # load.Rdata(filename = "IndInbreeding2_simulID.Rdata", "IndInbreeding2")
  
  for(i in 1:r) {
    # Founders' kinship matrix 
    rownames(fkin[[i]]) <- founders
    colnames(fkin[[i]]) <- founders
    fkin[[i]] <- as.matrix(fkin[[i]])
    
    # Simulation of pedigrees with founder genome (from source population) provided in the gene drop
    ped2[[i]] <- mysimul.pedigree(generations=g,
      ids=c(nfounders,m,m,m,m,
        m,m,m,m,m,
        m,m,m,m,m,
        m,m,m,m,m,
        m,m,m,m,m),
      genedrop=TRUE,
      foundergen = fg[[i]],
      comparegenedrop = TRUE,
      familySize = 0)
    names(ped2)[i] = i
    
    # Calculates pedigree based inbreeding
    IndInbreeding2[[i]] <- pedInbreeding(ped2[[i]])
    names(IndInbreeding2)[i] = i

    # Keeping generations 5, 10 and 25
    keep = c(ped2[[i]]$Indiv[ped2[[i]]$Born==5],
             ped2[[i]]$Indiv[ped2[[i]]$Born==10],
             ped2[[i]]$Indiv[ped2[[i]]$Born==25])
    
    # Calculates pedigree based kinship matrix (not considering kinship matrix for the founders)
    kin2[[i]]<-pedIBD(ped2[[i]],keep.only = keep)
    kin2[[i]]<-as(kin2[[i]], "matrix")
    names(kin2)[i] = i
    kin2[[i]] <- round(kin2[[i]], 3)
    
    # Calculates pedigree based kinship matrix (considering kinship matrix for the founders)
    kin_corrected[[i]]<-pedIBD(ped2[[i]], kinFounder = fkin[[i]],keep.only = keep)
    kin_corrected[[i]]<-as(kin_corrected[[i]], "matrix")
    names(kin_corrected)[i] = i
    kin_corrected[[i]] <- round(kin_corrected[[i]], 3)
  
    print(i)
   
  }
  
  T2 <- Sys.time()
  print(T2-T1)
 
# save(ped2,file="path/ped2_simulID.Rdata")
# save(IndInbreeding2,file="path/IndInbreeding2_simulID.Rdata")
# save(kin2,file="path/kin2_simulID.Rdata")
# save(kin_corrected,file="path/kin_corrected_simulID.Rdata")



