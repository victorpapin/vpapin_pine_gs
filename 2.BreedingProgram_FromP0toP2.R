library(AlphaSimR)
library(faux)
library(stringr)

path = "C:/Users/Victor/Documents/Simulations_maritime_pine/PIPELINE_CENTRAL/0.InputsPipeline"
load(paste(path,"founderBIG.RData", sep="/"))



################################################
## DEFINE A TRAIT AND A POPULATION P0
################################################

# TRAIT
SP              = SimParam$new(founderGenomes)
phenoMean       = 6.81    # Extract from data base - cf Laurent
heritability    = 0.33    # Extract from data base - cf Laurent (constante au cours des genrations)
CV              = 0.0505  # Extract from data base - cf Laurent
genVar          = (CV * phenoMean)^2  # sqrt(variance)/meanP = CV 
accuracy.BV.EBV = 0.97
SP$addTraitA(nQtlPerChr = nQTLPerChr, mean = phenoMean, var = genVar, name="Ht")
SP$setVarE(h2=heritability)

# POPULATION P0
popP0        = newPop(founderGenomes)


# CREATE AN ARRAY
QTL      = SP$invalidSnp
snpArray = vector("list", nChr)
for (chr in 1:nChr) {
    x           = pullSegSiteHaplo(founderGenomes, chr = chr)
    alleleFreq  = apply(X = x, MARGIN = 2, FUN = mean)
    alleleFreq  = alleleFreq[-QTL[[chr]]]
    allSegSites = 1:(10 * nSNPPerChr + nQTLPerChr) 
    
    alleleFreq[which(alleleFreq>=0.0&alleleFreq<0.12)] = alleleFreq[which(alleleFreq>=0.0&alleleFreq<0.12)] * 2
    alleleFreq[which(alleleFreq>=0.85)]                = alleleFreq[which(alleleFreq>=0.85)]  * 2
    
    sel = sample(allSegSites[-QTL[[chr]]], nSNPPerChr, prob = alleleFreq)
    snpArray[[chr]] = sort(sel) 
  }
  
snpArray = do.call("c", snpArray)
snpArray = new(Class = "LociMap", nLoci = as.integer(nChr) * as.integer(nSNPPerChr), lociPerChr = rep(as.integer(nSNPPerChr), as.integer(nChr)), lociLoc = snpArray)
SP$snpChips[[1]] = snpArray
  
id.SNP.array = NULL
for(i in 1:nChr){
    snp.array.chr = snpArray@lociLoc[((i-1)*nSNPPerChr+1):((i-1)*nSNPPerChr+nSNPPerChr)] 
    snp.array.chr = snp.array.chr + (10*nSNPPerChr+nQTLPerChr) * (i-1)
    id.SNP.array = c(id.SNP.array, snp.array.chr)
  }


###################################################
## STEP 2: SELECTION OF POPULATION G0
###################################################

if(T){
  
  EBV = rnorm_pre(bv(popP0), mu=0, sd=sd(bv(popP0)), r=accuracy.BV.EBV)
  popP0@ebv = as.matrix(tail(EBV, popP0@nInd))
  
  BV.P0.simul           = data.frame(popP0@id) 
  BV.P0.simul$BV        = as.numeric(ebv(popP0))
  colnames(BV.P0.simul) = c("Genotype_id","EBV")
  
  # ranking real individuals based on EBV
  BV.P0.real       = read.table(paste(path,"EBV_G0_real.csv",sep="/"), header=TRUE, sep=";")
  BV.P0.real       = BV.P0.real[order(BV.P0.real$EBV, decreasing=TRUE),]
  BV.P0.real$rank  = c(1:nrow(BV.P0.real))
  
  # number of real data missing
  nb.P0.NA         = nInd - length(which(!is.na(BV.P0.real$EBV)))
  
  # assign missing data in our simulated pop 
  BV.P0.simul$EBV[sample(1:nrow(BV.P0.simul),nb.P0.NA)] = NA
  
  # ranking simulated individuals based on BV
  BV.P0.simul      = BV.P0.simul[order(BV.P0.simul$EBV, decreasing=TRUE),]
  BV.P0.simul$rank = c(1:nrow(BV.P0.simul))
  
  # assign new id for simulated individuals based on the corresponding real individual
  for(i in 1:nrow(BV.P0.simul)){
    if(is.na(BV.P0.simul[i,"EBV"])){
      BV.P0.simul[i,"new.id"]=paste("old",BV.P0.simul[i,"Genotype_id"], sep="_")
    }else{
      BV.P0.simul[i,"new.id"]=BV.P0.real$Genotype_id[match(BV.P0.simul[i,"rank"], BV.P0.real$rank)]
    }
  }
  
  # assign random founder (not already used) for unknown parents
  CrossG0G1        = read.table(paste(path,"Cross_G0_G1.csv", sep="/"), header=TRUE, sep=";")
  nb.unknown       = length(which(CrossG0G1$P1_G0==0))+length(which(CrossG0G1$P2_G0==0))
  unknown.founders = paste("founder",c(1:nb.unknown),sep="")
  random.id        = sample(which(is.na(BV.P0.simul$EBV)), length(unknown.founders))
  BV.P0.simul[random.id, "new.id"] = unknown.founders
  
  # assign in popP0 object
  popP0@id = BV.P0.simul$new.id[match(popP0@id, BV.P0.simul$Genotype_id)]
  
  # create population G0
  parentsG0  = c(unique(c(CrossG0G1$P1_G0,CrossG0G1$P2_G0)), unknown.founders)
  parentsG0  = parentsG0[which(!parentsG0=="0")]
  popG0      = selectInd(popP0, nInd=length(parentsG0), use="rand", candidates=which(popP0@id%in%parentsG0))
  
}

###################################################
## STEP 3: CREATE POPULATION P1
###################################################

if(T){
  nbIndPerCross = 130 # Reference breeding program
  
  # create a founder for each unknown parent
  a = 1
  for(j in 1:ncol(CrossG0G1)){
    
    for(i in 1:nrow(CrossG0G1))
      
      if(CrossG0G1[i,j]==0){
        
        CrossG0G1[i,j] = paste("founder",a,sep="")
        a = a + 1 
      }
  }
  
  # cross parameters
  cross.names = paste(CrossG0G1$P1_G0, CrossG0G1$P2_G0, sep="x")
  cross.nb    = data.frame(table(cross.names))
  mothers     = CrossG0G1$P1_G0
  mothers     = mothers[which(!duplicated(cross.names))]
  fathers     = CrossG0G1$P2_G0
  fathers     = fathers[which(!duplicated(cross.names))]
  parentsP0   = unique(c(mothers, fathers))
  crossPlan   = matrix(c(rep(mothers, each=nbIndPerCross),rep(fathers, each=nbIndPerCross)), nrow=length(mothers)*nbIndPerCross, ncol=2)
  popP1       = makeCross(popP0, crossPlan, simParam=SP)
  popP1       = setPheno(popP1)
}


###################################################
## STEP 4: SELECT INDIVIDUAL FOR POPULATION G1
###################################################

if(T){
  selection.rate = 1.07 # see script Real_selection_familyG1.R
  
  EBV = rnorm_pre(bv(popP1), mu=0, sd=sd(bv(popP1)), r=accuracy.BV.EBV)
  popP1@ebv = as.matrix(tail(EBV, popP1@nInd))
  
  # select individuals in each cross
  IND.SELECTED = NULL
  for(i in unique(cross.names)){
    
    nb.select = cross.nb$Freq[match(i,cross.nb$cross.names)]
    
    mother = str_split(i,pattern="x")[[1]][1]
    father = str_split(i,pattern="x")[[1]][2]
    
    ind.table     = data.frame(matrix(ncol=2, nrow=nbIndPerCross))
    ind.table$id  = popP1@id[which(popP1@mother==mother & popP1@father==father)]
    ind.table$ebv = ebv(popP1)[which(popP1@mother==mother & popP1@father==father)]
    ind.table     = ind.table[order(ind.table$ebv, decreasing=TRUE),]
    
    mean.fam     = mean(ind.table$ebv, na.rm=TRUE)
    sd.fam       = sd(ind.table$ebv, na.rm=TRUE)
    value.select = selection.rate * sd.fam + mean.fam
    
    for(j in 1:nb.select){
      ind.selected = ind.table[which.min(abs(ind.table$ebv-value.select)),"id"]
      ind.table    = ind.table[-which.min(abs(ind.table$ebv-value.select)),]
      
      IND.SELECTED = c(IND.SELECTED, ind.selected)
    }
  }
  
  # create pop G1
  popG1 = selectInd(popP1, nInd=length(IND.SELECTED), use="rand", candidates=which(popP1@id%in%IND.SELECTED))
  
  # attribute new IDs to select individuals accroding real pedigree
  NEW.ID = NULL
  for(i in 1:length(popG1@id)){
    
    cross       = paste( popG1@mother[i] , popG1@father[i] , sep="x" )
    new.id      = CrossG0G1[match(cross, cross.names),"G1"]
    id.rm       = which(CrossG0G1$G1==new.id)[1]
    CrossG0G1   = CrossG0G1[-id.rm,]
    cross.names = cross.names[-id.rm]
    NEW.ID      = c(NEW.ID, new.id)
  }
  popG1@id = NEW.ID
  popG1    = setPheno(popG1)
  
}

save.image("BreedingProgramBeforeG2.RData")

###################################################
## STEP 5: CREATE POPULATION P2
###################################################

if(T){
  CrossP2G1 = read.table(paste(path,"Cross_P2_G1.csv", sep="/"), header=T, sep=";")
  
  mothers   = CrossP2G1$P1_G1
  fathers   = CrossP2G1$P2_G1
  
  crossPlan = matrix(c(mothers, fathers), nrow=length(mothers), ncol=2)
  popP2     = makeCross(popG1, crossPlan, simParam=SP)
  popP2     = setPheno(popP2, h2=0.13)
}

save.image("BreedingProgram.RData")
