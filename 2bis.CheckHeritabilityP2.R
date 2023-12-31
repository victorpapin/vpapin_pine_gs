library(AlphaSimR)
library(AGHmatrix)
library(breedR)

path = "C:/Users/Victor/Documents/Simulations_maritime_pine/PIPELINE_CENTRAL/0.InputsPipeline"
load(paste(path,"BreedingProgram.RData", sep="/"))

## REAL DATA
phenoP2.real     = read.table(paste(path,"phenoP2real.txt",sep="/"), header=T)
genoP2.real      = read.table(paste(path,"genoP2real.txt", sep="/"), header=T)
pedigreeP2.real  = read.table(paste(path,"pedigreeP2real.txt",sep="/"), header=TRUE)

## SIMULATED DATA
phenoP2.simul              = data.frame(pheno(popP2))
genoP2.simul               = data.frame(pullSegSiteGeno(popP2)[,id.SNP.array])
row.names(genoP2.simul)    = popP2@id
phenoP2.simul$ind          = popP2@id
pedigreeP2.simul           = rbind(data.frame(ind = c(popP2@id), mom = c(popP2@mother), dad = c(popP2@father)), data.frame(ind = c(popG1@id), mom = c(popG1@mother), dad = c(popG1@father)), data.frame(ind = c(popG0@id), mom = c(popG0@mother), dad = c(popG0@father)) )
pedigreeP2.simul           = pedigreeP2.simul[!duplicated(pedigreeP2.simul$ind),]

## RELATIONSHIP MATRICES
G.simul  = Gmatrix(as.matrix(genoP2.simul), method='VanRaden', missingValue=-9, maf=0)
G.real   = Gmatrix(as.matrix(genoP2.real), method='VanRaden', missingValue=-9, maf=0)
A.simul  = Amatrix(pedigreeP2.simul, ploidy=2)
A.real   = Amatrix(pedigreeP2.real, ploidy=2)

A.real   = A.real[phenoP2.real$ind,phenoP2.real$ind]
A.simul  = A.simul[phenoP2.simul$ind,phenoP2.simul$ind]
colnames(G.simul) = paste("ind",colnames(G.simul), sep="");colnames(G.real) = paste("ind",colnames(G.real), sep="");colnames(A.simul) = paste("ind",colnames(A.simul), sep="");colnames(A.real) = paste("ind",colnames(A.real), sep="");row.names(G.simul) = paste("ind",row.names(G.simul), sep="");row.names(G.real) = paste("ind",row.names(G.real), sep="");row.names(A.simul) = paste("ind",row.names(A.simul), sep="");row.names(A.real) = paste("ind",row.names(A.real), sep="")

pedigreeP2.simul$cross  = paste(pedigreeP2.simul$mom,pedigreeP2.simul$dad, sep="x")
pedigreeP2.simul$family = phenoP2.real$family[match(pedigreeP2.simul$cross, phenoP2.real$cross)]
phenoP2.simul$family    = pedigreeP2.simul$family[match(phenoP2.simul$ind, pedigreeP2.simul$ind)]

# Estimates h2 for accuracy

  # REAL DATA - GENOMIC
  phenoP2.real$ind         = as.factor(phenoP2.real$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.real)
  cov.mat                  = as.matrix(G.real[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.real$ind         = as.character(phenoP2.real$ind)
  gen.data.genomic.real    = gen.data
  model.G.real             = remlf90(fixed = Ht ~ 1, generic = gen.data.genomic.real, data = phenoP2.real, method = 'em')
  h2.G.real                = model.G.real$var["generic_genetic",1]/(model.G.real$var["generic_genetic",1]+model.G.real$var["Residual",1])
  
  # SIMULATED DATA - GENOMIC 
  phenoP2.simul$ind = as.factor(phenoP2.simul$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.simul)
  cov.mat                  = as.matrix(G.simul[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.simul$ind        = as.character(phenoP2.simul$ind)
  gen.data.genomic.simul   = gen.data
  model.G.simul            = remlf90(fixed = Ht ~ 1, generic = gen.data.genomic.simul, data = phenoP2.simul, method = 'em')
  h2.G.simul               = model.G.simul$var["generic_genetic",1]/(model.G.simul$var["generic_genetic",1]+model.G.simul$var["Residual",1])
  
  # REAL DATA - PEDIGREE 
  phenoP2.real$ind  = as.factor(phenoP2.real$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.real)
  cov.mat                  = as.matrix(A.real[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.real$ind         = as.character(phenoP2.real$ind)
  gen.data.pedigree.real   = gen.data
  model.A.real             = remlf90(fixed = Ht ~ 1, generic = gen.data.pedigree.real, data = phenoP2.real, method = 'em')
  h2.A.real                = model.A.real$var["generic_genetic",1]/(model.A.real$var["generic_genetic",1]+model.A.real$var["Residual",1])
  
  # SIMULATED DATA - PEDIGREE  
  phenoP2.simul$ind = as.factor(phenoP2.simul$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.simul)
  cov.mat                  = as.matrix(A.simul[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.simul$ind        = as.character(phenoP2.simul$ind)
  gen.data.pedigree.simul  = gen.data
  model.A.simul            = remlf90(fixed = Ht ~ 1, generic = gen.data.pedigree.simul, data = phenoP2.simul, method = 'em')
  h2.A.simul               = model.A.simul$var["generic_genetic",1]/(model.A.simul$var["generic_genetic",1]+model.A.simul$var["Residual",1])
  
  
print(h2.G.real)  
print(h2.G.simul)
print(h2.A.real)  
print(h2.A.simul)

