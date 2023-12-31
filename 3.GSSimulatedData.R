library(AlphaSimR)
library(AGHmatrix)
library(breedR)

path = "C:/Users/Victor/Documents/Simulations_maritime_pine/PIPELINE_CENTRAL/0.InputsPipeline"
load(paste(path,"BreedingProgram.RData", sep="/"))


###################################################
## GS MODELS WITH ONE SIMULATED DATA
###################################################

phenoP2.real     = read.table(paste(path, "phenoP2real.txt",sep="/"), header=T)

## SIMULATED DATA
phenoP2.simul              = data.frame(pheno(popP2))
gvP2.simul                 = data.frame(gv(popP2))
genoP2.simul               = data.frame(pullSegSiteGeno(popP2)[,id.SNP.array])
row.names(genoP2.simul)    = popP2@id
phenoP2.simul$ind          = popP2@id
gvP2.simul$ind             = popP2@id
pedigreeP2.simul           = rbind(data.frame(ind = c(popP2@id), mom = c(popP2@mother), dad = c(popP2@father)), data.frame(ind = c(popG1@id), mom = c(popG1@mother), dad = c(popG1@father)), data.frame(ind = c(popG0@id), mom = c(popG0@mother), dad = c(popG0@father)) )
pedigreeP2.simul           = pedigreeP2.simul[!duplicated(pedigreeP2.simul$ind),]

## RELATIONSHIP MATRICES
G.simul  = Gmatrix(as.matrix(genoP2.simul), method='VanRaden', missingValue=-9, maf=0)
A.simul  = Amatrix(pedigreeP2.simul, ploidy=2)
A.simul  = A.simul[phenoP2.simul$ind,phenoP2.simul$ind]

colnames(G.simul) = paste("ind",colnames(G.simul), sep="");colnames(A.simul) = paste("ind",colnames(A.simul), sep="");row.names(G.simul) = paste("ind",row.names(G.simul), sep="");row.names(A.simul) = paste("ind",row.names(A.simul), sep="")

pedigreeP2.simul$cross  = paste(pedigreeP2.simul$mom,pedigreeP2.simul$dad, sep="x")
pedigreeP2.simul$family = phenoP2.real$family[match(pedigreeP2.simul$cross, phenoP2.real$cross)]
phenoP2.simul$family    = pedigreeP2.simul$family[match(phenoP2.simul$ind, pedigreeP2.simul$ind)]

# Estimates h2 for accuracy

  # SIMULATED DATA - GENOMIC 
  phenoP2.simul$ind = as.factor(phenoP2.simul$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.simul)
  cov.mat                  = as.matrix(G.simul[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.simul$ind        = as.character(phenoP2.simul$ind)
  gen.data.genomic.simul   = gen.data
  model.G.simul            = remlf90(fixed = Ht ~ 1, generic = gen.data.genomic.simul, data = phenoP2.simul, method = 'em')
  h2.G.simul               = model.G.simul$var["generic_genetic",1]/(model.G.simul$var["generic_genetic",1]+model.G.simul$var["Residual",1])

  # SIMULATED DATA - PEDIGREE  
  phenoP2.simul$ind = as.factor(phenoP2.simul$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.simul)
  cov.mat                  = as.matrix(A.simul[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.simul$ind        = as.character(phenoP2.simul$ind)
  gen.data.pedigree.simul  = gen.data
  model.A.simul            = remlf90(fixed = Ht ~ 1, generic = gen.data.pedigree.simul, data = phenoP2.simul, method = 'em')
  h2.A.simul               = model.A.simul$var["generic_genetic",1]/(model.A.simul$var["generic_genetic",1]+model.A.simul$var["Residual",1])
  
print(h2.G.simul)
print(h2.A.simul)
  


## SCENARIO 1
###################

data          = c("simul")
type          = c("pedigree","genomic")
pourcent      = c(0.2,0.5)
nbIt          = 100
families      = as.character(unique(phenoP2.real$family))
RESULTATS     = data.frame(DATA=NA, TYPE=NA, POURCENT=NA, IT=NA, COR=NA, COR2=NA)
RESULTATS.FAM = data.frame(DATA=NA, TYPE=NA, POURCENT=NA, IT=NA, FAM=NA, COR=NA, COR2=NA, NB=NA)
a             = 1
b             = 1

time = Sys.time()
for(d in data){print(d)
  
  if(d=="real"){pheno.ref = phenoP2.real}else{pheno.ref = phenoP2.simul}
  
  for(t in type){print(t)
    
    if(t=="pedigree" & d=="real"){  GEN.MATRIX = A.real }
    if(t=="pedigree" & d=="simul"){ GEN.MATRIX = A.simul}
    if(t=="genomic"  & d=="real"){  GEN.MATRIX = G.real }
    if(t=="genomic"  & d=="simul"){ GEN.MATRIX = G.simul }
    
    for(p in pourcent){
      
      for(i in 1:nbIt){print(i)
        
        pheno.tmp      = pheno.ref
        individus.hide = NULL
        
        ## TRAINING SET
        ##############
        for(f in families){
          pheno.tmp.fam = pheno.tmp[pheno.tmp$family==f,]
          individus.hide.fam  = sample(as.character(unique(pheno.tmp.fam$ind)),round(nrow(pheno.tmp.fam)*p))
          individus.hide  = c(individus.hide, individus.hide.fam)
        }
        
        id.hide = as.numeric(which(pheno.tmp$ind%in%individus.hide))
        pheno.tmp[id.hide,"Ht"]  = NA
        
        ## GS MODEL
        ###########
        pheno.ref$ind = as.factor(pheno.ref$ind)
        inc.mat = model.matrix(~ 0 + ind, pheno.ref)
        cov.mat = as.matrix(GEN.MATRIX[colnames(inc.mat),colnames(inc.mat)])
        model.SG = remlf90(fixed = Ht ~ 1, generic = list(genetic=list(inc.mat, cov.mat)), data = pheno.tmp, method = 'em')
        
        ## RESULTS
        ##########
        PBV.full = ranef(model.SG)$generic
        PBV = model.matrix(model.SG)$generic %*% PBV.full
        cor1 = cor(pheno.ref[id.hide,"Ht"], PBV[id.hide])
        cor2 = cor(gvP2.simul[id.hide,"Ht"], PBV[id.hide])
        RESULTATS[a,c("DATA","TYPE","POURCENT","IT","COR","COR2")] = c(d, t, p, i, cor1, cor2)
        a = a + 1
        
        for(j in families){
          familles.hide = as.character(pheno.tmp$family[match(individus.hide, pheno.tmp$ind)])
          idx = which(familles.hide==j)
          idx = which(pheno.tmp$ind%in%individus.hide[idx])
          cor1.fam = cor(pheno.ref[idx,"Ht"], PBV[idx])
          cor2.fam = cor(gvP2.simul[idx,"Ht"], PBV[idx])
          RESULTATS.FAM[b,c("DATA","TYPE","POURCENT","FAM","IT","COR","COR2","NB")] = c(d, t, p, j, i, cor1.fam, cor2.fam, length(idx))
          b = b + 1
        }
        
      }
    }
  }
}
Sys.time()-time

## RESULTS
###########

for(i in 1:nrow(RESULTATS)){
  if(RESULTATS[i,"DATA"]=="real"&RESULTATS[i,"TYPE"]=="pedigree"){RESULTATS[i,"h2"] = h2.A.real}
  if(RESULTATS[i,"DATA"]=="real"&RESULTATS[i,"TYPE"]=="genomic"){RESULTATS[i,"h2"] = h2.G.real}
  if(RESULTATS[i,"DATA"]=="simul"&RESULTATS[i,"TYPE"]=="pedigree"){RESULTATS[i,"h2"] = h2.A.simul}
  if(RESULTATS[i,"DATA"]=="simul"&RESULTATS[i,"TYPE"]=="genomic"){RESULTATS[i,"h2"] = h2.G.simul}
}

RESULTATS$COR      = as.numeric(RESULTATS$COR)
RESULTATS$COR2     = as.numeric(RESULTATS$COR2)
RESULTATS$ACCURACY = RESULTATS$COR/sqrt(RESULTATS$h2)
RESULTATS$ACCURACY2= RESULTATS$COR2/sqrt(RESULTATS$h2)
RESULTATS.FAM$COR  = as.numeric(RESULTATS.FAM$COR)
RESULTATS.FAM$COR2 = as.numeric(RESULTATS.FAM$COR2)

write.table(RESULTATS, "RESULTATS.SIMUL.txt", row.names=F)
write.table(RESULTATS.FAM, "RESULTATS.SIMUL.FAM.txt", row.names=F)
save.image("Save.RData")
