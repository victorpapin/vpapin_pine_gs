library(AGHmatrix)
library(breedR)

path = "C:/Users/Victor/Documents/Simulations_maritime_pine/PIPELINE_CENTRAL/0.InputsPipeline"

###################################################
## GS MODELS WITH REAL DATA = REFERENCE
###################################################

## REAL DATA
phenoP2.real     = read.table(paste(path,"phenoP2real.txt", sep="/"), header=T)
genoP2.real      = read.table(paste(path, "genoP2real.txt", sep="/"),header=T)
pedigreeP2.real  = read.table(paste(path,"pedigreeP2real.txt", sep="/"), header=TRUE)

## RELATIONSHIP MATRICES
G.real   = Gmatrix(as.matrix(genoP2.real), method='VanRaden', missingValue=-9, maf=0)
A.real   = Amatrix(pedigreeP2.real, ploidy=2)
A.real   = A.real[phenoP2.real$ind,phenoP2.real$ind]

colnames(G.real) = paste("ind",colnames(G.real), sep="");colnames(A.real) = paste("ind",colnames(A.real), sep="");row.names(G.real) = paste("ind",row.names(G.real), sep="");row.names(A.real) = paste("ind",row.names(A.real), sep="")

# Estimates h2 for accuracy
if(T){
  
  # REAL DATA - GENOMIC
  phenoP2.real$ind         = as.factor(phenoP2.real$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.real)
  cov.mat                  = as.matrix(G.real[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.real$ind         = as.character(phenoP2.real$ind)
  gen.data.genomic.real    = gen.data
  model.G.real             = remlf90(fixed = Ht ~ 1, generic = gen.data.genomic.real, data = phenoP2.real, method = 'em')
  h2.G.real                = model.G.real$var["generic_genetic",1]/(model.G.real$var["generic_genetic",1]+model.G.real$var["Residual",1])
  
  # REAL DATA - PEDIGREE 
  phenoP2.real$ind  = as.factor(phenoP2.real$ind)
  inc.mat                  = model.matrix(~ 0 + ind, phenoP2.real)
  cov.mat                  = as.matrix(A.real[colnames(inc.mat),colnames(inc.mat)])
  gen.data                 = list(genetic=list(inc.mat, cov.mat))
  phenoP2.real$ind         = as.character(phenoP2.real$ind)
  gen.data.pedigree.real   = gen.data
  model.A.real             = remlf90(fixed = Ht ~ 1, generic = gen.data.pedigree.real, data = phenoP2.real, method = 'em')
  h2.A.real                = model.A.real$var["generic_genetic",1]/(model.A.real$var["generic_genetic",1]+model.A.real$var["Residual",1])

  print(h2.G.real)  
  print(h2.A.real)  
  
}



## SCENARIO 1
###################

data          = c("real")
type          = c("pedigree","genomic")
pourcent      = c(0.2,0.5)
nbIt          = 100
families      = as.character(unique(phenoP2.real$family))
RESULTATS     = data.frame(DATA=NA, TYPE=NA, POURCENT=NA, IT=NA, COR=NA)
RESULTATS.FAM = data.frame(DATA=NA, TYPE=NA, POURCENT=NA, IT=NA, FAM=NA, COR=NA, NB=NA)
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
        RESULTATS[a,c("DATA","TYPE","POURCENT","IT","COR")] = c(d, t, p, i, cor1)
        a = a + 1
        
        for(j in families){
          familles.hide = as.character(pheno.tmp$family[match(individus.hide, pheno.tmp$ind)])
          idx = which(familles.hide==j)
          idx = which(pheno.tmp$ind%in%individus.hide[idx])
          cor1.fam = cor(pheno.ref[idx,"Ht"], PBV[idx])
          RESULTATS.FAM[b,c("DATA","TYPE","POURCENT","FAM","IT","COR","NB")] = c(d, t, p, j, i, cor1.fam, length(idx))
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
RESULTATS$ACCURACY = RESULTATS$COR/sqrt(RESULTATS$h2)
RESULTATS.FAM$COR  = as.numeric(RESULTATS.FAM$COR)

write.table(RESULTATS, "RESULTATS.REAL.txt", row.names=F)
write.table(RESULTATS.FAM, "RESULTATS.REAL.FAM.txt", row.names=F)
save.image("Save.RData")