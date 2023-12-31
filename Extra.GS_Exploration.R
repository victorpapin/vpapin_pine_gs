library(breedR)
library(AlphaSimR)

for(index in 1:4800){
  
  if(index %in% c(001:100)  ){ h2=0.13 ;  NSNP = 1 ; NTSET = 400  }
  if(index %in% c(101:200)  ){ h2=0.13 ;  NSNP = 1 ; NTSET = 650  }
  if(index %in% c(201:300)  ){ h2=0.13 ;  NSNP = 1 ; NTSET = 1700 }
  if(index %in% c(301:400)  ){ h2=0.13 ;  NSNP = 1 ; NTSET = 2800 }
  if(index %in% c(401:500)  ){ h2=0.13 ;  NSNP = 2 ; NTSET = 400  }
  if(index %in% c(501:600)  ){ h2=0.13 ;  NSNP = 2 ; NTSET = 650  }
  if(index %in% c(601:700)  ){ h2=0.13 ;  NSNP = 2 ; NTSET = 1700 }
  if(index %in% c(701:800)  ){ h2=0.13 ;  NSNP = 2 ; NTSET = 2800 }
  if(index %in% c(801:900)  ){ h2=0.13 ;  NSNP = 3 ; NTSET = 400  }
  if(index %in% c(901:1000) ){ h2=0.13 ;  NSNP = 3 ; NTSET = 650  }
  if(index %in% c(1001:1100)){ h2=0.13 ;  NSNP = 3 ; NTSET = 1700 }
  if(index %in% c(1101:1200)){ h2=0.13 ;  NSNP = 3 ; NTSET = 2800 }
  if(index %in% c(1201:1300)){ h2=0.13 ;  NSNP = 4 ; NTSET = 400  }
  if(index %in% c(1301:1400)){ h2=0.13 ;  NSNP = 4 ; NTSET = 650  }
  if(index %in% c(1401:1500)){ h2=0.13 ;  NSNP = 4 ; NTSET = 1700 }
  if(index %in% c(1501:1600)){ h2=0.13 ;  NSNP = 4 ; NTSET = 2800 }
  
  if(index %in% c(1601:1700)){ h2=0.33 ;  NSNP = 1 ; NTSET = 400  }
  if(index %in% c(1701:1800)){ h2=0.33 ;  NSNP = 1 ; NTSET = 650  }
  if(index %in% c(1801:1900)){ h2=0.33 ;  NSNP = 1 ; NTSET = 1700 }
  if(index %in% c(1901:2000)){ h2=0.33 ;  NSNP = 1 ; NTSET = 2800 }
  if(index %in% c(2001:2100)){ h2=0.33 ;  NSNP = 2 ; NTSET = 400  }
  if(index %in% c(2101:2200)){ h2=0.33 ;  NSNP = 2 ; NTSET = 650  }
  if(index %in% c(2201:2300)){ h2=0.33 ;  NSNP = 2 ; NTSET = 1700 }
  if(index %in% c(2301:2400)){ h2=0.33 ;  NSNP = 2 ; NTSET = 2800 }
  if(index %in% c(2401:2500)){ h2=0.33 ;  NSNP = 3 ; NTSET = 400  }
  if(index %in% c(2501:2600)){ h2=0.33 ;  NSNP = 3 ; NTSET = 650  }
  if(index %in% c(2601:2700)){ h2=0.33 ;  NSNP = 3 ; NTSET = 1700 }
  if(index %in% c(2701:2800)){ h2=0.33 ;  NSNP = 3 ; NTSET = 2800 }
  if(index %in% c(2801:2900)){ h2=0.33 ;  NSNP = 4 ; NTSET = 400  }
  if(index %in% c(2901:3000)){ h2=0.33 ;  NSNP = 4 ; NTSET = 650  }
  if(index %in% c(3001:3100)){ h2=0.33 ;  NSNP = 4 ; NTSET = 1700 }
  if(index %in% c(3101:3200)){ h2=0.33 ;  NSNP = 4 ; NTSET = 2800 }
  
  if(index %in% c(3201:3300)){ h2=0.50 ;  NSNP = 1 ; NTSET = 400  }
  if(index %in% c(3301:3400)){ h2=0.50 ;  NSNP = 1 ; NTSET = 650  }
  if(index %in% c(3401:3500)){ h2=0.50 ;  NSNP = 1 ; NTSET = 1700 }
  if(index %in% c(3501:3600)){ h2=0.50 ;  NSNP = 1 ; NTSET = 2800 }
  if(index %in% c(3601:3700)){ h2=0.50 ;  NSNP = 2 ; NTSET = 400  }
  if(index %in% c(3701:3800)){ h2=0.50 ;  NSNP = 2 ; NTSET = 650  }
  if(index %in% c(3801:3900)){ h2=0.50 ;  NSNP = 2 ; NTSET = 1700 }
  if(index %in% c(3901:4000)){ h2=0.50 ;  NSNP = 2 ; NTSET = 2800 }
  if(index %in% c(4001:4100)){ h2=0.50 ;  NSNP = 3 ; NTSET = 400  }
  if(index %in% c(4101:4200)){ h2=0.50 ;  NSNP = 3 ; NTSET = 650  }
  if(index %in% c(4201:4300)){ h2=0.50 ;  NSNP = 3 ; NTSET = 1700 }
  if(index %in% c(4301:4400)){ h2=0.50 ;  NSNP = 3 ; NTSET = 2800 }
  if(index %in% c(4401:4500)){ h2=0.50 ;  NSNP = 4 ; NTSET = 400  }
  if(index %in% c(4501:4600)){ h2=0.50 ;  NSNP = 4 ; NTSET = 650  }
  if(index %in% c(4601:4700)){ h2=0.50 ;  NSNP = 4 ; NTSET = 1700 }
  if(index %in% c(4701:4800)){ h2=0.50 ;  NSNP = 4 ; NTSET = 2800 }
  
  
  if(h2==0.13){load("BreedingProgramG2_0.13.RData")}
  if(h2==0.33){load("BreedingProgramG2_0.33.RData")}
  if(h2==0.50){load("BreedingProgramG2_0.50.RData")}
  
  rm(.Random.seed)
  
  if(NSNP==1){G.simul = read.table("G.simul.1.txt",header=T)}
  if(NSNP==2){G.simul = read.table("G.simul.2.txt",header=T)}
  if(NSNP==3){G.simul = read.table("G.simul.3.txt",header=T)}
  if(NSNP==4){G.simul = read.table("A.simul.txt",header=T)}
  
  ## PREPARE
  ####################################################
  
  phenoP2.simul              = data.frame(pheno(popP2))
  gvP2.simul                 = data.frame(gv(popP2))
  phenoP2.simul$ind          = popP2@id
  gvP2.simul$ind             = popP2@id
  pedigreeP2.simul           = rbind(data.frame(ind = c(popP2@id), mom = c(popP2@mother), dad = c(popP2@father)), data.frame(ind = c(popG1@id), mom = c(popG1@mother), dad = c(popG1@father)), data.frame(ind = c(popG0@id), mom = c(popG0@mother), dad = c(popG0@father)) )
  pedigreeP2.simul           = pedigreeP2.simul[!duplicated(pedigreeP2.simul$ind),]
  phenoP2.real               = read.table("phenoP2real.txt", header=T)
  pedigreeP2.simul$cross     = paste(pedigreeP2.simul$mom,pedigreeP2.simul$dad, sep="x")
  pedigreeP2.simul$family    = phenoP2.real$family[match(pedigreeP2.simul$cross, phenoP2.real$cross)]
  phenoP2.simul$family       = pedigreeP2.simul$family[match(phenoP2.simul$ind, pedigreeP2.simul$ind)]
  
  
  ## LOOP FOR GS MODELS
  ###################################################
  
  NVSET         = 1200
  NSUPP         = 4200-NTSET-NVSET
  families      = as.character(unique(phenoP2.real$family))
  RESULTATS     = data.frame(H2=NA,NTSET=NA,NSNP=NA,COR1=NA,COR2=NA)
  RESULTATS.FAM = data.frame(H2=NA,NTSET=NA,NSNP=NA,FAM=NA,COR1.FAM=NA,COR2.FAM=NA,NB=NA)
  
  pheno.ref  = phenoP2.simul
  GEN.MATRIX = G.simul
  pheno.tmp  = pheno.ref
  
  ## TRAINING SET
  ###############
  
  pourcent.Tset  = NTSET/4200
  individus.Tset = NULL
  
  for(f in families){
    pheno.tmp.fam = pheno.tmp[pheno.tmp$family==f,]
    individus.Tset.fam = sample(as.character(unique(pheno.tmp.fam$ind)),round(nrow(pheno.tmp.fam)*pourcent.Tset))
    individus.Tset  = c(individus.Tset, individus.Tset.fam)
  }
  
  pheno.tmp = pheno.tmp[which(!pheno.tmp$ind %in% individus.Tset),]      
  
  ## VALIDATION SET
  #################
  
  pourcent.Vset  = 1200/nrow(pheno.tmp)
  individus.Vset = NULL
  
  for(f in families){
    pheno.tmp.fam = pheno.tmp[pheno.tmp$family==f,]
    individus.Vset.fam = sample(as.character(unique(pheno.tmp.fam$ind)),round(nrow(pheno.tmp.fam)*pourcent.Vset))
    individus.Vset  = c(individus.Vset, individus.Vset.fam)
  }
  
  pheno.tmp = pheno.tmp[which(!pheno.tmp$ind %in% individus.Vset),]      
  
  ## SUPPLEMENT
  #################
  
  individus.supp = as.character(pheno.tmp$ind)
  
  ## GS MODEL
  ###########
  
  print(h2) ; print(NTSET) ; print(NSNP)
  
  pheno.tmp = pheno.ref
  id.Vset   = as.numeric(which(pheno.tmp$ind%in%c(individus.Vset)))
  id.hide   = as.numeric(which(pheno.tmp$ind%in%c(individus.Vset, individus.supp)))
  pheno.tmp[id.hide,"Ht"]  = NA
  
  pheno.ref$ind = as.factor(pheno.ref$ind)
  inc.mat       = model.matrix(~ 0 + ind, pheno.ref)
  cov.mat       = as.matrix(GEN.MATRIX[colnames(inc.mat),colnames(inc.mat)])
  time = Sys.time()
  model.SG      = remlf90(fixed = Ht ~ 1, generic = list(genetic=list(inc.mat, cov.mat)), data = pheno.tmp, method = 'em')
  print(paste("time",Sys.time()-time, sep=":"))    
  
  ## RESULTS
  ##########
  
  PBV.full = ranef(model.SG)$generic
  PBV      = model.matrix(model.SG)$generic %*% PBV.full
  cor1     = cor(pheno.ref[id.Vset,"Ht"] , PBV[id.Vset])
  cor2     = cor(gvP2.simul[id.Vset,"Ht"], PBV[id.Vset])
  RESULTATS[1,c("H2","NTSET","NSNP","COR1","COR2")] = as.character(c(h2, NTSET, NSNP, cor1, cor2))
  print(paste("cor1",cor1,sep=":"))
  
  b=1
  for(j in families){
    familles.hide = as.character(pheno.tmp$family[match(individus.Vset, pheno.tmp$ind)])
    idx = which(familles.hide==j)
    idx = which(pheno.tmp$ind%in%individus.Vset[idx])
    cor1.fam = cor(pheno.ref[idx,"Ht"], PBV[idx])
    cor2.fam = cor(gvP2.simul[idx,"Ht"], PBV[idx])
    RESULTATS.FAM[b,c("H2","NTSET","NSNP","FAM","COR1.FAM","COR2.FAM","NB")] = as.character(c(h2, NTSET, NSNP, j,cor1.fam, cor2.fam, length(idx)))
    b = b + 1
    print(paste(j,cor1.fam, sep=":"))
  }
  
  
  if(file.exists(paste("ACC.GBL-",h2,"-",NSNP,"-",NTSET,".txt",sep=""))){
    
    ligne <- as.character(c(h2, NTSET, NSNP, cor1, cor2))
    chemin_fichier <- paste(getwd(),paste("ACC.GBL-",h2,"-",NSNP,"-",NTSET,".txt",sep=""), sep="/")
    write.table(as.data.frame(t(ligne)), file = chemin_fichier, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    
  }else{
    write.table(RESULTATS, paste("ACC.GBL-",h2,"-",NSNP,"-",NTSET,".txt",sep=""), row.names=F, sep=",")
  }
  
  if(file.exists(paste("ACC.FAM-",h2,"-",NSNP,"-",NTSET,".txt",sep=""))){
    
    write.table(RESULTATS.FAM, 
                file = paste(getwd(),paste("ACC.FAM-",h2,"-",NSNP,"-",NTSET,".txt",sep=""), sep="/"),
                append = TRUE, 
                sep=",",
                col.names = FALSE, 
                row.names = FALSE)
    
  }else{
    write.table(RESULTATS.FAM, paste("ACC.FAM-",h2,"-",NSNP,"-",NTSET,".txt",sep=""), row.names=F, sep=",")
  }
  
  
  
}

