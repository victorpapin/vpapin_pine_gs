library(AlphaSimR)
library(ggplot2)

path = "C:/Users/Victor/Documents/Simulations_maritime_pine/PIPELINE_CENTRAL/0.InputsPipeline"
load(paste(path,"founderBIG.RData", sep="/"))


##############################################
## DEFINE A POPULATION P0
##############################################

SP           = SimParam$new(founderGenomes)
popP0        = newPop(founderGenomes)

##############################################
## CHECK LD AND AF FOR SIMULATED DATA
##############################################

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


# CHECK FOR LD AND AF
genoP0.simul       = pullSegSiteGeno(popP0)
genoP0.simul.array = genoP0.simul[,id.SNP.array]
  
  
# LD and AF for simulated SNP array data
LD.R2.simul.array     = cor(genoP0.simul.array)
LD.R2.simul.array.val = LD.R2.simul.array[upper.tri(LD.R2.simul.array, diag=FALSE)]
AllFreq.simul.array   = colMeans(genoP0.simul.array)/2
  
# LD and AF for real snp array data
genoP0.real     = read.table(paste(path,"geno_P0.csv",sep="_"), header=TRUE)
LD.R2.real      = cor(genoP0.real)
LD.R2.real.val  = LD.R2.real[upper.tri(LD.R2.real, diag=FALSE)]
AllFreq.real    = colMeans(genoP0.real)/2
  
# FINAL PLOTS
pLD  = ggplot(data.frame(LD.R2.simul.array.val), aes(x=LD.R2.simul.array.val))+ 
  geom_density(alpha=.2, fill="red") +theme_bw()+xlab("LD with SIMULATED (red) and REAL (blue) SNP ARRAY DATA")+
  geom_density(data=data.frame(LD.R2.real.val), aes(x=LD.R2.real.val), fill="blue", alpha=0.2)
pLD
pLD + ylim(0,0.1)
pLD + ylim(0,0.01)
  
pAF  = ggplot(data.frame(AllFreq.simul.array), aes(x=AllFreq.simul.array))+ 
  geom_density(alpha=.2, fill="red") +theme_bw()+xlab("Allele frenquencies with SIMULATED (red) and REAL (blue) SNP ARRAY DATA")+
  geom_density(data=data.frame(AllFreq.real), aes(x=AllFreq.real), fill="blue", alpha=0.2)
pAF


ggsave("CheckLD.tiff", pLD, width= 20, height= 20, units= "cm", scale = 1, dpi= 600)
ggsave("CheckLDzoom1.tiff", pLD+ylim(0,0.1), width= 20, height= 20, units= "cm", scale = 1, dpi= 600)
ggsave("CheckLDzoom2.tiff", pLD+ylim(0,0.01), width= 20, height= 20, units= "cm", scale = 1, dpi= 600)
ggsave("CheckAF.tiff", pAF, width= 20, height= 20, units= "cm", scale = 1, dpi= 600)


