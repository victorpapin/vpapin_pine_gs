
library(AlphaSimR)

##############################################
## GENERATE FOUNDER GENOMES
##############################################

# INITIAL PARAMETERS
nInd       = 600     # Reference breeding program
nChr       = 12      # Chagne et al., 2002 
nSNPPerChr = 686     # 8235 SNPs 4TREE array (/12Chr)
nQTLPerChr = 8      # Grattapaglia et al., 2012, total of around 100 QTLs
Ne         = 100     # Reference breeding program

## GENERATE FOUNDER GENOMES
founderGenomes = runMacs2(nInd = nInd,      
                          nChr = nChr,       
                          segSites =  10*nSNPPerChr+nQTLPerChr, 
                          Ne = Ne,
                          bp = 2.14e+09,     # Chagne et al., 2002 
                          genLen = 1.2,      # Chancerel et al., 2013  
                          mutRate = 4.18e-8, # Reference from Santi: Jaramillo-Correa et al., 2020
                          histNe = c(50000), # Gonzalez-Martinez (not published)
                          histGen = c(1),    # Gonzalez-Martinez (not published)
                          inbred = FALSE,  
                          split = NULL,
                          ploidy = 2L,     
                          returnCommand = FALSE,
                          nThreads = NULL
)

save.image("founderBIG.RData")
