library(AlphaSimR)
library(ggplot2)

path = "C:/Users/Victor/Documents/Simulations_maritime_pine/PIPELINE_CENTRAL/0.InputsPipeline"
load(paste(path,"BreedingProgram.RData", sep="/"))

# compare simulated phenotypes with real data for P2
Pheno.P2.simul           = data.frame(popP2@id)
Pheno.P2.simul$pheno     = as.numeric(pheno(popP2))
colnames(Pheno.P2.simul) = c("ind","Ht")
Pheno.P2.simul$type      = "simul"
  
phenoP2.real     = read.table(paste(path,"phenoP2real.txt",sep="/"), header=T)
phenoP2.real     = phenoP2.real[,c("ind","Ht")]
phenoP2.real$type = "real"
  
Pheno.P2 = rbind(Pheno.P2.simul, phenoP2.real)
pPhenoP2 = ggplot(Pheno.P2, aes(x = Ht, fill=type)) + geom_density(alpha=0.2)+xlab("Trait Ht")+theme_bw()+ theme(legend.title = element_blank())
pPhenoP2

ggsave("CheckPhenoP2.tiff", pPhenoP2, width= 20, height= 20, units= "cm", scale = 1, dpi= 600)
