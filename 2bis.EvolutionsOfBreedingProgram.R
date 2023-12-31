library(AlphaSimR)
library(ggplot2)
library(tidyverse)
library(ggridges)

path = "C:/Users/Victor/Documents/Simulations_maritime_pine/PIPELINE_CENTRAL/0.InputsPipeline"
load(paste(path,"BreedingProgram.RData", sep="/"))


# evolutions across generations

values.gv       = c(gv(popP0), gv(popG0),  gv(popP1), gv(popG1), gv(popP2))
values.bv       = c(bv(popP0), bv(popG0),  bv(popP1), bv(popG1), bv(popP2))
values.pheno    = c(pheno(popP0), pheno(popG0),  pheno(popP1), pheno(popG1), pheno(popP2))
generation      = c(rep("P0",length(gv(popP0))),rep("G0",length(gv(popG0))), rep("P1",length(gv(popP1))),rep("G1",length(gv(popG1))),rep("P2",length(gv(popP2)))  )
data            = data.frame(matrix(ncol=3, nrow=length(values.gv)*3))
colnames(data)  = c("value","generation","variable")
data$value      = c(values.gv, values.bv, values.pheno)
data$generation = as.factor(rep(generation, 3))
data$variable   = c(rep(c("GV","BV","PHENO"), each=length(values.gv)))
  
plot =data %>%
    mutate(generation = fct_relevel(generation, 
                                    "P2", "G1", "P1", 
                                    "G0", "P0")) %>%
    ggplot(aes(x = value, y = generation)) +
    facet_grid(.~variable, scales = "free")+
    geom_density_ridges(aes(fill = generation), alpha = .4, rel_min_height = 0.01,quantile_lines=TRUE, quantile_fun=function(x,...)mean(x)) +
    ylab("Generations")+xlab("values")+theme(legend.position="none")

plot

ggsave("EvolutionsBP.tiff", plot, width= 40, height= 20, units= "cm", scale = 1, dpi= 600)
