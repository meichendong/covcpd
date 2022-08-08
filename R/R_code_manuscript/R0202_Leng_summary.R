# summarize A0009
source("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/FUNCTIONS.R")
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0011")

lf = list.files()
lf = lf[grep("order2",lf)]
list.ecp = list()
list.inspect = list() 
list.tcai = list()
pnames = NULL

nocp.factorcpt = nocp.inspect = nocp.my = nocp.ecp = 0

for(i in 1:length(lf)){
  dt1 = readRDS(lf[i])
  list.ecp[[i]] = dt1$ecp
  if(length(dt1$ecp) ==0){
    nocp.ecp = nocp.ecp + 1
  }
  tmpmy = dt1$TCai$cps_order[-c(1,2)]
  tmpmy = tmpmy[order(tmpmy)]
  list.tcai[[i]] = tmpmy
  if (length(tmpmy) ==0){
    nocp.my = nocp.my + 1
  }
  list.inspect[[i]] = dt1$inspect$changepoints[,1]
  if(length(dt1$inspect$changepoints[,1]) ==0){
    nocp.inspect = nocp.inspect +1
  }
  pnames = c(pnames, dt1$pathway)
}

pnames
nocp.ecp;nocp.inspect; nocp.factorcpt; nocp.my
list.tcai
# list.factorcpt

# pathways information
library(msigdbr) # this is for human 
# see A0011_cellcycle_cpd.R
# saveRDS(paths.all.genes.filter , "A0011_human_cellcycle_pathways.rds")
paths.all.genes.filter = readRDS("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0011_human_cellcycle_pathways.rds")


# ------------------- visualize cps ------------------------ #
library(SingleCellExperiment)
library(Biobase)
library(ggplot2)
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data")
Rcpp::sourceCpp('C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/FUNCTIONS_CPP_JGNSC.cpp')

eset = readRDS("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/cellcycledatasets_eset.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# feature plots
# marker genes from moussa et al.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene1 = "BUB3"
CCGENES = c("CCNA2","CCNA1","CCNB1","CCNB2","CCNB3","CCND1","CCND2","CCND3","CCNE1","CCNE2")

setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data")
sce = readRDS("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/cellcycledatasets_eset.rds")
ctvar = "cellcycle"
cellorder = c("G1","S","G2/M")
scename = "Leng et al."

ccmarkers = read.csv("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/data/CellCycleGenes_CellDiv2019paper_01202021.csv")

dt <- data.frame(ptime = seq(1:ncol(sce)), ct = sce$cellcycle)
rownames(dt) = colnames(sce)
# dt$paga = sce$paga
moussaorder = read.csv("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/Moussa_CellCycleOrder.csv")
rownames(moussaorder) = moussaorder$Cells
moussaorder$order2 = 248 - moussaorder$order
mo_g1 = moussaorder[1:90,]
mo_sg2 = moussaorder[91:247,]
mo_sg2$order3 = mo_sg2$order2 + 90
mo_g1$order3 = mo_g1$order2 - 157
mo_g1start = rbind(mo_g1, mo_sg2)

mo = mo_g1start[rownames(dt),]
dt$moussaorder = mo$order3


pfeature =ggplot(dt, aes(x = moussaorder ,y=factor(ct, levels = cellorder), 
                         colour = ct)) +
  geom_quasirandom(groupOnX = FALSE) + 
  coord_polar() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(size=10)) +
  xlab("SC1CC order") + ylab("Subgroups") +
  geom_vline(xintercept = 0, color = "red")+
  geom_vline(xintercept = c(list.tcai[[3]]), color = "blue", linetype="longdash")+
  xlim(c(0,250))

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0011_cellcycle_summary_g1start_moussaorder", Sys.Date(),".pdf", sep = ""), height =6, width = 6)
print(pfeature)
dev.off()


pline =ggplot(dt, aes(x = moussaorder ,y=factor(ct, levels = cellorder), 
                      colour = ct)) + geom_jitter() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(size=10)) +
  xlab("SC1CC order") + ylab("Subgroups") +
  geom_vline(xintercept = 0, color = "red")+
  geom_vline(xintercept = c(list.tcai[[3]]), color = "blue", linetype="longdash")+
  xlim(c(0,250))

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0011_cellcycle_summary_g1start_moussaorder", Sys.Date(),".pdf", sep = ""), height =4, width = 6)
print(pline)
dev.off()
