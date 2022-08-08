# from longleaf/A0009_v2

library(SingleCellExperiment)
library(ggbeeswarm)
library(ggthemes)
library(ggplot2)
# ============================
#  deng 2014
### DATA
library(R.utils) 
library(Seurat)
library(JGL)
library(pheatmap)
library(HDtest)
library(qusage)
library(huge)

deng <- readRDS("/pine/scr/m/e/meichen/hovestdta/data/H0605_dengSCE.rds")
# setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data")
# deng <- readRDS("H0605_dengSCE.rds")

library(gskb)
data("mm_pathway") 
mm_pathway_clean <- lapply(mm_pathway, function(x){
  y = x[-c(1,2)]
}) 

path.keywords <- c("Fgf","Hedgehog","Wnt","Notch","TGF","JAK_STAT","TYROSINE_KINASE","RHO","GPCR","ESTROGEN_RECEPTOR")

npath.select <- vector(length = length(path.keywords))
paths.all <- NULL
select.pathways2 <- NULL
for (i in 1:length(path.keywords)){
  select.pathways <- names(mm_pathway)[grep(path.keywords[i], names(mm_pathway), ignore.case = T)]
  select.pathways2 <- c(select.pathways2, select.pathways) 
}
select.pathways2

paths.all.genes <- mm_pathway_clean[select.pathways2]
# keep pathways with 20~100 genes
min.ngene <- 50
max.ngene <- 1000#150
keeppath <- lapply(paths.all.genes, function(x){
  length(x) >min.ngene & length(x) <max.ngene
})
paths.all.genes.filter <- paths.all.genes[unlist(keeppath)]
allgenes <- unique(unlist(paths.all.genes.filter))


paths.all.genes.filter$allgenes = rownames(assay(deng))
path0 <-  paths.all.genes.filter[[ipath]]

dengraw <- assay(deng)[toupper(rownames(deng)) %in% path0,]
dengraw <- dengraw[rowSums(dengraw > 0) > 3,] # filter out lowly expressed genes
y = dengraw

# ecp
library(ecp)
res.ecp1 = e.divisive(t(y), min.size = 20)


nsubgroup = length(res.ecp1$estimates) - 1
mcp = res.ecp1$estimates
colnorm_y = NULL
for(i in 1:nsubgroup){
  sgroupi = t(y)[mcp[i]:(mcp[i+1]-1),]
  colnorm_sgroupi = apply(sgroupi,2, function(x){
    (x-mean(x))
  })
  colnorm_y = rbind(colnorm_y, colnorm_sgroupi)
}

# covcpd
my2 = covcpd(X = t(colnorm_y), k=20, maxseg = 7, search_by = 2, nperm = 1000)


resall = list(ecp = res.ecp1$estimates[-c(1, length(res.ecp1$estimates))],
              TCai = my2[-c(1, length(my2))],
              pathway = names(paths.all.genes.filter)[ipath])

saveRDS(resall, "A0009_deng_pathway_ipath.rds")