#From A0011
# ----------------- CELL CYCLE DATASETS ------------------------- #
## setwd("/pine/scr/m/e/meichen/hovestdta/data/A0011")
## source("/pine/scr/m/e/meichen/hovestdta/scripts/FUNCTIONS.R")
## Rcpp::sourceCpp('/pine/scr/m/e/meichen/hovestdta/scripts/HDtest_CPP.cpp')
## Rcpp::sourceCpp('/pine/scr/m/e/meichen/hovestdta/scripts/FUNCTIONS_CPP_JGNSC.cpp')

library(SingleCellExperiment)
library(ggbeeswarm)
library(ggthemes)
library(ggplot2)
library(R.utils) 
library(Seurat)
library(JGL)
library(pheatmap)
library(HDtest)
library(qusage)
library(huge)

# =========================== load data

eset = readRDS("/pine/scr/m/e/meichen/hovestdta/data/cellcycledatasets_eset.rds")
seuratmarkers = readRDS("/pine/scr/m/e/meichen/hovestdta/data/Cellcycledatasets_seurat_clustermarkers.rds")

moussaorder = read.csv("/pine/scr/m/e/meichen/hovestdta/data/Moussa_CellCycleOrder.csv")
rownames(moussaorder) = moussaorder$Cells
moussaorder$order2 = 248 - moussaorder$order
mo_g1 = moussaorder[1:90,]
mo_sg2 = moussaorder[91:247,]
mo_sg2$order3 = mo_sg2$order2 + 90
mo_g1$order3 = mo_g1$order2 - 157
mo_g1start = rbind(mo_g1, mo_sg2)

mo = mo_g1start[colnames(eset),]
eset$moussaorder = mo$order3



# =========================== find other pathways related to cell cycles
library(msigdbr) # this is for human 
paths.all.genes.filter = readRDS("/pine/scr/m/e/meichen/hovestdta/data/A0011_human_cellcycle_pathways.rds")
path0 <-  paths.all.genes.filter[[ipath]]


# ========================================================
library(ecp)
countmat <- exprs(eset)[toupper(rownames(eset)) %in% path0,]
countmat <- countmat[rowSums(countmat > 0) > 3,] # filter out lowly expressed genes


# start from G2/M
if(estorder == 1){ 
  o <- order(eset$paga)
} else if(estorder ==2) { 
  o <- order(eset$moussaorder)
} else { 
  o <- order(eset$paga_g1)
}

countmat = countmat[,o]
y = countmat

# ecp
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

check.ecp = e.divisive(colnorm_y, alpha = 1)

my2 = covcpd(X = colnorm_y, k=20, maxseg = 7, search_by = 2, nperm = 1000)


resall = list(ecp = res.ecp1$estimates[-c(1, length(res.ecp1$estimates))],
              TCai = my2[-c(1, length(my2))],
              pathway = names(paths.all.genes.filter)[ipath])

saveRDS(resall, "A0011_cellcycle_orderestorder_pathway_ipath.rds")
