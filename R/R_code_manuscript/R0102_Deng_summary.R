# summarize A0009
source("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/FUNCTIONS.R")
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0009")
lf = list.files()
list.ecp = list()
list.inspect = list()
list.factorcpt = list()
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
  
  pnames = c(pnames, dt1$pathway)
}

# pathways information
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


paths.all.genes <- mm_pathway_clean[select.pathways2]
# keep pathways with 20~100 genes
min.ngene <- 20
max.ngene <- 1000#150
keeppath <- lapply(paths.all.genes, function(x){
  length(x) >min.ngene & length(x) <max.ngene
})
paths.all.genes.filter <- paths.all.genes[unlist(keeppath)]
res.paths = paths.all.genes.filter[pnames]


# ------------------- visualize cps ------------------------ #
library(SingleCellExperiment)
library(Biobase)
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data")
deng <- readRDS("H0605_dengSCE.rds")
dengraw <- assay(deng)[,order(deng$pseudotime_sling)]
dengid <- deng$cell_type2[order(deng$pseudotime_sling)]


dt1 = data.frame(cp = unlist(list.ecp))
dt4 = data.frame(cp = unlist(list.tcai))
p1 = ggplot(data = dt1, aes(cp)) + geom_histogram(bins = 54) + xlab("Chang point candidates") +
  xlim(c(0,270)) + ggtitle("ecp result, 29 pathways")
p4 = ggplot(data = dt4, aes(cp)) + geom_histogram(bins = 54) + xlab("Chang point candidates") + 
  xlim(c(0,270)) + ggtitle("covcpd result, 29 pathways") + ylim(c(0,8))


library(cowplot)
pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0009_changepoint_deng_summary_compareecp", Sys.Date(),".pdf", sep = ""), width = 6, height = 5)
plot_grid(p1, p4, ncol = 1)
dev.off()

pnames

# ----------------------------------
# example: INOH_MM_HEDGEHOG, #15
library(gridExtra)
library(grid)


sce = deng
ordervar = "pseudotime_sling"
ctvar = "cell_type2"
cellorder = c("zy", "early2cell", "mid2cell", "late2cell", "4cell", "8cell",
              "16cell", "earlyblast", "midblast", "lateblast")
scename = "Deng et al."

o <- order(colData(sce)[,ordervar])
sce <- sce[,o]
dt <- data.frame(ptime = seq(1:ncol(sce)), ct = colData(sce)[,ctvar])



pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0009_CHANGEPOINT_DENG_SUMMARY_allpaths", Sys.Date(),".pdf", sep = ""), height =5, width = 7)
for(i in 1:length(pnames)){
  if(length(list.tcai[[i]])>0){
    pdeng = ggplot(dt, aes(x = ptime, y = factor(ct, levels = cellorder),
                           colour = factor(ct, levels = cellorder))) +
      geom_quasirandom(groupOnX = FALSE) +
      scale_color_tableau() + theme_classic() +
      theme(legend.title = element_blank(),
            legend.position = "top",
            plot.title = element_text(size=10)) +
      xlab("Slingshot pseudotime order") + ylab("Subgroups") +
      ggtitle(paste("Cells ordered by Slingshot pseudotime\n",scename,"\n", pnames[i])) +
      geom_vline(xintercept = list.ecp[[i]], color = "red")+
      # geom_vline(xintercept = list.inspect[[i]], color = "blue", linetype="dotted")+
      geom_vline(xintercept = list.tcai[[i]], color = "blue", linetype="longdash")+ 
      labs(caption= "red: ecp, blue dash: TonyCai", x = "pseudotime order", hjust = 1, gp = gpar(fontface = 3L, fontsize = 9)) #caption="ggplot2 caption", 
    print(pdeng)
  }
} 
dev.off()


pnames
i = 6
pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/A0009_CHANGEPOINT_DENG_SUMMARY_INOH_MM_TGF_BETA", Sys.Date(),".pdf", sep = ""), height =5, width = 7)
pdeng = ggplot(dt, aes(x = ptime, y = factor(ct, levels = cellorder),
                       colour = factor(ct, levels = cellorder))) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(size=10)) +
  xlab("Slingshot pseudotime order") + ylab("Subgroups") +
  ggtitle(pnames[i]) + 
  geom_vline(xintercept = list.tcai[[i]], color = "blue", linetype="longdash")+ 
  labs(x = "pseudotime order") #caption="ggplot2 caption", 
print(pdeng)
dev.off()
