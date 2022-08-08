# s03
# SUMMARY FOR E0019 GAUSSIAN SIMULATIONS
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0019")
lf = list.files()
lf = lf[grep("txt", lf)]

tall = NULL
for (i in 1:length(lf)){
  tmp = read.table(lf[i])
  tall = rbind(tall, tmp)
}

library(plyr)
unique(tall$type)
tall$type2 = "Diff.Str, Diff.Mean"
tall$type2[tall$type =="diffstructure"] = "Diff.Str, Same.Mean"
tall$type2[tall$type =="samestructure"] = "Same.Str, Diff.Mean"
tall$package = factor(tall$package, levels = c("covcpd","ecp", "ccid -- euclidean","changepoint.geo: mean change",
                                               "changepoint.geo: cov change")) #,"factorcpt" "InspectChangepoint", 


table(tall$truecp)

nchange = c(1,2)


nx = "raw" 
nn = 1  ## for 1..3
ncpx = nchange[nn]

library(ggplot2)

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S03_Gaussian_SUMMARY_",nx,"_NCP_",ncpx,"_F1_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p2 = ggplot(tall[tall$ncp ==ncpx & tall$normalize == nx,], aes(x=package, y=F1score, color = truecp)) + 
  geom_boxplot() + facet_grid(type2 ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top")+ xlab("")
print(p2) 
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S03_Gaussian_SUMMARY_",nx,"_NCP_",ncpx,"_PREC_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p3 = ggplot(tall[tall$ncp ==ncpx & tall$normalize == nx,], aes(x=package, y=prec, color = truecp)) + 
  geom_boxplot() + facet_grid(type2 ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top")+ xlab("")
print(p3)
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S03_Gaussian_SUMMARY_",nx,"_NCP_",ncpx,"_RECALL_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p4 = ggplot(tall[tall$ncp ==ncpx & tall$normalize == nx,], aes(x=package, y=recall, color = truecp)) + 
  geom_boxplot() + facet_grid(type2 ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top")+ xlab("")
print(p4)
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S03_Gaussian_SUMMARY_",nx,"_NCP_",ncpx,"_ARI_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p5 = ggplot(tall[tall$ncp ==ncpx & tall$normalize == nx,], aes(x=package, y=ARI, color = truecp)) + 
  geom_boxplot() + facet_grid(type2 ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top")+ xlab("")
print(p5) 
dev.off()


