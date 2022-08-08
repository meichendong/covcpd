# S05 summary
# E0009V2 SUMMARY
# summary of E0009 longleaf results
# setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0009V2")
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0009ReduceMeanDiff")

lf = list.files()
scenarios = c("diffstructure", "samestructure", "diffstrdiffmean")
types =c("raw", "log", "imputeNPN")


tall = NULL
for(i in 1:length(lf)){
  tmp = read.table(lf[i])
  tall = rbind(tall, tmp)
}

tall = tall[!tall$package %in% c("factorcpt","TonyCai", "covcpd"),]

# setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/TCAI/E0014ReduceMeanDiff")
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/TCAI/E0022")

lf3 = list.files()
lf3 = lf3[grep("txt", lf3)]

tall3 = NULL
for(i in 1:length(lf3)){
  tmp = read.table(lf3[i])
  tall3 = rbind(tall3, tmp)
} 


tall3$TrueCP = factor(tall3$truecp, levels = c("200","500","200;400",
                                               "500;1000","200;400;600;800;1000"))

tall$TrueCP = factor(tall$truecp, levels = c("200","500","200;400",
                                             "500;1000","200;400;600;800;1000"))

colnames(tall3)
colnames(tall)
tall = tall[,-10]
tall3 = tall3[,-c(10,11)]
tall = rbind(tall, tall3)

# ---------------------------------
tall$prec[is.na(tall$prec)] <- 0
tall$recall[is.na(tall$recall)] <- 0
tall$F1score[is.na(tall$F1score)] <- 0
# saveRDS(tall, "C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0009_SUMMARY_TALL20211120.rds")
# saveRDS(tall, "C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_TALL_reduceMeandiff.rds")
# tall <- readRDS("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0009_SUMMARY_TALL20220226.rds")
# table(tall$package)
saveRDS(tall, "C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_TALL_reduceMeandiff20220321.rds")


# ----------------------- visualize ------------------------- #
library(ggplot2)
# ncpx = 10
nsample = c(200,500) #,200
nchange = c(1, 2) #, 10     
unique(tall$package)
tall2 = tall[tall$package %in% c("ecp", "ccid -- euclidean","changepoint.geo: mean change",
                                 "changepoint.geo: cov change",  "covcpd"),] # ,"factorcpt" "InspectChangepoint",
tall2$package = factor(tall2$package, levels = c("covcpd","ecp", "ccid -- euclidean","changepoint.geo: mean change",
                                                 "changepoint.geo: cov change")) #,"factorcpt" "InspectChangepoint", 
library(plyr)
unique(tall$type)
tall2$type2 = "Diff.Str, Diff.Mean"
tall2$type2[tall2$type =="diffstructure"] = "Diff.Str, Same.Mean"
tall2$type2[tall2$type =="samestructure"] = "Same.Str, Diff.Mean"

tall2 = tall2[tall2$type !="samestructure",]

table(tall2$diffpct, tall2$package)

tall2 = tall2[tall2$diffpct %in% c(0,5,10,20,50),]

nx = "raw" 
nn = 1  ## for 1..3
ncpx = nchange[nn]
# pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0009_SUMMARY_",nx,"_NCP_",ncpx,"_TIME_",Sys.Date(),".pdf", sep = ""),
#     height = 6, width = 6)
# p1 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=processtime, color = TrueCP)) + 
#   geom_boxplot() + facet_grid(type2 ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
#   theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + ylab("Process time (seconds)") + xlab("")
# 
# print(p1)
# dev.off()
pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_",nx,"_NCP_",ncpx,"_F1_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p2 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=F1score, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top")+ xlab("")
print(p2) 
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_",nx,"_NCP_",ncpx,"_PREC_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p3 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=prec, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + ylab("Precision") + xlab("")
print(p3)
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_",nx,"_NCP_",ncpx,"_RECALL_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p4 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=recall, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + ylab("Recall") + xlab("")
print(p4)
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_",nx,"_NCP_",ncpx,"_ARI_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p5 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=ARI, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + xlab("")
print(p5) 
dev.off()


pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_",nx,"_NCP_",ncpx,"_MSE_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p6 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=MSE, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + xlab("")
print(p6) 
dev.off()


# ==============================================
# SAME STRUCTURE
tall2 = tall[tall$package %in% c("ecp", "ccid -- euclidean","changepoint.geo: mean change",
                                 "changepoint.geo: cov change",  "covcpd"),] # ,"factorcpt" "InspectChangepoint",
tall2$package = factor(tall2$package, levels = c("covcpd","ecp", "ccid -- euclidean","changepoint.geo: mean change",
                                                 "changepoint.geo: cov change")) #,"factorcpt" "InspectChangepoint", 
library(plyr)
unique(tall$type)
tall2$type2 = "Diff.Str, Diff.Mean"
tall2$type2[tall2$type =="diffstructure"] = "Diff.Str, Same.Mean"
tall2$type2[tall2$type =="samestructure"] = "Same.Str, Diff.Mean"

tall2 = tall2[tall2$type =="samestructure",]

table(tall2$diffpct, tall2$package)
tall2 = tall2[tall2$diffpct %in% c(0,5,10,20,50),]


nx = "raw" 
nn = 1  ## for 1..3
ncpx = nchange[nn]
# 
# pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0009_SUMMARY_",nx,"_NCP_",ncpx,"_TIME_",Sys.Date(),".pdf", sep = ""),
#     height = 6, width = 6)
# p1 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=processtime, color = TrueCP)) + 
#   geom_boxplot() + facet_grid(type2 ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
#   theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + ylab("Process time (seconds)") + xlab("")
# 
# print(p1)
# dev.off()
pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_SAMESTR_",nx,"_NCP_",ncpx,"_F1_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p2 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=F1score, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top")+ xlab("")
print(p2) 
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_SAMESTR_",nx,"_NCP_",ncpx,"_PREC_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p3 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=prec, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + ylab("Precision") + xlab("")
print(p3)
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_SAMESTR_",nx,"_NCP_",ncpx,"_RECALL_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p4 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=recall, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + ylab("Recall") + xlab("")
print(p4)
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_SAMESTR_",nx,"_NCP_",ncpx,"_ARI_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p5 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=ARI, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + xlab("")
print(p5) 
dev.off()

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S05_SUMMARY_SAMESTR_",nx,"_NCP_",ncpx,"_MSE_",Sys.Date(),".pdf", sep = ""),
    height = 6, width = 6)
p6 = ggplot(tall2[tall2$ncp ==ncpx & tall2$normalize == nx,], aes(x=package, y=MSE, color = TrueCP)) + 
  geom_boxplot() + facet_grid(diffpct ~., scales = "free_y") + guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7), legend.position = "top") + xlab("")
print(p6) 
dev.off()