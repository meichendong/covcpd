

# ---------------------- null: benchmark all methods -------------------------
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0009NULL")
lfn = list.files()

tnull = NULL
for(i in 1:length(lfn)){
  tmp = read.table(lfn[i])
  # cat(i,"\t", nrow(tmp),"\n")
  if(nrow(tmp)>16){
    dup = table(tmp$package)
    dupn = names(dup[dup>1])
    duptmp = tmp[tmp$package %in% dupn,]
    nodup = tmp[!tmp$package %in% dupn,]
    nodup2 = duptmp[cumsum(dup[dupn])-1,]
    tmp = rbind(nodup, nodup2)
  }
  tmp$nsample = gsub("NSAM","",strsplit(lfn[i],"_")[[1]][5])
  tnull = rbind(tnull, tmp)
}
tnull$detect = tnull$cp > 0


# merge inspect, covcpd
# setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/TCAI/E0011NULL_20220124")

# setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/TCAI/E0011NULL") #20220214
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/TCAI/E0011NULLtest") #2022march 

lfn = list.files()
lfcp = lfn[grep("cp1", lfn)]
lfnull = setdiff(lfn, lfcp)

tnull.covcpd = NULL
tnull.v2 = NULL #data.frame(detect = NULL, package = NULL, simu = NULL, nsample = NULL)
for(i in 1:length(lfnull)){
  tmp = read.table(lfnull[i])
  if(nrow(tmp) ==26){
    tnull.covcpd = rbind(tnull.covcpd, t(tmp)) 
  } else {
    tnull.v2 = rbind(c(T, "covcpd", tmp[nrow(tmp)-3,], tmp[nrow(tmp)-2,]), tnull.v2)
  }
}
colnames(tnull.v2) =c("detect","package", "simu","nsample")
tnull.v2 = as.data.frame(tnull.v2)
tnull.v2$nsample = as.numeric(tnull.v2$nsample)
tnull.v2$simu = as.numeric(tnull.v2$simu)

table(tnull.covcpd[,24])
check500 = tnull.covcpd[tnull.covcpd[,24] == 500,]
table(check500[,1])
check500file = "E0011_distnull_8_500.txt"
tmp = read.table(check500file)
table(tnull.v2$nsample)


colnames(tnull)
tnull.covcpd.summary = as.data.frame(tnull.covcpd[,c(2,4,23,24,25)])
tnull.covcpd.summary$pget = as.numeric(tnull.covcpd.summary[,1])
tnull.covcpd.summary$pthresh = as.numeric(tnull.covcpd.summary[,2])
tnull.covcpd.summary$detect = tnull.covcpd.summary$pget <= tnull.covcpd.summary$pthresh
tnull.covcpd.summary$package = "covcpd"
tnull.covcpd.summary2 = tnull.covcpd.summary[,c(3,4,8,9)]
colnames(tnull.covcpd.summary2)[1:2] = c("simu","nsample")

table(tnull.covcpd.summary2$nsample)
table(tnull.covcpd.summary2$nsample, tnull.covcpd.summary2$detect)

tnull2 = tnull[!grepl("covcpd",tnull$package, ignore.case = T),]
tnull2 = tnull[!grepl("tonycai",tnull$package, ignore.case = T),c(1,4,7,8)]
tnull.all = rbind(tnull2, tnull.covcpd.summary2, tnull.v2)
tnull.all$detect = as.logical(tnull.all$detect)


table(tnull.all$package, tnull.all$nsample)

nns = unique(tnull.all$nsample)
j=1
sumtall = NULL

for(j in 1:length(nns)){
  dttmp = tnull.all[tnull.all$nsample == nns[j],]
  sumt = aggregate(dttmp$detect, list(package = dttmp$package), sum)
  sumt$t1err = sumt$x / length(unique(dttmp$simu))
  sumt$nsample = nns[j]