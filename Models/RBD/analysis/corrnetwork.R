
setwd(paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Models/RBD'))

library(dplyr)
library(ggplot2)


res<-as.tbl(read.csv('exploration/20190102_17_17_GRID.csv'))

# check counts
res %>% group_by(id) %>% summarise(count=n())

nBootstraps = 100
partitionDistancesNames = paste0('partitionDistances',0:(nBootstraps-1))
nullPartitionDistancesNames = paste0('nullPartitionDistances',0:(nBootstraps-1))

#hist(unlist(data.frame(res[which(res$id==89)[1],nullPartitionDistancesNames])),breaks=20)
#hist(unlist(data.frame(res[which(res$id==89)[5],nullPartitionDistancesNames])),breaks=20)
# distributions are a fucking mess
res$avgPartitionDistance = rowMeans(res[,partitionDistancesNames])
res$partitionDistanceIC = apply(res[,partitionDistancesNames],1,function(row){1.96*sd(row)/sqrt(length(row))})
res$avgNullPartitionDistance = rowMeans(res[,nullPartitionDistancesNames])
res$nullPartitionDistanceIC = apply(res[,nullPartitionDistancesNames],1,function(row){1.96*sd(row)/sqrt(length(row))})

# rq : there should exist an appropriate statistic for testing \neq nullmodel
res$significance = abs(res$avgPartitionDistance - res$avgNullPartitionDistance)/(res$partitionDistanceIC+res$nullPartitionDistanceIC)

res[res$significance>0.8,c("centerNumber","pCorrDist","paramMode")]
res[res$significance==max(res$significance),c("centerNumber","pCorrDist","paramMode")]





