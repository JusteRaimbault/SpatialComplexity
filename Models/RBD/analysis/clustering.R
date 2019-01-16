setwd(paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Models/RBD'))

library(dplyr)
library(ggplot2)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

experiment = '20190114_1101_GRID'
res<-as.tbl(read.csv(paste0('exploration/',experiment,'.csv')))
resdir = paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Results/RBD/',experiment,'/');dir.create(resdir)

nbootstraps=1000
params=c("id","replication","weightDensity","weightCenter","weightRoad","centerNumber","paramMode")
radiusNames=paste0("clustersRadius",0:(nbootstraps-1))
nullRadiusNames=paste0("nullClustersRadius",0:(nbootstraps-1))
distNames=paste0("partitionDistances",0:(nbootstraps-1))
nullDistNames=paste0("nullPartitionDistances",0:(nbootstraps-1))

sres = res[,params]
sres$dist = rowMeans(res[,distNames])
sres$nullDist = rowMeans(res[,nullDistNames])
sres$radius = rowMeans(res[,radiusNames])
sres$nullRadius = rowMeans(res[,nullRadiusNames])

g=ggplot(sres)
g+geom_smooth(aes(x=centerNumber,y=radius,color=interaction(weightDensity,weightCenter,weightRoad,paramMode)))+
  geom_smooth(aes(x=centerNumber,y=nullRadius,color=interaction(weightDensity,weightCenter,weightRoad,paramMode)),linetype=2)+
  scale_color_discrete(name=expression(w[D]*";"*w[C]*";"*w[R]*";mode"))+xlab("Number of centers")+ylab("Average cluster radius")+stdtheme
ggsave(file=paste0(resdir,'radius.png'),width=25,height=15,units='cm')

g=ggplot(sres)
g+geom_smooth(aes(x=centerNumber,y=dist,color=interaction(weightDensity,weightCenter,weightRoad,paramMode)))+
  geom_smooth(aes(x=centerNumber,y=nullDist,color=interaction(weightDensity,weightCenter,weightRoad,paramMode)),linetype=2)+
  scale_color_discrete(name=expression(w[D]*";"*w[C]*";"*w[R]*";mode"))+xlab("Number of centers")+ylab("Average distance to centers")+stdtheme
ggsave(file=paste0(resdir,'distance.png'),width=25,height=15,units='cm')







