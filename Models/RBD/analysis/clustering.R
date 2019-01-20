setwd(paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Models/RBD'))

library(dplyr)
library(ggplot2)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

#experiment = '20190114_1101_GRID'
experiment='20190117_1029_GRID'
res<-as.tbl(read.csv(paste0('exploration/',experiment,'.csv')))
resdir = paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Results/RBD/',experiment,'/');dir.create(resdir)

nbootstraps=1000
params=c("id","replication","weightDensity","weightCenter","weightRoad","centerNumber","paramMode")
radiusNames=paste0("clustersRadius",0:(nbootstraps-1))
nullRadiusNames=paste0("nullClustersRadius",0:(nbootstraps-1))
distNames=paste0("partitionDistances",0:(nbootstraps-1))
nullDistNames=paste0("nullPartitionDistances",0:(nbootstraps-1))
clustersWithinssNames=paste0("clustersWithinss",0:(nbootstraps-1))
nullClustersWithinssNames=paste0("nullClustersWithinss",0:(nbootstraps-1))
profileDistEuclNames=paste0("profileDistEucl",0:(nbootstraps-1))
nullProfileDistEuclNames=paste0("nullProfileDistEucl",0:(nbootstraps-1))
profileDistTaumaxNames=paste0("profileDistTaumax",0:(nbootstraps-1))
nullProfileDistTaumaxNames=paste0("nullProfileDistTaumax",0:(nbootstraps-1))
# forgot overlap ? anyway badly implemented

sres = res[,params]
sres$dist = rowMeans(res[,distNames])
sres$nullDist = rowMeans(res[,nullDistNames])
sres$radius = rowMeans(res[,radiusNames])
sres$nullRadius = rowMeans(res[,nullRadiusNames])
sres$withinss = rowMeans(res[,clustersWithinssNames])
sres$nullWithinss = rowMeans(res[,nullClustersWithinssNames])
sres$profileEucl = rowMeans(res[,profileDistEuclNames])
sres$nullProfileEucl = rowMeans(res[,nullProfileDistEuclNames])
sres$profileTaumax = rowMeans(res[,profileDistTaumaxNames])
sres$nullProfileTaumax = rowMeans(res[,nullProfileDistTaumaxNames])

g=ggplot(sres)
g+geom_smooth(aes(x=centerNumber,y=radius,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)))+
  geom_smooth(aes(x=centerNumber,y=nullRadius,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)),linetype=2)+
  scale_color_discrete(name=expression("mode;"*w[D]*";"*w[C]*";"*w[R]))+xlab("Number of centers")+ylab("Average cluster radius")+stdtheme
ggsave(file=paste0(resdir,'radius.png'),width=25,height=15,units='cm')

g+geom_smooth(aes(x=centerNumber,y=dist,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)))+
  geom_smooth(aes(x=centerNumber,y=nullDist,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)),linetype=2)+
  scale_color_discrete(name=expression("mode;"*w[D]*";"*w[C]*";"*w[R]))+xlab("Number of centers")+ylab("Average distance to centers")+stdtheme
ggsave(file=paste0(resdir,'distance.png'),width=25,height=15,units='cm')

g+geom_smooth(aes(x=centerNumber,y=withinss,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)))+
  geom_smooth(aes(x=centerNumber,y=nullWithinss,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)),linetype=2)+
  scale_color_discrete(name=expression("mode;"*w[D]*";"*w[C]*";"*w[R]))+xlab("Number of centers")+ylab("Within-cluster variance")+stdtheme
ggsave(file=paste0(resdir,'withinss.png'),width=25,height=15,units='cm')

g+geom_smooth(aes(x=centerNumber,y=profileEucl,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)))+
  geom_smooth(aes(x=centerNumber,y=nullProfileEucl,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)),linetype=2)+
  scale_color_discrete(name=expression("mode;"*w[D]*";"*w[C]*";"*w[R]))+xlab("Number of centers")+ylab("Average lagged profile distance")+stdtheme
ggsave(file=paste0(resdir,'profiledisteucl.png'),width=25,height=15,units='cm')

g+geom_smooth(aes(x=centerNumber,y=profileTaumax,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)))+
  geom_smooth(aes(x=centerNumber,y=nullProfileTaumax,color=paste0(paramMode,";",weightDensity,";",weightCenter,";",weightRoad)),linetype=2)+
  scale_color_discrete(name=expression("mode;"*w[D]*";"*w[C]*";"*w[R]))+xlab("Number of centers")+ylab("Average optimal lag distance")+stdtheme
ggsave(file=paste0(resdir,'profiledisttaumax.png'),width=25,height=15,units='cm')





