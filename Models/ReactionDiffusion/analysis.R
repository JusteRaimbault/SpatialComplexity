
setwd(paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Models/ReactionDiffusion'))

library(dplyr)
library(ggplot2)
library(reshape2)
library(segmented)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

###
## 1) Liapounov configurations

experiment='20190112_1620_LOCAL'

res<-as.tbl(read.csv(paste0('exploration/',experiment,'.csv')))

resdir=paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Results/ReactionDiffusion/',experiment,'/');dir.create(resdir)

# forgot id
#sres=res%>%group_by(growthrate,population,alphalocalization,diffusion,diffusionsteps)%>%summarize(count=n())
res$id=paste0(res$growthrate,res$population,res$alphalocalization,res$diffusion,res$diffusionsteps)
res$id=as.numeric(as.factor(res$id))
sres=res%>%group_by(id)%>%summarize(count=n())
#sres[sres$count>500,]
#sres[sres$count==1000,]

tres = melt(data=res,id.vars=c("growthrate","population","replication","alphalocalization","diffusion","diffusionsteps","id"))
tres$t=as.numeric(substring(as.character(tres$variable),first=10))

#id=643,694,686,1022,801,3784,3692,3642
id=3692
d=tres[tres$id==id&tres$value>0,]

tstep=1
#tstep=5
g=ggplot(d[d$t%%tstep==0,],aes(x=t,y=log(value),group=t))
#g+geom_point(pch='.')+stat_smooth(method = 'gam')
g+geom_boxplot(outlier.size = 0.2)+xlab("time")+ylab("log(d)")+
  ggtitle(bquote(N[G]/P[tot]*"="*.(round(d$growthrate[1]/d$population[1],digits=5))*" ; "*alpha*"="*.(round(d$alphalocalization[1],digits=2))*" ; "*beta*"="*.(round(d$diffusion[1],digits=2))*" ; "*n[d]*"="*.(round(d$diffusionsteps[1]))*" ; "*lambda[1]*"="*.(round(lambda1[ids==id],digits=2))*" ; "*lambda[2]*"="*.(round(lambda2[ids==id],digits=2))*" ; "*t[b]*"="*.(round(breaks[ids==id],digits=2))*" ; "*R^2*"="*.(round(adjr2[ids==id],digits=3))))+stdtheme
ggsave(file=paste0(resdir,'configdist_boxplot_id',id,'.png'),width=20,height = 15,units='cm')


# test piecewise linear fitting
#seg<-segmented(lm(data=d,log(value)~t))
#slope(seg)$t[1,1]
#summary(seg)$adj.r.squared

adjr2=c();lambda1=c();lambda2=c();ids=sres$id[sres$count==1000];breaks=c()
for(i in 1:length(ids)){
  # this should be parallelized when will have more runs

#library(doParallel)
#cl <- makeCluster(50)
#registerDoParallel(cl)
#res <- foreach(i=1:length(ids)) %dopar% {
#  library(segmented)
  id = ids[i]
  show(i)
  d=tres[tres$id==id&tres$value>0,]
  seg<-segmented(lm(data=d,log(value)~t))
  lambda1=append(lambda1,slope(seg)$t[1,1]);lambda2=append(lambda2,slope(seg)$t[2,1]);breaks=append(breaks,seg$psi[2]);adjr2=append(adjr2,summary(seg)$adj.r.squared)
  #return(c(lambda1=slope(seg)$t[1,1],lambda2=slope(seg)$t[2,1],breaks=seg$psi[2],adjr2=summary(seg)$adj.r.squared))
}

#stopCluster(cl)

save(adjr2,lambda1,lambda2,ids,breaks,file=paste0(resdir,'lyapounov_',experiment,'.RData'))
#save(res,file=paste0(resdir,'lyapounov_',experiment,'.RData'))

summary(lambda1)
summary(lambda2)

ids[lambda1==max(lambda1)]
adjr2[lambda1==max(lambda1)]
ids[adjr2==min(adjr2)]
ids[adjr2==max(adjr2)]

params=res%>%group_by(id)%>%summarize(alpha=mean(alphalocalization),beta=mean(diffusion),nd=mean(diffusionsteps),relgrowthrate=mean(growthrate)/mean(population),growthrate=mean(growthrate),population=mean(population))
rownames(params)=params$id

dd=data.frame(adjr2,lambda1,lambda2,ids,params[ids,])

summary(lm(data=dd,lambda1~alpha+beta+nd+growthrate+population+relgrowthrate))
summary(lm(data=dd,lambda2~alpha+beta+nd+growthrate+population+relgrowthrate))

#g=ggplot(dd,aes(x=alpha,y=beta*nd,color=lambda1))
#g+geom_point()+facet_wrap(~cut(log(relgrowthrate),4))

g=ggplot(dd,aes(x=alpha,y=lambda1,color=cut(beta,4)))
g+geom_point(alpha=0.5)+geom_smooth()+scale_color_discrete(name=expression(beta))+
  xlab(expression(alpha))+ylab(expression(lambda[1]))+stdtheme#+theme(legend.position = c(0.23, 0.85))
ggsave(file=paste0(resdir,'lambda1_alpha_colbeta.png'),width=25,height=15,units='cm')

g=ggplot(dd,aes(x=alpha,y=lambda2,color=cut(beta,4)))
g+geom_point(alpha=0.5)+geom_smooth()+scale_color_discrete(name=expression(beta))+
  xlab(expression(alpha))+ylab(expression(lambda[2]))+stdtheme#+theme(legend.position = c(0.23, 0.85))
ggsave(file=paste0(resdir,'lambda2_alpha_colbeta.png'),width=25,height=15,units='cm')

g=ggplot(dd,aes(x=alpha,y=lambda2,color=cut(relgrowthrate,breaks=c(0,quantile(relgrowthrate,c(0.25,0.5,0.75)),max(relgrowthrate)+0.0001),right=T)))
g+geom_point(alpha=0.5)+geom_smooth()+scale_color_discrete(name=expression(N[G]/P[max]))+
  xlab(expression(alpha))+ylab(expression(lambda[2]))+stdtheme#+theme(legend.position = c(0.23, 0.85))
ggsave(file=paste0(resdir,'lambda2_alpha_colrelgrowthrate.png'),width=25,height=15,units='cm')



###
## 2) Morpho trajectories

experiment='20190117_1029_MORPHO_GRID'
# pb with time steps not uniformized -> use a scan instead
#res<-as.tbl(read.csv(paste0('exploration/',experiment,'.csv')))
resdir=paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Results/ReactionDiffusion/',experiment,'/');dir.create(resdir)

tstep=5

#header<-scan(file=paste0('exploration/',experiment,'.csv'),nmax=32,what='character',sep=",")
res<-scan(file=paste0('exploration/',experiment,'.csv'),nmax=-1,skip = 1,what='character')
res<-sapply(res,function(s){
  splited = strsplit(s,",")[[1]]
  tsteps = (length(splited)-7)/5
  data.frame(alpha=as.numeric(rep(splited[1],tsteps)),
             beta=as.numeric(rep(splited[2],tsteps)),
             nd=as.numeric(rep(splited[3],tsteps)),
             distance=as.numeric(splited[4:(tsteps+3)]),
             entropy=as.numeric(splited[(tsteps+4):(2*tsteps+3)]),
             ng=as.numeric(rep(splited[2*tsteps+4],tsteps)),
             id=as.numeric(rep(splited[2*tsteps+5],tsteps)),
             moran=as.numeric(splited[(2*tsteps+6):(3*tsteps+5)]),
             pmax=as.numeric(rep(splited[3*tsteps+6],tsteps)),
             replication=as.numeric(rep(splited[3*tsteps+7],tsteps)),
             rsquared=as.numeric(splited[(3*tsteps+8):(4*tsteps+7)]),
             slope=as.numeric(splited[(4*tsteps+8):(5*tsteps+7)]),
             t=(0:(tsteps-1))*tstep
             )
})

sres=as.tbl(data.frame(alpha=unlist(res[seq(1,length(res),13)]),
                       beta=unlist(res[seq(2,length(res),13)]),
                       nd=unlist(res[seq(3,length(res),13)]),
                       distance=unlist(res[seq(4,length(res),13)]),
                       entropy=unlist(res[seq(5,length(res),13)]),
                       ng=unlist(res[seq(6,length(res),13)]),
                       id=unlist(res[seq(7,length(res),13)]),
                       moran=unlist(res[seq(8,length(res),13)]),
                       pmax=unlist(res[seq(9,length(res),13)]),
                       replication=unlist(res[seq(10,length(res),13)]),
                       rsquared=unlist(res[seq(11,length(res),13)]),
                       slope=unlist(res[seq(12,length(res),13)]),
                       t=unlist(res[seq(13,length(res),13)])
                       ))

#save(sres,file=paste0('exploration/',experiment,'.RData'))
load(paste0('exploration/',experiment,'.RData'))

ampl = sres%>%group_by(id,t)%>%summarize(dmoran=max(moran)-min(moran),dentropy=max(entropy)-min(entropy),dslope=max(slope)-min(slope),ddistance=max(distance)-min(distance))
ampl[ampl$dentropy==max(ampl$dentropy),]

# test some examples of trajectories

#id=unique(sres$id)[1]
#id=unique(sres$id)[1:5]
id=669
#replications=unique(sres$replication)[1:10]
replications=unique(sres$replication)


d=sres[sres$id%in%id,]
#replications=unique(d$replication)[1:10]
#d=d[d$replication%in%replications,]

g=ggplot(d,aes(x=moran,y=entropy,color=alpha,group=interaction(replication,id)))
g+geom_point(size=0.1)+geom_path(arrow=arrow())

#s=sample.int(10000,1);show(s)
s=8578
set.seed(s)
g=ggplot(d[d$replication%in%sample(d$replication,50),],aes(x=moran,y=distance,color=slope,group=interaction(replication,id)))
g+geom_point(size=0.1)+geom_path(arrow=arrow(length = unit(0.08,'cm')))+stdtheme













