
setwd(paste0(Sys.getenv('CS_HOME'),'/SpatialComplexity/Models/ReactionDiffusion'))

library(dplyr)
library(ggplot2)
library(reshape2)
library(segmented)

#source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

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
id=3642
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

res<-scan(file=paste0('exploration/',experiment,'.csv'))





