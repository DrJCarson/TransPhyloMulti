library(ape)
library(TransPhylo)
library(TransPhyloMulti)


tree=read.tree('pseudo.nwk')
h=as.matrix(read.csv('pseudo.csv',header = F))
host=h[,2]
names(host)=h[,1]
pt=ptreeFromPhyloM(tree,2008,host=host)
res<-inferTTreeM(pt,mcmc.length=1e5,mcmc.thinning=1,w.shape=2,w.scale=5,obs.end=2010)
