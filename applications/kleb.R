library(ape)
library(TransPhylo)
library(TransPhyloMulti)
set.seed(0)

td=read.nexus('kleb.nwk')
td$edge.length=td$edge.length/365
hosts=unlist(strsplit(td$tip.label,'-'))[seq(1,Ntip(td)*2,2)]
names(hosts)=td$tip.label
pt=ptreeFromPhyloM(td,lubridate::decimal_date(as.Date('2015/04/02')),hosts)
res=inferTTreeM(pt,mcmc.length=1e5,mcmc.thinning=1,w.shape=1,w.scale=0.5,obs.end=2016)

print(res)
plot(res[[length(res)]]$ctree)
mat=computeMatWIWM(res)
lattice::levelplot(mat,xlab='',ylab='')
