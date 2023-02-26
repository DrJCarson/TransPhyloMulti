#' Converts an ape phylo object into a phylogenetic tree with host information
#' @param tr phylo object with labels indicating host, eg 1.1 etc
#' @param dateLastSample date of the last sample
#' @return phylogenetic tree
#' @export
ptreeFromPhyloM <- function(tr,dateLastSample) {
  ptree=ptreeFromPhylo(tr,dateLastSample = dateLastSample)
  ptree$host=as.numeric(unlist(strsplit(ptree$nam,'\\.'))[seq(1,length(ptree$nam)*2,2)]) #Extract host information from leaf names
  n=length(ptree$host) # Number of leaves
  w=(n+1):nrow(ptree$ptree) # Indexes of internal nodes
  o=order(ptree$host*1e10+ptree$ptree[1:n,1]) #Order needed for leaves: by host first and then from oldest to most recent
  o2=order(-ptree$ptree[w,1]) #Order needed for internal nodes: from most recent to oldest (root)
  o=c(o,n+o2)
  ot=o;ot[o]=1:length(o)
  ptree$host=ptree$host[o[1:n]] #Reorder host
  ptree$nam=ptree$nam[o[1:n]] #Reorder nam
  ptree$ptree=ptree$ptree[o,] #Reorder ptree
  ptree$ptree[w,2]=ot[ptree$ptree[w,2]]
  ptree$ptree[w,3]=ot[ptree$ptree[w,3]]
  class(ptree)<-'ptreem'
  return(ptree)
}

#' Converts a resTransPhyloM into a resTransPhylo
#' @param r object of class resTransPhyloM which is the output of inferTTreeM
#' @return object of class resTransPhylo produced by removing non-first samples
#' @export
removeMulti <- function(r) {
  res=r
  ctree=res[[1]]$ctree$ctree
  remleaves=F
  for (i in 2:nrow(ctree)) {
    if (ctree[i,2]!=0 || ctree[i,3]!=0) break
    if (ctree[i,4]==ctree[i-1,4]) remleaves=c(remleaves,T) else remleaves=c(remleaves,F)
  }
  wrl=which(remleaves)
  for (i in 1:length(res)) {
    res[[i]]$neg=res[[i]]$lm_const
    nam=res[[i]]$ctree$nam[remleaves==F]
    nam=unlist(strsplit(nam,'\\.'))
    nam=nam[seq(1,length(nam),2)]
    res[[i]]$ctree$nam=nam
    ctree=res[[i]]$ctree$ctree
    rem=c(remleaves,rep(F,nrow(ctree)-length(remleaves)))
    parents=rep(NA,nrow(ctree)+1)
    parents[ctree[,2]+1]=1:nrow(ctree)
    parents[ctree[,3]+1]=1:nrow(ctree)
    parents=parents[-1]
    for (j in wrl) {
      cur=j
      while (rem[cur]==T || ctree[cur,3]==0) {
        rem[cur]=T
        cur=parents[cur]
      }
      rem[cur]=T
    }
    w=which(!rem)
    map=rep(NA,nrow(ctree))
    map[w]=1:length(w)
    for (j in 1:length(map)) {
      cur=j
      while (all(is.na(map[cur]))) {cur=ctree[cur,2:3];cur=cur[!is.na(cur)&cur>0];if (length(cur)==0) break;}
      if (length(cur)>0) {tmp=map[cur];map[j]=tmp[!is.na(tmp)]}
    }
    ctree=ctree[!rem,]
    w=which(ctree[,2]>0);ctree[w,2]=map[ctree[w,2]]
    w=which(ctree[,3]>0);ctree[w,3]=map[ctree[w,3]]
    ctree[,4]=TransPhylo:::.computeHost(ctree)
    res[[i]]$ctree$ctree=ctree
  }
  class(res)<-'resTransPhylo'
  return(res)
}

#' Return the medoid from a resTransPhyloM
#' @param record Output from inferTTreeM function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return The medoid
#' @export
medTTreeM = function(record,burnin=0.5)
{
  res2=removeMulti(record)
  med2=medTTree(res2,burnin)
  i=1
  while (nrow(res2[[i]]$ctree$ctree)!=nrow(med2$ctree) || any(res2[[i]]$ctree$ctree!=med2$ctree)) i=i+1
  return(res[[i]]$ctree)
}
