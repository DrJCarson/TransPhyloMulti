#' Converts an ape phylo object into a phylogenetic tree with host information
#' @param tr phylo object with labels indicating host, eg 1.1 etc
#' @param dateLastSample date of the last sample
#' @return phylogenetic tree
#' @export
ptreeFromPhyloM <- function(tr,dateLastSample) {
  ptree=ptreeFromPhylo(p,dateLastSample = dateLastSample)
  ptree$host=as.numeric(unlist(strsplit(ptree$nam,'\\.'))[seq(1,length(ptree$nam)*2,2)])
  o=order(ptree$host)
  o2=o;o2[o]=1:length(o)
  n=length(ptree$host)
  ptree$host=ptree$host[o]
  ptree$nam=ptree$nam[o]
  ptree$ptree[1:n,]=ptree$ptree[o,]
  w=which(ptree$ptree[,2]<=n & ptree$ptree[,2]>0)
  ptree$ptree[w,2]=o2[ptree$ptree[w,2]]
  w=which(ptree$ptree[,3]<=n & ptree$ptree[,3]>0)
  ptree$ptree[w,3]=o2[ptree$ptree[w,3]]
  class(ptree)<-'ptreem'
  return(ptree)
}
