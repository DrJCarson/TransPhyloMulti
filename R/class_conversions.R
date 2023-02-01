#' Converts an ape phylo object into a phylogenetic tree with host information
#' @param tr phylo object with labels indicating host, eg 1.1 etc
#' @param dateLastSample date of the last sample
#' @return phylogenetic tree
#' @export
ptreeFromPhyloM <- function(tr,dateLastSample) {
  ptree=ptreeFromPhylo(p,dateLastSample = dateLastSample)
  ptree$host=as.numeric(unlist(strsplit(ptree$nam,'\\.'))[seq(1,length(ptree$nam)*2,2)])
  class(ptree)<-'ptreem'
  return(ptree)
}
