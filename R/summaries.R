#' Build a matrix of probability of who infected whom from a MCMC output
#' @param record Output from inferTTreeM function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Matrix of probability of who infected whom
#' @export
computeMatWIWM = function(record,burnin=0.5)
{
  #Remove burnin
  if (burnin>0) record=record[round(length(record)*burnin):length(record)]
  m=length(record)
  t1=extractTTreeM(record[[1]]$ctree)
  n=length(unique(t1$obs[, 2])) #Number of sampled individuals
  mat=matrix(0,n,n)
  colnames(mat)<-1:n
  rownames(mat)<-1:n

  for (i in 1:length(record))
  {
    ttree=extractTTreeM(record[[i]]$ctree)$ttree
    infectors=ttree[1:n,3]
    infecteds=1:n
    w=which(infectors==0|infectors>n)
    infectors=infectors[setdiff(1:n,w)]
    infecteds=infecteds[setdiff(1:n,w)]
    mat[cbind(infectors,infecteds)]=mat[cbind(infectors,infecteds)]+1/length(record)
  }
  return(mat)
}

#' Convert resTransPhyloM object to coda mcmc format
#' @param record Output from inferTTreeM function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Object of class mcmc from coda package
#' @export
convertToCodaM = function(record,burnin=0.5) {
  record=record[max(1,round(length(record)*burnin)):length(record)]
  mat=cbind(
    sapply(record,function(x) x$pi),
    sapply(record,function(x) x$lambda),
    sapply(record,function(x) x$kappa),
    sapply(record,function(x) x$off.r))
  colnames(mat)<-c('pi','lambda','kappa','off.r')
  return(coda::as.mcmc(mat))
}

#' Print function for resTransPhyloM objects
#' @param x output from inferTTreeM
#' @param ... Additional parameters are passed on
#' @return Print out details of TransPhyloMulti results
#' @export
print.resTransPhyloM <- function(x, ...)
{
  stopifnot(inherits(x, "resTransPhyloM"))
  cat( 'Result from TransPhylo analysis\n')
  coda=convertToCodaM(x,0.5)
  for (nam in colnames(coda)) {
    v=coda[,nam]
    v=sort(v)
    vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
    cat(sprintf('%s=%.2e [%.2e;%.2e]\n',nam,vals[1],vals[2],vals[3]))
  }
  invisible(x)
}
