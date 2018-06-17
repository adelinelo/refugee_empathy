#function for ols, calculating regression outputs with possibility for robust, clustered, robust-clustered standard errors.
#from Mahmood Arai
clx <-function(fm, cluster=NULL){
  if(length(cluster)!=0){
    M <- length(unique(cluster))
    dfcw <- fm$df / (fm$df - (M -1))
    N <- length(cluster)
    dfc <- (M/(M-1))*((N-1)/(N-fm$rank))
    u <- apply(estfun(fm),2,
               function(x) tapply(x, cluster, sum))
    vcovCL <- dfc*sandwich(fm, meat=crossprod(u)/N)*dfcw
    coeftest(fm, vcovCL) }
  else{
    vcovR<-vcovHC(fm,sandwich=TRUE,type="HC")#heteroskedastic consistent cov matrix estimation, with White's estimator
    coeftest(fm,vcovR)
  }
}
#mahmood arai - just cluster
cl   <- function(fm, dat,cluster){
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  coeftest(fm, vcovCL) }