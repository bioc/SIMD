#' Compute P-values for Medip-seq and MRE-seq data.
#'
#' @description Compute P-values.
#' @param t The real value for random variable according to dataset.
#' @param size1 The sum of Medip-seq real reads of the each CpG site for
#' control and treatment sample.
#' @param size2 The sum of MRE-seq real reads of the each CpG site for control
#' and treatment sample.
#' @param c1 The scaling factor for MeDip-seq data.
#' @param c2 The scaling factor for MRE-seq data.
#' @return p The P-values for testing the methylation expression levels for each
#' CpG sites.
#' @examples 
#' set.seed(1234)
#' t<-0.1
#' size1<-sample(1:1000, 1, replace=TRUE)
#' size2<-sample(1:1000, 1, replace=TRUE)
#' c1<-1
#' c2<-2
#' result<-probBinom(t,size1,size2,c1,c2)
#' @export


probBinom<-function(t,size1,size2,c1,c2){
    prob1<-1/(1+c1)
    prob2<-1/(1+c2)
    i<-0:size1
    j<-size1:0
    pij<- dbinom(i,prob=prob1,size=size1)
    r<-0:size2
    k<-size2:0
    prk<- dbinom(r,prob=prob2,size=size2)
    ik<-outer(i,k)
    jr<-outer(j,r)
    ijrk<-c1*ik-c2*jr
    prob<-outer(pij,prk)
    ir<-which(abs(ijrk)>=t, arr.ind=TRUE)
    p<-sum(prob[ir])
    return(p)
}
