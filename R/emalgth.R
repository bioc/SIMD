#sum=1
#' Calculate the probability on condition that the sums equal to 1.
#' @description Calculate the probability on condition that only a single CpG 
#' contributes to a short read. 
#' @param X A matrix about X, the elements in X takes values on {0,1} and
#' satisfy the sums of each row equal to 1.
#' @return y1 The probability when sums equal to 1.
#' @examples 
#' set.seed(123)
#' d <- matrix(0, nrow=200, ncol=50)
#' random_num<-sample(1:50, 200, replace=TRUE)
#' for(i in 1:nrow(d)){
#'  d[i,random_num[i]]<-1
#' }
#' result<-emalgth(d)
#' head(result)
#' @export

emalgth<-function(X)
{
    N<-apply(X,2,sum)
    y<-X
    for(i in seq_len(length(N))){
    y[y[,i]!=0,i]<-N[i]
    }
    N1<-apply(y,1,sum)
    CpGnum<-rep(1,length(X[1,]))
    N11<-N1%o%CpGnum
    y1<-y/N11
    return(y1)
}
