#sum>=1
#' Calculate the probability on condition that the sums more than 1.
#' 
#' @description Calculate the probability on condition that at least a CpG 
#' contributes to a short read.
#' @param X A matrix about X, the elements in X takes values on {0,1} and
#' satisfy the sums of each row more than 1.
#' @return y1 The probability when sums more than 1.
#' @examples 
#' set.seed(123)
#' d <- matrix(0, nrow=200, ncol=50)
#' random_num <- sample(1:10, 200, replace=TRUE)
#' for(i in 1:nrow(d)){
#'     temp <- sample(1:50, random_num[i], replace=FALSE)
#'     d[i,temp] <- 1
#' }
#' result <- emalgth1(d)
#' head(result)
#' @export



emalgth1 <- function(X)
{
    N <- colSums(X)/length(X[, 1])
    p <- 1-N
    p1 <- 1-prod(p)
    N <- N*p1
    y <- X
    for(i in seq_len(length(N))){
        y[y[, i]!=0,i] <- N[i]
    }
    N1 <- 1-y
    N2 <- apply(N1, 1, prod)
    CpGnum <- rep(1, length(X[1, ]))
    N11 <- 1-N2 %o% CpGnum
    y1 <- y/N11
    return(y1)
}
