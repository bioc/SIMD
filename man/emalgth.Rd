\name{emalgth}
\alias{emalgth}
\title{Calculate the probability on condition that the sums equal to 1.}
\usage{
    emalgth(X)
}
\description{
    Calculate the probability on condition that only a single CpG contributes
    to a short read.
}
\arguments{
    \item{X}{A matrix about X, the elements in X takes values on {0,1} and
    satisfy the sums of each row equal to 1.}
}

\value{
    \item{y1}{The probability when sums equal to 1.}
}
\examples{
    set.seed(123)
    d <- matrix(0, nrow=200, ncol=50)
    random_num<-sample(1:50, 200, replace=TRUE)
    for(i in 1:nrow(d)){
        d[i,random_num[i]]<-1
    }
    result<-emalgth(d)
    head(result)
}
