\name{probBinom}
\alias{probBinom}
\title{Compute P-values for Medip-seq and MRE-seq data.}
\usage{
    probBinom(t,size1,size2,c1,c2)}
\description{
    Compute P-values.
}
\arguments{
    \item{t}{The real value for random variable according to dataset.}
    \item{size1}{The sum of Medip-seq real reads of the each CpG site for
    control and treatment sample.}
    \item{size2}{The sum of MRE-seq real reads of the each CpG site for control
    and treatment sample.}
    \item{c1}{The scaling factor for MeDip-seq data.}
    \item{c2}{The scaling factor for MRE-seq data.}
}
\value{
    \item{p}{The P-values for testing the methylation expression levels for each
    CpG sites.}
}
\examples{
    set.seed(1234)
    t<-0.1
    size1<-sample(1:1000, 1, replace=TRUE)
    size2<-sample(1:1000, 1, replace=TRUE)
    c1<-1
    c2<-2
    result<-probBinom(t,size1,size2,c1,c2)
}






