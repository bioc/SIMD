\name{EMtest}
\alias{EMtest}
\title{Inferring the methylation expression level of single sites.}
\usage{
    EMtest(datafile=NULL,chrstring=NULL,cpgfile,mrecpgfile=NULL,
    writefile=NULL,reportfile=NULL,mreratio=3/7, psd=2,mkadded=1,f=1)
}
\description{Using statistical framework and EM algorithm to infer the
methylation expression level of single sites.
}
\arguments{
    \item{datafile}{The files of sample. (datafile should be cbind(data1,data2,
    data3,data4), where data1 and data2 are Medip-seq data, data3 and data4 are
    MRE-seq data).}
    \item{chrstring}{The chromosome should be test.}
    \item{cpgfile}{The file of all CpG number.}
    \item{mrecpgfile}{The file of MRE-CpG number(If NULL, mrecpgfile will equal
    to cpgfile).}
    \item{writefile}{The path of file of output result. (If writefile=NULL,
    there will return the results back to main program)}
    \item{reportfile}{The path of output results of the number of bin, total
    reads before processing and total reads after processing.}
    \item{mreratio}{The ratio of total unmethylation level with total
    methylation level (Defaulted mreratio is 3/7).}
    \item{psd}{The parameters of pseudo count, which pseudo count added to
    Medip-seq and MRE-seq count.}
    \item{mkadded}{Added to all CpG and MRE CpG (We set psd=2 and mkadded=1 as
    defaulted for robust)}
    \item{f}{Adjustment weight, default to 1.}
}
\seealso{?EMtest}

\value{
    \item{values or file}{The output file "writefile" will own eleven columns,
that is, "chr", "chrSt", "chrEnd", "Medip1", "Medip2", "MRE1", "MRE2",
"cg","mrecg","pvalue" and "Ts". We also output a report file which will
include parameters "s1/s2", "s3/s4", "N1", "N2", "N3", "N4", "c1",
"c2", "Number of windows" and "Spend time".}
}
\examples{
    data(example_data)
    data1<-EM2_H1ESB1_MeDIP_sigleCpG
    data2<-EM2_H1ESB2_MeDIP_sigleCpG
    data3<-H1ESB1_MRE_sigleCpG
    data4<-H1ESB2_MRE_sigleCpG
    datafile<-cbind(data1,data2,data3,data4)
    allcpg<-all_CpGsite_bin_chr18
    mrecpg<-three_mre_cpg
    dirwrite<-paste(setwd(getwd()),"/",sep="")

    writefile<-paste(dirwrite,"pval_EM_H1ESB1_H1ESB21.bed",sep="")
    reportfile<-paste(dirwrite,"report_pvalH1ESB1_H1ESB21.bed",sep="")

    EMtest(datafile=datafile,chrstring=NULL,cpgfile=allcpg,
    mrecpgfile=mrecpg,writefile=writefile,reportfile=reportfile,
    mreratio=3/7, psd=2,mkadded=1,f=1)

}






