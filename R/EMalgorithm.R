#(2) EM algorithm
#' EM algorithm to infer CpG sites.
#' 
#' @description Using EM algorithm to infer the real number of CpG sites.
#' @param cpgsitefile The path of file to store CpG site.
#' @param allcpgfile The file to store CpG sites.
#' @param category Default to "1".
#' @param writefile The path of output results. (If writefile=NULL, there will
#' return the results back to main program.)
#' @param reportfile The path of output results.
#' @return values or file If writefile is NULL, then return the values of
#' results,otherwise output to write file.
#' @examples 
#' datafile<-system.file("extdata", package = "methylMnM")
#' data(example_data)
#' filepath<-datafile[1]
#' allcpgfile<-EM_H1ESB1_MeDIP_sigleCpG
#' dirwrite<-paste(setwd(getwd()),"/",sep="")
#' readshort<-paste(filepath,"/H1ESB1_Medip_18.extended.txt",sep="")
#' writefile<-paste(dirwrite,"EM2_H1ESB1_MeDIP_sigleCpG.bed",sep="")
#' reportfile<-paste(dirwrite,"EM2_H1ESB1_MeDIP_sigleCpG_report.bed",sep="")
#' f<-EMalgorithm(cpgsitefile=readshort,allcpgfile=allcpgfile,category="1",
#'               writefile=writefile,reportfile=reportfile)
#' @export

EMalgorithm<-function(cpgsitefile,allcpgfile,category="1",
writefile=NULL,reportfile=NULL)
{
    w=proc.time()
    if (is.null(allcpgfile)) {
        print("NO CpG sites file")
    }else {
    #allcpg <- read.table(allcpgfile, header = TRUE, sep = "\t",
    #                     as.is = TRUE)
        allcpg <- allcpgfile
    }
    if (is.null(cpgsitefile))
        print("Error, No cpgsitefile")
    nrec <- length(count.fields(cpgsitefile))
    tag_cpg1 <- read.table(cpgsitefile, header = FALSE,
    skip = 0,nrows = floor(nrec/4),as.is = TRUE)
    tag_cpg2 <- read.table(cpgsitefile, header = FALSE,
    skip = floor(nrec/4),nrows = floor(nrec/4), as.is = TRUE)
    tag_cpg3 <- read.table(cpgsitefile, header = FALSE,
    skip = 2 *floor(nrec/4),nrows = floor(nrec/4), as.is = TRUE)
    tag_cpg4 <- read.table(cpgsitefile, header = FALSE,
    skip = 3 *floor(nrec/4),nrows = nrec, as.is = TRUE)

    if (is.null(allcpg)) {
        chrtype1 <- unique(tag_cpg1[, 1])
        chrtype2 <- unique(tag_cpg2[, 1])
        chrtype3 <- unique(tag_cpg3[, 1])
        chrtype4 <- unique(tag_cpg4[, 1])
        chrstring <- unique(c(chrtype1, chrtype2, chrtype3, chrtype4))
    } else {
        chrstring <- unique(allcpg[, 1])
    }
    chrsizes <- length(chrstring)
    cpg <- matrix(0, 2, 4)
    for (i in seq_len(chrsizes)) {
        xx141 <- tag_cpg1[tag_cpg1[, 1] == chrstring[i], ]
        xx142 <- tag_cpg2[tag_cpg2[, 1] == chrstring[i], ]
        xx143 <- tag_cpg3[tag_cpg3[, 1] == chrstring[i], ]
        xx144 <- tag_cpg4[tag_cpg4[, 1] == chrstring[i], ]
        tag_cpgchr19 <- rbind(xx141, xx142, xx143, xx144)
        tag_cpgchr19 <- tag_cpgchr19[order(tag_cpgchr19[, 2]), ]
        cpg19 <- allcpg[allcpg[, 1] == chrstring[i], ]
        cpg19 <- cpg19[order(cpg19[, 2]), ]
        cpg20<-cpg19
        cpg20[,4]<-0
        zero<-cpg19[cpg19[,4]==0,]

        slid<-length(zero[,1])
        for(j in seq_len(slid)){
        if(j==1){
            region<-cpg19[cpg19[,2]<zero[j,2],]
        }else{
            if(j==slid){
            region<-cpg19[cpg19[,2]>zero[j,2],]
            }else{
            region<-cpg19[cpg19[,2]>zero[(j-1),2]&cpg19[,2]<zero[j,2],]
            }
        }
        num<-length(region[,1])
        if(num!=0){
        regiontag<-tag_cpgchr19[tag_cpgchr19[,2]<region[num,3]&
        tag_cpgchr19[,3]>region[1,2],seq_len(3)]
        tagm<-matrix(0,length(regiontag[,1]),num)
        for(r in seq_len(length(regiontag[,1]))){
            tagm[r,region[,2]>=regiontag[r,2]&region[,2]<regiontag[r,3]]<-1
        }
        tagm<-tagm[rowSums(tagm)!=0,]
        X1<-tagm+3
        Y<-tagm


        if(category=="1"){
            while(max(abs(X1-Y))>1e-6){
            X1<-Y
            if(is.matrix(tagm)){
                Y<-emalgth(tagm)
            }
            }
        }else{
            while(max(abs(X1-Y))>1e-6){
            X1<-Y
            if(is.matrix(tagm)){
                Y<-emalgth1(tagm)
            }
            }
        }


        if(is.matrix(tagm)){
            realnum<-apply(Y,2,sum)
        }else{
            if(num==1){
            realnum<-sum(Y)
            }else{
            realnum<-Y
            }
        }

        if(j==1){
            cpg20[cpg19[,2]<zero[j,2],4]<-realnum
        }else{
            if(j==slid){
            cpg20[cpg19[,2]>zero[j,2],4]<-realnum
            }else{
            cpg20[cpg19[,2]>zero[(j-1),2]&cpg19[,2]<zero[j,2],4]<-realnum
            }
        }
        }
        if(j %% 10000==0) print(j)
    }
    print(i)
    }

    d=proc.time()
    time=d-w
    if(!is.null(reportfile)) {
    write(paste("Number of CpG sites:", length(cpg20[,4])),
    reportfile, append = FALSE)
    write(paste("Total count of all CpG sites before EM
    algorithm:", sum(cpg19[, 4])),reportfile, append = TRUE)
    write(paste("Total count of all CpG sites after EM
    algorithm:", sum(cpg20[, 4])),reportfile, append = TRUE)
    write(paste("Spend time:", time), reportfile,append = TRUE)
    }

    if (is.null(writefile)) {
    return(cpg20)
    }else{
    write.table(cpg20, writefile, sep = "\t", quote = FALSE,
                row.names = FALSE)
    }

}
