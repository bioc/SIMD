#' Inferring the methylation expression level of single sites.
#' 
#' @description Using statistical framework and EM algorithm to infer the
#' methylation expression level of single sites.
#' @param datafile The files of sample. (datafile should be cbind(data1,data2,
#' data3,data4), where data1 and data2 are Medip-seq data, data3 and data4 are
#' MRE-seq data).
#' @param chrstring The chromosome should be test.
#' @param cpgfile The file of all CpG number.
#' @param mrecpgfile The file of MRE-CpG number(If NULL, mrecpgfile will equal
#' to cpgfile).
#' @param writefile The path of file of output result. (If writefile=NULL,
#' there will return the results back to main program)
#' @param reportfile The path of output results of the number of bin, total
#' reads before processing and total reads after processing.
#' @param mreratio The ratio of total unmethylation level with total
#' methylation level (Defaulted mreratio is 3/7).
#' @param psd The parameters of pseudo count, which pseudo count added to
#' Medip-seq and MRE-seq count.
#' @param mkadded Added to all CpG and MRE CpG (We set psd=2 and mkadded=1 as
#' defaulted for robust).
#' @param f Adjustment weight, default to 1.
#' @return values or file The output file "writefile" will own eleven columns,
#' that is, "chr", "chrSt", "chrEnd", "Medip1", "Medip2", "MRE1", "MRE2",
#' "cg","mrecg","pvalue" and "Ts". We also output a report file which will
#' include parameters "s1/s2", "s3/s4", "N1", "N2", "N3", "N4", "c1",
#' "c2", "Number of windows" and "Spend time".
#' @examples data(example_data)
#' data1 <- EM2_H1ESB1_MeDIP_sigleCpG
#' data2 <- EM2_H1ESB2_MeDIP_sigleCpG
#' data3 <- H1ESB1_MRE_sigleCpG
#' data4 <- H1ESB2_MRE_sigleCpG
#' datafile <- cbind(data1, data2, data3, data4)
#' allcpg <- all_CpGsite_bin_chr18
#' mrecpg <- three_mre_cpg
#' dirwrite <- paste(setwd(getwd()), "/", sep="")
#' writefile <- paste(dirwrite, "pval_EM_H1ESB1_H1ESB21.bed", sep="")
#' reportfile <- paste(dirwrite, "report_pvalH1ESB1_H1ESB21.bed", sep="")
#' EMtest(datafile=datafile, chrstring=NULL, cpgfile=allcpg,
#'        mrecpgfile=mrecpg, writefile=writefile, reportfile=reportfile,
#'        mreratio=3/7, psd=2, mkadded=1, f=1)
#' @export



EMtest <- function(datafile=NULL, chrstring=NULL, cpgfile, 
                   mrecpgfile=NULL,writefile=NULL, reportfile=NULL,
                   mreratio=3/7, psd=2, mkadded=1, f=1){
    
    starttime <- proc.time()
    message("reading data ....")
    if (length(datafile) == 2) {
        message("The datafile is only Medip-seq data")
        data1 <- datafile[1]
        data2 <- datafile[2]
        data3 <- data1
        data4 <- data2
        data3[, 4] <- 0
        data4[, 4] <- 0
    } else {
    if (length(datafile)!=8 & length(datafile)!=16) {
        stop("The input MeDIP and MRE file should be two or four files")
    }
    data1 <- datafile[, seq_len(4)]
    data2 <- datafile[, seq(5, 8)]
    data3 <- datafile[, seq(9, 12)]
    data4 <- datafile[, seq(13, 16)]
    }
    
    dataset <- cbind(data1, data2[, 4], data3[, 4], data4[, 4])
    if (is.null(f)) f <- 1
    message("adjustment weight f is :")
    print(f)
    if (f == 0) stop("f should larger than 0")
    message("Top line of dataset:")
    print(dataset[1, ])
    dataset[, 4:7] <-  dataset[,4:7]/f
    colnames(dataset) <- c("V1","V2","V3","V4","V5","V6","V7")
    if (is.null(cpgfile)) {
        stop(" No CpG file")
    } else {
        cpg <- cpgfile
        print(cpg[1, ])
    }
    if (is.null(mrecpgfile)) {
        mrecpg <- cpg
        message("warning: MRE-CpG file is same as CpG file")
    } else {
        mrecpg <- mrecpgfile
        print(mrecpg[1, ])
    }
    cpg[, 4] <- cpg[, 4]/f
    mrecpg[, 4] <- mrecpg[, 4]/f
    
    if (is.null(mkadded)) mkadded <- 1
    cpg[mrecpg[,4] != 0, 4] <- cpg[mrecpg[, 4] != 0, 4] + mkadded
    mrecpg[mrecpg[, 4] != 0, 4] <- mrecpg[mrecpg[, 4] != 0, 4] + mkadded
    newdata <- dataset
    newdata[mrecpg[,4] != 0, 6:7] <- newdata[mrecpg[,4] != 0, 6:7]*
    (cpg[mrecpg[, 4] != 0, 4]/mrecpg[mrecpg[, 4] != 0, 4])
    NN1 <- sum(dataset[, 4])
    NN2 <- sum(dataset[, 5])
    NN3 <- sum(dataset[, 6])*(sum(cpg[, 4])/sum(mrecpg[, 4]))
    NN4 <- sum(dataset[, 7])*(sum(cpg[, 4])/sum(mrecpg[, 4]))
    if (mreratio >= 0) {
        if (NN3 == 0){
            cc1 <- 1
        } else {
            cc1 <- (NN1/NN3)*mreratio}
        if(NN4 == 0){
            cc2 <- 1
        } else {
            cc2 <- (NN2/NN4)*mreratio}
        } else {
            cc1 <- 1
            cc2 <- 1
    }
    newdata[, 6] <- round(newdata[, 6]*cc1)
    newdata[, 7] <- round(newdata[, 7]*cc2)
    dataset <- newdata

    if (is.null(psd)) psd <- 1
    dataset[mrecpg[, 4] != 0, 4:7] <- round(dataset[mrecpg[, 4] != 0, 4:7]+psd)
    dataset[mrecpg[,4] == 0, 4:5] <- round(dataset[mrecpg[, 4]==0, 4:5]+psd)
    Medip <- dataset[, 4:5]
    Medip <- as.matrix(Medip)
    fws1s2 <- edgeR::calcNormFactors(Medip, logratioTrim=0.05, sumTrim=0.05)[2]
    mre <- dataset[rowSums(dataset[, 6:7]) != 2*round(psd) & mrecpg[, 4]!=0,6:7]
    if (length(mre[, 1]) == 0){
        fws3s4 <- 1
    } else {
        mre <- as.matrix(mre)
        mno2 <- cpg[rowSums(dataset[, 6:7]) != 2*round(psd) & mrecpg[, 4]!=0,4]
        kno2 <- mrecpg[rowSums(dataset[, 6:7]) != 2*round(psd)& mrecpg[,4]!=0,4]
        fws3s4 <- edgeR::calcNormFactors(mre,logratioTrim=0.05, sumTrim=0.05)[2]
    }
    
    if (is.null(chrstring)){
        sm <- as.matrix(dataset[, 4:7])
        cpg <- cpg
        kk <- mrecpg
    } else {
        sm <- as.matrix(dataset[dataset[, 1] == chrstring, 4:7])
        cpg <- cpg[cpg[, 1] == chrstring, ]
        kk <- mrecpg[mrecpg[, 1] == chrstring, ]
    if (length(sm[, 1]) == 0) print("No such chrstring data")
    }
    pall <- rep(0, length(sm[, 1]))
    if(length(datafile) == 2) {
        message("Compute all p-values with SAGE.test")
        kk[, 4] <- 0
    }

    medipct <- sm[kk[, 4] == 0, seq_len(2)]
    cpgS = cpg[kk[, 4] == 0, 4] 
    message(paste("Computing p-values for only MeDIP-seq windows:",
    length(medipct[, 1])))
    pmedip <- statmod::sage.test(medipct[, 1], medipct[, 2],
    n1=sum(medipct[, 1])/sqrt(fws1s2), n2=sum(medipct[, 2])*sqrt(fws1s2))
    pall[kk[, 4] == 0] <- pmedip
    ts <- cpg[, 4]
    ns1 <- round(sum(medipct[, 1])/sqrt(fws1s2))
    ns2 <- round(sum(medipct[, 2])*sqrt(fws1s2))
    probs <- ns1/(ns1+ns2)
    cs <- fws1s2*(sum(medipct[, 2])/sum(medipct[, 1]))
    tts1 <- cpg[kk[, 4] == 0, 4]
    nn <- rowSums(medipct)
    tts1[nn == 0] <- 0

    tts1[nn != 0] <- sqrt(nn[nn != 0])*((cs*medipct[nn != 0, 1]/
    nn[nn != 0]-medipct[nn != 0, 2]/nn[nn != 0])-((cs+1)*probs-1))/
    (sqrt(probs*(1-probs))*(1+cs))
    ts[kk[, 4] == 0] <- tts1

    if(length(cpg[kk[, 4] != 0, 1]) != 0){
        cpg4 <- cpg[kk[, 4]!=0, ]
        sm1 <- sm[kk[, 4] != 0, ]
        kk1 <- kk[kk[, 4] != 0, ]
        m <- cpg4[, 4]
        k <- kk1[, 4]
        message(paste("Computing p-values for MeDIP-seq
                    and MRE-seq windows:", length(m)))
        sm1 <- round(sm1)
        x1 <- sm1[, 1]  
        x2 <- sm1[, 2]  
        y1 <- sm1[, 3]  
        y2 <- sm1[, 4]  

        N1 <- sum(x1)
        N2 <- sum(x2)
        N3 <- sum(y1)
        N4 <- sum(y2)
        rmd21 <- fws1s2*(N2/N1)
        rmd43 <- fws3s4*(N4/N3)

        tt <- abs(rmd21*x1*y2 - rmd43*x2*y1)
        size <- round(rowSums(sm1))
        tts3 <- (rmd21*x1*y2 - rmd43*x2*y1)/(size^2+0.1)

        ts[kk[, 4] != 0] <- tts3

        pvalue <- rep(2, length(tt))
        smT <- TRUE
        smsize <- size
        smsm1 <- sm1
        typesm1 <- unique(smsm1)
        typermd21 <- rmd21
        typermd43 <- rmd43
        typeT <- abs(typermd21*typesm1[, 1]*typesm1[, 4] -
        typermd43*typesm1[, 2]*typesm1[, 3])
        typesize <- round(rowSums(typesm1))
        size1 <- typesm1[, 1]+typesm1[, 2]
        size2 <- typesm1[, 3]+typesm1[, 4]
        message("Exact compute p-values with product binomial distribution.")
        typepvalue <- rep(2, length(typesize))
        for(ii in seq_len(length(typesize))){
            typepvalue[ii] <- probBinom(typeT[ii], size1[ii], size2[ii],
                                        typermd21, typermd43)
    }

    smpvalue <- classifypvalue(typesm1[, 1], typesm1[, 2], typesm1[, 3], 
                               typesm1[, 4], smsm1[, 1], smsm1[, 2],
                               smsm1[, 3], smsm1[, 4],typepvalue, 
                               length(typesm1[, 1]),length(smsm1[, 1]), 
                               pvalue=rep(0, length(smsm1[, 1])))

    pvalue <- smpvalue
    pall[kk[, 4] != 0] <- pvalue
    }

    endtime <- proc.time()
    spendtime <- endtime-starttime
    
    if(!is.null(reportfile)) {
       str <- paste("s1/s2:", fws1s2, 
                    "\n",
                    "s3/s4:", fws3s4
      )
       write(str, reportfile)
       if (length(cpg[kk[, 4] != 0, 1]) != 0) {
           str <- paste("N1(sum of MRE-CpG!=0 ):", N1,
                         "\n",
                         "N2(sum of MRE-CpG!=0 ):", N2,
                         "\n",
                         "N3(sum of MRE-CpG!=0 ):", N3,
                         "\n",
                         "N4(sum of MRE-CpG!=0 ):", N4,
                         "\n",
                         "c1:", rmd21,
                         "\n",
                         "c2:", rmd43)
           write(str, reportfile, append = TRUE)
      }
      str <- paste("Number of windows:",length(pall),
                   "\n",
                   "Number of p<1e-8:",sum(pall<1e-8),
                   "\n",
                   "Number of p<1e-3:",sum(pall<1e-3),
                   "\n",
                   "Spend time:", spendtime)
      write(str, reportfile, append = TRUE)
    }

    cpg[, 2] <- as.integer(cpg[, 2])
    cpg[, 3] <- as.integer(cpg[, 3])
    cpgpq <- cbind(cpg[, seq_len(3)], sm, cpg[, 4], kk[, 4], pall, ts)
    colnames(cpgpq) <- c("chr", "chrSt", "chrEnd", "Medip1", "Medip2", "MRE1",
                         "MRE2", "cg", "mrecg", "pvalue", 'Ts')
    if (!is.null(writefile)) {
        message("Output to writefile.")
        write.table(cpgpq, writefile,sep="\t", quote=FALSE, row.names=FALSE)
    } else {
        message("Return value.")
        return(cpgpq)  
    }
    message("The End.")
    }

