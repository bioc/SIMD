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
