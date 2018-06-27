#sum>=1
emalgth1<-function(X)
{
    N<-apply(X,2,sum)/length(X[,1])
    p<-1-N
    p1<-1-prod(p)
    N<-N*p1
    y<-X
    for(i in seq_len(length(N))){
    y[y[,i]!=0,i]<-N[i]
    }
    N1<-1-y
    N2<-apply(N1,1,prod)
    CpGnum<-rep(1,length(X[1,]))
    N11<-1-N2%o%CpGnum
    y1<-y/N11
    return(y1)
}
