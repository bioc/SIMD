#sum=1
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
