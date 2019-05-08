levelplotSNP = function(Sigma){
  m = as.numeric(ncol(Sigma))
  a = m
  k = c(1:5,10,20,30,40,50,60,70,80,90,100,200,500,1000,5000,10000)
  i = 1
  while(a>10){
    i = i+1
    a = m%/%k[i]
  }
  SNPlabels = seq(0,m+k[i]-1,by=k[i])
  xat = SNPlabels
  if(length(SNPlabels<=10)){
    xat = xat[-1]
  }else{
    xat = xat[-c(1,length(xat))]
  }
  xlab = xat
  yat = (m:1)[xat]
  ylab = xlab
  pal = colorRampPalette(c("blue","chartreuse3","yellow","yellow","red","black"))
  levelplot((as.matrix(as.dist(Sigma))+diag(diag(Sigma)))[1:m,m:1],at=seq(-1,1,by=0.01),
            scales=list(x=list(at=xat,labels=xlab),y=list(at=yat,labels=ylab),tck=c(1,0)),
            xlab="",ylab="",col.regions=pal(200),
            colorkey=list(at=seq(-1,1,0.01),labels=list(at=seq(-1,1,0.5),labels=seq(-1,1,by=0.5))),
            useRaster=TRUE)
}
