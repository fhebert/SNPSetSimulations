PopulationSNPSet = function(n,Sigma=NULL,maf=NULL,marginal=NULL,X=NULL){
  if(!is.null(X)){
    m = ncol(X)
    X = as.matrix(X)
    Xvals = unique(X)
    if(!(all(Xvals%in%c(0,1,2)))){
      stop("Error: allowed values for X are 0/1/2")
    }
    sdX = apply(X,2,sd)
    if(any(sdX==0)){
      ind = which(sdX==0)
      X = X[,-ind]
      warning("Constant SNP variables were removed")
    }
    p0 = colMeans(X==0)
    p1 = colMeans(X==1)
    Mp = cbind(p0,p0+p1)
    marginal = lapply(seq_len(nrow(Mp)),function(i){
      if(Mp[i,2]==1){
        return(Mp[i,-2])
      }else{
        return(Mp[i,])
      }
    })
    Sigma = as.matrix(nearPD(cor(X),corr=TRUE)$mat)
    support = lapply(1:ncol(X),function(i){sort(unique(X[,i]))})
  }else{
    m = ncol(Sigma)
    if(is.null(marginal)){
      Mp = cbind((1-maf)^2,(1-maf)^2+2*maf*(1-maf))
      marginal = lapply(seq_len(nrow(Mp)),function(i){Mp[i,]})
    }
  }
  SigmaZ = ordcont2(marginal=marginal,Sigma=Sigma)$SigmaC
  C = chol(SigmaZ)
  matZ = matrix(rnorm(m*n),nrow=n)%*%C
  X = matrix(0,ncol=m,nrow=n)
  for(j in 1:ncol(matZ)){
    X[,j] = as.integer(cut(matZ[,j],breaks=c(min(matZ[,j])-1,qnorm(marginal[[j]]),max(matZ[,j])+1)))-1
  }
  return(X)
}
