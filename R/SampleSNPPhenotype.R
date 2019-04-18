SampleSNPPhenotype = function(X,beta0,beta,I,n0,n1,U=NULL,alpha=NULL,mod=rep("A",length(I))){
  Y = PopulationPhenotype(X,beta0,beta,I,U,alpha,mod)
  i = which(Y==0)
  ind = sample(c(sample(i,n0,FALSE),sample((1:length(Y))[-i],n1,FALSE)))
  XX = X[ind,]
  Y = matrix(Y[ind],ncol=1)
  if(!is.null(U)){
    UU = U[ind,]
  }else{
    UU = NULL
  }
  return(list(SNP=XX,Phenotype=Y,Covariates=UU))
}
