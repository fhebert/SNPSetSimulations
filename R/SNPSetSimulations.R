#The GenOrd package (Barbiero and Ferrari, 2015) provides tools for simulation of discrete/ordinal distributions
#with given marginal distributions and correlation matrix. These simulations are obtained by discretizing
#multivariate normal distributions. An intermediate correlation matrix is used to sample from the multivariate normal
#distribution so that after discretization, the desired correlation matrix is obtained for the discrete/ordinal
#distributions.
#References: Barbiero and Ferrari, 2015, GenOrd package (available on CRAN)
#            Barbiero and Ferrari, 2012, Simulating ordinal data, Multivariate Behavioral Research

#ordcont2 and contord2 are modified versions of the "ordcont" and "contord" functions from the GenOrd package
#(Barbiero and Ferrari 2015). These functions are used to compute the intermediate correlation matrix of the
#multivariate normal distribution used in the simulation. The most important modification is the following: if
#the correlation matrix obtained at a step of the iterative algorithm is not positive-definite, is it replaced by
#the nearest positive-definite correlation matrix (see the "nearPD" function from the Matrix package).
ordcont2 = function(marginal, Sigma, support = list(),epsilon = 1e-06, maxit = 100){
  if (!all(unlist(lapply(marginal, function(x) (sort(x) ==
                                                x & min(x) > 0 & max(x) < 1)))))
    stop("Error in assigning marginal distributions!")
  if (!isSymmetric(Sigma) | !all(diag(Sigma) == 1)){
    stop("Sigma is not a valid correlation matrix")
  }
  if(min(eigen(Sigma)$values) < 0 ){
    Sigma = as.matrix(nearPD(Sigma,corr=TRUE,maxit=500)$mat)
  }
  len <- length(support)
  k <- length(marginal)
  niter <- matrix(0, k, k)
  kj <- numeric(k)
  for (i in 1:k) {
    kj[i] <- length(marginal[[i]]) + 1
    if (len == 0) {
      support[[i]] <- 1:kj[i]
    }
  }
  for (g in 2:k) {
    if (det(Sigma[1:g, 1:g]) <= 0) {
      Sigma = as.matrix(nearPD(Sigma,corr=TRUE,maxit=500)$mat)
    }
  }
  mcmin <- corrcheck(marginal, support, FALSE)[[1]]
  mcmax <- corrcheck(marginal, support, FALSE)[[2]]
  if (sum(mcmin <= Sigma & Sigma <= mcmax) != k^2) {
    stop("Some correlation coefficients are not feasible!",
         "\n", "Please use function corrcheck to get lower and upper bounds!",
         "\n")
  }
  Sigma0 <- Sigma
  Sigmaord <- Sigma
  Sigmaord <- contord2(marginal, Sigma, support)
  Sigmaold <- Sigma
  Sigmaordold <- Sigmaord
  for (q in 1:(k - 1)) {
    for (r in (q + 1):k) {
      if (Sigma0[q, r] == 0) {
        Sigma[q, r] <- 0
      }
      else {
        it <- 0
        while (max(abs(Sigmaord[q, r] - Sigma0[q, r]) >
                   epsilon) & it < maxit) {
          if (Sigma0[q, r] * (Sigma0[q, r]/Sigmaordold[q,
                                                       r]) >= 1) {
            Sigma[q, r] <- Sigmaold[q, r] * (1 + 0.1 *
                                               (1 - Sigmaold[q, r]) * sign(Sigma0[q,
                                                                                  r] - Sigmaord[q, r]))
          }
          else {
            Sigma[q, r] <- Sigmaold[q, r] * (Sigma0[q,
                                                    r]/Sigmaord[q, r])
          }
          Sigma[r, q] <- Sigma[q, r]
          Sigmaord[r, q] <- contord2(list(marginal[[q]],marginal[[r]]),
                                     matrix(c(1, Sigma[q, r],Sigma[q, r], 1), 2, 2),
                                     list(support[[q]],support[[r]]))[2]
          Sigmaord[q, r] <- Sigmaord[r, q]
          Sigmaold[q, r] <- Sigma[q, r]
          Sigmaold[r, q] <- Sigmaold[q, r]
          it <- it + 1
        }
        niter[q, r] <- it
        niter[r, q] <- it
      }
    }
  }
  if (eigen(Sigma)$values[k] <= 0) {
    Sigma = as.matrix(nearPD(Sigma,corr=TRUE,maxit=500)$mat)
  }
  emax <- max(abs(Sigmaord - Sigma0))
  list(SigmaC = Sigma, SigmaO = Sigmaord, Sigma = Sigma0,
       niter = niter, maxerr = emax)
}
contord2 = function(marginal, Sigma, support = list()){
  if (!all(unlist(lapply(marginal, function(x) (sort(x) ==
                                                x & min(x) > 0 & max(x) < 1)))))
    stop("Error in assigning marginal distributions!")
  if (!isSymmetric(Sigma) | min(eigen(Sigma)$values) < 0 |
      !all(diag(Sigma) == 1))
    Sigma = as.matrix(nearPD(Sigma,corr=TRUE,maxit=500)$mat)
  len <- length(support)
  k <- length(marginal)
  kj <- numeric(k)
  Sigmaord <- diag(k)
  for (i in 1:k) {
    kj[i] <- length(marginal[[i]]) + 1
    if (len == 0) {
      support[[i]] <- 1:kj[i]
    }
    # if (Spearman) {
    #   s1 <- c(marginal[[i]], 1)
    #   s2 <- c(0, marginal[[i]])
    #   support[[i]] <- (s1 + s2)/2
    # }
  }
  L <- vector("list", k)
  for (i in 1:k) {
    L[[i]] <- qnorm(marginal[[i]])
    L[[i]] <- c(-Inf, L[[i]], +Inf)
  }
  for (q in 1:(k - 1)) {
    for (r in (q + 1):k) {
      pij <- matrix(0, kj[q], kj[r])
      for (i in 1:kj[q]) {
        for (j in 1:kj[r]) {
          low <- rep(-Inf, k)
          upp <- rep(Inf, k)
          low[q] <- L[[q]][i]
          low[r] <- L[[r]][j]
          upp[q] <- L[[q]][i + 1]
          upp[r] <- L[[r]][j + 1]
          pij[i, j] <- pmvnorm(low, upp, rep(0, k),
                               corr = Sigma)
          low <- rep(-Inf, k)
          upp <- rep(Inf, k)
        }
      }
      my <- sum(apply(pij, 2, sum) * support[[r]])
      sigmay <- sqrt(sum(apply(pij, 2, sum) * support[[r]]^2) -
                       my^2)
      mx <- sum(apply(pij, 1, sum) * support[[q]])
      sigmax <- sqrt(sum(apply(pij, 1, sum) * support[[q]]^2) -
                       mx^2)
      mij <- support[[q]] %*% t(support[[r]])
      muij <- sum(mij * pij)
      covxy <- muij - mx * my
      corxy <- covxy/(sigmax * sigmay)
      Sigmaord[q, r] <- corxy
    }
  }
  as.matrix(forceSymmetric(Sigmaord))
}

#Generates a binary phenotype conditionally to genotype values of a population of n individuals using
#a logistic model. The matrix of genotype values is given in the "X" argument. For a vector x of genotype values,
#the from a Bernoulli distribution with probability P[Y=1|X=x] = 1/(1+exp(-beta0-x%*%beta)), beta0 being an
#intercept term and beta the vector of association paramters. In the "beta" argument, only the non-zero coefficients
#must be given. The indices of the causative SNPs (corresponding to non-zero coefficients) are given in the "I"
#argument. The "mod" argument specifies an association model for each causative SNP ("A" for an additive association,
#"R" for a recessive association and "D" for a dominant association).
PopulationPhenotype = function(X,beta0,beta,I,U=NULL,alpha=NULL,mod=rep("A",length(I))){
  if(!is.null(U)){
    Ua = U%*%alpha
  }else{
    Ua = 0
  }
  if(any(mod=="R")){
    indR = I[which(mod=="R")]
    XR = X[,indR]
    ind0 = which(XR<2)
    XR[ind0] = 0
    XR[-ind0] = 1
    X[,indR] = XR
  }
  if(any(mod=="D")){
    indD = I[which(mod=="D")]
    XD = X[,indD]
    ind0 = which(XD<1)
    XD[ind0] = 0
    XD[-ind0] = 1
    X[,indD] = XD
  }
  Y = beta0+Ua+X[,I,drop=FALSE]%*%beta
  Y = 1/(1+exp(-Y))
  Y = rbinom(length(Y),1,Y)
  return(Y)
}




