\name{SNPSetSimulations-package}
\alias{SNPSetSimulations-package}
\alias{SNPSetSimulations}
\docType{package}
\encoding{UTF-8}
\title{
  Simulation of genotypic profiles and binary phenotypes for GWAS
}
\description{
  The SNPSetSimulations package provides functions to generate genotypic profiles according to given correlation
  structures, and to generate binary phenotypes conditionally to genotypic profiles using a logistic model.
}
\author{
Florian Hébert, Mathieu Emily, David Causeur

Maintainer: Florian Hébert <florian.hebert@agrocampus-ouest.fr>
}
\keyword{Genome-wide association studies}
\seealso{
  \code{\link{PopulationSNPSet}}
  \code{\link{SampleSNPPhenotype}}
  \code{\link{levelplotSNP}}
}
\examples{
#Simulation of 100,000 genotypic profiles of 10 autocorrelated SNPs,
#with correlation coefficient 0.8. All SNPs have a minor allele frequency equal to 0.4.
X = PopulationSNPSet(100000,Sigma=0.8^abs(outer(1:10,1:10,"-")),maf=rep(0.4,10))

\dontrun{
#Assuming XX is an observed matrix of genotypic profiles, simulation of 100,000
#genotypic profiles according to the correlation matrix and marginal distributions observed in XX:
X = PopulationSNPSet(100000,X=XX)
}

#Generation of a sample of 1,000 cases and 1,000 controls conditionally to the genotypic
#profiles contained in X. The intercept of the logistic model is set to -3, the phenotype
#is associated to SNPs 3 and 7 (parameters 0.2 and 0.4, respectively) under an additive
#and a recessive models, respectively.
tmp = SampleSNPPhenotype(X,-3,c(0.2,0.4),c(3,7),1000,1000,mod=c("A","R"))

#Correlation matrix of the SNPs, estimated on the obtained sample:
levelplotSNP(cor(tmp$SNP))

#Chi-square test between the simulated phenotype and each simulated SNP variable:
apply(tmp$SNP,2,function(x){chisq.test(x,tmp$Phenotype,correct=FALSE)$p.value})
}
