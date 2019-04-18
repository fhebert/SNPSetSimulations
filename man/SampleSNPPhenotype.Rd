\name{SampleSNPPhenotype}
\alias{SampleSNPPhenotype}
\title{Simulation of a sample of genotype and binary phenotype values}
\usage{
SampleSNPPhenotype(X, beta0, beta, I, n0, n1, U = NULL, alpha = NULL,
  mod = rep("A", length(I)))
}
\arguments{
\item{X}{a N x p matrix of genotype values (coded as 0/1/2) for a given gene for a population}

\item{beta0}{the intercept of the logistic model}

\item{beta}{the vector of effect parameters for the logistic model}

\item{I}{the vector of coordinates corresponding to the coefficients in \code{beta}}

\item{n0}{the desired number of controls}

\item{n1}{the desired number of cases}

\item{U}{an optional matrix of covariates profiles}

\item{alpha}{an optional vector of covariates effects parameters}

\item{mod}{an optional vector of characters specifying the genetic association model for each SNP position
given in \code{I}. Possibles values: \code{"A"} (additive model), \code{"D"} (dominant) and \code{"R"}
(recessive). Models can be combined (causal SNPs can be associated to the phenotype with different models).}
}
\value{
A list containing 3 elements:
\item{SNP}{the sample matrix of genotypic profiles}
\item{Y}{the vector of phenotype values}
\item{Covariates}{the sample matrix of covariates profiles (if any)}
}
\description{
\code{SampleSNPPhenotype} is used to generate a sample of genotype values for a gene and of phenotype values.
The phenotype is generated conditionally to the gene using a logistic model.
}
