% Generated by roxygen2 (4.0.1): do not edit by hand
\name{logLikelihoodCalculator1s}
\alias{logLikelihoodCalculator1s}
\title{Function that given observed count data along a known haplotype returns a function that can calculate the likelihood of observing that data for a supplied underlying haplotype frequency.}
\usage{
logLikelihoodCalculator1s(lociHapACounts, lociTotalCounts, lociHapANoASEProbs,
  lociRhos, checkArgs = FALSE)
}
\arguments{
\item{lociHapACounts}{counts of haplotype A-supporting reads at individual loci. Must be a vector of non-negative integers.}

\item{lociTotalCounts}{total read counts of at individual loci. Must be a vector of positive integers.}

\item{lociHapANoASEProbs}{probabilities of observing haplotype A-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus). Must be a vector with entries >0 and <1.}

\item{lociRhos}{dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial). Must be a numeric vector with entries >=0 and <1.}

\item{checkArgs}{single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE}
}
\value{
a function of a single argument pHapA that calculates log likelihood of the observed data if the true underlying haplotype A frequency is pHapA.
}
\description{
Function that given observed count data along a known haplotype returns a function that can calculate the likelihood of observing that data for a supplied underlying haplotype frequency.
}
\details{
Given observed read counts supporting hapltoype A at a collection of loci, the total read counts at those loci, the probablities of observing haplotype A-supporting reads under conditions of no ASE and the dispersion parameters, this function returns a function of a single argument, pHapA, that calculates the likelihood of observing the given haplotype A-supporting counts under the assumption that the true underlying frequency of haplotype A is pHapA.
}
\examples{
LLC <- MBASED:::logLikelihoodCalculator1s(lociHapACounts=c(5, 12), lociTotalCounts=c(10, 24), lociHapANoASEProbs=c(0.5, 0.5), lociRhos=c(0,0))
LLC(0.5) ## the MLE estimate of hapA frequency
LLC(0.1) ## highly implausible value of pHapA
LLC (0.51)
}

