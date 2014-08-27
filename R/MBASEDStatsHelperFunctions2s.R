#' Helper function to adjust proportions for pre-existing allelic bias and also to obtain estimate of proportion variance based on attenuated read counts (adding pseudocount of 0.5 to each allele in each sample).
#'
#' @param countsMat matrix of observed major allele counts. Each row represents a specific genomic locus, while each column represents a set of observed major allele counts across loci (in practice, multiple columns represent different outcomes of count simulations).
#' @param totalsMat matrix of total read counts across both alleles. The interpretation of rows and columns is the same as for countsMat.
#' @param probsMat matrix of underlying probabilites of observing the major allele. The interpretation of rows and columns is the same as for countsMat.
#' @param rhosMat matrix of dispersion parameters of beta distributions for each locus. The interpretation of rows and columns is the same as for countsMat.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with 2 elements:
#' \item{propsShifted}{a 1-row marix of shifted major allele frequencies}
#' \item{propsShiftedVars}{a 1-row matrix of estimated variances of obtained MAF estimates}
#'
#' @examples
#' SNVCoverageTumor=sample(10:100,10) ## 2 genes with 5 loci each
#' SNVAllele1CountsTumor=rbinom(length(SNVCoverageTumor), SNVCoverageTumor, 0.5)
#' MBASED:::shiftAndAttenuateProportions(countsMat=matrix(SNVAllele1CountsTumor, ncol=2), totalsMat=matrix(SNVCoverageTumor, ncol=2), probsMat=matrix(rep(0.5, length(SNVCoverageTumor)), ncol=2),  rhosMat=matrix(rep(0, length(SNVCoverageTumor)), ncol=2)) 
#' 
#'
shiftAndAttenuateProportions <- function (
    countsMat, 
    totalsMat, 
    probsMat, 
    rhosMat,
    checkArgs=FALSE
) {
    if (checkArgs) {
        if (
            !is.matrix(totalsMat) || 
            any(is.na(totalsMat)) || 
            !is.numeric(totalsMat) || 
            !isTRUE(all.equal(totalsMat, round(totalsMat))) 
            || any(totalsMat<1)
        ) {
            stop('MBASED:shiftAndAttenuateProportions: argument totalsMat must be a matrix of positive integers')
        }
        if (
            !is.matrix(countsMat) || 
            any(is.na(countsMat)) || 
            !is.numeric(countsMat) || 
            !isTRUE(all.equal(countsMat, round(countsMat))) || 
            any(countsMat<0)
        ) {
            stop('MBASED:shiftAndAttenuateProportions: argument countsMat must be a matrix of non-negative integers')
        }
        if (
            !is.matrix(probsMat) || 
            any(is.na(probsMat)) || 
            !is.numeric(probsMat) || 
            any(probsMat<=0) || 
            any(probsMat>=1)
        ) {
            stop('MBASED:shiftAndAttenuateProportions: argument probsMat must be a matrix with entries >0 and <1')
        }
        if (
            !is.matrix(rhosMat) || 
            any(is.na(rhosMat)) || 
            !is.numeric(rhosMat) || 
            any(rhosMat<0) || 
            any(rhosMat>=1)
        ) {
            stop('MBASED:shiftAndAttenuateProportions: argument rhosMat must be a matrix with entries >=0 and <1')
        }
        if (
            !isTRUE(all.equal(dim(countsMat), dim(totalsMat))) ||  
            !isTRUE(all.equal(dim(probsMat), dim(totalsMat))) || 
            !isTRUE(all.equal(dim(rhosMat), dim(totalsMat)))    
        ) {
            stop('MBASED:shiftAndAttenuateProportions: arguments countsMat, totalsMat, probsMat, and rhosMat must have same dimension')
        }
        if (any(countsMat>totalsMat)) {
            stop('MBASED:shiftAndAttenuateProportions: entries of argument countsMat must be <= corresponding entries in  argument totalsMatSample1')
        }
    }
    propsShifted <- unFT(
        z=FTAdjust(
            x=countsMat, 
            n=totalsMat, 
            p=probsMat
        ), 
        n=totalsMat
    )
    ## calculate variances
    totalsMatAttenuated <- (totalsMat+1)
    propsShiftedAttenuated <- (
        0.5+round(propsShifted*totalsMat)
     )/totalsMatAttenuated
    propsShiftedVars <- propsShiftedAttenuated*
        (1-propsShiftedAttenuated)*
        (rhosMat+(1-rhosMat)/totalsMatAttenuated)
    return(
        list(
            propsShifted=propsShifted,
            propsShiftedVars=propsShiftedVars
        )
    )
}


#' Vectorized wrapper around a test for difference of 2 proportions.
#'
#' @details MBASEDVectorizedPropDiffTest implements meta-analysis-like apporoach using proportion differences at each locus as variables to be aggregated. Input matrices countsMatSample1, totalsMatSample1, countsMatSample2, totalsMatSample2, probsMatSample1, probsMatSample2, rhosMatSample1, and rhosMatSample2 have 1 column for each set of loci ('independent studies') to be combined, with each row corresponding to an individual locus. MBASEDVectorizedPropDiffTest uses meta analysis approach by transforming counts at each locus into proportions and combininig the proportion differences (between sample1 and sample2) using the inverse-variance weighted schema. The function reports proportion difference estimates, corresponding standard errors, z-values (based on expected value of 0 under the null hypothesis of overall difference of 0) , and corresponding p-values based on normal distribution assumption of z-values, where alternative hypothesis of 'two.sided', 'greater', and 'less' can be specified, with the latter two specified w.r.t. 0. Adjustment for pre-existing allelic bias is performed by taking observed proportion in each sample, transforming it with FT transformation, adjusting for allelic bias as in 1-sample case and back-transforming to get a shifted proportion. The shifted proportion is then used to estimate its variance. The function is used to cacluate p-values in ASE settings, where countsMatSample1 and countsMatSample2 represent major allele counts in sample1 and sample2, respectively, and totalsMatSample1 and totalsMatSample2 represent total allele counts. Matrices probsMatSample1 and probsMatSample2 capture the pre-existing allelic bias by supplying the underlying probabilities of observing alleles currently specified as major in absence of any allele-specific expression, and rhosMatSample1 and rhosMatSample2 provide values of dispersion parameter for beta-binomial counts (0, in case of binomial model) for individual loci within each sample. 
#'
#' @param countsMatSample1,countsMatSample2 matrices of observed major allele counts in sample1 and sample2, respectively. Each row represents a specific genomic locus, while each column represents a set of observed major allele counts across loci (in practice, multiple columns represent different outcomes of count simulations).
#' @param totalsMatSample1,totalsMatSample2 matrices of total read counts across both alleles in sample1 and sample2, respectively. The interpretation of rows and columns is the same as for countsMatSample1.
#' @param probsMatSample1,probsMatSample2 matrices of underlying probabilites of observing the major allele in sample1 and sample2, respectively. The interpretation of rows and columns is the same as for countsMatSample1.
#' @param rhosMatSample1,rhosMatSample2 matrices of dispersion parameters of beta distributions for each locus in sample1 and sample2, respectively. The interpretation of rows and columns is the same as for countsMatSample1.
#' @param alternative one of 'two.sided', 'greater', 'less'.  DEFAULT: 'two.sided'
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with 7 elements:
#' \item{hetPval}{a 1-row marix of heterogeneity P-values}
#' \item{hetQ}{a 1-row matrix of heterogeneity statistics}
#' \item{TEFinal}{a 1-row matrix of estimated proportion differences}
#' \item{seTEFinal}{a 1-row matrix of estimated SEs of prop differences estimates }
#' \item{propDifferenceFinal}{a 1-row matrix of estimated proportion differences}
#' \item{pValue}{a 1-row matrix of corresponding p-values.}
#' \item{propDifferenceLoci}{a matrix of same dimension as original input matrices giving estimated proportion differences on transformed scale at each individual locus.}
#'
#' @family vectorizedTests
#'
#' @examples
#' SNVCoverageTumor=sample(10:100,10) ## 2 genes with 5 loci each
#' SNVCoverageNormal=sample(10:100,10) ## 2 genes with 5 loci each
#' SNVAllele1CountsTumor=rbinom(length(SNVCoverageTumor), SNVCoverageTumor, 0.5)
#' SNVAllele1CountsNormal=rbinom(length(SNVCoverageNormal), SNVCoverageNormal, 0.5)
#' MBASED:::MBASEDVectorizedPropDiffTest(countsMatSample1=matrix(SNVAllele1CountsTumor, ncol=2), totalsMatSample1=matrix(SNVCoverageTumor, ncol=2), countsMatSample2=matrix(SNVAllele1CountsNormal, ncol=2), totalsMatSample2=matrix(SNVCoverageNormal, ncol=2), probsMatSample1=matrix(rep(0.5, length(SNVCoverageTumor)), ncol=2), probsMatSample2=matrix(rep(0.5, length(SNVCoverageNormal)), ncol=2), rhosMatSample1=matrix(rep(0, length(SNVCoverageTumor)), ncol=2), rhosMatSample2=matrix(rep(0, length(SNVCoverageNormal)), ncol=2), alternative='two.sided') 
#'
MBASEDVectorizedPropDiffTest <- function (
    countsMatSample1, 
    totalsMatSample1, 
    countsMatSample2, 
    totalsMatSample2, 
    probsMatSample1, 
    probsMatSample2, 
    rhosMatSample1, 
    rhosMatSample2, 
    alternative='two.sided', 
    checkArgs=FALSE
) { 
	## rows are loci, columns are simulations 
	## (in case of no simulations, one column each matrix)
    if (checkArgs) {
        if (
            !is.matrix(totalsMatSample1) || 
            any(is.na(totalsMatSample1)) || 
            !is.numeric(totalsMatSample1) || 
            !isTRUE(all.equal(totalsMatSample1, round(totalsMatSample1))) 
            || any(totalsMatSample1<1)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument totalsMatSample1 must be a matrix of positive integers')
        }
        if (
            !is.matrix(totalsMatSample2) || 
            any(is.na(totalsMatSample2)) || 
            !is.numeric(totalsMatSample2) || 
            !isTRUE(all.equal(totalsMatSample2, round(totalsMatSample2))) || 
            any(totalsMatSample2<1)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument totalsMatSample2 must be a matrix of positive integers')
        }
        if (
            !is.matrix(countsMatSample1) || 
            any(is.na(countsMatSample1)) || 
            !is.numeric(countsMatSample1) || 
            !isTRUE(all.equal(countsMatSample1, round(countsMatSample1))) || 
            any(countsMatSample1<0)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument countsMatSample1 must be a matrix of non-negative integers')
        }
        if (
            !is.matrix(countsMatSample2) || 
            any(is.na(countsMatSample2)) || 
            !is.numeric(countsMatSample2) || 
            !isTRUE(all.equal(countsMatSample2, round(countsMatSample2))) || 
            any(countsMatSample2<0)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument countsMatSample2 must be a matrix of non-negative integers')
        }
        if (
            !is.matrix(probsMatSample1) || 
            any(is.na(probsMatSample1)) || 
            !is.numeric(probsMatSample1) || 
            any(probsMatSample1<=0) || 
            any(probsMatSample1>=1)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument probsMatSample1 must be a matrix with entries >0 and <1')
        }
        if (
            !is.matrix(probsMatSample2) || 
            any(is.na(probsMatSample2)) || 
            !is.numeric(probsMatSample2) || 
            any(probsMatSample2<=0) || 
            any(probsMatSample2>=1)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument probsMatSample2 must be a matrix with entries >0 and <1')
        }
        if (
            !is.matrix(rhosMatSample1) || 
            any(is.na(rhosMatSample1)) || 
            !is.numeric(rhosMatSample1) || 
            any(rhosMatSample1<0) || 
            any(rhosMatSample1>=1)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument rhosMatSample1 must be a matrix with entries >=0 and <1')
        }
        if (
            !is.matrix(rhosMatSample2) || 
            any(is.na(rhosMatSample2)) || 
            !is.numeric(rhosMatSample2) || 
            any(rhosMatSample2<0) || 
            any(rhosMatSample2>=1)
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: argument rhosMatSample2 must be a matrix with entries >=0 and <1')
        }
        if (
            !isTRUE(all.equal(dim(countsMatSample1), dim(totalsMatSample1))) ||  
            !isTRUE(all.equal(dim(probsMatSample1), dim(totalsMatSample1))) || 
            !isTRUE(all.equal(dim(rhosMatSample1), dim(totalsMatSample1))) || 
            !isTRUE(all.equal(dim(countsMatSample2), dim(totalsMatSample1))) || 
            !isTRUE(all.equal(dim(totalsMatSample2), dim(totalsMatSample1))) || 
            !isTRUE(all.equal(dim(probsMatSample2), dim(totalsMatSample1))) || 
            !isTRUE(all.equal(dim(rhosMatSample2), dim(totalsMatSample1)))     
        ) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: arguments countsMatSample1, countsMatSample2, totalsMatSample1, totalsMatSample2, probsMatSample1, probsMatSample2, rhosMatSample1, and rhosMatSample2 must have same dimension')
        }
        if (any(countsMatSample1>totalsMatSample1)) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: entries of argument countsMatSample1 must be <= corresponding entries in  argument totalsMatSample1')
        }
        if (any(countsMatSample2>totalsMatSample2)) {
            stop('MBASED:MBASEDVectorizedPropDiffTest: entries of argument countsMatSample2 must be <= corresponding entries in  argument totalsMatSample2')
        }
        if ( 
            !(is.vector(alternative)) || 
            length(alternative)!=1 || 
            !(alternative%in%c('greater','less', 'two.sided')) 
        ) {
            stop("MBASED:MBASEDVectorizedPropDiffTest: argument alternative must be a single character string with a value of 'greater', 'less', or 'equal' ")
        }
    }
    ## get proportions after adjusting for prior allelic bias
    propsShiftedAndAttenuatedSample1 <- shiftAndAttenuateProportions(
        countsMat=countsMatSample1, 
        totalsMat=totalsMatSample1, 
        probsMat=probsMatSample1, 
        rhosMat=rhosMatSample1
    ) 
    propsShiftedAndAttenuatedSample2 <- shiftAndAttenuateProportions(
        countsMat=countsMatSample2, 
        totalsMat=totalsMatSample2, 
        probsMat=probsMatSample2, 
        rhosMat=rhosMatSample2
    ) 

    propDifferences <- (
        propsShiftedAndAttenuatedSample1$propsShifted-
        propsShiftedAndAttenuatedSample2$propsShifted
    )
    propDifferencesVars <- (
        propsShiftedAndAttenuatedSample1$propsShiftedVars+
        propsShiftedAndAttenuatedSample2$propsShiftedVars
    )
    metaAnalysisResults <- MBASEDMetaAnalysis(
        zValuesMat=propDifferences, 
        zVariancesMat=propDifferencesVars, 
        alternative=alternative
    ) 
    return(
        list(
            hetPVal=metaAnalysisResults$hetPVal, 
            hetQ=metaAnalysisResults$hetQ, 
            TEFinal=metaAnalysisResults$fixedEffectsMeans, 
            seTEFinal=metaAnalysisResults$fixedEffectsSEs, 
            ## propDifference is the z-statistic 
            propDifferenceFinal=metaAnalysisResults$fixedEffectsMeans, 
            pValue=metaAnalysisResults$pValsFixed,
            propDifferenceLoci=propDifferences
        )
    )
} 

#' Function that adjusts true underlying allele frequency for pre-existing allelic bias to produce actual generating probability of observing allele-supporting read
#'
#' @details Given true underlying allele frequency AF and probability of observing reads supporting that allele under conditiosn of no ASE (P(allele, noASE)), it calculates the generating probability for observed allele-supporting reads as P(allele-supporting read)=AF*P(allele, noASE)/(AF*P(allele, noASE) + (1-AF)*(1-P(allele, noASE))). 
#'
#' @param trueAF true underlying allele frequency. Must be a single number >=0 and <=1.
#' @param noASEAF probability of observing allele-supporting read under conditions of no ASE. Must be a vector of numbers >0 and <1.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a vector of generating probabilities of the same length as noASEAF
#'
#' @examples
#' MBASED:::getPFinal(trueAF=1, noASEAF=seq(0.1, 0.9, by=0.1)) ## is always 1
#' MBASED:::getPFinal(trueAF=0, noASEAF=seq(0.1, 0.9, by=0.1)) ## is always 0
#' MBASED:::getPFinal(trueAF=0.3, noASEAF=0.5) ## no pre-existing allelic bias
#' c(MBASED:::getPFinal(trueAF=0.3, noASEAF=0.9), MBASED:::getPFinal(trueAF=1-0.3, noASEAF=1-0.9)) ## strong pre-existing allelic bias
#'
getPFinal <- function (
    trueAF, 
    noASEAF, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if (
            !is.vector(trueAF) || 
            length(trueAF)!=1 || 
            is.na(trueAF) || 
            !is.numeric(trueAF) || 
            trueAF<0 || 
            trueAF>1
        ) {
            stop('MBASED:getPFinal : argument trueAF must be a single number >=0 and <=1')
        }
        if (
            !is.vector(noASEAF) || 
            any(is.na(noASEAF)) || 
            !is.numeric(noASEAF) || 
            any(noASEAF<=0) || 
            any(noASEAF>=1)
        ) {
            stop('MBASED:getPFinal : argument noASEAF must be a vector with entries >0 and <1')
        }
    }
    return(
        trueAF*noASEAF/(trueAF*noASEAF+(1-trueAF)*(1-noASEAF))
    )
}

#' Function that given observed count data along a known haplotype returns a function that can calculate the likelihood of observing that data for a supplied underlying haplotype frequency.
#'
#' @details Given observed read counts supporting hapltoype A at a collection of loci, the total read counts at those loci, the probablities of observing haplotype A-supporting reads under conditions of no ASE and the dispersion parameters, this function returns a function of a single argument, pHapA, that calculates the likelihood of observing the given haplotype A-supporting counts under the assumption that the true underlying frequency of haplotype A is pHapA.
#'
#' @param lociHapACounts counts of haplotype A-supporting reads at individual loci. Must be a vector of non-negative integers.
#' @param lociTotalCounts total read counts of at individual loci. Must be a vector of positive integers.
#' @param lociHapANoASEProbs probabilities of observing haplotype A-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus). Must be a vector with entries >0 and <1.
#' @param lociRhos dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial). Must be a numeric vector with entries >=0 and <1.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a function of a single argument pHapA that calculates log likelihood of the observed data if the true underlying haplotype A frequency is pHapA.
#'
#' @examples
#' LLC <- MBASED:::logLikelihoodCalculator1s(lociHapACounts=c(5, 12), lociTotalCounts=c(10, 24), lociHapANoASEProbs=c(0.5, 0.5), lociRhos=c(0,0)) 
#' LLC(0.5) ## the MLE estimate of hapA frequency
#' LLC(0.1) ## highly implausible value of pHapA
#' LLC (0.51)
#'
logLikelihoodCalculator1s <- function (
    lociHapACounts, 
    lociTotalCounts, 
    lociHapANoASEProbs, 
    lociRhos, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(lociHapACounts) || 
            any(is.na(lociHapACounts)) || 
            !is.numeric(lociHapACounts) || 
            !isTRUE(all.equal(lociHapACounts, round(lociHapACounts))) || 
            any(lociHapACounts<0) 
        ) {
            stop('MBASED:logLikelihoodCalculator1s: argument lociHapACounts must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociTotalCounts) || 
            any(is.na(lociTotalCounts)) || 
            !is.numeric(lociTotalCounts) || 
            !isTRUE(all.equal(lociTotalCounts, round(lociTotalCounts))) || 
            any(lociTotalCounts<=0) 
        ) {
            stop('MBASED:logLikelihoodCalculator1s: argument lociTotalCounts must be a vector of positive integers')
        }
        if (
            !is.vector(lociHapANoASEProbs) || 
            any(is.na(lociHapANoASEProbs)) || 
            !is.numeric(lociHapANoASEProbs) || 
            any(lociHapANoASEProbs<=0) || 
            any(lociHapANoASEProbs>=1)
        ) {
            stop('MBASED:logLikelihoodCalculator1s: argument lociHapANoASEProbs must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhos) || 
            any(is.na(lociRhos)) || 
            !is.numeric(lociRhos) || 
            any(lociRhos<0) || 
            any(lociRhos>=1)
        ) {
            stop('MBASED:logLikelihoodCalculator1s: argument lociRhos must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociHapACounts), 
                length(lociTotalCounts), 
                length(lociHapANoASEProbs), 
                length(lociRhos)
            )))>1
        ) {
            stop('MBASED:logLikelihoodCalculator1s: arguments lociHapACounts, lociTotalCounts, lociHapANoASEProbs, and lociRhos must be of same length')
        }
        if (any(lociHapACounts>lociTotalCounts)) {
            stop('MBASED:logLikelihoodCalculator1s: entries of argument lociHapACounts must be less than or equal to corresponding entries in argument lociTotalCounts')
        }
    }
    lociHapBCounts <- (lociTotalCounts-lociHapACounts)
    returnFunc <- function(pHapA) {
        finalPHapA <- getPFinal(
            trueAF=pHapA, 
            noASEAF=lociHapANoASEProbs
        )
        finalPHapB <- (1-finalPHapA)
        ABhapA <- getAB(
            mu=finalPHapA, 
            rho=lociRhos
        )
        logLikes <- ifelse(
            lociRhos==0,
            ## binomial
            log(finalPHapA)*lociHapACounts+
                log(finalPHapB)*(lociHapBCounts), 
            ## beta-binomial
            lbeta(lociHapACounts+ABhapA$a, lociHapBCounts+ABhapA$b)-
                lbeta(ABhapA$a, ABhapA$b) 
        )
        return(sum(logLikes))
    }
    return(returnFunc)
}

#' Function that given observed count data along a known haplotype returns a maximum likelihood estimate of the underlying haplotype frequency.
#'
#' @details Given observed read counts supporting hapltoype A at a collection of loci, the total read counts at those loci, the probablities of observing haplotype A-supporting reads under conditions of no ASE and the dispersion parameters, this function returns a maximum likelihood estimate of the true underlying frequency of haplotype A as well as corresponding value of log-likelihood.
#'
#' @param lociHapACounts counts of haplotype A-supporting reads at individual loci. Must be a vector of non-negative integers.
#' @param lociTotalCounts total read counts of at individual loci. Must be a vector of positive integers.
#' @param lociHapANoASEProbs probabilities of observing haplotype A-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus). Must be a vector with entries >0 and <1.
#' @param lociRhos dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial). Must be a numeric vector with entries >=0 and <1.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with two elements: maximum (MLE of haplotype A frequency) and objective (loglikelihood at MLE). These are the two elements that are output by the optimize() function, which is used internally by the maxLogLikelihoodCalculator1s.
#'
#' @examples
#' MBASED:::maxLogLikelihoodCalculator1s(lociHapACounts=c(5, 12), lociTotalCounts=c(10, 24), lociHapANoASEProbs=c(0.5, 0.5), lociRhos=c(0,0)) 
#'    
maxLogLikelihoodCalculator1s <- function (
    lociHapACounts, 
    lociTotalCounts, 
    lociHapANoASEProbs, 
    lociRhos, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(lociHapACounts) || 
            any(is.na(lociHapACounts)) || 
            !is.numeric(lociHapACounts) || 
            !isTRUE(all.equal(lociHapACounts, round(lociHapACounts))) || 
            any(lociHapACounts<0) 
        ) {
            stop('MBASED:maxLogLikelihoodCalculator1s: argument lociHapACounts must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociTotalCounts) || 
            any(is.na(lociTotalCounts)) || 
            !is.numeric(lociTotalCounts) || 
            !isTRUE(all.equal(lociTotalCounts, round(lociTotalCounts))) || 
            any(lociTotalCounts<=0) 
        ) {
            stop('MBASED:maxLogLikelihoodCalculator1s: argument lociTotalCounts must be a vector of positive integers')
        }
        if (
            !is.vector(lociHapANoASEProbs) || 
            any(is.na(lociHapANoASEProbs)) || 
            !is.numeric(lociHapANoASEProbs) || 
            any(lociHapANoASEProbs<=0) || 
            any(lociHapANoASEProbs>=1)
        ) {
            stop('MBASED:maxLogLikelihoodCalculator1s: argument lociHapANoASEProbs must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhos) || 
            any(is.na(lociRhos)) || 
            !is.numeric(lociRhos) || 
            any(lociRhos<0) || 
            any(lociRhos>=1)
        ) {
            stop('MBASED:maxLogLikelihoodCalculator1s: argument lociRhos must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociHapACounts), 
                length(lociTotalCounts), 
                length(lociHapANoASEProbs), 
                length(lociRhos)
            )))>1
        ) {
            stop('MBASED:maxLogLikelihoodCalculator1s: arguments lociHapACounts, lociTotalCounts, lociHapANoASEProbs, and lociRhos must be of same length')
        }
        if (any(lociHapACounts>lociTotalCounts)) {
            stop('MBASED:maxLogLikelihoodCalculator1s: entries of argument lociHapACounts must be less than or equal to corresponding entries in argument lociTotalCounts')
        }
    }
    logLikeFunc <- logLikelihoodCalculator1s(
        lociHapACounts=lociHapACounts, 
        lociTotalCounts=lociTotalCounts, 
        lociHapANoASEProbs=lociHapANoASEProbs, 
        lociRhos=lociRhos
    )
    return(optimize(logLikeFunc, c(0,1), maximum=TRUE))
}
 
#' Function that given observed count data returns a maximum likelihood estimate of the underlying haplotype frequency. Both situations where the haplotype are known and unknown are handled. In the latter case, likelihood is further maximized over all possible assignments of alleles to haplotypes.
#'
#' @details Given observed read counts supporting allele1 at a collection of loci, the total read counts at those loci, the probablities of observing allele1-supporting reads under conditions of no ASE and the dispersion parameters, this function returns a maximum likelihood estimate of the major haplotype frequency as well as corresponding assignment of alleles to haplotypes.
#'
#' @param lociAllele1Counts counts of allele1-supporting reads at individual loci. Must be a vector of non-negative integers.
#' @param lociTotalCounts total read counts of at individual loci. Must be a vector of positive integers.
#' @param lociAllele1NoASEProbs probabilities of observing allele1-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus). Must be a vector with entries >0 and <1.
#' @param lociRhos dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial). Must be a numeric vector with entries >=0 and <1.
#' @param isPhased single boolean specifying whether the phasing has already been performed, in which case the lociAllele1Counts represent the same haplotype. If FALSE (DEFAULT), likelihood is further maximized over all possible assignments of alleles to haplotypes.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with two elements: MAF (MLE of major allele frequency) and allele1IsMajor (whether allele1 is assigned to haplotype corresponding to maximum likleihood MAF). 
#'
#' @examples
#' MBASED:::estimateMAF1s(lociAllele1Counts=c(5, 24), lociTotalCounts=c(15, 36), lociAllele1NoASEProbs=c(0.5, 0.5), lociRhos=c(0,0), isPhased=TRUE) 
#' MBASED:::estimateMAF1s(lociAllele1Counts=c(5, 24), lociTotalCounts=c(15, 36), lociAllele1NoASEProbs=c(0.5, 0.5), lociRhos=c(0,0), isPhased=FALSE) 
#'    
estimateMAF1s <- function (
    lociAllele1Counts, 
    lociTotalCounts, 
    lociAllele1NoASEProbs, 
    lociRhos, 
    isPhased=FALSE, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(lociAllele1Counts) || 
            any(is.na(lociAllele1Counts)) || 
            !is.numeric(lociAllele1Counts) || 
            !isTRUE(all.equal(lociAllele1Counts, round(lociAllele1Counts))) || 
            any(lociAllele1Counts<0) 
        ) {
            stop('MBASED:estimateMAF1s: argument lociAllele1Counts must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociTotalCounts) || 
            any(is.na(lociTotalCounts)) || 
            !is.numeric(lociTotalCounts) || 
            !isTRUE(all.equal(lociTotalCounts, round(lociTotalCounts))) || 
            any(lociTotalCounts<=0) 
        ) {
            stop('MBASED:estimateMAF1s: argument lociTotalCounts must be a vector of positive integers')
        }
        if (
            !is.vector(lociAllele1NoASEProbs) || 
            any(is.na(lociAllele1NoASEProbs)) || 
            !is.numeric(lociAllele1NoASEProbs) || 
            any(lociAllele1NoASEProbs<=0) || 
            any(lociAllele1NoASEProbs>=1)
        ) {
            stop('MBASED:estimateMAF1s: argument lociAllele1NoASEProbs must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhos) || 
            any(is.na(lociRhos)) || 
            !is.numeric(lociRhos) || 
            any(lociRhos<0) || 
            any(lociRhos>=1)
        ) {
            stop('MBASED:estimateMAF1s: argument lociRhos must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociAllele1Counts), 
                length(lociTotalCounts), 
                length(lociAllele1NoASEProbs), 
                length(lociRhos)
            )))>1
        ) {
            stop('MBASED:estimateMAF1s: arguments lociAllele1Counts, lociTotalCounts, lociAllele1NoASEProbs, and lociRhos must be of same length')
        }
        if (any(lociAllele1Counts>lociTotalCounts)) {
            stop('MBASED:estimateMAF1s: entries of argument lociAllele1Counts must be less than or equal to corresponding entries in argument lociTotalCounts')
        }
        if ( 
            !(is.vector(isPhased)) || 
            !is.logical(isPhased) || 
            length(isPhased)!=1
        ) {
            stop('MBASED:estimateMAF1s: argument isPhased must be a single TRUE or FALSE value')
        }
    }
    lociAllele2Counts <- (lociTotalCounts-lociAllele1Counts)
    lociAllele2NoASEProbs <- (1-lociAllele1NoASEProbs)
    ## create matrix of all possible assignments of reference and 
    ## alternative alleles to haplotype A. Only 1 possible assignment 
    ## if the data is phased (actually, 2, bu they are equivalent). 
    ## Number of columns is equal to number of loci, while number of rows
    ## is equal to number of possible (non-equivalent) allele assignments 
    ## to haplotypes
    allele1IsHapAMat <- expand.grid(as.data.frame(matrix(
        rep(c(TRUE, FALSE), length(lociAllele1Counts)), 
        nrow=2
    ))) 
    if (!isPhased) {
        ## since each assignment is basically duplicated 
        ## (each row is equivalent to !row), 
        ## require that allele1 at first locus is haplotype A.
        allele1IsHapAMat <- allele1IsHapAMat[allele1IsHapAMat[,1],] 
    } else {
        ## only one non-redundant possible assignment
        allele1IsHapAMat <- allele1IsHapAMat[
            rowSums(!allele1IsHapAMat)==0,
        ]
    }
    maxLogLikelihoodByHapAssignment <- lapply(
        1:nrow(allele1IsHapAMat), 
        function(ind) {
            allele1IsHapA <- as.logical(allele1IsHapAMat[ind,])
            return(
                maxLogLikelihoodCalculator1s(
                    lociHapACounts=ifelse(
                        allele1IsHapA, 
                        lociAllele1Counts, 
                        lociAllele2Counts
                    ), 
                    lociTotalCounts=lociTotalCounts , 
                    lociHapANoASEProbs=ifelse(
                        allele1IsHapA, 
                        lociAllele1NoASEProbs, 
                        lociAllele2NoASEProbs
                    ), 
                    lociRhos=lociRhos
                )  
            )
        }
    )
    if (!isPhased) { ## pick the best assignment
        maxLogLikelihoodInd <- which.max(
            vapply(
                maxLogLikelihoodByHapAssignment, 
                function(el) {
                    el$objective
                }, 
                numeric(1)
            )
        )
    } else {
        maxLogLikelihoodInd <- 1
    }
    estimatedMAF <- 
        maxLogLikelihoodByHapAssignment[[maxLogLikelihoodInd]]$maximum
    allele1IsMajor <- as.logical(allele1IsHapAMat[maxLogLikelihoodInd,])
    if (estimatedMAF < 0.5) {
        estimatedMAF <- 1-estimatedMAF
        allele1IsMajor <- !(allele1IsMajor)
    }
    return(
        list(
            MAF=estimatedMAF,
            allele1IsMajor=allele1IsMajor
        )
    )
}

#' Function that given observed count data along a known haplotype returns a function that can calculate the likelihood of observing that data for a supplied underlying haplotype frequency.
#'
#' @details Given observed read counts supporting hapltoype A at a collection of loci in two samples, the total read counts at those loci, the probablities of observing haplotype A-supporting reads under conditions of no ASE and the dispersion parameters, this function returns a function of a single argument, pHapA, that calculates the likelihood of observing the given haplotype A-supporting counts under the assumption that the true underlying frequency of haplotype A is pHapA.
#'
#' @param lociHapACountsSample1,lociHapACountsSample2 counts of haplotype A-supporting reads at individual loci in sample1 and sample2, respectively. Both arguments must be vectors of non-negative integers.
#' @param lociTotalCountsSample1,lociTotalCountsSample2 total read counts of at individual loci in sample1 and sample2, respectively. Both arguments must be vectors of non-negative integers.
#' @param lociHapANoASEProbsSample1,lociHapANoASEProbsSample2 probabilities of observing haplotype A-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus) in sample1 and sample2, respectively. Both arguments must be vectors with entries >0 and <1.
#' @param lociRhosSample1,lociRhosSample2 dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial) in sample1 and sample2, respectively. Both arguments must be vectors with entries >=0 and <1.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a function of a single argument pHapA that calculates log likelihood of the observed data if the true underlying haplotype A frequency is pHapA.
#'
#' @examples
#' LLC <- MBASED:::logLikelihoodCalculator2s(lociHapACountsSample1=c(5, 12), lociTotalCountsSample1=c(15, 36), lociHapACountsSample2=c(15, 22), lociTotalCountsSample2=c(45, 66), lociHapANoASEProbsSample1=c(0.5, 0.5), lociHapANoASEProbsSample2=c(0.5, 0.5), lociRhosSample1=c(0,0), lociRhosSample2=c(0,0)) 
#' LLC(1/3) ## the MLE estimate of hapA frequency
#' LLC(0.9) ## highly implausible value of pHapA
#' LLC (0.334)
#' 
logLikelihoodCalculator2s <- function (
    lociHapACountsSample1, 
    lociTotalCountsSample1, 
    lociHapACountsSample2, 
    lociTotalCountsSample2, 
    lociHapANoASEProbsSample1, 
    lociHapANoASEProbsSample2, 
    lociRhosSample1, 
    lociRhosSample2, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(lociHapACountsSample1) || 
            any(is.na(lociHapACountsSample1)) || 
            !is.numeric(lociHapACountsSample1) || 
            !isTRUE(all.equal(lociHapACountsSample1, round(lociHapACountsSample1))) || 
            any(lociHapACountsSample1<0) 
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociHapACountsSample1 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociHapACountsSample2) || 
            any(is.na(lociHapACountsSample2)) || 
            !is.numeric(lociHapACountsSample2) || 
            !isTRUE(
                all.equal(
                    lociHapACountsSample2, 
                    round(lociHapACountsSample2)
                )
            ) || 
            any(lociHapACountsSample2<0) 
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociHapACountsSample2 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociTotalCountsSample1) || 
            any(is.na(lociTotalCountsSample1)) || 
            !is.numeric(lociTotalCountsSample1) || 
            !isTRUE(
                all.equal(
                    lociTotalCountsSample1, 
                    round(lociTotalCountsSample1)
                )
            ) || 
            any(lociTotalCountsSample1<=0) 
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociTotalCountsSample1 must be a vector of positive integers')
        }
        if ( 
            !is.vector(lociTotalCountsSample2) || 
            any(is.na(lociTotalCountsSample2)) || 
            !is.numeric(lociTotalCountsSample2) || 
            !isTRUE(
                all.equal(
                    lociTotalCountsSample2, 
                    round(lociTotalCountsSample2)
                )
            ) || 
            any(lociTotalCountsSample2<=0) 
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociTotalCountsSample2 must be a vector of positive integers')
        }
        if (
            !is.vector(lociHapANoASEProbsSample1) || 
            any(is.na(lociHapANoASEProbsSample1)) || 
            !is.numeric(lociHapANoASEProbsSample1) || 
            any(lociHapANoASEProbsSample1<=0) || 
            any(lociHapANoASEProbsSample1>=1)
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociHapANoASEProbsSample1 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociHapANoASEProbsSample2) || 
            any(is.na(lociHapANoASEProbsSample2)) || 
            !is.numeric(lociHapANoASEProbsSample2) || 
            any(lociHapANoASEProbsSample2<=0) || 
            any(lociHapANoASEProbsSample2>=1)
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociHapANoASEProbsSample2 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhosSample1) || 
            any(is.na(lociRhosSample1)) || 
            !is.numeric(lociRhosSample1) || 
            any(lociRhosSample1<0) || 
            any(lociRhosSample1>=1)
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociRhosSample1 must be a vector with entries >=0 and <1')
        }
        if (
            !is.vector(lociRhosSample2) || 
            any(is.na(lociRhosSample2)) || 
            !is.numeric(lociRhosSample2) || 
            any(lociRhosSample2<0) || 
            any(lociRhosSample2>=1)
        ) {
            stop('MBASED:logLikelihoodCalculator2s: argument lociRhosSample2 must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociHapACountsSample1), 
                length(lociTotalCountsSample1), 
                length(lociHapANoASEProbsSample1), 
                length(lociRhosSample1), 
                length(lociHapACountsSample2), 
                length(lociTotalCountsSample2), 
                length(lociHapANoASEProbsSample2), 
                length(lociRhosSample2)
            )))>1
        ) {
            stop('MBASED:logLikelihoodCalculator2s: arguments lociHapACountsSample1, lociTotalCountsSample1, lociHapANoASEProbsSample1, lociRhosSample1, lociHapACountsSample2, lociTotalCountsSample2, lociHapANoASEProbsSample2, and lociRhosSample2 must be of same length')
        }
        if (
            any(lociHapACountsSample1>lociTotalCountsSample1) || 
            any(lociHapACountsSample2>lociTotalCountsSample2)
        ) {
            stop('MBASED:logLikelihoodCalculator2s: entries of arguments lociHapACountsSample1 and lociHapACountsSample2 must be less than or equal to corresponding entries in arguments lociTotalCountsSample1 and lociTotalCountsSample2, respectively')
        }
    }
    logLikelihoodCalculatorSample1 <- logLikelihoodCalculator1s(
        lociHapACounts=lociHapACountsSample1, 
        lociTotalCounts=lociTotalCountsSample1, 
        lociHapANoASEProbs=lociHapANoASEProbsSample1, 
        lociRhos=lociRhosSample1
    ) 
    logLikelihoodCalculatorSample2 <- logLikelihoodCalculator1s(
        lociHapACounts=lociHapACountsSample2, 
        lociTotalCounts=lociTotalCountsSample2, 
        lociHapANoASEProbs=lociHapANoASEProbsSample2, 
        lociRhos=lociRhosSample2
    ) 
    returnFunc=function(pHapA) {
        logLikelihoodCalculatorSample1(pHapA)+
            logLikelihoodCalculatorSample2(pHapA)
    }
    return(returnFunc)
}
    
#' Function that given observed count data along a known haplotype returns a maximum likelihood estimate of the underlying haplotype frequency.
#'
#' @details Given observed read counts supporting hapltoype A at a collection of loci in two samples, the total read counts at those loci, the probablities of observing haplotype A-supporting reads under conditions of no ASE and the dispersion parameters, this function returns a maximum likelihood estimate of the true underlying frequency of haplotype A as well as corresponding value of log-likelihood.
#'
#' @param lociHapACountsSample1,lociHapACountsSample2 counts of haplotype A-supporting reads at individual loci in sample1 and sample2, respectively. Both arguments must be vectors of non-negative integers.
#' @param lociTotalCountsSample1,lociTotalCountsSample2 total read counts of at individual loci in sample1 and sample2, respectively. Both arguments must be vectors of non-negative integers.
#' @param lociHapANoASEProbsSample1,lociHapANoASEProbsSample2 probabilities of observing haplotype A-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus) in sample1 and sample2, respectively. Both arguments must be vectors with entries >0 and <1.
#' @param lociRhosSample1,lociRhosSample2 dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial) in sample1 and sample2, respectively. Both arguments must be vectors with entries >=0 and <1.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with two elements: maximum (MLE of haplotype A frequency) and objective (loglikelihood at MLE). These are the two elements that are output by the optimize() function, which is used internally by the maxLogLikelihoodCalculator2s.
#'
#' @examples
#' MBASED:::maxLogLikelihoodCalculator2s(lociHapACountsSample1=c(5, 12), lociTotalCountsSample1=c(15, 36), lociHapACountsSample2=c(15, 22), lociTotalCountsSample2=c(45, 66), lociHapANoASEProbsSample1=c(0.5, 0.5), lociHapANoASEProbsSample2=c(0.5, 0.5), lociRhosSample1=c(0,0), lociRhosSample2=c(0,0)) 
#'            
maxLogLikelihoodCalculator2s <- function (
    lociHapACountsSample1, 
    lociTotalCountsSample1, 
    lociHapACountsSample2, 
    lociTotalCountsSample2, 
    lociHapANoASEProbsSample1, 
    lociHapANoASEProbsSample2, 
    lociRhosSample1, 
    lociRhosSample2, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(lociHapACountsSample1) || 
            any(is.na(lociHapACountsSample1)) || 
            !is.numeric(lociHapACountsSample1) || 
            !isTRUE(
                all.equal(
                    lociHapACountsSample1, 
                    round(lociHapACountsSample1)
                )
            ) || 
            any(lociHapACountsSample1<0) 
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociHapACountsSample1 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociHapACountsSample2) || 
            any(is.na(lociHapACountsSample2)) || 
            !is.numeric(lociHapACountsSample2) || 
            !isTRUE(
                all.equal(
                    lociHapACountsSample2, 
                    round(lociHapACountsSample2)
                )
            ) || 
            any(lociHapACountsSample2<0) 
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociHapACountsSample2 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociTotalCountsSample1) || 
            any(is.na(lociTotalCountsSample1)) || 
            !is.numeric(lociTotalCountsSample1) || 
            !isTRUE(
                all.equal(
                    lociTotalCountsSample1, 
                    round(lociTotalCountsSample1)
                )
            ) || 
            any(lociTotalCountsSample1<=0) 
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociTotalCountsSample1 must be a vector of positive integers')
        }
        if ( 
            !is.vector(lociTotalCountsSample2) || 
            any(is.na(lociTotalCountsSample2)) || 
            !is.numeric(lociTotalCountsSample2) || 
            !isTRUE(
                all.equal(
                    lociTotalCountsSample2, 
                    round(lociTotalCountsSample2)
                )
            ) || 
            any(lociTotalCountsSample2<=0) 
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociTotalCountsSample2 must be a vector of positive integers')
        }
        if (
            !is.vector(lociHapANoASEProbsSample1) || 
            any(is.na(lociHapANoASEProbsSample1)) || 
            !is.numeric(lociHapANoASEProbsSample1) || 
            any(lociHapANoASEProbsSample1<=0) || 
            any(lociHapANoASEProbsSample1>=1)
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociHapANoASEProbsSample1 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociHapANoASEProbsSample2) || 
            any(is.na(lociHapANoASEProbsSample2)) || 
            !is.numeric(lociHapANoASEProbsSample2) || 
            any(lociHapANoASEProbsSample2<=0) || 
            any(lociHapANoASEProbsSample2>=1)
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociHapANoASEProbsSample2 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhosSample1) || 
            any(is.na(lociRhosSample1)) || 
            !is.numeric(lociRhosSample1) || 
            any(lociRhosSample1<0) || 
            any(lociRhosSample1>=1)
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociRhosSample1 must be a vector with entries >=0 and <1')
        }
        if (
            !is.vector(lociRhosSample2) || 
            any(is.na(lociRhosSample2)) || 
            !is.numeric(lociRhosSample2) || 
            any(lociRhosSample2<0) || 
            any(lociRhosSample2>=1)
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: argument lociRhosSample2 must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociHapACountsSample1), 
                length(lociTotalCountsSample1), 
                length(lociHapANoASEProbsSample1), 
                length(lociRhosSample1), 
                length(lociHapACountsSample2), 
                length(lociTotalCountsSample2), 
                length(lociHapANoASEProbsSample2), 
                length(lociRhosSample2)
            )))>1
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: arguments lociHapACountsSample1, lociTotalCountsSample1, lociHapANoASEProbsSample1, lociRhosSample1, lociHapACountsSample2, lociTotalCountsSample2, lociHapANoASEProbsSample2, and lociRhosSample2 must be of same length')
        }
        if (
            any(lociHapACountsSample1>lociTotalCountsSample1) || 
            any(lociHapACountsSample2>lociTotalCountsSample2)
        ) {
            stop('MBASED:maxLogLikelihoodCalculator2s: entries of arguments lociHapACountsSample1 and lociHapACountsSample2 must be less than or equal to corresponding entries in arguments lociTotalCountsSample1 and lociTotalCountsSample2, respectively')
        }
    }
    logLikeFunc <- logLikelihoodCalculator2s(
        lociHapACountsSample1=lociHapACountsSample1, 
        lociTotalCountsSample1=lociTotalCountsSample1, 
        lociHapACountsSample2=lociHapACountsSample2, 
        lociTotalCountsSample2=lociTotalCountsSample2, 
        lociHapANoASEProbsSample1=lociHapANoASEProbsSample1, 
        lociHapANoASEProbsSample2=lociHapANoASEProbsSample2, 
        lociRhosSample1=lociRhosSample1, 
        lociRhosSample2=lociRhosSample2
    )
    return(optimize(logLikeFunc, c(0,1), maximum=TRUE))
}
 
#' Function that given observed count data returns a maximum likelihood estimate of the underlying haplotype frequency. Both situations where the haplotype are known and unknown are handled. In the latter case, likelihood is further maximized over all possible assignments of alleles to haplotypes.
#'
#' @details Given observed read counts supporting allele1 at a collection of loci in two samples, the total read counts at those loci, the probablities of observing allele1-supporting reads under conditions of no ASE and the dispersion parameters, this function returns a maximum likelihood estimate of the major haplotype frequency as well as corresponding assignment of alleles to haplotypes.
#'
#' @param lociAllele1CountsSample1,lociAllele1CountsSample2 counts of allele1-supporting reads at individual loci in sample1 and sample2, respectively. Both arguments must be vectors of non-negative integers.
#' @param lociTotalCountsSample1,lociTotalCountsSample2 total read counts of at individual loci in sample1 and sample2, respectively. Both arguments must be vectors of non-negative integers.
#' @param lociAllele1NoASEProbsSample1,lociAllele1NoASEProbsSample2 probabilities of observing haplotype A-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus) in sample1 and sample2, respectively. Both arguments must be vectors with entries >0 and <1.
#' @param lociRhosSample1,lociRhosSample2 dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial) in sample1 and sample2, respectively. Both arguments must be vectors with entries >=0 and <1.
#' @param isPhased single boolean specifying whether the phasing has already been performed, in which case the lociAllele1CountsSample1 (and, therefore, lociAllele1CountsSample2) represent the same haplotype. If FALSE (DEFAULT), likelihood is further maximized over all possible assignments of alleles to haplotypes.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with two elements: MAF (MLE of major allele frequency) and allele1IsMajor (whether allele1 is assigned to haplotype corresponding to maximum likleihood MAF). 
#'
#' @examples
#' MBASED:::estimateMAF2s(lociAllele1CountsSample1=c(5, 24), lociTotalCountsSample1=c(15, 36), lociAllele1CountsSample2=c(15, 44), lociTotalCountsSample2=c(45, 66), lociAllele1NoASEProbsSample1=c(0.5, 0.5), lociAllele1NoASEProbsSample2=c(0.5, 0.5), lociRhosSample1=c(0,0), lociRhosSample2=c(0,0), isPhased=TRUE) 
#' MBASED:::estimateMAF2s(lociAllele1CountsSample1=c(5, 12), lociTotalCountsSample1=c(15, 36), lociAllele1CountsSample2=c(15, 22), lociTotalCountsSample2=c(45, 66), lociAllele1NoASEProbsSample1=c(0.5, 0.5), lociAllele1NoASEProbsSample2=c(0.5, 0.5), lociRhosSample1=c(0,0), lociRhosSample2=c(0,0), isPhased=FALSE) 
#'    
estimateMAF2s <- function (
    lociAllele1CountsSample1, 
    lociTotalCountsSample1, 
    lociAllele1CountsSample2, 
    lociTotalCountsSample2, 
    lociAllele1NoASEProbsSample1, 
    lociAllele1NoASEProbsSample2, 
    lociRhosSample1, 
    lociRhosSample2, 
    isPhased=FALSE, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(lociAllele1CountsSample1) || 
            any(is.na(lociAllele1CountsSample1)) || 
            !is.numeric(lociAllele1CountsSample1) || 
            !isTRUE(
                all.equal(
                    lociAllele1CountsSample1, 
                    round(lociAllele1CountsSample1)
                )
            ) || 
            any(lociAllele1CountsSample1<0) 
        ) {
            stop('MBASED:estimateMAF2s: argument lociAllele1CountsSample1 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociAllele1CountsSample2) || 
            any(is.na(lociAllele1CountsSample2)) || 
            !is.numeric(lociAllele1CountsSample2) || 
            !isTRUE(
                all.equal(
                    lociAllele1CountsSample2, 
                    round(lociAllele1CountsSample2)
                )
            ) || 
            any(lociAllele1CountsSample2<0) 
        ) {
            stop('MBASED:estimateMAF2s: argument lociAllele1CountsSample2 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociTotalCountsSample1) || 
            any(is.na(lociTotalCountsSample1)) || 
            !is.numeric(lociTotalCountsSample1) || 
            !isTRUE(
                all.equal(
                    lociTotalCountsSample1, 
                    round(lociTotalCountsSample1)
                )
            ) || 
            any(lociTotalCountsSample1<=0) 
        ) {
            stop('MBASED:estimateMAF2s: argument lociTotalCountsSample1 must be a vector of positive integers')
        }
        if ( 
            !is.vector(lociTotalCountsSample2) || 
            any(is.na(lociTotalCountsSample2)) || 
            !is.numeric(lociTotalCountsSample2) || 
            !isTRUE(
                all.equal(
                    lociTotalCountsSample2,
                    round(lociTotalCountsSample2)
                )
            ) || 
            any(lociTotalCountsSample2<=0) 
        ) {
            stop('MBASED:estimateMAF2s: argument lociTotalCountsSample2 must be a vector of positive integers')
        }
        if (
            !is.vector(lociAllele1NoASEProbsSample1) || 
            any(is.na(lociAllele1NoASEProbsSample1)) || 
            !is.numeric(lociAllele1NoASEProbsSample1) || 
            any(lociAllele1NoASEProbsSample1<=0) || 
            any(lociAllele1NoASEProbsSample1>=1)
        ) {
            stop('MBASED:estimateMAF2s: argument lociAllele1NoASEProbsSample1 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociAllele1NoASEProbsSample2) || 
            any(is.na(lociAllele1NoASEProbsSample2)) || 
            !is.numeric(lociAllele1NoASEProbsSample2) || 
            any(lociAllele1NoASEProbsSample2<=0) || 
            any(lociAllele1NoASEProbsSample2>=1)
        ) {
            stop('MBASED:estimateMAF2s: argument lociAllele1NoASEProbsSample2 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhosSample1) || 
            any(is.na(lociRhosSample1)) || 
            !is.numeric(lociRhosSample1) || 
            any(lociRhosSample1<0) || 
            any(lociRhosSample1>=1)
        ) {
            stop('MBASED:estimateMAF2s: argument lociRhosSample1 must be a vector with entries >=0 and <1')
        }
        if (
            !is.vector(lociRhosSample2) || 
            any(is.na(lociRhosSample2)) || 
            !is.numeric(lociRhosSample2) || 
            any(lociRhosSample2<0) || 
            any(lociRhosSample2>=1)
        ) {
            stop('MBASED:v: argument lociRhosSample2 must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociAllele1CountsSample1), 
                length(lociTotalCountsSample1), 
                length(lociAllele1NoASEProbsSample1), 
                length(lociRhosSample1), 
                length(lociAllele1CountsSample2), 
                length(lociTotalCountsSample2), 
                length(lociAllele1NoASEProbsSample2), 
                length(lociRhosSample2)
            )))>1
        ) {
            stop('MBASED:estimateMAF2s: arguments lociAllele1CountsSample1, lociTotalCountsSample1, lociAllele1NoASEProbsSample1, lociRhosSample1, lociAllele1CountsSample2, lociTotalCountsSample2, lociAllele1NoASEProbsSample2, and lociRhosSample2 must be of same length')
        }
        if (
            any(lociAllele1CountsSample1>lociTotalCountsSample1) || 
            any(lociAllele1CountsSample2>lociTotalCountsSample2)
        ) {
            stop('MBASED:estimateMAF2s: entries of arguments lociAllele1CountsSample1 and lociAllele1CountsSample2 must be less than or equal to corresponding entries in arguments lociTotalCountsSample1 and lociTotalCountsSample2, respectively')
        }
        if ( 
            !(is.vector(isPhased)) || 
            !is.logical(isPhased) || 
            length(isPhased)!=1
        ) {
            stop('MBASED:estimateMAF2s: argument isPhased must be a single TRUE or FALSE value')
        }
    }
    lociAllele2CountsSample1 <- (
        lociTotalCountsSample1-lociAllele1CountsSample1
    )
    lociAllele2CountsSample2 <- (
        lociTotalCountsSample2-lociAllele1CountsSample2
    )
    lociAllele2NoASEProbsSample1 <- (1-lociAllele1NoASEProbsSample1)
    lociAllele2NoASEProbsSample2 <- (1-lociAllele1NoASEProbsSample2)
    ## create matrix of all possible assignments of reference and 
    ## alternative alleles to haplotype A. Only 1 possible assignment 
    ## if the data is phased (actually, 2, bu they are equivalent). 
    ## Number of columns is equal to number of loci, while 
    ## number of rows is equal to number of possible (non-equivalent) 
    ## allele assignments to haplotypes
    allele1IsHapAMat <- expand.grid(as.data.frame(matrix(
        rep(c(TRUE, FALSE), length(lociAllele1CountsSample1)), 
        nrow=2
    ))) 
    if (!isPhased) {
        ## since each assignment is basically duplicated 
        ## (each row is equivalent to !row), require that allele1 at first locus
        ##  is haplotype A.
        allele1IsHapAMat <- allele1IsHapAMat[
            allele1IsHapAMat[,1], , drop=FALSE
        ] 
    } else {
        ## only one non-redundant possible assignment
        allele1IsHapAMat <- allele1IsHapAMat[
            rowSums(!allele1IsHapAMat)==0, , drop=FALSE
        ]
    }
    maxLogLikelihoodByHapAssignment <- lapply(
        1:nrow(allele1IsHapAMat), 
        function(ind) {
            allele1IsHapA <- as.logical(allele1IsHapAMat[ind,])
            return(
                maxLogLikelihoodCalculator2s(
                    lociHapACountsSample1=ifelse(
                        allele1IsHapA, 
                        lociAllele1CountsSample1, 
                        lociAllele2CountsSample1
                    ), 
                    lociTotalCountsSample1=lociTotalCountsSample1, 
                    lociHapACountsSample2=ifelse(
                        allele1IsHapA, 
                        lociAllele1CountsSample2, 
                        lociAllele2CountsSample2
                    ), 
                    lociTotalCountsSample2=lociTotalCountsSample2, 
                    lociHapANoASEProbsSample1=ifelse(
                        allele1IsHapA, 
                        lociAllele1NoASEProbsSample1, 
                        lociAllele2NoASEProbsSample1
                    ), 
                    lociHapANoASEProbsSample2=ifelse(
                        allele1IsHapA, 
                        lociAllele1NoASEProbsSample2, 
                        lociAllele2NoASEProbsSample2
                    ),
                    lociRhosSample1=lociRhosSample1,
                    lociRhosSample2=lociRhosSample2
                )
            )
        }
    )
    if (!isPhased) { ## pick the best assignment
        maxLogLikelihoodInd <- which.max(
            vapply(
                maxLogLikelihoodByHapAssignment, 
                function(el) {
                    el$objective
                }, 
                numeric(1)
            )
        )
    } else {
        maxLogLikelihoodInd <- 1
    }
    estimatedMAF <- 
        maxLogLikelihoodByHapAssignment[[maxLogLikelihoodInd]]$maximum
    allele1IsMajor <- as.logical(allele1IsHapAMat[maxLogLikelihoodInd,])
    if (estimatedMAF < 0.5) {
        estimatedMAF <- 1-estimatedMAF
        allele1IsMajor <- !(allele1IsMajor)
    }
    return(
        list(
            MAF=estimatedMAF,
            allele1IsMajor=allele1IsMajor
        )
    )
}
    
