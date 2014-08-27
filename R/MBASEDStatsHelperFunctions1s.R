#' Freeman-Tukey transformation functions.
#'
#' @details FT takes integers x and n, where x is observed Bin(n, p) random variable, and performs Freeman-Tukey transformation. Arguments x and n are vectorized and must be of the same length (if vectors) or dimension (if matrices).
#'
#' @param n number of trials, (vector/matrix of) positive number(s).
#' @param x number of successes, (vector/matrix of) non-negative number(s) <=n.
#' @param p probability of success on each trial, (vector/matrix of) value(s) between 0 and 1.
#' @param tieBreakRandom if FALSE, a backransformed value of 0.5 in isCountMajorFT() will be called major; if TRUE, it will be called major with probability 0.5 and minor with probability 0.5. DEFAULT: FALSE
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return FT returns (vector of) transformed proportion(s) of successes.
#' 
#' @rdname FreemanTukeyFunctions
#'
#' @aliases FT unFT FTAdjust
#'
#' @examples
#' isTRUE(all.equal(MBASED:::FT(x=5,n=10), pi/2))
#'
FT <- function (
    x, 
    n, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            (!is.vector(n) && !is.matrix(n)) || 
            any(is.na(n)) || 
            !is.numeric(n) || 
            any(n<=0) 
        ) {
            stop('MBASED:FT: argument n must be (vector/matrix of) positive number(s)')
        }
        if (
            (!is.vector(x) && !is.matrix(x)) || 
            any(is.na(x)) || 
            !is.numeric(x) || 
            any(x<0) 
        ) {
            stop('MBASED:FT: argument x must be (vector/matrix of) non-negative number(s)')
        }
        if (
            !isTRUE(all.equal(length(x),length(n))) || 
            !isTRUE(all.equal(dim(x), dim(n)))
        ) {
            stop('MBASED:FT: arguments x and n must be of same length/dimension')
        }
        if (any(x>n)) {
            stop('MBASED:FT: (entries of) argument x must be <= (corresponding entries in ) argument n')
        }
    }
    return(
        asin(sqrt(x/(n+1)))+asin(sqrt((x+1)/(n+1)))
    )
} 

#'
#' @details unFT takes transformed proportion and original total count and untransforms it, using the same approach as metaprop() function from R package "meta", with one correction: to avoid situations that arise in practice when z takes a value that cannot result from the supplied value of n (e.g. z corresponding to a count of < 0 out of n or > n out of n), we assign z to be the smallest/largest allowed value. Arguments z, and n are vectorized and must be of same length (if vectors) or dimension (if matrices).
#'
#' @param z (vector/matrix of) transformed proportion(s).
#'
#' @return unFT returns (vector/matrix of) backtransformed proportion(s).
#' 
#' @rdname FreemanTukeyFunctions
#'
#' @examples
#' MBASED:::unFT(z=MBASED:::FT(x=5, n=10), n=10)
#' MBASED:::unFT(z=MBASED:::FT(x=7, n=10), n=10)
#' isTRUE(all.equal(MBASED:::unFT(z=MBASED:::FT(x=7, n=10), n=10), 0.7))
#'
unFT <- function (
    z, 
    n, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            (!is.vector(n) && !is.matrix(n)) || 
            any(is.na(n)) || 
            !is.numeric(n) || 
            any(n<=0) 
        ) {
            stop('MBASED:unFT: argument n must be (vector/matrix of) positive number(s)')
        }
        if ( 
            (!is.vector(z) && !is.matrix(z)) || 
            any(is.na(z)) || 
            !is.numeric(z) 
        ) {
            stop('MBASED:unFT: argument z must be (vector/matrix of)  number(s)')
        }
        if (
            !isTRUE(all.equal(length(z),length(n))) || 
            !isTRUE(all.equal(dim(z),dim(n)))
        ) {
            stop('MBASED:unFT: arguments z and n must be of same length/dimension')
        }
    }
    minVals <- FT(
    ## can't use x=0 because I don't know if 
    ## n is vector/matrix of what length/dimension
        x=n-n, 
        n=n
    )
    maxVals <- FT(
        x=n, 
        n=n
    )
    z <- pmax(z, minVals)
    z <- pmin(z, maxVals)
    return(
        0.5*(1-sign(cos(z))*sqrt(1-(sin(z)+(sin(z)-1/sin(z))/n)^2))
    )
}

#'
#' @details FTAdjust takes integers x and n, and probability p, where x is observed Bin(n, p) random variable and performs Freeman-Tukey transformation, followed by shifting the transformed variable so that its mean is 2*arcsin(sqrt(0.5)) instead of 2*arcsin(sqrt(p)). Arguments x, n and p are vectorized and must be of the same length (if vectors) or dimension (if matrices).
#'
#' @return FTAdjust returns (vector/matrix of) shifted transformed proportion(s).
#' 
#' @rdname FreemanTukeyFunctions
#'
#' @examples
#' MBASED:::FT(x=50, n=100) 
#' MBASED:::FTAdjust(x=50, n=100, p=0.5) ## transformation is trivial if underlying probability of success is 0.5
#' MBASED:::FT(x=80, n=100) 
#' MBASED:::FTAdjust(x=80, n=100, p=0.8) ## if underlying probability of success is 0.8, the shift adjusts transformed proportion to have mean close to pi/2
#'
FTAdjust <- function (
    x, 
    n, 
    p, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            (!is.vector(n) && !is.matrix(n)) || 
            any(is.na(n)) || 
            !is.numeric(n) || 
            any(n<=0) 
        ) {
            stop('MBASED:FTAdjust: argument n must be (vector/matrix of) positive number(s)')
        }
        if ( 
            (!is.vector(x) && !is.matrix(x)) || 
            any(is.na(x)) || 
            !is.numeric(x) || 
            any(x<0) 
        ) {
            stop('MBASED:FTAdjust: argument x must be (vector/matrix of) non-negative number(s)')
        }
        if ( 
            (!is.vector(p) && !is.matrix(p)) || 
            any(is.na(p)) || 
            !is.numeric(p) || 
            any(p<=0) || 
            any(p>=1) 
        ) {
            stop('MBASED:FTAdjust: argument p must be (vector/matrix of) value(s) between 0 and 1')
        }
        if (
            !isTRUE(all.equal(length(x),length(n))) || 
            !isTRUE(all.equal(dim(x), dim(n))) || 
            !isTRUE(all.equal(length(p),length(n))) || 
            !isTRUE(all.equal(dim(p), dim(n)))
        ) {
            stop('MBASED:FTAdjust: arguments x, n and p must be of same length/dimension')
        }
        if (any(x>n)) {
            stop('MBASED:FTAdjust: (entries of) argument x must be <= (corresponding entries in ) argument n')
        }
    }
    return(
        FT(x=x, n=n)-2*asin(sqrt(p))+2*asin(sqrt(0.5))
    )
} 
    
#'
#' @details isCountMajorFT takes original observed count and total count, transforms, adjusting for underlying probability of success and returns TRUE or FALSE depending on whether the count is major (back-transformed proportion >=0.5 or not). Arguments x, n and p are vectorized and must be of same length (if vectors) or dimension (if matrices).
#'
#' @return isCountMajorFT returns  (vector/matrix of) TRUE or FALSE, depending on whether count is judged to be from 'major' allele or not.
#' 
#' @rdname FreemanTukeyFunctions
#'
#' @examples
#' MBASED:::isCountMajorFT(x=6, n=10, p=0.5, tieBreakRandom=FALSE)
#' MBASED:::isCountMajorFT(x=6, n=10, p=0.8, tieBreakRandom=FALSE)
#' MBASED:::isCountMajorFT(x=4, n=10, p=0.2, tieBreakRandom=FALSE) 
#' table(replicate(1000, MBASED:::isCountMajorFT(x=5, n=10, p=0.5, tieBreakRandom=FALSE)))
#' table(replicate(1000, MBASED:::isCountMajorFT(x=5, n=10, p=0.5, tieBreakRandom=TRUE)))
#'
isCountMajorFT <- function (
    x, 
    n, 
    p, 
    tieBreakRandom=FALSE, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            (!is.vector(n) && !is.matrix(n)) || 
            any(is.na(n)) || 
            !is.numeric(n) || 
            any(n<=0) 
        ) {
            stop('MBASED:isCountMajorFT: argument n must be (vector/matrix of) positive number(s)')
        }
        if ( 
            (!is.vector(x) && !is.matrix(x)) || 
            any(is.na(x)) || 
            !is.numeric(x) || 
            any(x<0) 
        ) {
            stop('MBASED:isCountMajorFT: argument x must be (vector/matrix of) non-negative number(s)')
        }
        if (
            (!is.vector(p) && !is.matrix(p)) || 
            any(is.na(p)) || 
            !is.numeric(p) ||  
            any(p<=0) || 
            any(p>=1) 
        ) {
            stop('MBASED:isCountMajorFT: argument p must be (vector/matrix of) value(s) between 0 and 1')
        }
        if (
            !isTRUE(all.equal(length(x),length(n))) || 
            !isTRUE(all.equal(dim(x), dim(n))) || 
            !isTRUE(all.equal(length(p),length(n))) || 
            !isTRUE(all.equal(dim(p), dim(n)))
        ) {
            stop('MBASED:isCountMajorFT: arguments x, n and p must be of same length/dimension')
        }
        if (any(x>n)) {
            stop('MBASED:isCountMajorFT: (entries of) argument x must be <= (corresponding entries in ) argument n')
        }
        if ( 
            !(is.vector(tieBreakRandom)) || 
            !is.logical(tieBreakRandom) || 
            length(tieBreakRandom)!=1 
        ) {
            stop('MBASED:isCountMajorFT: argument tieBreakRandom must be a single TRUE or FALSE value')
        }
    }
    unFTVals <- unFT(
        z=FTAdjust(
            x=x, 
            n=n, 
            p=p
        ), 
        n=n
    )
    returnVals <- (unFTVals>=0.5)
    if (tieBreakRandom) {
        tiedValsSubv <- (unFTVals==0.5)
        returnVals[tiedValsSubv] <- sample(
            c(TRUE,FALSE), 
            sum(tiedValsSubv), 
            replace=TRUE
        )
    }
    return(returnVals)
}

#' Helper function to obtain estimate of underlying mean and the standard error of the estimate in meta analysis framework.
#'
#' @details MBASEDMetaAnalysisGetMeansAndSEs is a helper function employed by MBASEDMetaAnalysis(). For each column of input matrices, it calculates the inverse-variance weighted column average and provides an estimate of the standard error of this mean estimator. Input matrices zValuesMat and zVariancesMat have one column for each set of loci ('independent studies') to be combined, with each row corresponding to an individual locus.
#'
#' @param zValuesMat matrix of z-values, on standard normal scale. Each row represents a specific genomic locus, while each column represents a set of observed values across loci (in practice, multiple columns represent different outcomes of simulations).
#' @param zVariancesMat matrix of (estimated) variances of each z-value in zValuesMat. The interpretation of rows and columns is the same as for zValuesMat.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with 4 elements:
#' \item{weightsMat}{a matrix of same dimension as zValuesMat, giving the assigned weight for each observation}
#' \item{totalWeights}{a vector of length equal to number of rows in zValuesMat, giving the column sum of assigned weights}
#' \item{hetQ}{a vector of length equal to number of rows in zValuesMat, giving the estimated standard error for the corresponding entries in meanValues}
#' \item{meanValues}{a vector of length equal to number of rows in zValuesMat, giving for each column the estimated average value.}
#' \item{hetQ}{a vector of length equal to number of rows in zValuesMat, giving the estimated standard error for the corresponding entries in meanValues}
#' 
#'
#' @examples
#' set.seed(127000)
#' zVals1=rnorm(5, mean=rep(2,5), sd=sqrt(1:5))
#' zVals2=rnorm(5, mean=0, sd=1)+c(0,0,5,0,0) ## one outlier
#' MBASED:::MBASEDMetaAnalysisGetMeansAndSEs(zValuesMat=matrix(c(zVals1, zVals2), ncol=2), zVariancesMat=matrix(c(1:5, rep(1,5)), ncol=2)) 
#' 
MBASEDMetaAnalysisGetMeansAndSEs <- function (
    zValuesMat, 
    zVariancesMat, 
    checkArgs=FALSE
) { 
	## rows are loci, columns are simulations 
	## (in case of no simulations, one column each matrix)
    if (checkArgs) {
        if (
            !is.matrix(zValuesMat) || 
            any(is.na(zValuesMat)) || 
            !is.numeric(zValuesMat)
        ) {
            stop('MBASED:MBASEDMetaAnalysisGetMeansAndSEs: argument zValuesMat must be a matrix of numbers')
        }
        if (
            !is.matrix(zVariancesMat) || 
            any(is.na(zVariancesMat)) || 
            !is.numeric(zVariancesMat) || 
            any(zVariancesMat<=0)
        ) {
            stop('MBASED:MBASEDMetaAnalysisGetMeansAndSEs: argument zVariancesMat must be a matrix of positive numbers')
        }
        if (
            !isTRUE(all.equal(dim(zValuesMat), dim(zVariancesMat))) 
        ) {
            stop('MBASED:MBASEDMetaAnalysisGetMeansAndSEs: arguments zValuesMat, and zVariancesMat must have same dimension')
        }
    }
    weightsMat <- 1/zVariancesMat
    totalWeights <- colSums(weightsMat)
    returnMeans <- colSums(weightsMat*zValuesMat)/totalWeights
    returnSEs <- sqrt(1/totalWeights)
    return(
        list(
            weightsMat=weightsMat,
            totalWeights=totalWeights,
            meanValues=returnMeans,
            meanSEs=returnSEs
        )
    )
}

#' Generic function to perform standard meta analysis.
#'
#' @details MBASEDMetaAnalysis performs meta analysis calculations in a vectorized fashion. Input matrices zValuesMat and zVariancesMat have one column for each set of loci ('independent studies') to be combined, with each row corresponding to an individual locus. MBASEDMetaAnalysis uses meta analysis approach to combine values in each column of zValuesMat into a single column-specific value of z (using corresponding supplied variances to appropriately weight contributions of each individual z). The function reports the resulting averaged z values, together with corresponding standard deviations (standard errors), for fixed-effects setting (note: random effects are not meaningful in the context of SNVs in ASE). If the supplied matrices have a single row (only one locus), no meta-analysis is possible, and the original value and corresponding standard deviations are returned.
#'
#' @param zValuesMat matrix of z-values, on standard normal scale. Each row represents a specific genomic locus, while each column represents a set of observed values across loci (in practice, multiple columns represent different outcomes of simulations).
#' @param zVariancesMat matrix of (estimated) variances of each z-value in zValuesMat. The interpretation of rows and columns is the same as for zValuesMat.
#' @param alternative one of 'two.sided', 'greater', 'less'. DEFAULT: 'two.sided'.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with 5 elements:
#' \item{hetPVal}{a 1-row marix of heterogeneity p-values.}
#' \item{hetQ}{a 1-row matrix of heterogeneity statistics.}
#' \item{fixedEffectsMeans}{a 1-row matrix of column-specific fixed-effects meta analysis restults.}
#' \item{fixedEffectsSEs}{a 1-row matrix of estimated SEs of fixed-effects meta analysis results.}
#' \item{pvalueFixed}{a 1-row matrix of p-values for fixed-effects analysis.}
#' 
#' @examples
#' set.seed(127000)
#' zVals1=rnorm(5, mean=rep(2,5), sd=sqrt(1:5))
#' zVals2=rnorm(5, mean=0, sd=1)+c(0,0,5,0,0) ## one outlier
#' MBASED:::MBASEDMetaAnalysis(zValuesMat=matrix(c(zVals1, zVals2), ncol=2), zVariancesMat=matrix(c(1:5, rep(1,5)), ncol=2), alternative='two.sided') 
#' 
MBASEDMetaAnalysis <- function (
    zValuesMat, 
    zVariancesMat, 
    alternative='two.sided', 
    checkArgs=FALSE
) { 
	## rows are loci, columns are simulations 
	## (in case of no simulations, one column each matrix)
    if (checkArgs) {
        if (
            !is.matrix(zValuesMat) || 
            any(is.na(zValuesMat)) || 
            !is.numeric(zValuesMat)
        ) {
            stop('MBASED:MBASEDMetaAnalysis: argument zValuesMat must be a matrix of numbers')
        }
        if (
            !is.matrix(zVariancesMat) || 
            any(is.na(zVariancesMat)) || 
            !is.numeric(zVariancesMat) || 
            any(zVariancesMat<=0)
        ) {
            stop('MBASED:MBASEDMetaAnalysis: argument zVariancesMat must be a matrix of positive numbers')
        }
        if (
            !isTRUE(all.equal(dim(zValuesMat), dim(zVariancesMat))) 
        ) {
            stop('MBASED:MBASEDMetaAnalysis: arguments zValuesMat, and zVariancesMat must have same dimension')
        }
        if ( 
            !(is.vector(alternative)) || 
            length(alternative)!=1 || 
            !(alternative%in%c('greater','less', 'two.sided')) 
        ) {
            stop("MBASED:MBASEDMetaAnalysis: argument alternative must be a single character string with a value of 'greater', 'less', or 'equal' ")
        }
    }
    numLoci <- nrow(zValuesMat)
    numRepetitions <- ncol(zValuesMat)
    fixedEffectsResults <- MBASEDMetaAnalysisGetMeansAndSEs(
        zValuesMat=zValuesMat, 
        zVariancesMat=zVariancesMat
    )
    fixedEffectsMeans <- fixedEffectsResults$meanValues
    fixedEffectsSEs <- fixedEffectsResults$meanSEs
    if (numLoci==1) { ## single locus, no inter-loci variability
        extraVars <- 0
        hetPVals <- rep(NA, numRepetitions)
        hetQs <- rep(NA, numRepetitions)
    } else {
        fixedZsMat <- matrix(
            rep(fixedEffectsMeans, each=numLoci), 
            nrow=numLoci, 
            ncol=numRepetitions
        ) ## matrix to be used for computing hetQ
        fixedWeights <- fixedEffectsResults$weightsMat
        totalFixedWeights <- fixedEffectsResults$totalWeights
        hetQs <- colSums(fixedWeights*((zValuesMat-fixedZsMat)^2))
        hetPVals <- 1-pchisq(hetQs,numLoci-1)
    }
    ## calculate p-values:
    zScoresFixedEffects <- fixedEffectsMeans/fixedEffectsSEs
    if (alternative=='two.sided') {
        pValsFixed <- 2*(1-pnorm(abs(zScoresFixedEffects)))
    } else if (alternative=='greater') { 
        pValsFixed <- (1-pnorm(zScoresFixedEffects))
    } else if (alternative=='less') { 
        pValsFixed <- pnorm(zScoresFixedEffects)
    }
    return(
        list(
            hetPVal=matrix(
                hetPVals, 
                nrow=1
            ), 
            hetQ=matrix(
                hetQs, 
                nrow=1
            ), 
            fixedEffectsMeans=matrix(
                fixedEffectsMeans,
                nrow=1
            ), 
            fixedEffectsSEs=matrix(
                fixedEffectsSEs,
                nrow=1
            ), 
            pValsFixed=matrix(
                pValsFixed,
                nrow=1
            )
        )
    )
}

#' Vectorized wrapper around metaprop() function from R package "meta" with some modifications and extensions to beta-binomial count models.
#'
#' @details MBASEDVectorizedMetaprop performs computations similar to metaprop() with default options (fixed-effects only), but with less overhang, in a vectorized fashion, and accomodating extensions to beta-binomial distribution. It also allows the input counts to come from loci with different underlying binomial probabilities (means of beta distribution, in cases of beta-binomial extensions). One technical difference is the way the value of 'n' is calculated for Freeman-Tukey back-transformation of average 'z' into a proportion. While metaprop() uses harmonic mean of n's at individual loci (which puts more weight toward loci with small read counts), we use the weighted mean of n's with weights proportional to n's, by analogy with how the value of average 'z' is calculated from individual z's. Input matrices countsMat and totalsMat have one column for each set of loci ('independent studies') to be combined, with each row corresponding to an individual locus. Matrix probsMat provides the underlying binomial probabilities (means of beta distributions, in case of beta-binomial extensions), while matrix rhosMat gives the values of the dispersion at the loci. MBASEDVectorizedMetaprop uses meta analysis approach with Freeman-Tukey transformation to report for each set of loci (each column) its estimated overall proportion on both transformed and untransformed scale, corresponding standard errors (on transformed scale), z-values (based on expected value of 2*asin(sqrt(0.5)) on transformed scale under the null hypothesis of common underlying proportion (binomial probability or mean of beta distribution) of 0.5), and corresponding p-values based on imposing normal distribution assumption on z-values, where alternative hypothesis of 'two.sided', 'greater', and 'less' can be specified, with the latter two specified w.r.t. 2*asin(sqrt(0.5)). If some of the supplied entries in probsMat are different from 0.5, then the corresponding transformed proportions are shifted, so that the new mean for the resulting z-values is still approximately 2*asin(sqrt(0.5)). Extensions to beta-binomial counts are accomodated by increasing the variance of each individual z from 1/(n+0.5) to rho+1/(n+0.5), where rho is the dispersion parameter of the beta distribution. The function is used to cacluate p-values in ASE settings, where countsMat represents major allele counts, totalsMat represents total allele counts, probsMat represents the underlying binomial probabilities of observing major allele-supporting read (means of beta distributions in case of beta-binomial extensions), which may be different, e.g.,  for major allele counts coming from reference and alternative alleles in case of pre-existing allelic bias, and rhosMat provides values of dispersion parameter for beta-binomial counts (0, in case of binomial model) for individual loci. 
#'
#' @param countsMat matrix of observed major allele counts. Each row represents a specific genomic locus, while each column represents a set of observed major allele counts across loci (in practice, multiple columns represent different outcomes of count simulations).
#' @param totalsMat matrix of total read counts across both alleles. The interpretation of rows and columns is the same as for countsMat.
#' @param probsMat matrix of probabilities of success (means of beta distributions in case of beta-binomial extensions). The interpretation of rows and columns is the same as for countsMat.
#' @param rhosMat matrix of dispersion parameters of beta distribution in case of beta-binomial counts. The interpretation of rows and columns is the same as for countsMat.
#' @param alternative one of 'two.sided', 'greater', 'less'. DEFAULT: 'two.sided'.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a list with 7 elements:
#' \item{hetPVal}{a 1-row marix of heterogeneity p-values.}
#' \item{hetQ}{a 1-row matrix of heterogeneity statistics.}
#' \item{TEFinal}{a 1-row matrix of estimated proportions on transformed scale.}
#' \item{seTEFinal}{a 1-row matrix of estimated SEs of proportion estimates on transformed scale.}
#' \item{propFinal}{a 1-row matrix of estimated proportions on 0-1 scale.}
#' \item{pValue}{a 1-row matrix of corresponding p-values.}
#' \item{propLoci}{a matrix of same dimension as original input matrices giving estimated proportions on transformed scale at each individual locus.}
#'
#' @examples
#' SNVCoverage=rep(sample(10:100,5),2) ## 2 genes with 5 loci each
#' SNVAllele1Counts=rbinom(length(SNVCoverage), SNVCoverage, 0.5)
#' SNVMajorAlleleCounts=pmax(SNVAllele1Counts, SNVCoverage-SNVAllele1Counts)
#' MBASED:::MBASEDVectorizedMetaprop(countsMat=matrix(SNVAllele1Counts, ncol=2), totalsMat=matrix(SNVCoverage, ncol=2), probsMat=matrix(rep(0.5, length(SNVCoverage)), ncol=2), rhosMat=matrix(rep(0, length(SNVCoverage)), ncol=2), alternative='two.sided') ## ideal situation when phasing is known
#' MBASED:::MBASEDVectorizedMetaprop(countsMat=matrix(SNVMajorAlleleCounts, ncol=2), totalsMat=matrix(SNVCoverage, ncol=2), probsMat=matrix(rep(0.5, length(SNVCoverage)), ncol=2), rhosMat=matrix(rep(0, length(SNVCoverage)), ncol=2), alternative='two.sided') ## what happens if we put all major alleles together into a single haplotype and obtain nominal p-value
#' 
MBASEDVectorizedMetaprop <- function (
    countsMat, 
    totalsMat, 
    probsMat, 
    rhosMat, 
    alternative='two.sided', 
    checkArgs=FALSE
) { 
	## rows are loci, columns are simulations 
	## (in case of no simulations, one column each matrix)
    if (checkArgs) {
        if (
            !is.matrix(totalsMat) || 
            any(is.na(totalsMat)) || 
            !is.numeric(totalsMat) || 
            !isTRUE(all.equal(totalsMat, round(totalsMat))) || 
            any(totalsMat<1)
        ) {
            stop('MBASED:MBASEDVectorizedMetaprop: argument totalsMat must be a matrix of positive integers')
        }
        if (
            !is.matrix(countsMat) ||
            any(is.na(countsMat)) || 
            !is.numeric(countsMat) || 
            !isTRUE(all.equal(countsMat, round(countsMat))) || 
            any(countsMat<0)
        ) {
            stop('MBASED:MBASEDVectorizedMetaprop: argument countsMat must be a matrix of non-negative integers')
        }
        if (
            !is.matrix(probsMat) || 
            any(is.na(probsMat)) || 
            !is.numeric(probsMat) || 
            any(probsMat<=0) || 
            any(probsMat>=1)
        ) {
            stop('MBASED:MBASEDVectorizedMetaprop: argument probsMat must be a matrix with entries >0 and <1')
        }
        if (
            !is.matrix(rhosMat) || 
            any(is.na(rhosMat)) || 
            !is.numeric(rhosMat) || 
            any(rhosMat<0) || 
            any(rhosMat>=1)
        ) {
            stop('MBASED:MBASEDVectorizedMetaprop: argument rhosMat must be a matrix with entries >=0 and <1')
        }
        if (
            !isTRUE(all.equal(dim(countsMat), dim(totalsMat))) ||  
            !isTRUE(all.equal(dim(probsMat), dim(totalsMat))) || 
            !isTRUE(all.equal(dim(rhosMat), dim(totalsMat)))
        ) {
            stop('MBASED:MBASEDVectorizedMetaprop: arguments countsMat, totalsMat, probsMat and rhosMat must have same dimension')
        }
        if (any(countsMat>totalsMat)) {
            stop('MBASED:MBASEDVectorizedMetaprop: entries of argument countsMat must be <= corresponding entries in  argument totalsMat')
        }
        if ( 
            !(is.vector(alternative)) || 
            length(alternative)!=1 || 
            !(alternative%in%c('greater','less', 'two.sided')) 
        ) {
            stop("MBASED:MBASEDVectorizedMetaprop: argument alternative must be a single character string with a value of 'greater', 'less', or 'equal' ")
        }
    }
    ## get props after adjusting for prior allelic bias
    FTs <- FTAdjust(
        x=countsMat, 
        n=totalsMat, 
        p=probsMat
    )
    lociFreqs <- unFT(
        z=FTs, 
        n=totalsMat
    )
    FTvars <- (rhosMat+1/(totalsMat+0.5))
    metaAnalysisResults <- MBASEDMetaAnalysis(
        zValuesMat=FTs-pi/2, 
        zVariancesMat=FTvars, 
        alternative=alternative
    ) 
    finalFTs <- metaAnalysisResults$fixedEffectsMeans+pi/2
    finalSEs <- metaAnalysisResults$fixedEffectsSEs
    pVals <- metaAnalysisResults$pValsFixed
    pooledFreqs <- unFT(
        z=finalFTs, 
        ## n for back-transformation is weighted sum of locus-specific n's, 
        ## where the weight at each nLocus is proportional to nLocus itself.
        n=matrix(
            colSums(totalsMat^2)/colSums(totalsMat), 
            nrow=1
        )
    ) 
    return(
        list(
            hetPVal=metaAnalysisResults$hetPVal, 
            hetQ=metaAnalysisResults$hetQ, 
            TEFinal=finalFTs, 
            seTEFinal=finalSEs, 
            propFinal=pooledFreqs, 
            pValue=pVals,
            propLoci=lociFreqs
        )
    )
}


#' Function to calculate simulations-based p-values
#'
#' @details this function calculates fraction of simulated values (statistics from null distribution) that are >= (direction='greater') or <= (direction='less') than the observed statistic. The choice of direction depends on the nature of the statistic (i.e., direction is 'greater' if large values of statistic indicate departure from null hypothesis, and direction is 'less' if the opposite is the case)
#'
#' @param observedVal observed statistic (single number)
#' @param simulatedVals statistics observed in simulations of the outcomes based on assumed null distribution.
#' @param direction one of 'greater' or 'less', depending on the nature of statistic.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a fraction of simulated statistics that are as or more extreme as the observed one
#'
#' @examples
#' MBASED:::getSimulationPvalue(observedVal=2, simulatedVals=1:10, direction='greater')
#' MBASED:::getSimulationPvalue(observedVal=2, simulatedVals=1:10, direction='less')
#'
getSimulationPvalue <- function (
    observedVal, 
    simulatedVals, 
    direction='greater', 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(observedVal) || 
            is.na(observedVal) || 
            !is.numeric(observedVal) || 
            length(observedVal)!=1
        ) {
            stop('MBASED:getEmpiricalPvalue: argument observedVal must be a single numerical value')
        }
        if (
            !is.vector(simulatedVals) || 
            any(is.na(simulatedVals)) || 
            !is.numeric(simulatedVals) || 
            length(observedVal)==0
        ) {
            stop('MBASED:getEmpiricalPvalue: argument simulatedVals must be a vector of numerical values')
        }
        if ( 
            !(is.vector(direction)) || 
            length(direction)!=1 || 
            !(direction%in%c('greater','less')) 
        ) {
            stop("MBASED:getEmpiricalPvalue: argument direction must be a single character string with a value of 'greater' or 'less' ")
        }
    }
    ## need this to deal with issue of round-off errors
    tolerance=(.Machine$double.eps^0.5) 
    if (direction=='greater') {
        empPval <- sum(
            simulatedVals>=(observedVal-tolerance)
        )/length(simulatedVals)
    } else if (direction=='less') {
        empPval <- sum(
            simulatedVals<=(observedVal+tolerance)
        )/length(simulatedVals)
    } else {
        stop("MBASED:getEmpiricalPvalue: argument direction must be a single character string with a value of 'greater' or 'less' ")
    }
    return(empPval)
} 

