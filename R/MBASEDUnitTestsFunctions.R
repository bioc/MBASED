 #' Function to test quantile equality for theoretical and observed distributions
#'
#' @details For some random variable X, observed sample x1, x2, .., xN, and attainable value x, we compare theoretical P(X<=x) to observed Num(xi <= x)/N.
#'
#' @param theoreticalCumDist for (unspecified) value of x, P(X<=x)
#' @param observedCumDist for (unspecified) value of x, observed Fraction(values<=x) = Num(values<=x)/Num(total values). Actual values of x must be the same as those for corresponding entries in theoreticalCumDist
#' @param numTotalCounts Num(total values) (see argument observedCumDist)
#' @param numSEsToCheck number of standard errors to go in each direction from theoretical quantity to see if the estimate falls into the confidence interval
#' @param errorMessage error message to return if observed fraction falls outside of confidence interval
#'
#' @return TRUE (all tests were passed, otherwise exits with error message).
#'
#' @family unitTestsFunctions
#'
testQuantiles <- function (
    theoreticalCumDist, 
    observedCumDist, 
    numTotalCounts, 
    numSEsToCheck, 
    errorMessage
) {
    for (quantileInd in 1:length(theoreticalCumDist)) {
        obsF <- observedCumDist[quantileInd]
        predF <- theoreticalCumDist[quantileInd]
        ## standard error (based on theory)
        predFVar <- predF*(1-predF)/numTotalCounts 
        if (isTRUE(all.equal(predF,1))) {
            checkEqualsNumeric(
                predF, 
                obsF, 
                msg=errorMessage
            )
        } else {
            predFSD <- sqrt(predFVar)
            checkTrue(
                (obsF>=(predF-numSEsToCheck*predFSD)) && 
                (obsF<=(predF+numSEsToCheck*predFSD)), 
                msg=errorMessage
            )
        }
    }
    
    return(TRUE)
}


#' Function that checks to see if the difference between 2 number is small enough.
#'
#' @details for 2 numbers a and b, the function checks to see if |a-b|/min(a,b) <= cutoff.
#'
#' @param queryVals,targetVals vectors of values to be compared (pairwise comparison will be performed)
#' @param cutoffFraction the value of cutoff to be used to declare if the two numbers are close enough.
#'
#' @return vector of same length as input vectors queryVals and targetVals, recording for each pair of numbers whether they pas the cutoff (TRUE) or not (FALSE).
#'
#' @family unitTestsFunctions
#'
testNumericDiff <- function(
    queryVals, 
    targetVals, 
    cutoffFraction
) {
    myPrecision <- .Machine$double.eps^0.5
    diffBetweenEstimates <- abs(queryVals-targetVals)
    minOfTwoEstimates <- pmin(abs(queryVals), abs(targetVals))
    return(
        ifelse(
            minOfTwoEstimates==0,
            (queryVals+myPrecision)>=targetVals & 
                (queryVals-myPrecision)<=targetVals,
            (diffBetweenEstimates/minOfTwoEstimates)<=cutoffFraction
        )
    )
}
    
