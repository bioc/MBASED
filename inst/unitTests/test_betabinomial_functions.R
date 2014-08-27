## size of round-off error I'm tolerating 
##(quantities with differences of this much or smaller are considered the same)
myPrecision <- .Machine$double.eps^0.5 
FULLTESTING <- FALSE ## whether everything should be tested
if (FULLTESTING) {
    myBPPARAM <- MulticoreParam(workers=24)
} else {
    myBPPARAM <- SerialParam()
}

## test getMuRho() function
test_getMuRho <- function() {
    ## check1: we get expected values for known case: uniform distribution
    checkEqualsNumeric(
        MBASED:::getMuRho(a=1, b=1)$mu, 
        1/2, 
        msg='test_getMuRho: fail getMuRho uniform distribution test'
    )
    checkEqualsNumeric(
        MBASED:::getMuRho(a=1, b=1)$rho, 
        1/3, 
        msg='test_getMuRho: fail getMuRho uniform distribution test'
    )
    
    ## check2: vectorization works:
    testCasesMat <- expand.grid(
        test_a=10^(-4:4),
        test_b=10^(-4:4)
    )
    testCasesMat$test_mu=with(testCasesMat, 
        test_a/(test_a+test_b)
    )
    testCasesMat$test_rho=with(testCasesMat, 
        1/(1+test_a+test_b)
    )    
    testCasesMat_getMuRho <- with(testCasesMat, 
        MBASED:::getMuRho(a=test_a, b=test_b)
    )
    checkEqualsNumeric(
        testCasesMat$test_mu, 
        testCasesMat_getMuRho$mu, 
        msg='test_getMuRho: fail getMuRho vectorization test'
    )
    checkEqualsNumeric(
        testCasesMat$test_rho, 
        testCasesMat_getMuRho$rho, 
        msg='test_getMuRho: fail getMuRho vectorization test'
    )
    
    return(TRUE)
}

## test getAB() function
test_getAB <- function() {
    ## check1: we get expected values for known case: uniform distribution
    checkEqualsNumeric(
        MBASED:::getAB(mu=1/2, rho=1/3)$a, 
        1, 
        msg='test_getAB: fail getAB uniform distribution test'
    )
    checkEqualsNumeric(
        MBASED:::getAB(mu=1/2, rho=1/3)$b, 
        1, 
        msg='test_getAB: fail getAB uniform distribution test'
    )
    
    ## check2: vectorization works:
    testCasesMat <- expand.grid(
        test_a=10^(-4:4),
        test_b=10^(-4:4)
    )
    testCasesMat$test_mu=with(testCasesMat, 
        test_a/(test_a+test_b)
    )
    testCasesMat$test_rho=with(testCasesMat, 
        1/(1+test_a+test_b)
    )    
    testCasesMat_getAB <- with(testCasesMat, 
        MBASED:::getAB(mu=test_mu, rho=test_rho)
    )
    checkEqualsNumeric(
        testCasesMat$test_a, 
        testCasesMat_getAB$a, 
        msg='test_getAB: fail getAB vectorization test'
    )
    checkEqualsNumeric(
        testCasesMat$test_b, 
        testCasesMat_getAB$b, 
        msg='test_getAB: fail getAB vectorization test'
    )
    
    return(TRUE)
}
    
## test wheter getMuRho and getAB are inverses of each other.
test_getMuRho_getAB <- function() {
    ## check1: functions are inverses of each other:
    testCasesMat <- expand.grid(
        test_a=10^(-4:4),
        test_b=10^(-4:4)
    )
    testCasesMat$test_mu=with(testCasesMat, 
        test_a/(test_a+test_b)
    )
    testCasesMat$test_rho=with(testCasesMat, 
        1/(1+test_a+test_b)
    )
    testResults <- bplapply(1:nrow(testCasesMat), function(testInd) { 
        testCase <- testCasesMat[testInd,]
        testCase_getMuRho <- with(testCase, 
            MBASED:::getMuRho(a=test_a, b=test_b)
        )
        testCase_getAB <- with(testCase, 
            MBASED:::getAB(mu=test_mu, rho=test_rho)
        )
        checkEqualsNumeric(
            testCase$test_mu, 
            testCase_getMuRho$mu, 
            msg=paste('test_getMuRho_getAB: fail getMuRho is inverse of getAB test case #', testInd, sep='')
        )                                                                                                                 
        checkEqualsNumeric(
            testCase$test_rho, 
            testCase_getMuRho$rho, 
            msg=paste('test_getMuRho_getAB: fail getMuRho is inverse of getAB test case #',testInd, sep='')
        )
        checkEqualsNumeric(
            testCase$test_a, 
            testCase_getAB$a, 
            msg=paste('test_getMuRho_getAB: fail getMuRho is inverse of getAB test case #', testInd, sep='')
        )
        checkEqualsNumeric(
            testCase$test_b, 
            testCase_getAB$b, 
            msg=paste('test_getMuRho_getAB: fail getMuRho is inverse of getAB test case #', testInd, sep='')
        )
    }, BPPARAM=myBPPARAM)
    
    return(TRUE)
}

## don't want to require package VGAM, unless full testing should be done
if (FULLTESTING) { 
    library(VGAM)
    ## test vectorizedRbetabinomAB() function
    test_vectorizedRbetabinomAB <- function() {
        ## check1: generated counts have quantiles close to what we expect    
        set.seed(15775)
        quantilesToCheck <- seq(0.1, 0.9, by=0.1)
        numCountsToGenerate <- 10^6
        ## will compare to see if observed quantile is within 
        ## this many standard errors of the true value
        numSEsToCheck=5 
        testCasesMat <- expand.grid(
            test_mu=c(0.001, 0.3, 0.5, 0.7, 0.999),
            test_rho=c(0.001, 0.1, 1/3, 0.7, 0.999),
            test_size=c(1,2,10,100)
        )
        testCasesMat$test_a <- with(testCasesMat, 
            MBASED:::getAB(
                mu=test_mu, 
                rho=test_rho
            )$a
        )
        testCasesMat$test_b <- with(testCasesMat, 
            MBASED:::getAB(
                mu=test_mu, 
                rho=test_rho
            )$b
        )
        testResults <- bplapply(1:nrow(testCasesMat), function(testInd) { 
            testCase <- testCasesMat[testInd,]
            countsAB <- with(testCase, 
                MBASED:::vectorizedRbetabinomAB(
                    n=numCountsToGenerate, 
                    size=rep(test_size, numCountsToGenerate),
                    a=rep(test_a, numCountsToGenerate), 
                    b=rep(test_b, numCountsToGenerate)
                )
            )        
            ## distribution may be chunky, 
            ##so some quantiles might correspond to same number.
            countsABQuantileCounts <- unique(
                unname(quantile(countsAB, quantilesToCheck))
            ) 
            obsFs <- sapply(
                countsABQuantileCounts, 
                function(countsABQuantile) {
                    mean(countsAB<=countsABQuantile)
                }
            )
            predFs <- with(testCase, 
                pbetabinom.ab( 
                    q=countsABQuantileCounts, 
                    size=test_size, 
                    shape1=test_a, 
                    shape2=test_b
                )
            )
            testResults <- MBASED:::testQuantiles(
                theoreticalCumDist=predFs, 
                observedCumDist=obsFs, 
                numTotalCounts=numCountsToGenerate, 
                numSEsToCheck=numSEsToCheck, 
                errorMessage=paste('test_vectorizedRbetabinomAB: fail vectorizedRbetabinomAB quantile test case #',testInd, sep='')
            )
        }, BPPARAM=myBPPARAM)
        
        ## check2: vectorization works 
        set.seed(631844)
        countsABSizeScrambleVec <- sample(
            1:(numCountsToGenerate*nrow(testCasesMat)), 
            replace=FALSE
        )
        countsABSizeUnscrambleVec <- order(countsABSizeScrambleVec)
        countsABSizeScrambled <- with(testCasesMat, 
            MBASED:::vectorizedRbetabinomAB(
                n=numCountsToGenerate*nrow(testCasesMat), 
                size=rep(
                    test_size, 
                    each=numCountsToGenerate
                )[countsABSizeScrambleVec],
                a=rep(
                    test_a, 
                    each=numCountsToGenerate
                )[countsABSizeScrambleVec], 
                b=rep(
                    test_b, 
                    each=numCountsToGenerate
                )[countsABSizeScrambleVec]
            )
        )
        testResults <- bplapply(1:nrow(testCasesMat), function(testInd) { 
            testCase <- testCasesMat[testInd,]
            ## check AB counts generation first
            vecStartInd <- (testInd-1)*numCountsToGenerate+1
            vecEndInd <- testInd*numCountsToGenerate
            countsAB <- countsABSizeScrambled[countsABSizeUnscrambleVec]
            countsAB <- countsAB[vecStartInd:vecEndInd]
            ## distribution may be chunky, 
            ## so some quantiles might correspond to same number.
            countsABQuantileCounts <- unique(
                unname(quantile(countsAB, quantilesToCheck))
            ) 
            obsFs <- sapply(
                countsABQuantileCounts, 
                function(countsABQuantile) {
                    mean(countsAB<=countsABQuantile)
                }
            )
            predFs <- with(testCase, 
                pbetabinom.ab( 
                    q=countsABQuantileCounts, 
                    size=test_size, 
                    shape1=test_a, 
                    shape2=test_b
                )
            )
            testResults <- MBASED:::testQuantiles(
                theoreticalCumDist=predFs, 
                observedCumDist=obsFs, 
                numTotalCounts=numCountsToGenerate, 
                numSEsToCheck=numSEsToCheck, 
                errorMessage=paste('test_vectorizedRbetabinomAB: fail vectorizedRbetabinomAB vectorization test case #',testInd, sep='')
            )
        }, BPPARAM=myBPPARAM)
    
        return(TRUE)
    }
}    

if (FULLTESTING) { ## don't want to require package VGAM, unless full testing should be done
    library(VGAM)
    ## test vectorizedRbetabinomMR() function
    test_vectorizedRbetabinomMR <- function() {
        ## check1:  generated counts have quantiles close to what we expect   
        set.seed(913943)
        quantilesToCheck <- seq(0.1, 0.9, by=0.1)
        numCountsToGenerate <- 10^6
        ## will compare to see if observed quantile is within 
        ## this many standard errors of the true value
        numSEsToCheck=5 
        testCasesMat <- expand.grid(
            test_mu=c(0.001, 0.3, 0.5, 0.7, 0.999),
            test_rho=c(0.001, 0.1, 1/3, 0.7, 0.999),
            test_size=c(1,2,10,100)
        )
        testResults <- bplapply(1:nrow(testCasesMat), function(testInd) { 
            testCase <- testCasesMat[testInd,]
            countsMR <- with(testCase, 
                MBASED:::vectorizedRbetabinomMR(
                    n=numCountsToGenerate, 
                    size=rep(test_size, numCountsToGenerate), 
                    mu=rep(test_mu, numCountsToGenerate), 
                    rho=rep(test_rho, numCountsToGenerate)
                )
            )
            countsMRQuantileCounts <- unique(
                unname(quantile(countsMR, quantilesToCheck))
            )
            obsFs <- sapply(
                countsMRQuantileCounts, 
                function(countsMRQuantile) {
                    mean(countsMR<=countsMRQuantile)
                }
            )
            predFs <- with(testCase, 
                pbetabinom( 
                    q=countsMRQuantileCounts, 
                    size=test_size, 
                    prob=test_mu, 
                    rho=test_rho
                )
            )
            testResults <- MBASED:::testQuantiles(
                theoreticalCumDist=predFs, 
                observedCumDist=obsFs, 
                numTotalCounts=numCountsToGenerate, 
                numSEsToCheck=numSEsToCheck, 
                errorMessage=paste('test_vectorizedRbetabinomMR: fail vectorizedRbetabinomMR quantile test case #',testInd, sep='')
            )
        }, BPPARAM=myBPPARAM)
        
        ## check2: vectorization works for beta-binomial mode
        set.seed(497556)
        countsMuRhoSizeScrambleVec <- sample(
            1:(numCountsToGenerate*nrow(testCasesMat)), 
            replace=FALSE
        )
        countsMuRhoSizeUnscrambleVec <- order(countsMuRhoSizeScrambleVec)
        countsMuRhoSizeScrambled <- with(testCasesMat, 
            MBASED:::vectorizedRbetabinomMR(
                n=numCountsToGenerate*nrow(testCasesMat), 
                size=rep(
                    test_size, 
                    each=numCountsToGenerate
                )[countsMuRhoSizeScrambleVec],
                mu=rep(
                    test_mu, 
                    each=numCountsToGenerate
                )[countsMuRhoSizeScrambleVec], 
                rho=rep(
                    test_rho, 
                    each=numCountsToGenerate
                )[countsMuRhoSizeScrambleVec]
            )
        )
        testResults <- bplapply(1:nrow(testCasesMat), function(testInd) { 
            testCase <- testCasesMat[testInd,]
            vecStartInd <- (testInd-1)*numCountsToGenerate+1
            vecEndInd <- testInd*numCountsToGenerate
            countsMR <- countsMuRhoSizeScrambled[countsMuRhoSizeUnscrambleVec]
            countsMR <- countsMR[vecStartInd:vecEndInd]
            countsMRQuantileCounts <- unique(
                unname(quantile(countsMR, quantilesToCheck))
            )
            obsFs <- sapply(
                countsMRQuantileCounts, 
                function(countsMRQuantile) {
                    mean(countsMR<=countsMRQuantile)
                }
            )
            predFs <- with(testCase, 
                pbetabinom( 
                    q=countsMRQuantileCounts, 
                    size=test_size, 
                    prob=test_mu, 
                    rho=test_rho
                )
            )
            testResults <- MBASED:::testQuantiles(
                theoreticalCumDist=predFs, 
                observedCumDist=obsFs, 
                numTotalCounts=numCountsToGenerate, 
                numSEsToCheck=numSEsToCheck, 
                errorMessage=paste('test_vectorizedRbetabinomMR: fail vectorizedRbetabinomMR vectorization test case #',testInd, sep='')
            )
        }, BPPARAM=myBPPARAM)
    
        ## check3: vectorization works for binomial and beta-binomial modes 
        set.seed(327526)
        testCasesMat <- expand.grid(
            test_mu=c(0.001, 0.3, 0.5, 0.7, 0.999),
            test_rho=c(0, 0.1, 1/3,  0.999),
            test_size=c(1,2,100)
        )
        countsMuRhoSizeScrambleVec <- sample(
            1:(numCountsToGenerate*nrow(testCasesMat)), 
            replace=FALSE
        )
        countsMuRhoSizeUnscrambleVec <- order(countsMuRhoSizeScrambleVec)
        countsMuRhoSizeScrambled <- with(testCasesMat, 
            MBASED:::vectorizedRbetabinomMR(
                n=numCountsToGenerate*nrow(testCasesMat), 
                size=rep(
                    test_size, 
                    each=numCountsToGenerate
                )[countsMuRhoSizeScrambleVec],
                mu=rep(
                    test_mu, 
                    each=numCountsToGenerate
                )[countsMuRhoSizeScrambleVec], 
                rho=rep(
                    test_rho, 
                    each=numCountsToGenerate
                )[countsMuRhoSizeScrambleVec]
            )
        )
        testResults <- bplapply(1:nrow(testCasesMat), function(testInd) { 
            testCase <- testCasesMat[testInd,]
            countsMR <- countsMuRhoSizeScrambled[countsMuRhoSizeUnscrambleVec]
            countsMRStartPos <- ((testInd-1)*numCountsToGenerate+1)
            countsMREndPos <- (testInd*numCountsToGenerate)
            countsMR <- countsMR[countsMRStartPos:countsMREndPos]
            countsMRQuantileCounts <- unique(
                unname(quantile(countsMR, quantilesToCheck))
            )
            obsFs <- sapply(
                countsMRQuantileCounts, 
                function(countsMRQuantile) {
                    mean(countsMR<=countsMRQuantile)
                }
            )
            predFs <- with(testCase, 
                if (isTRUE(all.equal(test_rho, 0))) {
                    pbinom(
                        q=countsMRQuantileCounts, 
                        size=test_size, 
                        prob=test_mu, 
                    )
                } else {
                    pbetabinom( 
                        q=countsMRQuantileCounts, 
                        size=test_size, 
                        prob=test_mu, 
                        rho=test_rho
                    )
                }
            )
            testResults <- MBASED:::testQuantiles(
                theoreticalCumDist=predFs, 
                observedCumDist=obsFs, 
                numTotalCounts=numCountsToGenerate, 
                numSEsToCheck=numSEsToCheck, 
                errorMessage=paste('test_vectorizedRbetabinomMR: fail vectorizedRbetabinomMR vectorization test case #',testInd, sep='')
            )
        }, BPPARAM=myBPPARAM)
    
        return(TRUE)
    }
}    


    