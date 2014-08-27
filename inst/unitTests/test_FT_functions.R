## size of round-off error I'm tolerating 
## (quantities with differences of this much or smaller are considered the same)
myPrecision <- .Machine$double.eps^0.5 
FULLTESTING <- FALSE
if (FULLTESTING) {
    myBPPARAM <- MulticoreParam(workers=24)
} else {
    myBPPARAM <- SerialParam()
}

## test FT() function
test_FT <- function() {
    ## check1: FT is monotone increasing function of x
    test_n <- 100
    test_FT_values <- MBASED:::FT(
        x=0:test_n, 
        n=rep(test_n, length(0:test_n))
    )
    checkEqualsNumeric(
        1:(test_n+1),
        order(test_FT_values), 
        msg='test_FT: fail FT monotonic test'
    )
    checkTrue(
        MBASED:::FT(x=0,n=1)<MBASED:::FT(x=1,n=1), 
        msg='test_FT: fail FT n=1 test'
    ) ## also tests border case of n=1
    
    ## check2: FT(n/2,n)=pi/2
    test_n <- 100
    checkEqualsNumeric(
        MBASED:::FT(x=(test_n/2), n=test_n), 
        pi/2,
        msg='test_FT: fail FT pi/2 test'
    )
    
    ## check3: vectorization works
    ## test for matrices 
    set.seed(777891)
    numCases <- 100
    numCols <- 5
    numRows <- numCases/numCols
    test_xs <- sample(0:10, numCases, replace=T)
    test_ns <- sample(10:20, numCases, replace=T)
    test_xs_mat <- matrix(
        test_xs, 
        ncol=numCols
    )
    test_ns_mat <- matrix(
        test_ns, 
        ncol=numCols
    )
    individualResults <- matrix(
        NA,
        nrow=numRows, 
        ncol=numCols
    )
    for (i in 1:nrow(individualResults)) {
        for (j in 1:ncol(individualResults)) {
            individualResults[i,j] <- MBASED:::FT(
                x=test_xs_mat[i,j], 
                n=test_ns_mat[i,j]
            )
        }
    }
    jointResults <- MBASED:::FT(
        x=test_xs_mat, 
        n=test_ns_mat
    )
    checkEqualsNumeric(
        individualResults, 
        jointResults,
        msg='test_FT: fail FT vectorization test'
    )
    ## also test for vectors
    individualResults <- sapply(1:numCases, function(ind) {
        MBASED:::FT(
            x=test_xs[ind], 
            n=test_ns[ind]
        )
    })
    jointResults <- MBASED:::FT(
        x=test_xs, 
        n=test_ns
    )
    checkEqualsNumeric(
        jointResults, 
        individualResults, 
        msg='test_FT: fail FT vectorization test'
    )
    
    
    ## check4: if n is large and x~Bin(n,p), 
    ## then FT(x,n) ~ N(2*arcsin(sqrt(p)), 1/(n+0.5))
    test_n <- 10^5
    test_ps <- seq(0.1, 0.9, by=0.1)
    if (FULLTESTING) {
        numCountsToGenerate <- 10^6
    } else {
    	## need this to be high to make quantiles close to truth
        numCountsToGenerate <- 10^6 
    }
    set.seed(443678)
    testResults <- bplapply(1:length(test_ps), function(testInd) {
        test_p <- test_ps[testInd]
        test_counts <- rbinom(
            n=numCountsToGenerate, 
            size=test_n, 
            prob=test_p
        )
        test_FT_values <- MBASED:::FT(
            x=test_counts, 
            n=rep(test_n, length(test_counts))
        )
        ## check mean
        obsMean <- mean(test_FT_values)
        predMean <- 2*asin(sqrt(test_p))
        checkTrue(
            MBASED:::testNumericDiff(
                queryVals=obsMean, 
                targetVals=predMean, 
                cutoffFraction=0.01 ## within 1%
            ),
            msg=paste('test_FT: fail predicted mean test for p=', test_p, sep='')
        ) 
        ## check variance
        obsVar <- var(test_FT_values)
        predVar <- 1/(test_n+0.5)
        checkTrue(
            MBASED:::testNumericDiff(
                queryVals=obsVar, 
                targetVals=predVar, 
                cutoffFraction=0.01 ## within 1%
            ),
            msg=paste('test_FT: fail predicted mean test for p=', test_p, sep='')
        ) 
        ## check normality
        ## exclude 0.5 to avoid quantile value of 0
        quantilesToCheck <- setdiff(seq(0.1, 0.9, by=0.1), 0.5) 
        test_z_values <- (test_FT_values-predMean)/sqrt(predVar)
        obsQuantiles <- unname(quantile(test_z_values, quantilesToCheck))
        predQuantiles <- qnorm(quantilesToCheck)
        checkTrue(
            all(
                MBASED:::testNumericDiff(
                    queryVals=obsQuantiles, 
                    targetVals=predQuantiles, 
                    cutoffFraction=0.05 ## within 5%
                )
            ),
            msg=paste('test_FT: fail normality test for p=', test_p, sep='')
        )
    }, BPPARAM=myBPPARAM)
    
    return(TRUE)
}    
    
## test unFT() function
test_unFT <- function() {    
    ## check1: unFT is the reverse of FT
    test_n <- 100
    testResults <- bplapply(0:test_n, function(testInd) {
        test_x <- testInd
        trueProps <- test_x/test_n
        FTProps <- MBASED:::unFT(
            z=MBASED:::FT(
                x=test_x, 
                n=test_n
            ), 
            n=test_n
        )
        checkEqualsNumeric(
            FTProps, 
            trueProps, 
            msg='test_unFT: fail unFT is reverse of FT test'
        )
    }, BPPARAM=myBPPARAM)
    
    ## check2: unFT handles situations with extreme z values
    test_n <- 100
    maxZ <- MBASED:::FT(
        x=test_n,
        n=test_n
    )
    minZ <- MBASED:::FT(
        x=0,
        n=test_n
    )
    checkEqualsNumeric(
        MBASED:::unFT(
            z=maxZ, 
            n=test_n
        ), 
        MBASED:::unFT(
            z=maxZ+10, 
            n=test_n
        ), 
        msg='test_unFT: fail extreme value test'
    )
    checkEqualsNumeric(
        MBASED:::unFT(
            z=minZ, 
            n=test_n
        ), 
        MBASED:::unFT(
            z=minZ-10, 
            n=test_n
        ), 
        msg='test_unFT: fail extreme value test'
    )

    ## check3: vectorization works
    set.seed(863249)
    ## check for matrices
    numCases <- 1000
    numCols <- 20
    numRows <- numCases/numCols
    test_zs <- runif(numCases, min=(-pi), max=pi)
    test_ns <- sample(10:20, numCases, replace=T)
    test_zs_mat <- matrix(
        test_zs, 
        ncol=numCols
    )
    test_ns_mat <- matrix(
        test_ns, 
        ncol=numCols
    )
    individualResults <- matrix(
        NA,
        nrow=numRows, 
        ncol=numCols
    )
    for (i in 1:nrow(individualResults)) {
        for (j in 1:ncol(individualResults)) {
            individualResults[i,j] <- MBASED:::unFT(
                z=test_zs_mat[i,j], 
                n=test_ns_mat[i,j]
            )
        }
    }
    jointResults <- MBASED:::unFT(
        z=test_zs_mat, 
        n=test_ns_mat
    )
    checkEqualsNumeric(
        jointResults,
        individualResults,
        msg='test_unFT: fail unFT vectorization test'
    ) 
    ## check vectorization for vectors
    individualResults <- unlist(bplapply(1:numCases, function(testInd) {
        MBASED:::unFT(
            z=test_zs[testInd], 
            n=test_ns[testInd]
        )
    }, BPPARAM=myBPPARAM))
    jointResults <- MBASED:::unFT(
        z=test_zs, 
        n=test_ns
    )
    checkEqualsNumeric(
        jointResults,
        individualResults,
        msg='test_unFT: fail unFT vectorization test'
    ) 
    
    return(TRUE)
}
    
## test FTAdjust() function
test_FTAdjust <- function() {        
    ## check1: vectorization works
    ## check for matrices
    set.seed(324515)
    numCases <- 100
    numCols <- 20
    numRows <- numCases/numCols
    test_xs <- sample(0:10, numCases, replace=T)
    test_ns <- sample(10:20, numCases, replace=T)
    test_ps <- sample(seq(0.1, 0.9, by=0.1), numCases, replace=T)
    test_xs_mat <- matrix(
        test_xs, 
        ncol=numCols
    )
    test_ns_mat <- matrix(
        test_ns, 
        ncol=numCols
    )
    test_ps_mat <- matrix(
        test_ps, 
        ncol=numCols
    )
    individualResults <- matrix(
        NA,
        nrow=numRows, 
        ncol=numCols
    )
    for (i in 1:nrow(individualResults)) {
        for (j in 1:ncol(individualResults)) {
            individualResults[i,j] <- MBASED:::FTAdjust(
                x=test_xs_mat[i,j], 
                n=test_ns_mat[i,j], 
                p=test_ps_mat[i,j]
            )
        }
    }
    jointResults <- MBASED:::FTAdjust(
        x=test_xs_mat, 
        n=test_ns_mat, 
        p=test_ps_mat
    )
    checkEqualsNumeric(
        jointResults,
        individualResults,
        msg='test_FTAdjust: fail FTAdjust vectorization test'
    ) 
    ## also check for vectors
    individualResults <- sapply(1:numCases, function(ind) {
        MBASED:::FTAdjust(
            x=test_xs[ind], 
            n=test_ns[ind], 
            p=test_ps[ind]
        )
    })
    jointResults <- MBASED:::FTAdjust(
        x=test_xs, 
        n=test_ns, 
        p=test_ps
    )
    checkEqualsNumeric(
        jointResults,
        individualResults,
        msg='test_FTAdjust: fail FTAdjust vectorization test'
    ) 
    
    ## check2: adjustment is exact when p=0.5
    test_n <- 100
    unadjustedFTs <- MBASED:::FT(
        x=0:test_n, 
        n=rep(test_n, test_n+1)
    )
    adjustedFTs <- MBASED:::FTAdjust(
        x=0:test_n, 
        n=rep(test_n, test_n+1), 
        p=rep(0.5, test_n+1)
    )
    checkEqualsNumeric(
        unadjustedFTs,
        adjustedFTs,
        msg='test_FTAdjust: fail FTAdjust p=0.5 test'
    ) 
    
    ## check3: adjustment actually does what it's supposed to do:
    test_n <- 100
    test_ps <- seq(0.1, 0.9, by=0.1)
    for (test_p in test_ps) {
        for (test_x in 0:test_n) {
            automatedResult <- MBASED:::FTAdjust(
                x=test_x, 
                n=test_n, 
                p=test_p
            )
            manualResult <- MBASED:::FT(
                x=test_x, 
                n=test_n
            )-2*asin(sqrt(test_p))+2*asin(sqrt(0.5))
            checkEqualsNumeric(
                automatedResult, 
                manualResult, 
                msg='test_FTAdjust: fail FTAdjust p!=0.5 test'
            )
        }
    }
    
    ## check4: adjustment works ok for mild levels of 
    ## distorting ps and for high counts
    test_n=10^4
    test_ps=seq(0.4, 0.6, by=0.05)
    test_props=seq(0.1, 0.9, by=0.1)
    for (test_p in test_ps) {
        for (test_prop in test_props) {
        	test_generating_p_num <- test_p*test_prop
        	test_generating_p_denom <- (test_p*test_prop+(1-test_p)*(1-test_prop))
            test_generating_p <- test_generating_p_num/test_generating_p_denom
            predProp <- test_prop
            obsProp <- MBASED:::unFT(
                z=MBASED:::FTAdjust(
                    x=test_n*test_generating_p, 
                    n=test_n, 
                    p=test_p
                ), 
                n=test_n
            )
            checkTrue(
                (abs(obsProp-predProp)<=0.05) ||
                MBASED:::testNumericDiff(
                    queryVals=obsProp, 
                    targetVals=predProp, 
                    cutoffFraction=0.05 ## within 5%
                ),
                msg='test_FTAdjust: fail FTAdjust check if works for reasonable values of distortion'
            ) 
        }
    }
        
    return(TRUE)
}
    
## test isCountMajor() function
test_isCountMajor<- function() {    
    ## check1: vectorization works (no breaking ties randomly)
    set.seed(804225)
    ## check for matrices
    numCases <- 100
    numCols <- 20
    numRows <- numCases/numCols
    test_xs <- sample(0:10, numCases, replace=T)
    test_ns <- sample(10:20, numCases, replace=T)
    test_ps <- sample(seq(0.1, 0.9, by=0.1), numCases, replace=T)
    test_xs_mat <- matrix(
        test_xs, 
        ncol=numCols
    )
    test_ns_mat <- matrix(
        test_ns, 
        ncol=numCols
    )
    test_ps_mat <- matrix(
        test_ps, 
        ncol=numCols
    )
    individualResults <- matrix(
        NA,
        nrow=numRows, 
        ncol=numCols
    )
    for (i in 1:nrow(individualResults)) {
        for (j in 1:ncol(individualResults)) {
            individualResults[i,j] <- MBASED:::isCountMajorFT(
                x=test_xs_mat[i,j], 
                n=test_ns_mat[i,j], 
                p=test_ps_mat[i,j], 
                tieBreakRandom=FALSE
            )
        }
    }
    jointResults <- MBASED:::isCountMajorFT(
        x=test_xs_mat, 
        n=test_ns_mat, 
        p=test_ps_mat, 
        tieBreakRandom=FALSE
    )    
    checkIdentical(
        jointResults,
        individualResults,
        msg='test_isCountMajor: fail isCountMajorFT vectorization test'
    ) 
    ## also check for vectors
    individualResults <- sapply(1:numCases, function(ind) {
        MBASED:::isCountMajorFT(
            x=test_xs[ind], 
            n=test_ns[ind], 
            p=test_ps[ind], 
            tieBreakRandom=FALSE
        )
    })
    jointResults <- MBASED:::isCountMajorFT(
        x=test_xs, 
        n=test_ns, 
        p=test_ps, 
        tieBreakRandom=FALSE
    )
    checkIdentical(
        jointResults,
        individualResults,
        msg='test_isCountMajor: fail isCountMajorFT vectorization test'
    ) 

    ## check2: isCountMajor does what it's supposed to do 
    ## (with possible exception of borderline cases):
    test_ns <- 1:100
    test_ps <- seq(0.1, 0.9, by=0.1)
    for (test_p in test_ps) {
        for (test_n in test_ns) {
        	## the  cases at the breakpoint can go either way 
        	## due to roundoff during adjustment
            test_xs <- setdiff(0:test_n, c(floor(test_p*test_n), ceiling(test_p*test_n))) 
            test_xs_major <- test_xs[test_xs>( test_p*test_n)]
            test_xs_minor <- test_xs[test_xs<( test_p*test_n)]
            test_xs_major_results <- MBASED:::isCountMajorFT(
                x=test_xs_major, 
                n=rep(test_n, length(test_xs_major)), 
                p=rep(test_p, length(test_xs_major)), 
                tieBreakRandom=FALSE
            )
            test_xs_minor_results <- MBASED:::isCountMajorFT(
                x=test_xs_minor, 
                n=rep(test_n, length(test_xs_minor)), 
                p=rep(test_p, length(test_xs_minor)), 
                tieBreakRandom=FALSE
            )
            checkTrue(
                all(test_xs_major_results), 
                msg='test_isCountMajor: fail isCountMajorFT basic functionality test'
            )
            checkTrue(
                all(!test_xs_minor_results), 
                msg='test_isCountMajor: fail isCountMajorFT basic functionality test'
            )
        }
    }
    
    ## check3: the breaking of the ties works
    numCases <- 10^6
    midVal <- 5
    test_results_no_tie_break <- MBASED:::isCountMajorFT(
        x=rep(midVal, numCases), 
        n=rep(2*midVal, numCases), 
        p=rep(0.5, numCases), 
        tieBreakRandom=FALSE
    )
    checkTrue(
        ## all have same assignment
        length(unique(test_results_no_tie_break))==1, 
        msg='test_isCountMajor: fail isCountMajorFT tie-breaking test'
    ) 
    test_results_tie_break <- MBASED:::isCountMajorFT(
        x=rep(midVal, numCases), 
        n=rep(2*midVal, numCases), 
        p=rep(0.5, numCases), 
        tieBreakRandom=TRUE
    ) ## should be about 50% true, 50% false. 
    MBASED:::testQuantiles(
        theoreticalCumDist=0.5, 
        observedCumDist=mean(test_results_tie_break), 
        numTotalCounts=numCases, 
        numSEsToCheck=5, 
        errorMessage='test_isCountMajor: fail isCountMajorFT tie-breaking test'
    ) 
    
    return(TRUE)
}    
    
