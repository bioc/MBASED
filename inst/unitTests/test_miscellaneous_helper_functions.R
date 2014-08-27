## size of round-off error I'm tolerating 
## (quantities with differences of this much or smaller 
## are considered the same)

myPrecision <- .Machine$double.eps^0.5 
FULLTESTING <- FALSE
if (FULLTESTING) {
    myBPPARAM <- MulticoreParam(workers=24)
} else {
    myBPPARAM <- SerialParam()
}

## test getSimulationPvalue() function
test_getSimulationPvalue <- function() {
    testVals <- c(1:10, 2)
    manualPvalsLess <- sapply(0:11, function(tVal) {
        if (tVal<2) {
            max(0,tVal/11)
        } else {
            min(1,(tVal+1)/11)
        }
    })
    manualPvalsGreater <- sapply(0:11, function(tVal) {
        if (tVal<=2) {
            min(1,(11-tVal+1)/11)
        } else {
            max(0,(11-tVal)/11)
        }
    })
    automatedPvalsLess <- sapply(0:11, function(tVal) {
        MBASED:::getSimulationPvalue(
            observedVal=tVal, 
            simulatedVals=testVals, 
            direction='less'
        )
    })
    automatedPvalsGreater <- sapply(0:11, function(tVal) {
        MBASED:::getSimulationPvalue(
            observedVal=tVal, 
            simulatedVals=testVals, 
            direction='greater'
        )
    })
    checkEqualsNumeric(
        manualPvalsLess, 
        automatedPvalsLess,
        msg='test_getSimulationPvalue: fail p-values less test'
    )
    checkEqualsNumeric(
        manualPvalsGreater, 
        automatedPvalsGreater,
        msg='test_getSimulationPvalue: fail p-values greater test'
    )
    
    return(TRUE)
}

## test shiftAndAttenuateProportions() function
test_shiftAndAttenuateProportions <- function() {
    ## check1: I will simply check that it does what it's meant to 
    ## for the case of no bias, no over-dispersion
    set.seed(312003)
    numCases <- 100
    for (testInd in numCases) {
        for (numLoci in c(1,2)) {
            for (numCols in 1:2) {
                totalCounts <- sample(10:100, numLoci*numCols) 
                allele1Counts <- rbinom(
                    n=length(totalCounts), 
                    size=totalCounts, 
                    prob=0.5
                )
                countsMat=matrix(
                    allele1Counts, 
                    ncol=numCols
                )
                totalsMat=matrix(
                    totalCounts, 
                    ncol=numCols
                )
                probsMat=matrix(
                    rep(0.5, length(totalCounts)), 
                    ncol=numCols
                )
                rhosMat=matrix(
                    rep(0, length(totalCounts)), 
                    ncol=numCols
                )
                output <- MBASED:::shiftAndAttenuateProportions(
                    countsMat=countsMat, 
                    totalsMat=totalsMat, 
                    probsMat=probsMat,  
                    rhosMat=rhosMat
                ) 
                propsShiftedManual <- countsMat/totalsMat
                propsAttenuatedManual <- (countsMat+0.5)/(totalsMat+1)
                propsShiftedVarsManual <- 
                    propsAttenuatedManual*(1-propsAttenuatedManual)/(totalsMat+1)
                checkEqualsNumeric(
                    propsShiftedManual,
                    output$propsShifted,
                    msg='test shiftAndAttenuateProportions: fail manual comparison test'
                )
                checkEqualsNumeric(
                    dim(propsShiftedManual),
                    dim(output$propsShifted),
                    msg='test shiftAndAttenuateProportions: fail manual comparison test'
                )
                checkEqualsNumeric(
                    propsShiftedVarsManual,
                    output$propsShiftedVars,
                    msg='test shiftAndAttenuateProportions: fail manual comparison test'
                )
                checkEqualsNumeric(
                    dim(propsShiftedVarsManual),
                    dim(output$propsShiftedVars),
                    msg='test shiftAndAttenuateProportions: fail manual comparison test'
                )
            }
        }
    }
    
    return(TRUE)
}
    
## test getPFinal() function
test_getPFinal <- function() {
    ## check1: if true AF is 1, pFinal should also be 1
    always1 <- MBASED:::getPFinal(
        trueAF=1, 
        noASEAF=seq(0.1, 0.9, by=0.1)
    )
    checkTrue(
        all(always1==1),
        msg='test_getPFinal: fail all equal to 1 test'
    )
    
    ## check2: if true AF is 0, pFinal should also be 0
    always0 <- MBASED:::getPFinal(
        trueAF=0, 
        noASEAF=seq(0.1, 0.9, by=0.1)
    )
    checkTrue(
        all(always0==0),
        msg='test_getPFinal: fail all equal to 0 test'
    )
    
    ## check3: if no pre-existing alleic bias, pFinal should be trueAF
    noBiasInitialProbs <- seq(0.1, 0.9, by=0.1)
    noBiasFinalProbs <- sapply(noBiasInitialProbs, function(prob) {
        MBASED:::getPFinal(
            trueAF=prob, 
            noASEAF=0.5
        )
    })
    checkEqualsNumeric(
        noBiasInitialProbs,
        noBiasFinalProbs,
        msg='test_getPFinal: fail no bias test'
    )
    
    ## check4: one manual check for pre-existing allelic bias
    test_trueAF <- 0.3
    test_noASEAFs  <- c(0.2, 0.9)
    biasProbsManual <- sapply(
        test_noASEAFs, 
        function(test_noASEAF) {
            test_noASEAF*test_trueAF/
                (test_noASEAF*test_trueAF + (1-test_noASEAF)*(1-test_trueAF))
        }
    )
    biasProbsAutomated <- MBASED:::getPFinal(
        trueAF=test_trueAF, 
        noASEAF=test_noASEAFs
    )
    checkEqualsNumeric(
        biasProbsManual,
        biasProbsAutomated,
        msg='test_getPFinal: fail manual test'
    )
    
    return(TRUE)
}


## test logLikelihoodCalculator1s() function
test_logLikelihoodCalculator1s <- function() {
    ## rudimentary manual check
    test_noASEProbs <- c(0.2, 0.5)
    test_rhos <- c(0.01, 0.1)
    test_trueHapA <- 0.8
    test_generatingProbs <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbs
    ) ## this particular combination gives 0.5 and 0.8 as return probs
    totalCounts <- c(200, 1000) ## high counts lead to greater confidence
    meanCounts <- test_generatingProbs*totalCounts
    LLC <- MBASED:::logLikelihoodCalculator1s(
        lociHapACounts=meanCounts, 
        lociTotalCounts=totalCounts, 
        lociHapANoASEProbs=test_noASEProbs, 
        lociRhos=test_rhos
    )
    neighborhoodRange <- 0.01
    probsToTest <- test_trueHapA + 
        c(
            -neighborhoodRange, 
            0, 
            neighborhoodRange
        )
    neighborhoodLLCVals <- sapply(probsToTest, LLC)
    checkTrue(
        which.max(neighborhoodLLCVals)==2,
        msg='test_logLikelihoodCalculator1s: fail manual test'
    )
    
    return(TRUE)
}

    

## test maxLogLikelihoodCalculator1s() function
test_maxLogLikelihoodCalculator1s <- function() {
    ## rudimentary manual check
    test_noASEProbs <- c(0.2, 0.5)
    test_rhos <- c(0.01, 0.1)
    test_trueHapA <- 0.8
    test_generatingProbs <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbs
    ) ## this particular combination gives 0.5 and 0.8 as return probs
    ## high counts lead to greater confidence
    totalCounts <- c(200, 1000) 
    meanCounts <- test_generatingProbs*totalCounts
    MLL <- MBASED:::maxLogLikelihoodCalculator1s(
        lociHapACounts=meanCounts, 
        lociTotalCounts=totalCounts, 
        lociHapANoASEProbs=test_noASEProbs, 
        lociRhos=test_rhos
    )
    neighborhoodRange <- 0.01
    probsToTest <- test_trueHapA + 
        c(
            -neighborhoodRange, 
            0, 
            neighborhoodRange
        )
    distMLEEstimateToProbs <- abs(probsToTest-MLL$maximum)
    checkTrue(
        which.min(distMLEEstimateToProbs)==2,
        msg='maxLogLikelihoodCalculator1s: fail manual test'
    )
    MBASED:::testNumericDiff(
        MLL$maximum, 
        test_trueHapA, 
        cutoff=0.01 ## within 1% of true hapA frequency
    )
    
    return(TRUE)
}

    
## test estimateMAF1s() function
test_estimateMAF1s <- function() {
    ## generate data: both phased and non-phased analysis 
    ## should be performed on my toy example.
    test_trueHapA <- 0.8
    test_trueHapB <- 1-test_trueHapA
    test_noASEProbsHapA <- c(0.25, 0.5)
    test_noASEProbsHapB <- 1-test_noASEProbsHapA
    test_rhos <- c(0.01, 0.1)
    test_generatingProbsHapA <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbsHapA
    ) 
    ## high counts lead to greater confidence
    totalCounts <- c(200, 1000) 
    meanCountsHapA <- round(test_generatingProbsHapA*totalCounts)
    ## define phased and non-phased counts
    hapACountsPhased <- meanCountsHapA
    hapBCountsPhased <- totalCounts-meanCountsHapA
    allele1CountsUnphased <- c(
        hapACountsPhased[1], 
        hapBCountsPhased[2]
    )
    allele2CountsUnphased <- totalCounts-allele1CountsUnphased
    test_noASEProbsAllele1 <- c(
        test_noASEProbsHapA[1], 
        test_noASEProbsHapB[2]
    )
    test_noASEProbsAllele2 <- 1-test_noASEProbsAllele1
    outputEstimates <- vector('list')
    outputEstimates[['hapA']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF1s(
                lociAllele1Counts=hapACountsPhased, 
                lociTotalCounts=totalCounts, 
                lociAllele1NoASEProbs=test_noASEProbsHapA, 
                lociRhos=test_rhos,
                isPhased=test_isPhased
            )
        }
    )
    outputEstimates[['hapB']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF1s(
                lociAllele1Counts=hapBCountsPhased, 
                lociTotalCounts=totalCounts, 
                lociAllele1NoASEProbs=test_noASEProbsHapB, 
                lociRhos=test_rhos,
                isPhased=test_isPhased
            )
        }
    )
    outputEstimates[['allele1']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF1s(
                lociAllele1Counts=allele1CountsUnphased, 
                lociTotalCounts=totalCounts, 
                lociAllele1NoASEProbs=test_noASEProbsAllele1, 
                lociRhos=test_rhos,
                isPhased=test_isPhased
           )
        }
    )
    outputEstimates[['allele2']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF1s(
                lociAllele1Counts=allele2CountsUnphased, 
                lociTotalCounts=totalCounts, 
                lociAllele1NoASEProbs=test_noASEProbsAllele2, 
                lociRhos=test_rhos,
                isPhased=test_isPhased
            )
        }
    )
    for (outputInd in 1:length(outputEstimates)) {
        names(outputEstimates[[outputInd]]) <- c('phased', 'nonPhased')
    }

    ## check1: when data is unphased, all MAF estimates 
    ## should be identical, and major haplotype should be haplotypeA
    for (outputInd in 2:length(outputEstimates)) {
        checkEqualsNumeric(
            outputEstimates[[1]]$nonPhased$MAF,
            outputEstimates[[outputInd]]$nonPhased$MAF,
            msg='test_estimateMAF1s: fail unphased MAF test'
        )
    }
    checkIdentical(
        outputEstimates$hapA$nonPhased$allele1IsMajor,
        c(TRUE, TRUE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )
    checkIdentical(
        outputEstimates$hapB$nonPhased$allele1IsMajor,
        c(FALSE, FALSE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )
    checkIdentical(
        outputEstimates$allele1$nonPhased$allele1IsMajor,
        c(TRUE, FALSE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )
    checkIdentical(
        outputEstimates$allele2$nonPhased$allele1IsMajor,
        c(FALSE, TRUE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )

    ## check2: when data is phased, hapA and hapB results 
    #3 are the same as when it's unphased
    for (outputName in c('hapA', 'hapB')) {
        checkEquals(
            outputEstimates[[outputName]]$nonPhased,
            outputEstimates[[outputName]]$phased,
            msg='test_estimateMAF1s: fail unphased equals to phased for known haplotypes test'
        )
    }
    
    ## check3: when data is phased, allele1 and allele2 results 
    ## are different from before
    for (outputName in c('allele1', 'allele2')) {
        checkTrue(
            (outputEstimates[[outputName]]$nonPhased$MAF) != 
            (outputEstimates[[outputName]]$phased$MAF),
            msg='test_estimateMAF1s: fail unphased different from phased for unknown haplotypes test'
        )
    }
    
    ## check4: when data is phased, allele1 and allele2 results are 
    ## opposites of each other in terms of haplotype assignment 
    ## but show same MAF
    checkEqualsNumeric(
        outputEstimates$allele1$nonPhased$MAF,
        outputEstimates$allele2$nonPhased$MAF,
        msg='test_estimateMAF1s: fail unphased test'
    )
    checkIdentical(
        outputEstimates$allele1$nonPhased$allele1IsMajor,
        !(outputEstimates$allele2$nonPhased$allele1IsMajor),
        msg='test_estimateMAF1s: fail unphased test'
    )
    
    return(TRUE)
}    
    
    
## test logLikelihoodCalculator2s() function
test_logLikelihoodCalculator2s <- function() {
    ## rudimentary manual check
    test_noASEProbs1s <- c(0.2, 0.5)
    test_noASEProbs2s <- c(0.3, 0.7)
    test_rhos1s <- c(0.01, 0.1)
    test_rhos2s <- c(0.005, 0.05)
    test_trueHapA <- 0.8
    test_generatingProbs1s <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbs1s
    ) 
    test_generatingProbs2s <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbs2s
    ) 
    ## high counts lead to greater confidence
    totalCounts1s <- c(200, 1000) 
    totalCounts2s <- c(900, 600) 
    meanCounts1s <- round(test_generatingProbs1s*totalCounts1s)
    meanCounts2s <- round(test_generatingProbs2s*totalCounts2s)
    LLC <- MBASED:::logLikelihoodCalculator2s(
        lociHapACountsSample1=meanCounts1s, 
        lociTotalCountsSample1=totalCounts1s, 
        lociHapANoASEProbsSample1=test_noASEProbs1s, 
        lociRhosSample1=test_rhos1s,
        lociHapACountsSample2=meanCounts2s, 
        lociTotalCountsSample2=totalCounts2s, 
        lociHapANoASEProbsSample2=test_noASEProbs2s, 
        lociRhosSample2=test_rhos2s
    )
    neighborhoodRange <- 0.01
    probsToTest <- test_trueHapA + 
        c(
            -neighborhoodRange, 
            0, 
            neighborhoodRange
        )
    neighborhoodLLCVals <- sapply(probsToTest, LLC)
    checkTrue(
        which.max(neighborhoodLLCVals)==2,
        msg='test_logLikelihoodCalculator2s: fail manual test'
    )
    
    return(TRUE)
}
    
## test maxLogLikelihoodCalculator2s() function
test_maxLogLikelihoodCalculator2s <- function() {
    ## rudimentary manual check
    test_noASEProbs1s <- c(0.2, 0.5)
    test_noASEProbs2s <- c(0.3, 0.7)
    test_rhos1s <- c(0.01, 0.1)
    test_rhos2s <- c(0.005, 0.05)
    test_trueHapA <- 0.8
    test_generatingProbs1s <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbs1s
    ) 
    test_generatingProbs2s <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbs2s
    ) 
    ## high counts lead to greater confidence
    totalCounts1s <- c(200, 1000) 
    totalCounts2s <- c(900, 600) 
    meanCounts1s <- round(test_generatingProbs1s*totalCounts1s)
    meanCounts2s <- round(test_generatingProbs2s*totalCounts2s)
    MLL <- MBASED:::maxLogLikelihoodCalculator2s(
        lociHapACountsSample1=meanCounts1s, 
        lociTotalCountsSample1=totalCounts1s, 
        lociHapANoASEProbsSample1=test_noASEProbs1s, 
        lociRhosSample1=test_rhos1s,
        lociHapACountsSample2=meanCounts2s, 
        lociTotalCountsSample2=totalCounts2s, 
        lociHapANoASEProbsSample2=test_noASEProbs2s, 
        lociRhosSample2=test_rhos2s
    )
    neighborhoodRange <- 0.01
    probsToTest <- test_trueHapA + 
        c(
            -neighborhoodRange, 
            0, 
            neighborhoodRange
        )
    distMLEEstimateToProbs <- abs(probsToTest-MLL$maximum)
    checkTrue(
        which.min(distMLEEstimateToProbs)==2,
        msg='maxLogLikelihoodCalculator2s: fail manual test'
    )
    MBASED:::testNumericDiff(
        MLL$maximum, 
        test_trueHapA, 
        cutoff=0.01 ## within 1% of true hapA frequency
    )
    
    return(TRUE)
}





## test estimateMAF2s() function
test_estimateMAF2s <- function() {
    ## generate data: both phased and non-phased analysis 
    ## should be performed on my toy example.
    test_trueHapA <- 0.8
    test_trueHapB <- 1-test_trueHapA
    test_noASEProbsHapA1s <- c(0.25, 0.5)
    test_noASEProbsHapA2s <- c(0.3, 0.7)
    test_noASEProbsHapB1s <- 1-test_noASEProbsHapA1s
    test_noASEProbsHapB2s <- 1-test_noASEProbsHapA2s
    test_rhos1s <- c(0.01, 0.1)
    test_rhos2s <- c(0.005, 0.05)
    test_generatingProbsHapA1s <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbsHapA1s
    ) 
    test_generatingProbsHapA2s <- MBASED:::getPFinal(
        trueAF=test_trueHapA,
        noASEAF=test_noASEProbsHapA2s
    ) 
    ## high counts lead to greater confidence
    totalCounts1s <- c(200, 1000) 
    totalCounts2s <- c(900, 600) 
    meanCountsHapA1s <- round(
        test_generatingProbsHapA1s*totalCounts1s
    )
    meanCountsHapA2s <- round(
        test_generatingProbsHapA2s*totalCounts2s
    )
    ## define phased and non-phased counts
    hapACountsPhased1s <- meanCountsHapA1s
    hapACountsPhased2s <- meanCountsHapA2s
    hapBCountsPhased1s <- totalCounts1s-meanCountsHapA1s
    hapBCountsPhased2s <- totalCounts2s-meanCountsHapA2s
    allele1CountsUnphased1s <- c(
        hapACountsPhased1s[1], 
        hapBCountsPhased1s[2]
    )
    allele1CountsUnphased2s <- c(
        hapACountsPhased2s[1], 
        hapBCountsPhased2s[2]
    )
    allele2CountsUnphased1s <- totalCounts1s-allele1CountsUnphased1s
    allele2CountsUnphased2s <- totalCounts2s-allele1CountsUnphased2s
    test_noASEProbsAllele11s <- c(
        test_noASEProbsHapA1s[1], 
        test_noASEProbsHapB1s[2]
    )
    test_noASEProbsAllele12s <- c(
        test_noASEProbsHapA2s[1], 
        test_noASEProbsHapB2s[2]
    )
    test_noASEProbsAllele21s <- 1-test_noASEProbsAllele11s
    test_noASEProbsAllele22s <- 1-test_noASEProbsAllele12s
    outputEstimates <- vector('list')
    outputEstimates[['hapA']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF2s(
                lociAllele1CountsSample1=hapACountsPhased1s, 
                lociTotalCountsSample1=totalCounts1s, 
                lociAllele1NoASEProbsSample1=test_noASEProbsHapA1s, 
                lociRhosSample1=test_rhos1s,
                lociAllele1CountsSample2=hapACountsPhased2s, 
                lociTotalCountsSample2=totalCounts2s, 
                lociAllele1NoASEProbsSample2=test_noASEProbsHapA2s, 
                lociRhosSample2=test_rhos2s,
                isPhased=test_isPhased
            )
        }
    )
    outputEstimates[['hapB']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF2s(
                lociAllele1CountsSample1=hapBCountsPhased1s, 
                lociTotalCountsSample1=totalCounts1s, 
                lociAllele1NoASEProbsSample1=test_noASEProbsHapB1s, 
                lociRhosSample1=test_rhos1s,
                lociAllele1CountsSample2=hapBCountsPhased2s, 
                lociTotalCountsSample2=totalCounts2s, 
                lociAllele1NoASEProbsSample2=test_noASEProbsHapB2s, 
                lociRhosSample2=test_rhos2s,
                isPhased=test_isPhased
            )
        }
    )
    outputEstimates[['allele1']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF2s(
                lociAllele1CountsSample1=allele1CountsUnphased1s, 
                lociTotalCountsSample1=totalCounts1s, 
                lociAllele1NoASEProbsSample1=test_noASEProbsAllele11s, 
                lociRhosSample1=test_rhos1s,
                lociAllele1CountsSample2=allele1CountsUnphased2s, 
                lociTotalCountsSample2=totalCounts2s, 
                lociAllele1NoASEProbsSample2=test_noASEProbsAllele12s, 
                lociRhosSample2=test_rhos2s,
                isPhased=test_isPhased
            )
        }
    )
    outputEstimates[['allele2']] <- lapply(
        c(TRUE, FALSE), 
        function(test_isPhased) {
            MBASED:::estimateMAF2s(
                lociAllele1CountsSample1=allele2CountsUnphased1s, 
                lociTotalCountsSample1=totalCounts1s, 
                lociAllele1NoASEProbsSample1=test_noASEProbsAllele21s, 
                lociRhosSample1=test_rhos1s,
                lociAllele1CountsSample2=allele2CountsUnphased2s, 
                lociTotalCountsSample2=totalCounts2s, 
                lociAllele1NoASEProbsSample2=test_noASEProbsAllele22s, 
                lociRhosSample2=test_rhos2s,
                isPhased=test_isPhased
            )
        }
    )
    for (outputInd in 1:length(outputEstimates)) {
        names(outputEstimates[[outputInd]]) <- c('phased', 'nonPhased')
    }

    ## check1: when data is unphased, all MAF estimates 
    ## should be identical, and major haplotype should be haplotypeA
    for (outputInd in 2:length(outputEstimates)) {
        checkEqualsNumeric(
            outputEstimates[[1]]$nonPhased$MAF,
            outputEstimates[[outputInd]]$nonPhased$MAF,
            msg='test_estimateMAF1s: fail unphased MAF test'
        )
    }
    checkIdentical(
        outputEstimates$hapA$nonPhased$allele1IsMajor,
        c(TRUE, TRUE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )
    checkIdentical(
        outputEstimates$hapB$nonPhased$allele1IsMajor,
        c(FALSE, FALSE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )
    checkIdentical(
        outputEstimates$allele1$nonPhased$allele1IsMajor,
        c(TRUE, FALSE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )
    checkIdentical(
        outputEstimates$allele2$nonPhased$allele1IsMajor,
        c(FALSE, TRUE),
        msg='test_estimateMAF1s: fail unphased MAF test'
    )

    ## check2: when data is phased, hapA and hapB results are 
    ## the same as when it's unphased
    for (outputName in c('hapA', 'hapB')) {
        checkEquals(
            outputEstimates[[outputName]]$nonPhased,
            outputEstimates[[outputName]]$phased,
            msg='test_estimateMAF1s: fail unphased equals to phased for known haplotypes test'
        )
    }
    
    ## check3: when data is phased, allele1 and allele2 results are 
    ## different from before
    for (outputName in c('allele1', 'allele2')) {
        checkTrue(
            (outputEstimates[[outputName]]$nonPhased$MAF) != 
            (outputEstimates[[outputName]]$phased$MAF),
            msg='test_estimateMAF1s: fail unphased different from phased for unknown haplotypes test'
        )
    }
    
    ## check4: when data is phased, allele1 and allele2 results are 
    ## opposites of each other in terms of haplotype assignment 
    ## but show same MAF
    checkEqualsNumeric(
        outputEstimates$allele1$nonPhased$MAF,
        outputEstimates$allele2$nonPhased$MAF,
        msg='test_estimateMAF1s: fail unphased test'
    )
    checkIdentical(
        outputEstimates$allele1$nonPhased$allele1IsMajor,
        !(outputEstimates$allele2$nonPhased$allele1IsMajor),
        msg='test_estimateMAF1s: fail unphased test'
    )
    
    return(TRUE)
}    



