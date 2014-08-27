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

test_MBASEDVectorizedPropDiffTest <- function() { 
    ## check1: vectorization works for both a single locus gene 
    ## and multi-loci genes
    set.seed(176907)
    numCases <- 100
    for (numLoci in c(1, 4)) {
        test_ns_mats <- lapply(1:2, function(sampleInd) {
            matrix(
                sample(10:20, numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
        })
        test_xs_mats <- lapply(1:2, function(sampleInd) {
            matrix(
                sample(0:10, numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
        })
        test_mus_mats <- lapply(1:2, function(sampleInd) {
            matrix(
                sample(seq(0.1, 0.9, by=0.1), numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
        })
        test_rhos_mats <- lapply(1:2, function(sampleInd) {
            matrix(
                sample(c(0,0.01, 0.3, 0.9), numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
        })
        for (test_AH in c('greater', 'less', 'two.sided')) {
            individualResults <- bplapply(1:numCases, function(colInd) { 
                MBASED:::MBASEDVectorizedPropDiffTest(
                    countsMatSample1=matrix(
                        test_xs_mats[[1]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ), 
                    countsMatSample2=matrix(
                        test_xs_mats[[2]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ),
                    totalsMatSample1=matrix(
                        test_ns_mats[[1]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ), 
                    totalsMatSample2=matrix(
                        test_ns_mats[[2]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ), 
                    probsMatSample1=matrix(
                        test_mus_mats[[1]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ), 
                    probsMatSample2=matrix(
                        test_mus_mats[[2]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ), 
                    rhosMatSample1=matrix(
                        test_rhos_mats[[1]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ), 
                    rhosMatSample2=matrix(
                        test_rhos_mats[[2]][,colInd], 
                        nrow=numLoci, 
                        ncol=1
                    ), 
                    alternative=test_AH
                )
            }, BPPARAM=myBPPARAM)
            jointResults <- MBASED:::MBASEDVectorizedPropDiffTest(
                countsMatSample1=test_xs_mats[[1]], 
                countsMatSample2=test_xs_mats[[2]], 
                totalsMatSample1=test_ns_mats[[1]], 
                totalsMatSample2=test_ns_mats[[2]], 
                probsMatSample1=test_mus_mats[[1]], 
                probsMatSample2=test_mus_mats[[2]], 
                rhosMatSample1=test_rhos_mats[[1]], 
                rhosMatSample2=test_rhos_mats[[2]],
                alternative=test_AH
            )
            for (outputName in names(jointResults)) {
                checkEqualsNumeric(
                    jointResults[[outputName]], 
                    do.call(
                        cbind,
                        lapply(
                            individualResults, 
                            function(el) {
                                el[[outputName]]
                            }
                        )
                    ), 
                    msg='test_MBASEDVectorizedPropDiffTest: fail vectorization test'
                )
            }
        }
    }
            
    ## check2: for large n, reasonable p and rho=0 (binomial setting), 
    ## single-locus approach should produce pvalues similar to 
    ## proportion difference test (apart from possibly very small p-values). 
    ## Further, the estimate should be (possibly adjusted) proportion 
    ## difference (we previously checked that the adjustment of 
    ## proportions is meaningful). Here I need to keep test_mu the same 
    ## for two samples since prop.test assumes the underlying 
    ## probabilities of success are the same for both samples.
    set.seed(27567)
    testCasesMat <- expand.grid(
        test_mu=seq(0.2, 0.8, by=0.1),
        test_rho=0,
        test_AH=c('greater', 'less', 'two.sided')
    )
    if (FULLTESTING) {
        numCountsToGenerate <- 10^6
    } else {
        numCountsToGenerate <- 10^3  
    }
    testResults <- bplapply(1:nrow(testCasesMat), function(testInd) {
        testCase <- testCasesMat[testInd,]
        testNs <- lapply(1:2, function(sampleInd) {
            sample((10^4):(10^5), 1, replace=TRUE)
        })
        testCounts <- lapply(1:2, function(sampleInd) {
            MBASED:::vectorizedRbetabinomMR(
                n=numCountsToGenerate,
                size=rep(testNs[[sampleInd]], numCountsToGenerate),
                mu=rep(testCase$test_mu, numCountsToGenerate),
                rho=rep(testCase$test_rho, numCountsToGenerate)
            )
        })
        predEstimatesList <- lapply(1:2, function(sampleInd) {
            MBASED:::shiftAndAttenuateProportions(
                countsMat=matrix(
                    testCounts[[sampleInd]], 
                    nrow=1
                ), 
                totalsMat=matrix(
                    rep(testNs[[sampleInd]], numCountsToGenerate),
                    nrow=1
                ),
                probsMat=matrix(
                    rep(testCase$test_mu, numCountsToGenerate),
                    nrow=1
                ),
                rhosMat=matrix(
                    rep(testCase$test_rho, numCountsToGenerate),
                    nrow=1
                )
            )
        })
        predEstimates <- as.vector(
            predEstimatesList[[1]]$propsShifted- 
            predEstimatesList[[2]]$propsShifted
        )
        predPvals <- sapply(1:numCountsToGenerate, function(countInd) {
            prop.test(
                x=c(
                    testCounts[[1]][countInd], 
                    testCounts[[2]][countInd]
                ),
                n=c(
                    testNs[[1]], 
                    testNs[[2]]
                ),
                alternative=as.character(testCase$test_AH)
            )$p.value
        })
        MBASEDVectorizedPropDiffTestOutput <- 
            MBASED:::MBASEDVectorizedPropDiffTest(
                countsMatSample1=matrix(
                    testCounts[[1]], 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                countsMatSample2=matrix(
                    testCounts[[2]], 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                totalsMatSample1=matrix(
                    rep(testNs[[1]], numCountsToGenerate), 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                totalsMatSample2=matrix(
                    rep(testNs[[2]], numCountsToGenerate), 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                probsMatSample1=matrix(
                    rep(testCase$test_mu, numCountsToGenerate), 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                probsMatSample2=matrix(
                    rep(testCase$test_mu, numCountsToGenerate), 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                rhosMatSample1=matrix(
                    rep(testCase$test_rho, numCountsToGenerate), 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                rhosMatSample2=matrix(
                    rep(testCase$test_rho, numCountsToGenerate), 
                    nrow=1, 
                    ncol=numCountsToGenerate
                ), 
                alternative=as.character(testCase$test_AH)
            )
        obsEstimates <- as.vector(
            MBASEDVectorizedPropDiffTestOutput$propDifferenceFinal
        )
        checkEqualsNumeric(
            obsEstimates, 
            predEstimates, 
            msg='test_MBASEDVectorizedPropDiffTestOutput: fail comparison to prop.test'
        )
        obsPvals <- as.vector(
            MBASEDVectorizedPropDiffTestOutput$pValue
        )
        checkTrue(
            all(
                (pmax(obsPvals, predPvals)<=0.005) |
                    MBASED:::testNumericDiff(
                    queryVals=obsPvals, 
                    targetVals=predPvals, 
                    cutoffFraction=0.2 ## within 20%
                )
            ),
            msg='test_MBASEDVectorizedPropDiffTestOutput: fail comparison to binom/betabinom test'
        ) 
    }, BPPARAM=myBPPARAM)
    
    
    return(TRUE)
}    
    
    

## GENERAL FRAMEWORK FOR TEST CASES FOR
##  test_runMBASED2s1aseID() and test_runMBASED2s(): 
## Try out n=10,11,12,..,20, 50, 51, 100, 101, 1000, 1001. 
## For each n, test (when possible) x={0,1,2}(small), 
## {n, n-1, n-2}(large), {floor(n/2-1), floor(n/2), floor(n/2+1)} (mean), 
## and 10 random values.
## For each combination of n and x, try mu={0.4, 0.45, 0.5, 0.55, 0.6}
## For each combination of n, x, mu, try rho={0, 0.01, 0.1}
set.seed(554475)
my_ns <- c(10:20, 50, 51, 100, 101, 1000, 1001)
my_xn_combos <- unlist(lapply(my_ns, function(n) {
    candidate_xs <- c(
        0:2, 
        n-(0:2), 
        floor(n/2-1), 
        floor(n/2), 
        floor(n/2+1), 
        sample(0:n, 10, replace=TRUE)
    )
    candidate_xs <- candidate_xs[
        candidate_xs>=0 & candidate_xs<=n
    ]
    candidate_xs <- unique(
        sort(
            candidate_xs, 
            decreasing=FALSE
        )
    )
    return(paste(n, candidate_xs, sep='_'))
}))
AllTestCasesMat <- expand.grid(
    test_xn_combo=my_xn_combos,
    test_mu=seq(0.4, 0.6, by=0.05),
    test_rho=c(0, 0.01, 0.1)
)
AllTestCasesMat$test_n <- as.numeric(
    matrix(
        unlist(
            strsplit(
                as.character(
                    AllTestCasesMat$test_xn_combo
                ), 
                '_'
            )
        ), 
        ncol=2, 
        byrow=TRUE
    )[,1]
)
AllTestCasesMat$test_x <- as.numeric(
    matrix(
        unlist(
            strsplit(
                as.character(
                    AllTestCasesMat$test_xn_combo
                ), 
                '_'
            )
        ), 
        ncol=2, 
        byrow=TRUE
    )[,2]
)
rm(
    my_ns,
    my_xn_combos
)

## helper function to get test cases from the matrix of all test cases
## returns a matrix with test cases of interest 
## spread over a desired matrix
getTestCases <- function(numCases, numLoci) {
    testCasesInds <- sample(
        1:nrow(AllTestCasesMat), 
        numCases*numLoci, 
        replace=TRUE
    )
    test_ns_mat <- matrix(
        AllTestCasesMat$test_n[testCasesInds], 
        nrow=numLoci, 
        ncol=numCases
    )
    test_xs_mat <- matrix(
        AllTestCasesMat$test_x[testCasesInds], 
        nrow=numLoci, 
        ncol=numCases
    )
    test_mus_mat <- matrix(
        AllTestCasesMat$test_mu[testCasesInds], 
        nrow=numLoci, 
        ncol=numCases
    )
    test_rhos_mat <- matrix(
        AllTestCasesMat$test_rho[testCasesInds], 
        nrow=numLoci, 
        ncol=numCases
    )
    return(
        list(
            n=test_ns_mat,
            x=test_xs_mat,
            mu=test_mus_mat,
            rho=test_rhos_mat
        )
    )
}

## helper function to get individual test case (column) from 
## output of getTestCases()
getTestCase <- function(testCasesMats, testInd) {
    return(
        list(
            n=testCasesMats$n[,testInd],
            x=testCasesMats$x[,testInd],
            mu=testCasesMats$mu[,testInd],
            rho=testCasesMats$rho[,testInd]
        )
    )
}



## test behavior of runMBASED2s1aseID() function 
test_runMBASED2s1aseID <- function() {
    ## check 1:  with no simulations return the same result as 
    ## MBASEDVectorizedPropDiffTest with alternative hypothesis 
    ## set to 'two.sided' and with additional step of 
    ## picking the larger of allelic frequencies in sample1 to be 'major'
    ## use tieBreakRandom=FALSE for simplicity. 
    ## Later will check that tie-breaking works as intended.
    if (FULLTESTING) {
    	## number of test cases
        numCases <- 1000 
    } else {
    	## number of test cases
        numCases <- 3 
    }
    ## number of loci
    for (numLoci in c(1,5)) { 
        testCasesMatsSample1 <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        testCasesMatsSample2 <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        testResults <- bplapply(1:numCases, function(testInd) {
            testCaseSample1 <- getTestCase(
                testCasesMats=testCasesMatsSample1, 
                testInd=testInd
            )
            testCaseSample2 <- getTestCase(
                testCasesMats=testCasesMatsSample2, 
                testInd=testInd
            )
            for (test_isPhased in c(TRUE, FALSE)) {
            	## data will be phased by MBASED based on sample1 only
                test_x_is_major <- MBASED:::runMBASED1s1aseID(
                    lociAllele1Counts=testCaseSample1$x, 
                    lociAllele2Counts=testCaseSample1$n-
                        testCaseSample1$x, 
                    lociAllele1NoASEProbs=testCaseSample1$mu, 
                    lociRhos=testCaseSample1$rho, 
                     ## no need to do simulations: 
                     ## we only require information about which allele is major
                    numSim=0,
                    isPhased=test_isPhased,
                    tieBreakRandom=FALSE
                )$lociAllele1IsMajor
                test_maj_sample1 <- ifelse(
                    test_x_is_major, 
                    testCaseSample1$x, 
                    testCaseSample1$n-testCaseSample1$x
                )
                test_maj_sample2 <- ifelse(
                    test_x_is_major, 
                    testCaseSample2$x, 
                    testCaseSample2$n-testCaseSample2$x
                )
                test_prob_maj_sample1 <- ifelse(
                    test_x_is_major, 
                    testCaseSample1$mu, 
                    1-testCaseSample1$mu
                )
                test_prob_maj_sample2 <- ifelse(
                    test_x_is_major, 
                    testCaseSample2$mu, 
                    1-testCaseSample2$mu
                )
                MBASEDVectorizedPropDiffTestOutput <- 
                    MBASED:::MBASEDVectorizedPropDiffTest(
                        countsMatSample1=matrix(
                            test_maj_sample1, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        countsMatSample2=matrix(
                            test_maj_sample2, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        totalsMatSample1=matrix(
                            testCaseSample1$n, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        totalsMatSample2=matrix(
                            testCaseSample2$n, 
                            nrow=numLoci, 
                            ncol=1
                        ),
                        probsMatSample1=matrix(
                            test_prob_maj_sample1,
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        probsMatSample2=matrix(
                            test_prob_maj_sample2,
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        rhosMatSample1=matrix(
                            testCaseSample1$rho, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        rhosMatSample2=matrix(
                            testCaseSample2$rho, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        alternative='two.sided'
                    )    
                runMBASED2s1aseIDOutput <- 
                    MBASED:::runMBASED2s1aseID(
                        lociAllele1CountsSample1=testCaseSample1$x, 
                        lociAllele2CountsSample1=testCaseSample1$n-
                            testCaseSample1$x, 
                        lociAllele1CountsSample2=testCaseSample2$x, 
                        lociAllele2CountsSample2=testCaseSample2$n-
                            testCaseSample2$x, 
                        lociAllele1NoASEProbsSample1=testCaseSample1$mu, 
                        lociAllele1NoASEProbsSample2=testCaseSample2$mu, 
                        lociRhosSample1=testCaseSample1$rho, 
                        lociRhosSample2=testCaseSample2$rho,
                        numSim=0, 
                        isPhased=test_isPhased,
                        tieBreakRandom=FALSE
                    )
                checkEqualsNumeric(
                    as.vector(
                        MBASEDVectorizedPropDiffTestOutput$propDifferenceFinal
                    ), 
                    runMBASED2s1aseIDOutput$majorAlleleFrequencyDifference, 
                    msg=paste('test_runMBASED2s1aseID: failed comaprison with MBASEDVectorizedPropDiffTest test')
                )
                checkEqualsNumeric(
                    as.vector(
                        MBASEDVectorizedPropDiffTestOutput$pValue
                    ), 
                    runMBASED2s1aseIDOutput$pValueASE, 
                    msg=paste('test_runMBASED2s1aseID: failed comaprison with MBASEDVectorizedPropDiffTest test')
                )
                checkEqualsNumeric(
                    as.vector(
                        MBASEDVectorizedPropDiffTestOutput$hetQ
                    ), 
                    runMBASED2s1aseIDOutput$heterogeneityQ, 
                    msg=paste('test_runMBASED2s1aseID: failed comaprison with MBASEDVectorizedPropDiffTest test')
                )
                checkTrue(
                    ifelse(
                        numLoci==1,
                        is.na(runMBASED2s1aseIDOutput$heterogeneityQ), 
                        !is.na(runMBASED2s1aseIDOutput$heterogeneityQ)
                    ),
                    msg=paste('test_runMBASED2s1aseID: failed comaprison with MBASEDVectorizedPropDiffTest  test')
                )
                checkEqualsNumeric(
                    as.vector(
                        MBASEDVectorizedPropDiffTestOutput$hetPVal
                    ), 
                    runMBASED2s1aseIDOutput$pValueHeterogeneity, 
                    msg=paste('test_runMBASED2s1aseID: failed comaprison with MBASEDVectorizedPropDiffTest test')
                )
                checkTrue(
                    ifelse(
                        numLoci==1,
                        is.na(runMBASED2s1aseIDOutput$pValueHeterogeneity), 
                        !is.na(runMBASED2s1aseIDOutput$pValueHeterogeneity)
                    ),
                    msg=paste('test_runMBASED2s1aseID: failed comaprison with MBASEDVectorizedPropDiffTest test')
                )
                checkIdentical(
                    runMBASED2s1aseIDOutput$lociAllele1IsMajor,
                    test_x_is_major, 
                     msg=paste('test_runMBASED2s1aseID: failed comaprison with MBASEDVectorizedPropDiffTest test')
                )
            }
        }, BPPARAM=myBPPARAM)
    }
    rm(
        numCases, 
        numLoci, 
        testCasesMatsSample1,
        testCasesMatsSample2,
        testResults
    )
    
    ## check2: breaking ties works as intended .
    if (FULLTESTING) {
    	 ## number of test cases
        numCases <- 1000
    } else {
    	## number of test cases
        numCases <- 100 
    }
    tiedVal <- 10
    for (numLoci in c(1, 5)) {
        for (test_numSim in c(0, 3)) {
            noTieBreaking <- unlist(bplapply(1:numCases, function(testInd) {
                MBASED:::runMBASED2s1aseID(
                    lociAllele1CountsSample1=rep(tiedVal, numLoci), 
                    lociAllele2CountsSample1=rep(tiedVal, numLoci), 
                    lociAllele1CountsSample2=rep(tiedVal, numLoci), 
                    lociAllele2CountsSample2=rep(tiedVal, numLoci),
                    lociAllele1NoASEProbsSample1=rep(0.5, numLoci), 
                    lociAllele1NoASEProbsSample2=rep(0.5, numLoci),
                    lociRhosSample1=rep(0, numLoci), 
                    lociRhosSample2=rep(0, numLoci), 
                    numSim=test_numSim, 
                    isPhased=sample(c(TRUE,FALSE),1), 
                    tieBreakRandom=FALSE
                )$lociAllele1IsMajor
            }, BPPARAM=myBPPARAM))
            checkTrue(
                length(unique(noTieBreaking))==1, 
                msg='test_runMBASED2s1aseID: failed tieBreak test: no simulations, no tie breaks'
            )
            yesTieBreakingNoNeed <- unlist(bplapply(
                1:numCases, 
                function(testInd) {
                    MBASED:::runMBASED2s1aseID(
                        lociAllele1CountsSample1=rep(tiedVal, numLoci), 
                        lociAllele2CountsSample1=rep(2*tiedVal, numLoci), 
                        lociAllele1CountsSample2=rep(tiedVal, numLoci), 
                        lociAllele2CountsSample2=rep(tiedVal, numLoci),
                        lociAllele1NoASEProbsSample1=rep(0.5, numLoci), 
                        lociAllele1NoASEProbsSample2=rep(0.5, numLoci),
                        lociRhosSample1=rep(0, numLoci), 
                        lociRhosSample2=rep(0, numLoci), 
                        numSim=test_numSim, 
                        isPhased=sample(c(TRUE,FALSE),1), 
                        tieBreakRandom=TRUE
                    )$lociAllele1IsMajor
                }, BPPARAM=myBPPARAM
            ))
            checkTrue(
                length(unique(yesTieBreakingNoNeed))==1, 
                msg='test_runMBASED2s1aseID: failed tieBreak test: no simulations, no tie breaks'
            )
            yesTieBreaking <- unlist(bplapply(1:numCases, function(testInd) {
                MBASED:::runMBASED2s1aseID(
                    lociAllele1CountsSample1=rep(tiedVal, numLoci), 
                    lociAllele2CountsSample1=rep(tiedVal, numLoci), 
                    lociAllele1CountsSample2=rep(tiedVal, numLoci), 
                    lociAllele2CountsSample2=rep(tiedVal, numLoci),
                    lociAllele1NoASEProbsSample1=rep(0.5, numLoci), 
                    lociAllele1NoASEProbsSample2=rep(0.5, numLoci),
                    lociRhosSample1=rep(0, numLoci), 
                    lociRhosSample2=rep(0, numLoci), 
                    numSim=test_numSim, 
                    isPhased=sample(c(TRUE,FALSE),1), 
                    tieBreakRandom=TRUE
                )$lociAllele1IsMajor
            }, BPPARAM=myBPPARAM))
            MBASED:::testQuantiles(
                theoreticalCumDist=0.5, 
                observedCumDist=mean(yesTieBreaking), 
                numTotalCounts=numCases, 
                numSEsToCheck=5, 
                errorMessage='test_runMBASED2s1aseID: failed tieBreak test: no simulations, yes tie breaks'
            )
        }
    }
    rm(
        numLoci,
        numCases,
        test_numSim,
        tiedVal,
        noTieBreaking,
        yesTieBreakingNoNeed,
        yesTieBreaking
    )
    
    ## check3: with simulations, MAF difference estimate and 
    ## assignment of alleles to haplotypes is the same 
    ## as without simulations
    if (FULLTESTING) {
    	## number of test cases
        numCases <- 1000 
    } else {
    	## number of test cases
        numCases  <- 3 
    }
     ## number of loci
    for (numLoci in c(1, 5)) {
        testCasesMatsSample1 <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        testCasesMatsSample2 <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        ## check for all possible values of tieBreakRandom
        for (test_tieBreakRandom in c(TRUE,FALSE)) { 
            for (test_isPhased in c(TRUE, FALSE)) {
                testIndSeeds <- sample(
                    1:(10^6), 
                    numCases, 
                    replace=TRUE
                )
                testResults <- bplapply(1:numCases, function(testInd) {
                    testCaseSample1 <- getTestCase(
                        testCasesMats=testCasesMatsSample1, 
                        testInd=testInd
                    )
                    testCaseSample2 <- getTestCase(
                        testCasesMats=testCasesMatsSample2, 
                        testInd=testInd
                    )
                    test_seed <- testIndSeeds[testInd]
                    testOutputs <- lapply(c(0, 10), function(test_numSim) {
                    	## make sure ties get broken in same way
                        set.seed(test_seed) 
                        MBASED:::runMBASED2s1aseID(
                            lociAllele1CountsSample1=testCaseSample1$x, 
                            lociAllele2CountsSample1=testCaseSample1$n-
                                testCaseSample1$x, 
                            lociAllele1CountsSample2=testCaseSample2$x, 
                            lociAllele2CountsSample2=testCaseSample2$n-
                                testCaseSample2$x, 
                            lociAllele1NoASEProbsSample1=testCaseSample1$mu, 
                            lociAllele1NoASEProbsSample2=testCaseSample2$mu, 
                            lociRhosSample1=testCaseSample1$rho, 
                            lociRhosSample2=testCaseSample2$rho,
                            numSim=test_numSim, 
                            isPhased=test_isPhased, 
                            tieBreakRandom=test_tieBreakRandom
                        )
                    })
                    checkEqualsNumeric(
                        testOutputs[[1]]$majorAlleleFrequencyDifference, 
                        testOutputs[[2]]$majorAlleleFrequencyDifference, 
                        msg='test_runMBASED2s1aseID: failed test checking that MAF difference estimate is the same for sim and no sim'
                    )
                    checkEquals(
                        testOutputs[[1]]$lociAllele1IsMajor, 
                        testOutputs[[2]]$lociAllele1IsMajor, 
                        msg='test_runMBASED2s1aseID: failed test checking that assignment of alleles to haplotypes is the same for sim and no sim'
                    )
                }, BPPARAM=myBPPARAM)
            }
        }
    }
    rm(
        numCases, 
        numLoci, 
        testCasesMatsSample1,
        testCasesMatsSample2,
        test_tieBreakRandom,
        test_isPhased,
        testIndSeeds,
        testResults
    )
    
    ## check 4: under null hypothesis and for reasonable rho, 
    ## the ASE and heterogeneity p-values are approximately uniform. 
    ## The actual distribution will NOT be uniform because 
    ## of the step where MAF gets estimated, resulting in p-values 
    ## being generated based on incorrect (but, hopefully, 
    ## good approximation to ) null distribution. 
    
     ## currently this test is not run: there are cases where test fails 
     ## due to estimate of joint MAF being sufficiently distinct from 
     ## true MAF (here: 0.5), where 'sufficiently' means different things 
     ## for different coverage levels.
    if (FALSE) {
        set.seed(29921)
        ## number of test cases
        numCases <- 5 
        ## number of simulations to perform by MBASED
        test_numSim <- 10^4 
        ## number of null p-values to be calculated to assess if 
        ## distribution is close to uniform or not 
        ## (~ 20 min for each set of settings) 
        numTests <- 10^4 
        ## number of loci
        for (numLoci in c(1,3)) { 
            testCasesMatsSample1 <- getTestCases(
                numCases=numCases, 
                numLoci=numLoci
            )
            testCasesMatsSample2 <- getTestCases(
                numCases=numCases, 
                numLoci=numLoci
            )
            quantilesToCheck <- seq(0.05, 0.95, by=0.05)
            ## check for all possible values of test_isPhased
            for (test_isPhased in c(TRUE,FALSE)) { 
            	## check for all possible values of tieBreakRandom
                for (test_tieBreakRandom in c(TRUE,FALSE)) { 
                	## no parallelizing here, because 
                	## parallelizing is done over numTests
                    for (testInd in 1:numCases) { 
                        testCaseSample1 <- getTestCase(
                            testCasesMats=testCasesMatsSample1, 
                            testInd=testInd
                        )
                        testCaseSample2 <- getTestCase(
                            testCasesMats=testCasesMatsSample2, 
                            testInd=testInd
                        )
                        ## make sure counts are high enough
                        testCaseSample1$n <- sample(200:1000, 1)
                        testCaseSample2$n <- sample(200:1000, 1)
                        test_nullResults <- bplapply(1:numTests, function(runInd) {
                            test_x_null_sample1 <- 
                                MBASED:::vectorizedRbetabinomMR(
                                    n=length(testCaseSample1$n), 
                                    size=testCaseSample1$n, 
                                    mu=testCaseSample1$mu, 
                                    rho=testCaseSample1$rho
                                )
                            test_x_null_sample2 <- 
                                MBASED:::vectorizedRbetabinomMR(
                                    n=length(testCaseSample2$n), 
                                    size=testCaseSample2$n, 
                                    mu=testCaseSample2$mu, 
                                    rho=testCaseSample2$rho
                                )
                            MBASED:::runMBASED2s1aseID(
                                lociAllele1CountsSample1=test_x_null_sample1, 
                                lociAllele2CountsSample1=testCaseSample1$n-
                                    test_x_null_sample1, 
                                lociAllele1CountsSample2=test_x_null_sample2, 
                                lociAllele2CountsSample2=testCaseSample2$n-
                                    test_x_null_sample2, 
                                lociAllele1NoASEProbsSample1=testCaseSample1$mu, 
                                lociAllele1NoASEProbsSample2=testCaseSample2$mu, 
                                lociRhosSample1=testCaseSample1$rho, 
                                lociRhosSample2=testCaseSample2$rho, 
                                numSim=test_numSim, 
                                isPhased=test_isPhased,
                                tieBreakRandom=test_tieBreakRandom
                            )
                        }, BPPARAM=myBPPARAM)
                        test_null_ASEPvals <- sapply(test_nullResults, function(el) {
                            el$pValueASE
                        })
                        obs_ASEPvals <- sapply(
                            quantilesToCheck, 
                            function(quantileToCheck) {
                                mean(test_null_ASEPvals<=quantileToCheck)
                            }
                        )
                        checkTrue(
                            all(
                                abs(quantilesToCheck-obs_ASEPvals)<0.05 | 
                                MBASED:::testNumericDiff(
                                    queryVals=quantilesToCheck, 
                                    targetVals=obs_ASEPvals, 
                                    ## within 20% of each other (e.g, 0.12 vs. 0.1 p-value)
                                    cutoffFraction=0.2 
                                )
                            ),
                            msg='test_runMBASED2s1aseID: fail test for relative uniformity of ASE p-values under null'
                        ) 
                        ## het p-values only defined for multi-locus genes
                        if (numLoci >1) { 
                            test_null_hetPvals <- sapply(
                                test_nullResults, 
                                function(el) {
                                    el$pValueHeterogeneity
                                }
                            )
                            obs_hetPvals <- sapply(
                                quantilesToCheck, 
                                function(quantileToCheck) {
                                    mean(test_null_hetPvals<=quantileToCheck)
                                }
                            )
                            checkTrue(
                                all(
                                    abs(quantilesToCheck-obs_hetPvals)<0.03 |
                                    MBASED:::testNumericDiff(
                                        queryVals=quantilesToCheck, 
                                        targetVals=obs_hetPvals, 
                                        ## within 20% of each other (e.g, 0.12 vs. 0.1 p-value)
                                        cutoffFraction=0.2 
                                    )    
                                ),
                                msg='test_runMBASED2s1aseID: fail test for relative uniformity of heterogeneity p-values under null'
                            )
                        }     
                    }
                }
            }
        }
        rm(
            numCases, 
            numLoci, 
            test_numSim,
            numTests,
            testCasesMats,
            quantilesToCheck,
            numSEsToCheck,
            test_tieBreakRandom,
            test_isPhased,
            testInd,
            testCase,
            test_nullResults
        )
    }
        
    ## Check 5: check that we can detect true signal: 
    ## p-value should be small. 
    ## only do for full testing, otherwise takes too long: ~ 20-25 min
    if (FULLTESTING) { 
        set.seed(433938)
        ## number of test cases
        numCases <- 100 
        ## number of loci
        for (numLoci in c(1, 5))  {
        	## number of simulations to perform by MBASED
            test_numSim <- 10^6 
            trueMAF=0.9
            testCasesMatsSample1 <- getTestCases(
                numCases=numCases, 
                numLoci=numLoci
            )
            testCasesMatsSample2 <- getTestCases(
                numCases=numCases, 
                numLoci=numLoci
            )
            ## check for all possible values of isPhased
            for (test_isPhased in c(TRUE,FALSE)) { 
            	## check for all possible values of tieBreakRandom
                for (test_tieBreakRandom in c(TRUE,FALSE)) { 
                    testResults <- bplapply(1:numCases, function(testInd) {
                    #for (testInd in 1:numCases) {
                        testCaseSample1 <- getTestCase(
                            testCasesMats=testCasesMatsSample1, 
                            testInd=testInd
                        )
                        testCaseSample2 <- getTestCase(
                            testCasesMats=testCasesMatsSample2, 
                            testInd=testInd
                        )
                        ## make sure we have power to detect ASE
                        testCaseSample1$n <- 1000*testCaseSample1$n 
                        ## make sure we have power to detect ASE
                        testCaseSample2$n <- 1000*testCaseSample2$n 
                        ## trueAF must be a single number
                        test_generating_prob_sample1 <- MBASED:::getPFinal( 
                            trueAF=trueMAF, 
                            noASEAF=testCaseSample1$mu
                        )
                        ## trueAF must be a single number
                        test_generating_prob_sample2 <- MBASED:::getPFinal( 
                            trueAF=0.5, 
                            noASEAF=testCaseSample2$mu
                        )
                        test_x_sample1 <- round(
                            test_generating_prob_sample1*testCaseSample1$n
                        )
                        test_x_sample2 <- round(
                            test_generating_prob_sample2*testCaseSample2$n
                        )
                        runMBASED2s1aseIDOutput <- 
                            MBASED:::runMBASED2s1aseID(
                                lociAllele1CountsSample1=test_x_sample1, 
                                lociAllele2CountsSample1=testCaseSample1$n-
                                    test_x_sample1, 
                                lociAllele1CountsSample2=test_x_sample2, 
                                lociAllele2CountsSample2=testCaseSample2$n-
                                    test_x_sample2, 
                                lociAllele1NoASEProbsSample1=testCaseSample1$mu, 
                                lociAllele1NoASEProbsSample2=testCaseSample2$mu, 
                                lociRhosSample1=testCaseSample1$rho, 
                                lociRhosSample2=testCaseSample2$rho,
                                numSim=test_numSim, 
                                isPhased=test_isPhased,
                                tieBreakRandom=test_tieBreakRandom
                            )
                        checkTrue(
                            runMBASED2s1aseIDOutput$pValueASE <= 0.1, 
                            msg='test_runMBASED2s1aseID: 2-sample does not detect the signal'
                        )
                    #}    
                    }, BPPARAM=myBPPARAM)
                }
            }
        }
        rm(
            numCases, 
            numLoci, 
            test_numSim,
            trueMAF,
            testCasesMatsSample1,
            testCasesMatsSample2,
            test_isPhased,
            test_tieBreakRandom,
            testResults
        )
    }
    
    ## Check 6: For multi-SNV genes only: check that we can 
    ## detect true SNV-specific ASE: heterogeneity p-value should be 
    ## small (for largish n and smallish rho, otherwise not enough evidence). 
    ## only do for full testing, otherwise takes too long
    if (FULLTESTING) { 
        set.seed(708872)
        ## number of test cases
        numCases <- 100 
        ## number of loci
        numLoci <- 3 
        ## number of simulations to perform by MBASED
        test_numSim <- 10^6 
        trueMAF=0.99
        testCasesMatsSample1 <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        testCasesMatsSample2 <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        ## check for all possible values of isPhased
        for (test_isPhased in c(TRUE,FALSE)) { 
        	## check for all possible values of tieBreakRandom
            for (test_tieBreakRandom in c(TRUE,FALSE)) { 
                testResults <- bplapply(1:numCases, function(testInd) {
                    testCaseSample1 <- getTestCase(
                        testCasesMats=testCasesMatsSample1, 
                        testInd=testInd
                    )
                    testCaseSample2 <- getTestCase(
                        testCasesMats=testCasesMatsSample2, 
                        testInd=testInd
                    )
                    whichSNVisASE <- sample(1:length(testCaseSample1$n), 1)
                    trueMAFsSample1 <- rep(0.5, length(testCaseSample1$n))
                    trueMAFsSample2 <- rep(0.5, length(testCaseSample2$n))
                    trueMAFsSample1[whichSNVisASE] <- trueMAF
                    ## make sure we have power to detect ASE
                    testCaseSample1$n <- 1000*testCaseSample1$n 
                    ## make sure we have power to detect ASE
                    testCaseSample2$n <- 1000*testCaseSample2$n 
                    ## make sure we have power to detect ASE
                    testCaseSample1$rho <- testCaseSample1$rho/10 
                    ## make sure we have power to detect ASE
                    testCaseSample2$rho <- testCaseSample2$rho/10  
                    test_generating_prob_sample1 <- sapply(
                        1:length(trueMAFsSample1), 
                        function(MAFInd) {
                        	## trueAF must be a single number
                            MBASED:::getPFinal(
                                trueAF=trueMAFsSample1[MAFInd], 
                                noASEAF=testCaseSample1$mu[MAFInd]
                            )
                        }
                    )
                    test_generating_prob_sample2 <- sapply(
                        1:length(trueMAFsSample2), 
                        function(MAFInd) {
                        	## trueAF must be a single number
                            MBASED:::getPFinal(
                                trueAF=trueMAFsSample2[MAFInd], 
                                noASEAF=testCaseSample2$mu[MAFInd]
                            )
                        }
                    )
                    test_x_sample1 <- round(
                        test_generating_prob_sample1*testCaseSample1$n
                    )
                    test_x_sample2 <- round(
                        test_generating_prob_sample2*testCaseSample2$n
                    )
                    runMBASED2s1aseIDOutput <- 
                        MBASED:::runMBASED2s1aseID(
                            lociAllele1CountsSample1=test_x_sample1, 
                            lociAllele2CountsSample1=testCaseSample1$n-
                                test_x_sample1, 
                            lociAllele1CountsSample2=test_x_sample2, 
                            lociAllele2CountsSample2=testCaseSample2$n-
                                test_x_sample2, 
                            lociAllele1NoASEProbsSample1=testCaseSample1$mu, 
                            lociAllele1NoASEProbsSample2=testCaseSample2$mu, 
                            lociRhosSample1=testCaseSample1$rho, 
                            lociRhosSample2=testCaseSample2$rho,
                            numSim=test_numSim, 
                            isPhased=test_isPhased,
                            tieBreakRandom=test_tieBreakRandom
                       )
                    checkTrue(
                        runMBASED2s1aseIDOutput$pValueHeterogeneity <= 0.25, 
                        msg='test_runMBASED2s1aseID: 1-sample does not detect the heterogeneity signal '
                    )
                }, BPPARAM=myBPPARAM)
            }
        }
        rm(
            numCases, 
            numLoci, 
            test_numSim,
            trueMAF,
            testCasesMatsSample1,
            testCasesMatsSample2,
            test_isPhased,
            test_tieBreakRandom,
            testResults
        )
    }
    
    return(TRUE)
}

## test behavior of runMBASED2s1() function 
test_runMBASED2s <- function() {
    ## check1: runMBASED2s gives same results as 
    ## performing analysis on each aseID individually.
    set.seed(353948)
    ## Complication due to random sampling I employ: bplapply 
    ## messes with random seed, so I need to write slightly modified 
    ## version of functions to directly specificy seeds.
    runMBASED2s1aseIDTest <- function (
        lociAllele1CountsSample1, 
        lociAllele2CountsSample1, 
        lociAllele1CountsSample2, 
        lociAllele2CountsSample2, 
        lociAllele1NoASEProbsSample1, 
        lociAllele1NoASEProbsSample2, 
        lociRhosSample1, 
        lociRhosSample2, 
        numSim=0, 
        isPhased=FALSE, 
        tieBreakRandom=FALSE, 
        testSeed
    ) {
        set.seed(testSeed)
        res <- MBASED:::runMBASED2s1aseID(
            lociAllele1CountsSample1=lociAllele1CountsSample1, 
            lociAllele2CountsSample1=lociAllele2CountsSample1, 
            lociAllele1CountsSample2=lociAllele1CountsSample2, 
            lociAllele2CountsSample2=lociAllele2CountsSample2, 
            lociAllele1NoASEProbsSample1=lociAllele1NoASEProbsSample1, 
            lociAllele1NoASEProbsSample2=lociAllele1NoASEProbsSample2, 
            lociRhosSample1=lociRhosSample1, 
            lociRhosSample2=lociRhosSample2, 
            numSim=numSim, 
            isPhased=isPhased, 
            tieBreakRandom=tieBreakRandom
        )
        return(res)
    }
    runMBASED2sTest <- function (
        lociAllele1CountsSample1, 
        lociAllele2CountsSample1, 
        lociAllele1CountsSample2, 
        lociAllele2CountsSample2, 
        lociAllele1NoASEProbsSample1, 
        lociAllele1NoASEProbsSample2, 
        lociRhosSample1, 
        lociRhosSample2, 
        aseIDs, 
        numSim=0, 
        BPPARAM=SerialParam(), 
        isPhased=FALSE, 
        tieBreakRandom=FALSE, 
        testSeeds ## one seed for each aseID
    ) {
        ASEResults <- bplapply(unique(aseIDs), function(aseID) {
            aseIDsubv <- (aseIDs==aseID)
            aseIDInd <- which(unique(aseIDs)==aseID)
            res <- runMBASED2s1aseIDTest(
                lociAllele1CountsSample1=lociAllele1CountsSample1[
                    aseIDsubv
                ], 
                lociAllele2CountsSample1=lociAllele2CountsSample1[
                    aseIDsubv
                ], 
                lociAllele1CountsSample2=lociAllele1CountsSample2[
                    aseIDsubv
                ], 
                lociAllele2CountsSample2=lociAllele2CountsSample2[
                    aseIDsubv
                ], 
                lociAllele1NoASEProbsSample1=lociAllele1NoASEProbsSample1[
                    aseIDsubv
                ], 
                lociAllele1NoASEProbsSample2=lociAllele1NoASEProbsSample2[
                    aseIDsubv
                ], 
                lociRhosSample1=lociRhosSample1[
                    aseIDsubv
                ], 
                lociRhosSample2=lociRhosSample2[
                    aseIDsubv
                ], 
                numSim=numSim, 
                isPhased=isPhased, 
                tieBreakRandom=tieBreakRandom, 
                testSeed=testSeeds[
                    aseIDInd
                ]
            ) 
            return(res)
        }, BPPARAM=BPPARAM)
        names(ASEResults) <- unique(aseIDs)
        allele1IsMajor <- unname(
            unlist(
                lapply(
                    ASEResults, 
                    function(el) {
                        el$lociAllele1IsMajor
                    }
                )
            )
        )
        lociMAFDifference <- unname(
            unlist(
                lapply(
                    ASEResults, 
                    function(el) {
                        el$lociMAFDifference
                    }
                )
            )
        )
        ASEResults <- lapply(ASEResults, function(el) {
            el$lociAllele1IsMajor <- NULL
            el$lociMAFDifference <- NULL
            return(el)
        })
        ASEResults <- as.data.frame(
            do.call(
                rbind,
                lapply(ASEResults, function(el) {
                    unlist(el)
                })
            )
        )
        return(
            list(
                ASEResults=ASEResults, 
                allele1IsMajor=allele1IsMajor,
                lociMAFDifference=lociMAFDifference
            )
        )
    }
    ## actual testing
    if (FULLTESTING) {
    	## number of test cases
        numCases <- 100 
    } else {
    	## number of test cases
        numCases <- 3 
    }
    ## 3 aseIDs, with 1, 3, and 4 loci each
    numLoci_for_aseIDs <- c(1, 3, 4) 
    ## number of total loci
    totalNumLoci <- sum(numLoci_for_aseIDs) 
    testCasesMatsSample1 <- getTestCases(
        numCases=numCases, 
        numLoci=totalNumLoci
    )
    testCasesMatsSample2 <- getTestCases(
        numCases=numCases, 
        numLoci=totalNumLoci
    )
    ## check for all possible values of tieBreakRandom
    for (test_tieBreakRandom in c(TRUE,FALSE)) { 
    	## check for all possible values of test_isPhased
        for (test_isPhased in c(TRUE, FALSE)) { 
        	## try sims and no-sims
            for (test_numSim in c(0, 3)) { 
                testResults <- bplapply(1:numCases, function(testInd) {
                    testCaseSample1 <- getTestCase(
                        testCasesMats=testCasesMatsSample1, 
                        testInd=testInd
                    )
                    testCaseSample2 <- getTestCase(
                        testCasesMats=testCasesMatsSample2, 
                        testInd=testInd
                    )
                    test_aseIDs <- rep(
                        paste(
                            'gene',
                            1:length(numLoci_for_aseIDs), 
                            sep=''
                        ), 
                        numLoci_for_aseIDs
                    )
                    test_seeds=sample(
                        1:(10^6), 
                        length(unique(test_aseIDs))
                    )
                    ## get joint results:
                    jointResults <- runMBASED2sTest(
                        lociAllele1CountsSample1=testCaseSample1$x, 
                        lociAllele2CountsSample1=testCaseSample1$n-
                            testCaseSample1$x, 
                        lociAllele1CountsSample2=testCaseSample2$x, 
                        lociAllele2CountsSample2=testCaseSample2$n-
                            testCaseSample2$x, 
                        lociAllele1NoASEProbsSample1=testCaseSample1$mu,
                        lociAllele1NoASEProbsSample2=testCaseSample2$mu, 
                        lociRhosSample1=testCaseSample1$rho, 
                        lociRhosSample2=testCaseSample2$rho,
                        aseIDs=test_aseIDs,
                        numSim=test_numSim, 
                        BPPARAM=SerialParam(),
                        isPhased=test_isPhased,
                        tieBreakRandom=test_tieBreakRandom,
                        testSeeds=test_seeds
                    )
                    ## get individual results:
                    runMBASED2s1aseIDOutputs <- lapply(
                        unique(test_aseIDs), 
                        function(test_aseID) {
                            test_aseID_subv <- (test_aseIDs==test_aseID)
                            test_aseID_ind <- which(
                                unique(test_aseIDs)==test_aseID
                            )
                            runMBASED2s1aseIDOutput <- 
                                runMBASED2s1aseIDTest(
                                    lociAllele1CountsSample1=testCaseSample1$x[
                                        test_aseID_subv
                                    ], 
                                    lociAllele2CountsSample1=(
                                        testCaseSample1$n-testCaseSample1$x
                                    )[test_aseID_subv], 
                                    lociAllele1CountsSample2=testCaseSample2$x[
                                        test_aseID_subv
                                    ], 
                                    lociAllele2CountsSample2=(
                                        testCaseSample2$n-testCaseSample2$x
                                    )[test_aseID_subv],
                                    lociAllele1NoASEProbsSample1=testCaseSample1$mu[
                                        test_aseID_subv
                                    ], 
                                    lociAllele1NoASEProbsSample2=testCaseSample2$mu[
                                        test_aseID_subv
                                    ], 
                            lociRhosSample1=testCaseSample1$rho[
                                test_aseID_subv
                            ], 
                            lociRhosSample2=testCaseSample2$rho[
                                test_aseID_subv
                            ],
                            numSim=test_numSim, 
                            isPhased=test_isPhased,
                            tieBreakRandom=test_tieBreakRandom,
                            testSeed=test_seeds[
                                test_aseID_ind
                            ]
                        )
                        return(runMBASED2s1aseIDOutput)
                    })
                    individualResults <- list(
                        ASEResults=data.frame(
                            majorAlleleFrequencyDifference=sapply(
                                runMBASED2s1aseIDOutputs, 
                                function(el) {
                                    el$majorAlleleFrequencyDifference
                                }
                            ),
                            pValueASE=sapply(
                                runMBASED2s1aseIDOutputs, 
                                function(el) {
                                    el$pValueASE
                                }
                            ),
                            heterogeneityQ=sapply(
                                runMBASED2s1aseIDOutputs, 
                                function(el) {
                                    el$heterogeneityQ
                                }
                            ),
                            pValueHeterogeneity=sapply(
                                runMBASED2s1aseIDOutputs, 
                                function(el) {
                                    el$pValueHeterogeneity
                                }
                            ),
                            nullHypothesisMAF=sapply(
                                runMBASED2s1aseIDOutputs, 
                                function(el) {
                                    el$nullHypothesisMAF
                                }
                            ),
                            row.names=unique(test_aseIDs)
                        ),
                        allele1IsMajor=unlist(
                            lapply(
                                runMBASED2s1aseIDOutputs, 
                                function(el) {
                                    el$lociAllele1IsMajor
                                }
                            )
                        ),
                        lociMAFDifference=unlist(
                            lapply(
                                runMBASED2s1aseIDOutputs, 
                                function(el) {
                                    el$lociMAFDifference
                                }
                            )
                        )
                    )
                    checkEquals(
                        jointResults, 
                        individualResults, 
                        msg='test_runMBASED2s: failed test checking that runMBASED1s does the same job as individual calls to runMBASED2s1aseID'
                    )
                }, BPPARAM=myBPPARAM)
            }
        }
    }
    
    return(TRUE)
}    
    

    
    

