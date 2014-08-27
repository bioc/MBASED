## size of round-off error I'm tolerating 
##(quantities with differences of this much or smaller are considered the same)
myPrecision <- .Machine$double.eps^0.5 
FULLTESTING <- FALSE
if (FULLTESTING) {
    myBPPARAM <- MulticoreParam(workers=24)
} else {
    myBPPARAM <- SerialParam()
}

## test MBASEDVectorizedMetaprop() function
## do not want to load packages VGAM, meta, 
## unless full testing mode is on
if (FULLTESTING) { 
    library(VGAM)
    library(meta)
    test_MBASEDVectorizedMetaprop <- function() { 
        ## check1: vectorization works for both a single locus gene 
        ## and multi-loci genes
        set.seed(125450)
        numCases <- 100
        for (numLoci in c(1, 4)) {
            test_ns_mat <- matrix(
                sample(10:20, numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
            test_xs_mat <- matrix(
                sample(0:10, numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
            test_mus_mat <- matrix(
                sample(seq(0.1, 0.9, by=0.1), numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
            test_rhos_mat <- matrix(
                sample(c(0,0.01, 0.3, 0.9), numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
            for (test_AH in c('greater', 'less', 'two.sided')) {
                individualResults <- bplapply(1:numCases, function(colInd) { 
                    MBASED:::MBASEDVectorizedMetaprop(
                        countsMat=matrix(
                            test_xs_mat[,colInd], 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        totalsMat=matrix(
                            test_ns_mat[,colInd], 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        probsMat=matrix(
                            test_mus_mat[,colInd], 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        rhosMat=matrix(
                            test_rhos_mat[,colInd], 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        alternative=test_AH
                    )
                }, BPPARAM=myBPPARAM)
                jointResults <- MBASED:::MBASEDVectorizedMetaprop(
                    countsMat=test_xs_mat, 
                    totalsMat=test_ns_mat, 
                    probsMat=test_mus_mat, 
                    rhosMat=test_rhos_mat, 
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
                        msg='test_MBASEDVectorizedMetaprop: fail vectorization test'
                    )
                }
            }
        }
            
        ## check2: for large n, reasonable mu/rho (or p, in case of binomial), 
        ## single-locus approach should produce pvalues similar to 
        ## (beta)binomial tests (apart from possibly very small p-values). 
        ## Further, the estimate should be (possibly adjusted) 
        ## back-transformed FT estimate 
        ## (we previously checked the behavior of adjustment 
        ## and back-transformation).
        set.seed(72286)
        test_n <- 10^6
        testCasesMat <- expand.grid(
            test_mu=seq(0.2, 0.8, by=0.1),
            test_rho=c(0, 0.01),
            test_AH=c('greater', 'less', 'two.sided')
        )
        numCountsToGenerate <- 10^6
        testResults <- bplapply(1:nrow(testCasesMat), function(testInd) {
            testCase <- testCasesMat[testInd,]
            testCounts <- with(testCase, 
                MBASED:::vectorizedRbetabinomMR(
                    n=numCountsToGenerate,
                    size=rep(test_n, numCountsToGenerate),
                    mu=rep(test_mu, numCountsToGenerate),
                    rho=rep(test_rho, numCountsToGenerate)
                )
            )
            predEstimates <- MBASED:::unFT(
                z=MBASED:::FTAdjust(
                    x=testCounts, 
                    n=rep(test_n, numCountsToGenerate), 
                    p=rep(testCase$test_mu, numCountsToGenerate)
                ), 
                n=rep(test_n, numCountsToGenerate)
            )
            if (isTRUE(all.equal(testCase$test_rho,0))) { ## binomial setting
                cumProbs <- pbinom(
                    q=testCounts, 
                    size=test_n, 
                    prob=testCase$test_mu
                )
            } else { ## beta-binomial setting
                cumProbs <- pbetabinom(
                    q=testCounts, 
                    size=test_n, 
                    prob=testCase$test_mu, 
                    rho=testCase$test_rho
                )
            }            
            if (testCase$test_AH=='greater') {
                predPvals <- 1-cumProbs
            } else if (testCase$test_AH=='less') {
                predPvals <- cumProbs
            } else {
                predPvals <- pmin(1,2*pmin(cumProbs, 1-cumProbs))
            }
            MBASEDVectorizedMetapropOutput <- with(testCase, 
                MBASED:::MBASEDVectorizedMetaprop(
                    countsMat=matrix(
                        testCounts, 
                        nrow=1, 
                        ncol=numCountsToGenerate
                    ), 
                    totalsMat=matrix(
                        rep(test_n, numCountsToGenerate), 
                        nrow=1, 
                        ncol=numCountsToGenerate
                    ), 
                    probsMat=matrix(
                        rep(test_mu, numCountsToGenerate), 
                        nrow=1, 
                        ncol=numCountsToGenerate
                    ), 
                    rhosMat=matrix(
                        rep(test_rho, numCountsToGenerate), 
                        nrow=1, 
                        ncol=numCountsToGenerate
                    ), 
                    alternative=as.character(test_AH)
                )
            )
            obsEstimates <- as.vector(MBASEDVectorizedMetapropOutput$propFinal)
            checkEqualsNumeric(
                obsEstimates, 
                predEstimates, 
                msg='test_MBASEDVectorizedMetaprop: fail comparison to binom/betabinom test'
            )
            obsPvals <- as.vector(MBASEDVectorizedMetapropOutput$pValue)
            checkTrue(
                all(
                    (pmax(obsPvals, predPvals)<=0.01) |
                    MBASED:::testNumericDiff(
                        queryVals=obsPvals, 
                        targetVals=predPvals, 
                        cutoffFraction=0.2 ## within 20%
                    )
                ),
                msg='test_MBASEDVectorizedMetaprop: fail comparison to binom/betabinom test'
            ) 
        }, BPPARAM=myBPPARAM)
    
    ## check3: should produce same results as package meta 
    ## (except possibly for MAF estimate, since that depends 
    ## on slightly different backtransformation) 
    ## for binomial counts and two-sided alternative.
        set.seed(336924)
        numCases <- 100 ## number of samples to test
        test_AHs <- c('greater', 'less', 'two.sided')
        for (numLoci in c(2,10)) { 
            test_ns_mat <- matrix(
                sample(10:20,  numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
            test_xs_mat <- matrix(
                sample(0:10,  numCases*numLoci, replace=T), 
                nrow=numLoci, 
                ncol=numCases
            )
            test_ps_mat <- matrix(
                0.5, 
                nrow=numLoci, 
                ncol=numCases
            )
            test_rhos_mat <- matrix(
                0, 
                nrow=numLoci, 
                ncol=numCases
            )
            for (test_AH in c('greater', 'less', 'two.sided')) {
                testResults <- bplapply(1:numCases, function(colInd) { 
                    MBASEDVectorizedMetapropOutputFixed <-
                        MBASED:::MBASEDVectorizedMetaprop(
                            countsMat=test_xs_mat[,colInd, drop=F], 
                            totalsMat=test_ns_mat[,colInd, drop=F],  
                            probsMat=test_ps_mat[,colInd, drop=F],  
                            rhosMat=test_rhos_mat[,colInd, drop=F],  
                            alternative=test_AH
                        )
                    metaOutput <- metaprop(
                        test_xs_mat[,colInd], 
                        test_ns_mat[,colInd], 
                        sm="PFT"
                    )
                    checkEqualsNumeric(
                        MBASEDVectorizedMetapropOutputFixed$hetQ, 
                        metaOutput$Q, 
                        msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                    )
                    checkEqualsNumeric(
                        MBASEDVectorizedMetapropOutputFixed$hetPVal, 
                        1-pchisq(metaOutput$Q, metaOutput$df.Q), 
                        msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                    )
                    checkEqualsNumeric(
                        MBASEDVectorizedMetapropOutputFixed$TEFinal, 
                        metaOutput$TE.fixed, 
                        msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                    )
                    checkEqualsNumeric(
                        MBASEDVectorizedMetapropOutputFixed$seTEFinal, 
                        metaOutput$seTE.fixed, 
                        msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                    )
                    MBASEDVectorizedMetaprop_backtransformN <- 
                        weighted.mean(
                            x=test_ns_mat[,colInd], 
                            w=test_ns_mat[,colInd]
                        )
                    metaOutput_backtransformN <- 
                        exp(mean(log(test_ns_mat[,colInd])))
                    checkEqualsNumeric(
                        MBASEDVectorizedMetapropOutputFixed$propFinal, 
                        meta:::asin2p(
                            metaOutput$TE.fixed,
                            MBASEDVectorizedMetaprop_backtransformN
                        ), 
                        msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                    )
                    metaOutputPropFixed <- suppressWarnings(
                        meta:::asin2p(
                            metaOutput$TE.fixed, 
                            metaOutput_backtransformN
                        )
                    )
                    diffBetweenEstimates <- abs(
                        MBASEDVectorizedMetapropOutputFixed$propFinal-
                        metaOutputPropFixed
                    )
                    minOfTwoEstimates <- min(c(
                        MBASEDVectorizedMetapropOutputFixed$propFinal,
                        metaOutputPropFixed
                    ))
                    checkTrue( 
                        (abs(
                            MBASEDVectorizedMetapropOutputFixed$propFinal-
                            metaOutputPropFixed
                        )<=0.01) || 
                        MBASED:::testNumericDiff(
                            queryVals=obsPvals, 
                            targetVals=predPvals, 
                            cutoffFraction=0.05 ## within 5%
                        ),
                        msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                    )
                    metaOutputZvalFixed <- 
                        (metaOutput$TE.fixed-pi/2)/(metaOutput$seTE.fixed)
                    if (test_AH=='greater') {
                        checkEqualsNumeric(
                            MBASEDVectorizedMetapropOutputFixed$pValue, 
                            1-pnorm(metaOutputZvalFixed), 
                            msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                        )
                    } else if (test_AH=='less') {
                        checkEqualsNumeric(
                            MBASEDVectorizedMetapropOutputFixed$pValue, 
                            pnorm(metaOutputZvalFixed), 
                            msg='test_MBASEDVectorizedMetaprop: fail meta comparison test'
                        )
                    } else if (test_AH=='two.sided') {
                        checkEqualsNumeric(
                            MBASEDVectorizedMetapropOutputFixed$pValue, 
                            min(c(
                                1,
                                2*min(c(
                                    pnorm(metaOutputZvalFixed), 
                                    1-pnorm(metaOutputZvalFixed)
                                ))
                            )), 
                            msg='test_MBASEDVectorizedMetaprop: fail meta comparison test')
                    } else {
                        stop('test_MBASEDVectorizedMetaprop: fail meta comparison test: unknown alternative hypothesis')
                    }
                }, BPPARAM=myBPPARAM)    
            }
        }        
        return(TRUE)
    }
}    
    
    

## GENERAL FRAMEWORK FOR TEST CASES FOR
## test_runMBASED1s1aseID() and test_runMBASED1s(): 
## Try out n=10,11,12,..,20, 50, 51, 100, 101, 1000, 1001. 
## For each n, test (when possible) x={0,1,2}(small), 
## {n, n-1, n-2}(large), {floor(n/2-1), floor(n/2), floor(n/2+1)} (mean), 
## and 10 random values.
## For each combination of n and x, try mu={0.4, 0.45, 0.5, 0.55, 0.6}
## For each combination of n, x, mu, try rho={0, 0.01, 0.1}
set.seed(897134)
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
    candidate_xs <- candidate_xs[candidate_xs>=0 & candidate_xs<=n]
    candidate_xs <- unique(sort(candidate_xs, decreasing=FALSE))
    return(paste(n, candidate_xs, sep='_'))
}))
AllTestCasesMat <- expand.grid(
    test_xn_combo=my_xn_combos,
    test_mu=seq(0.4, 0.6, by=0.05),
    test_rho=c(0, 0.01, 0.1)
)
AllTestCasesMat$test_n <- as.numeric(
    matrix(
        unlist(strsplit(as.character(AllTestCasesMat$test_xn_combo), '_')), 
        ncol=2, 
        byrow=TRUE
    )[,1]
)
AllTestCasesMat$test_x <- as.numeric(
    matrix(
        unlist(strsplit(as.character(AllTestCasesMat$test_xn_combo), '_')), 
        ncol=2, 
        byrow=TRUE
    )[,2]
)
rm(
    my_ns,
    my_xn_combos
)

## helper function to get test cases from the matrix of all test cases
## returns a matrix with test cases of interest spread over a desired matrix
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

## helper function to get individual test case (column) 
## from output of getTestCases()
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


## test behavior of runMBASED1s1aseID() function 
test_runMBASED1s1aseID <- function() { 
    ## check1: for single-locus gene argument isPhased 
    ## has no effect on the output.
    set.seed(898639)
    if (FULLTESTING) {
        numCases <- 1000 ## number of test cases
    } else {
        numCases <- 3 ## number of test cases
    }
    numLoci <- 1 ## number of loci
    testCasesMats <- getTestCases(
        numCases=numCases, 
        numLoci=numLoci
    )
    ## check for all possible values of tieBreakRandom
    for (test_tieBreakRandom in c(TRUE,FALSE)) { 
    	## check for both simulations and no-simulations cases
        for (test_numSim in c(0, 3)) { 
        	## this is to ensure that the same simulations are used 
        	## to produce p-values in case of simulated data
            testIndSeeds <- sample(1:(10^6), numCases, replace=TRUE) 
            testResults <- bplapply(1:numCases, function(testInd) {
                testCase <- getTestCase(
                    testCasesMats=testCasesMats, 
                    testInd=testInd
                )
                test_seed <- testIndSeeds[testInd]
                testOutputs <- lapply(c(TRUE, FALSE), function(test_isPhased) {
                    set.seed(test_seed) 
                    MBASED:::runMBASED1s1aseID(
                        lociAllele1Counts=testCase$x, 
                        lociAllele2Counts=testCase$n-testCase$x, 
                        lociAllele1NoASEProbs=testCase$mu, 
                        lociRhos=testCase$rho, 
                        numSim=test_numSim, 
                        isPhased=test_isPhased, 
                        tieBreakRandom=test_tieBreakRandom
                    )
                })
                checkEquals(
                    testOutputs[[1]], 
                    testOutputs[[2]], 
                    msg='test_runMBASED1s1aseID: failed test to see that isPhased does not affect output'    
                )            
            }, BPPARAM=myBPPARAM)
        }
    }
    rm(
        numCases, 
        numLoci, 
        testCasesMats,
        test_tieBreakRandom,
        test_numSim,
        testIndSeeds,
        testResults
    )

    ## check 2: argument tieBreakRandom does not affect allele frequencies 
    ## or heterogeneity statistic. It does not affect p-values if no simulations 
    ## are performed (some effect will be observed with simulations 
    ## due to random variation)
    set.seed(643129)
    if (FULLTESTING) {
        numCases <- 1000 ## number of test cases
    } else {
        numCases <- 3 ## number of test cases
    }
    for (numLoci in c(1,5)) {## number of loci 
        testCasesMats <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        ## check for both simulations and no-simulations cases
        for (test_numSim in c(0, 3)) { 
        	## check for both settings of isPhased
            for (test_isPhased in c(TRUE, FALSE)) { 
                testResults <- bplapply(1: numCases, function(testInd) {
                    testCase <- getTestCase(
                        testCasesMats=testCasesMats, 
                        testInd=testInd
                    )
                    testOutputs <- lapply(
                        c(TRUE, FALSE), 
                        function(test_tieBreakRandom) {
                            MBASED:::runMBASED1s1aseID(
                                lociAllele1Counts=pmax(1,testCase$x), 
                                ## make sure we have ties 
                                lociAllele2Counts=pmax(1,testCase$x), 
                                lociAllele1NoASEProbs=testCase$mu, 
                                lociRhos=testCase$rho, 
                                numSim=test_numSim, 
                                isPhased=test_isPhased, 
                                tieBreakRandom=test_tieBreakRandom
                            )    
                        }
                    )
                    outputNamesToCheck <- setdiff(
                        names(testOutputs[[1]]), 
                        'lociAllele1IsMajor'
                    )
                    if (test_numSim>0) {
                        outputNamesToCheck <- setdiff(
                            outputNamesToCheck, 
                            c('pValueASE', 'pValueHeterogeneity')
                        )
                    }
                    for (outputName in outputNamesToCheck) {
                        checkEqualsNumeric(
                            testOutputs[[1]][[outputName]],
                            testOutputs[[2]][[outputName]],
                            msg='test_runMBASED1s1aseID: failed test to see that tieBreakRandom does not affect output apart from assigning alleles to haplotypes'
                        )
                    }
                }, BPPARAM=myBPPARAM)
            }
        }
    }
    rm(
        numCases, 
        numLoci, 
        testCasesMats,
        test_numSim,
        test_isPhased,
        testResults
    )
        
    ## check 3:  with no simulations return the same result as 
    ## MBASEDVectorizedMetaprop with alternative hypothesis set to
    ## 'two.sided' and with additional step of picking the larger 
    ## of allelic frequencies to be 'major'
    ## use tieBreakRandom=FALSE for simplicity. 
    ## Later will check that tie-breaking works as intended.
    set.seed(136565)
    if (FULLTESTING) {
        numCases <- 1000 ## number of test cases
    } else {
        numCases <- 3 ## number of test cases
    }
    for (numLoci in c(1,5)) { ## number of loci
        testCasesMats <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        testResults <- bplapply(1:numCases, function(testInd) {
            testCase <- getTestCase(
                testCasesMats=testCasesMats, 
                testInd=testInd
            )
            for (test_isPhased in c(TRUE, FALSE)) {
                if (!test_isPhased) { ## data will be phased by MBASED
                    test_x_is_major <- MBASED:::isCountMajorFT(
                        x=testCase$x, 
                        n=testCase$n, 
                        p=testCase$mu, 
                        tieBreakRandom=FALSE
                    )
                    test_maj <- ifelse(
                        test_x_is_major, 
                        testCase$x, 
                        testCase$n-testCase$x
                    )
                    test_prob_maj <- ifelse(
                        test_x_is_major, 
                        testCase$mu, 
                        1-testCase$mu
                    )
                } else { ## data will not be phased by MBASED
                	## whether it's really major or not, will be determined later
                    test_maj <- testCase$x 
                    test_prob_maj <- testCase$mu
                }
                MBASEDVectorizedMetapropOutput <- 
                    MBASED:::MBASEDVectorizedMetaprop(
                        countsMat=matrix(
                            test_maj, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        totalsMat=matrix(
                            testCase$n, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        probsMat=matrix(
                            test_prob_maj,
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        rhosMat=matrix(
                            testCase$rho, 
                            nrow=numLoci, 
                            ncol=1
                        ), 
                        alternative='two.sided'
                    )    
               ## now we can check to see if the supplied allele was major or not
                if (test_isPhased) { 
                    if (MBASEDVectorizedMetapropOutput$propFinal<0.5) {
                        test_x_is_major=rep(FALSE, length(testCase$x))
                    } else {
                        test_x_is_major=rep(TRUE, length(testCase$x))
                    }
                }
                runMBASED1s1aseIDOutput <- 
                    MBASED:::runMBASED1s1aseID(
                        lociAllele1Counts=testCase$x, 
                        lociAllele2Counts=testCase$n-testCase$x, 
                        lociAllele1NoASEProbs=testCase$mu, 
                        lociRhos=testCase$rho, 
                        numSim=0, 
                        isPhased=test_isPhased,
                        tieBreakRandom=FALSE
                    )
                checkEqualsNumeric(
                    as.vector(pmax(
                        MBASEDVectorizedMetapropOutput$propFinal, 
                        1-MBASEDVectorizedMetapropOutput$propFinal)
                    ), 
                    runMBASED1s1aseIDOutput$majorAlleleFrequency, 
                    msg=paste('test_runMBASED1s1aseID: failed comaprison with MBASEDVectorizedMetaprop test')
                )
                checkEqualsNumeric(
                    as.vector(MBASEDVectorizedMetapropOutput$pValue), 
                    runMBASED1s1aseIDOutput$pValueASE, 
                    msg=paste('test_runMBASED1s1aseID: failed comaprison with MBASEDVectorizedMetaprop test')
                )
                checkEqualsNumeric(
                    as.vector(MBASEDVectorizedMetapropOutput$hetQ), 
                    runMBASED1s1aseIDOutput$heterogeneityQ, 
                    msg=paste('test_runMBASED1s1aseID: failed comaprison with MBASEDVectorizedMetaprop test')
                )
                checkTrue(
                    ifelse(
                        numLoci==1,
                        is.na(runMBASED1s1aseIDOutput$heterogeneityQ), 
                        !is.na(runMBASED1s1aseIDOutput$heterogeneityQ)
                    ),
                    msg=paste('test_runMBASED1s1aseID: failed comaprison with MBASEDVectorizedMetaprop test')
                )
                checkEqualsNumeric(
                    as.vector(MBASEDVectorizedMetapropOutput$hetPVal), 
                    runMBASED1s1aseIDOutput$pValueHeterogeneity, 
                    msg=paste('test_runMBASED1s1aseID: failed comaprison with MBASEDVectorizedMetaprop test')
                )
                checkTrue(
                    ifelse(
                        numLoci==1,
                        is.na(runMBASED1s1aseIDOutput$pValueHeterogeneity), 
                        !is.na(runMBASED1s1aseIDOutput$pValueHeterogeneity)
                    ),
                    msg=paste('test_runMBASED1s1aseID: failed comaprison with MBASEDVectorizedMetaprop test')
                )
                checkIdentical(
                    runMBASED1s1aseIDOutput$lociAllele1IsMajor,
                     test_x_is_major, 
                     msg=paste('test_runMBASED1s1aseID: failed comaprison with MBASEDVectorizedMetaprop test')
                )
            }
        }, BPPARAM=myBPPARAM)
    }
    rm(
        numCases, 
        numLoci, 
        testCasesMats,
        testResults
    )
    
    ## check4: breaking ties works as intended.
    set.seed(378473)
    numCases <- 10^3
    tiedVal <- 10
    for (numLoci in c(1, 5)) {
        for (test_numSim in c(0, 3)) {
            noTieBreaking <- unlist(bplapply(1:numCases, function(testInd) {
                MBASED:::runMBASED1s1aseID(
                    lociAllele1Counts=rep(tiedVal, numLoci), 
                    lociAllele2Counts=rep(tiedVal, numLoci), 
                    lociAllele1NoASEProbs=rep(0.5, numLoci), 
                    lociRhos=rep(0, numLoci), 
                    numSim=test_numSim, 
                    isPhased=sample(c(TRUE,FALSE),1), 
                    tieBreakRandom=FALSE
                )$lociAllele1IsMajor
            }, BPPARAM=myBPPARAM))
            checkTrue(
                length(unique(noTieBreaking))==1, 
                msg='test_runMBASED1s1aseID: failed tieBreak test: no simulations, no tie breaks'
            )
            yesTieBreaking <- unlist(bplapply(
                1:numCases, 
                function(testInd) {
                    MBASED:::runMBASED1s1aseID(
                        lociAllele1Counts=rep(tiedVal, numLoci), 
                        lociAllele2Counts=rep(tiedVal, numLoci), 
                        lociAllele1NoASEProbs=rep(0.5, numLoci), 
                        lociRhos=rep(0, numLoci), 
                        numSim=0, 
                        isPhased=sample(c(TRUE,FALSE),1), 
                        tieBreakRandom=TRUE
                    )$lociAllele1IsMajor
                }, BPPARAM=myBPPARAM
            ))
            MBASED:::testQuantiles(
                theoreticalCumDist=0.5, 
                observedCumDist=mean(yesTieBreaking), 
                numTotalCounts=numCases, 
                numSEsToCheck=5, 
                errorMessage='test_runMBASED1s1aseID: failed tieBreak test: no simulations, yes tie breaks'
            )
        }
    }
    rm(
        numLoci,
        numCases,
        test_numSim,
        tiedVal,
        noTieBreaking,
        yesTieBreaking
    )
    
    ## check5: with simulations, MAF estimate and assignment of alleles to haplotypes is the same as without simulations
    set.seed(784069)
    if (FULLTESTING) {
        numCases <- 1000 ## number of test cases
    } else {
        numCases  <- 3 ## number of test cases
    }
    for (numLoci in c(1, 5)) { ## number of loci
        testCasesMats <- getTestCases(
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
                    testCase <- getTestCase(
                        testCasesMats=testCasesMats, 
                        testInd=testInd
                    )
                    test_seed <- testIndSeeds[testInd]
                    testOutputs <- lapply(c(0, 10), function(test_numSim) {
                    	## make sure ties get broken in same way
                        set.seed(test_seed) 
                        MBASED:::runMBASED1s1aseID(
                            lociAllele1Counts=testCase$x, 
                            lociAllele2Counts=testCase$n-testCase$x, 
                            lociAllele1NoASEProbs=testCase$mu, 
                            lociRhos=testCase$rho, 
                            numSim=test_numSim, 
                            isPhased=test_isPhased, 
                            tieBreakRandom=test_tieBreakRandom
                        )
                    })
                    checkEqualsNumeric(
                        testOutputs[[1]]$majorAlleleFrequency, 
                        testOutputs[[2]]$majorAlleleFrequency, 
                        msg='test_runMBASED1s1aseID: failed test checking that MAF estimate is the same for sim and no sim'
                    )
                    checkEquals(
                        testOutputs[[1]]$lociAllele1IsMajor, 
                        testOutputs[[2]]$lociAllele1IsMajor, 
                        msg='test_runMBASED1s1aseID: failed test checking that assignment of alleles to haplotypes is the same for sim and no sim'
                    )
                }, BPPARAM=myBPPARAM)
            }
        }
    }
    rm(
        numCases, 
        numLoci, 
        testCasesMats,
        test_tieBreakRandom,
        test_isPhased,
        testIndSeeds,
        testResults
    )
    
    ## check6: for single-locus gene with simulations, 
    ## results should be similar to those obtained from binomial 
    ## or beta-binomial tests (as long as number of simulations 
    ## is large enough)
    ## no looping over test_isPhased, since this argument 
    ## has no effect on output for single-locus genes
    set.seed(988358)
    ## only want to load package VGAM if full testing mode is on
    if (FULLTESTING) { 
        library(VGAM)
        numCases <- 1000 ## number of test cases
        numLoci <- 1 ## number of loci
        test_numSim <- 10^6 ## number of simulations
        testCasesMats <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        ## the following logic is used: if the true p-values are very small, 
        ## we may be prevented from properly assessing them 
        ## due to granularity of empirical p-value estimates 
        ## (smallest p-value we can achieve is 1/test_numSim). 
        ## Such cases are identified in the code below and 
        ## for such instances, p-value comparison is NOT performed.
        
        ## if CI for the predicted p-value is not at least this many SDs 
        ## away from 0 or 1, we will not test that p-value
        numSDsToCheck <- 5 
        testResults <- bplapply(1:numCases, function(testInd) {
            testCase <- getTestCase(
                testCasesMats=testCasesMats, 
                testInd=testInd
            )
            ## NOTE: I am looking at likelihood of observing MAFs, 
            ## and NOT actual count outcomes. In extreme cases 
            ## (e.g., n=17, x=7, mu=0.6, rho=0.3), there may be counts 
            ## yielding higher MAF (after phasing) but also having 
            ## higher likelihood of being observed (the density of counts 
            ## is not modal). In such situations the standard two-sided tests 
            ## (that look at likelihoods of individual counts) 
            ## will give different results.
            runMBASED1s1aseIDOutputMAF <- 
            ## know that MAF is the same under simultions and no-simulations
                MBASED:::runMBASED1s1aseID(
                    lociAllele1Counts=testCase$x, 
                    lociAllele2Counts=testCase$n-testCase$x, 
                    lociAllele1NoASEProbs=testCase$mu, 
                    lociRhos=testCase$rho, 
                    numSim=0, 
                    isPhased=sample(c(TRUE,FALSE),1), 
                    tieBreakRandom=sample(c(TRUE,FALSE),1)
                )$majorAlleleFrequency 
            test_all_AFs <- MBASED:::unFT(
                z=MBASED:::FTAdjust(
                    x=0:(testCase$n), 
                    n=rep(testCase$n, testCase$n+1), 
                    p=rep(testCase$mu, testCase$n+1)
                ), 
                n=rep(testCase$n, testCase$n+1)
            )
            test_all_MAFs <- pmax(
                test_all_AFs, 
                1-test_all_AFs
            )
            predLikelihoods <- dbetabinom(
                x=0:(testCase$n), 
                size=testCase$n, 
                prob=testCase$mu, 
                rho=testCase$rho
            )
            predPval <- sum(
                predLikelihoods[
                    test_all_MAFs>=(runMBASED1s1aseIDOutputMAF-myPrecision)
                 ]
            )
            if (
                isTRUE(all.equal(predPval, 1)) || 
                isTRUE(all.equal(predPval, 0))
            ) {
                obsPval <- MBASED:::runMBASED1s1aseID(
                    lociAllele1Counts=testCase$x, 
                    lociAllele2Counts=testCase$n-testCase$x, 
                    lociAllele1NoASEProbs=testCase$mu, 
                    lociRhos=testCase$rho, 
                    numSim=test_numSim, 
                    isPhased=sample(c(TRUE,FALSE),1), 
                    tieBreakRandom=sample(c(TRUE,FALSE),1)
                )$pValueASE
                checkEqualsNumeric(
                    predPval, 
                    obsPval, 
                    tolerance=myPrecision, 
                    msg='test_runMBASED1s1aseID_singleLocus: 1-sample does not approximate (beta)-binomial p-values'
                ) 
            } else {
                predPvalSD <- sqrt(predPval*(1-predPval)/test_numSim)
                if (
                    (predPval-numSDsToCheck*predPvalSD)>0 && 
                    (predPval+numSDsToCheck*predPvalSD)<1 
                ) {
                    obsPval <- MBASED:::runMBASED1s1aseID(
                        lociAllele1Counts=testCase$x, 
                        lociAllele2Counts=testCase$n-testCase$x, 
                        lociAllele1NoASEProbs=testCase$mu, 
                        lociRhos=testCase$rho, 
                        numSim=test_numSim, 
                        isPhased=sample(c(TRUE,FALSE),1), 
                        tieBreakRandom=sample(c(TRUE,FALSE),1)
                    )$pValueASE
                    checkTrue(
                        (abs(obsPval-predPval)/predPvalSD)<5, 
                        msg='test_runMBASED1s1aseID_singleLocus: 1-sample does not approximate (beta)-binomial p-values'
                    )
                }
            }
        }, BPPARAM=myBPPARAM)
        rm(
            numCases, 
            numLoci, 
            testCasesMats,
            test_numSim,
            numSDsToCheck,
            testResults
        )
    }

    ## check 7: for multi-locus genes, under null hypothesis 
    ## and for reasonable rho, the ASE and heterogeneity p-values 
    ## are approximately uniform. No need to check for single-locus gene: 
    ## there's no heterogeneity and already saw that single-locus 
    ## p-values are close to those from (beta)binomial tests. 
    set.seed(296092)
    ## only do this in full testing mode, otherwise it takes way too long
    if (FULLTESTING) { 
    	## number of test cases
        numCases <- 5
        ## number of loci
        numLoci <- 3 
        ## number of simulations to perform by MBASED
        test_numSim <- 10^4 
        ## number of null p-values to be calculated to assess 
        ## if distribution is close to uniform or not 
        ## (~ 2 min for each set of settings)
        numTests <- 10^4  
        testCasesMats <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
        quantilesToCheck <- seq(0.1, 0.9, by=0.1)
        ## will compare to see if observed quantile is within 
        ## this many standard errors of the true value. 
        ## Need to have this number relatively high, because of 
        ## a) many tests, 
        ## b) actual chunkiness in observed p-value for smallish n. 
        ## Given current numTests of 10k, the biggest SE I allow 
        ## is for the case of the median (F=0.5, predFSD=0.03), 
        ## which is not super huge
        numSEsToCheck=6 
        ## check for all possible values of test_isPhased
        for (test_isPhased in c(TRUE,FALSE)) { 
        	## check for all possible values of tieBreakRandom
            for (test_tieBreakRandom in c(TRUE,FALSE)) { 
            	## no parallelizing here, because 
            	## parallelizing is done over numTests
                for (testInd in 1:numCases) { 
                    testCase <- getTestCase(
                        testCasesMats=testCasesMats, 
                        testInd=testInd
                    )
                    test_nullResults <- bplapply(1:numTests, function(runInd) {
                        test_x_null <- MBASED:::vectorizedRbetabinomMR(
                            n=length(testCase$n), 
                            size=testCase$n, 
                            mu=testCase$mu, 
                            rho=testCase$rho
                        )
                        MBASED:::runMBASED1s1aseID(
                            lociAllele1Counts=test_x_null, 
                            lociAllele2Counts=testCase$n-test_x_null, 
                            lociAllele1NoASEProbs=testCase$mu, 
                            lociRhos=testCase$rho, 
                            numSim=test_numSim, 
                            isPhased=test_isPhased,
                            tieBreakRandom=test_tieBreakRandom
                        )
                    }, BPPARAM=myBPPARAM)
                    test_null_ASEPvals <- sapply(test_nullResults, function(el) {
                        el$pValueASE
                    })
                    test_null_hetPvals <- sapply(test_nullResults, function(el) {
                        el$pValueHeterogeneity
                    })
                    MBASED:::testQuantiles(
                        theoreticalCumDist=quantilesToCheck, 
                        observedCumDist=sapply(quantilesToCheck, function(el) {
                            mean(test_null_ASEPvals<=(el+myPrecision))
                        }), 
                        numTotalCounts=numTests, 
                        numSEsToCheck=numSEsToCheck, 
                        errorMessage='test_runMBASED1s1aseID_allCases: failed uniformity of ASE p-values test for multi-locus genes'
                    )
                    MBASED:::testQuantiles(
                        theoreticalCumDist=quantilesToCheck, 
                        observedCumDist=sapply(quantilesToCheck, function(el) {
                            mean(test_null_hetPvals<=(el+myPrecision))
                        }), 
                        numTotalCounts=numTests, 
                        numSEsToCheck=numSEsToCheck, 
                        errorMessage='test_runMBASED1s1aseID_allCases: failed uniformity of heterogeneity p-values test for multi-locus genes'
                    )

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
        
    ## Check 8: check that we can detect true signal: 
    ## p-value should be small. 
    set.seed(588586)
    ## only do for full testing, otherwise takes too long
    if (FULLTESTING) { 
    	## number of test cases
        numCases <- 100 
        ## number of loci
        for (numLoci in c(1, 5))  {
        	## number of simulations to perform by MBASED
            test_numSim <- 10^6 
            trueMAF=0.9
            testCasesMats <- getTestCases(
                numCases=numCases, 
                numLoci=numLoci
            )
            ## check for all possible values of isPhased
            for (test_isPhased in c(TRUE,FALSE)) { 
            	## check for all possible values of tieBreakRandom
                for (test_tieBreakRandom in c(TRUE,FALSE)) { 
                    testResults <- bplapply(1:numCases, function(testInd) {
                        testCase <- getTestCase(
                            testCasesMats=testCasesMats, 
                            testInd=testInd
                        )
                        ## make sure we have power to detect ASE
                        testCase$n <- 100*testCase$n 
                        ## trueAF must be a single number
                        test_generating_prob <- MBASED:::getPFinal( 
                            trueAF=trueMAF, 
                            noASEAF=testCase$mu
                        )
                        test_x <- round(test_generating_prob*testCase$n)
                        runMBASED1s1aseIDOutput <- 
                            MBASED:::runMBASED1s1aseID(
                                lociAllele1Counts=test_x, 
                                lociAllele2Counts=testCase$n-test_x, 
                                lociAllele1NoASEProbs=testCase$mu, 
                                lociRhos=testCase$rho, 
                                numSim=test_numSim, 
                                isPhased=test_isPhased,
                                tieBreakRandom=test_tieBreakRandom
                            )
                        checkTrue(
                            runMBASED1s1aseIDOutput$pValueASE <= 0.05, 
                            msg='test_runMBASED1s1aseID: 1-sample does not detect the signal'
                        )
                    }, BPPARAM=myBPPARAM)
                }
            }
        }
        rm(
            numCases, 
            numLoci, 
            test_numSim,
            trueMAF,
            testCasesMats,
            test_isPhased,
            test_tieBreakRandom,
            testResults
        )
    }
    
    ## Check 9: For multi-SNV genes only: check that we can detect 
    ## true SNV-specific ASE: heterogeneity p-value should be small 
    ## (for largish n and smallish rho, otherwise not enough evidence). 
    ## This takes about 10 mins
    set.seed(660289)
    ## only do for full testing, otherwise takes too long
    if (FULLTESTING) { 
    	## number of test cases
        numCases <- 100 
        ## number of loci
        numLoci <- 3 
        ## number of simulations to perform by MBASED
        test_numSim <- 10^6
        trueMAF=0.9
        testCasesMats <- getTestCases(
            numCases=numCases, 
            numLoci=numLoci
        )
         ## check for all possible values of isPhased
        for (test_isPhased in c(TRUE,FALSE)) {
        	## check for all possible values of tieBreakRandom
            for (test_tieBreakRandom in c(TRUE,FALSE)) { 
                testResults <- bplapply(1:numCases, function(testInd) {
                    testCase <- getTestCase(
                        testCasesMats=testCasesMats, 
                        testInd=testInd
                    )
                    whichSNVisASE <- sample(1:length(testCase$n), 1)
                    trueMAFs <- rep(0.5, length(testCase$n))
                    trueMAFs[whichSNVisASE] <- trueMAF
                    ## make sure we have power to detect ASE
                    testCase$n <- 100*testCase$n 
                    test_generating_prob <- sapply(
                        1:length(trueMAFs), 
                        function(MAFInd) {
                        	## trueAF must be a single number
                            MBASED:::getPFinal(
                                trueAF=trueMAFs[MAFInd], 
                                noASEAF=testCase$mu[MAFInd]
                            )
                        }
                    )
                    test_x <- round(test_generating_prob*testCase$n)
                    runMBASED1s1aseIDOutput <- 
                        MBASED:::runMBASED1s1aseID(
                            lociAllele1Counts=test_x, 
                            lociAllele2Counts=testCase$n-test_x, 
                            lociAllele1NoASEProbs=testCase$mu, 
                            lociRhos=testCase$rho, 
                            numSim=test_numSim, 
                            isPhased=test_isPhased,
                            tieBreakRandom=test_tieBreakRandom
                        )
                    checkTrue(
                        runMBASED1s1aseIDOutput$pValueHeterogeneity <= 0.2, 
                        msg='test_runMBASED1s1aseID: 1-sample does not detect the heterogeneity signal '
                    )
                }, BPPARAM=myBPPARAM)
            }
        }
        rm(
            numCases, 
            numLoci, 
            test_numSim,
            trueMAF,
            testCasesMats,
            test_isPhased,
            test_tieBreakRandom,
            testResults
        )
    }
    
    return(TRUE)
}

## test behavior of runMBASED1s() function 
test_runMBASED1s <- function() {
    ## check1: runMBASED1s gives same results as 
    ## performing analysis on each aseID individually.
    ## Complication due to random sampling I employ: 
    ## bplapply messes with random seed, so I need to write 
    ## slightly modified version of functions to directly specificy seeds.
    set.seed(817815)
    runMBASED1s1aseIDTest <- function(
        lociAllele1Counts, 
        lociAllele2Counts,
        lociAllele1NoASEProbs, 
        lociRhos, 
        numSim=0, 
        isPhased=FALSE, 
        tieBreakRandom=FALSE, 
        testSeed
    ) {
        set.seed(testSeed)
        res <- MBASED:::runMBASED1s1aseID(
            lociAllele1Counts=lociAllele1Counts, 
            lociAllele2Counts=lociAllele2Counts, 
            lociAllele1NoASEProbs=lociAllele1NoASEProbs, 
            lociRhos=lociRhos, 
            numSim=numSim, 
            isPhased=isPhased, 
            tieBreakRandom=tieBreakRandom
        )
        return(res)
    }
    runMBASED1sTest <- function(
        lociAllele1Counts, 
        lociAllele2Counts, 
        lociAllele1NoASEProbs, 
        lociRhos, 
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
            res <- runMBASED1s1aseIDTest(
                lociAllele1Counts=lociAllele1Counts[aseIDsubv],
                lociAllele2Counts=lociAllele2Counts[aseIDsubv],
                lociAllele1NoASEProbs=lociAllele1NoASEProbs[aseIDsubv],
                lociRhos=lociRhos[aseIDsubv],
                numSim=numSim,
                isPhased=isPhased,
                tieBreakRandom=tieBreakRandom,
                testSeed=testSeeds[aseIDInd]
            ) 
            return(res)
        }, BPPARAM=BPPARAM)
        names(ASEResults) <- unique(aseIDs)
        allele1IsMajor <- as.vector(
            unlist(
                lapply(
                    ASEResults, 
                    function(el) el$lociAllele1IsMajor
                )
            )
        )
        lociMAF <- as.vector(
            unlist(
                lapply(
                    ASEResults, 
                    function(el) el$lociMAF
                )
            )
        )
        ASEResults <- lapply(ASEResults, function(el) {
            el$lociAllele1IsMajor <- NULL
            el$lociMAF <- NULL
            return(el)
        })
        ASEResults <- as.data.frame(
            do.call(
                rbind,
                lapply(ASEResults, function(el) unlist(el))
            )
        )
        return(
            list(
                ASEResults=ASEResults, 
                allele1IsMajor=allele1IsMajor,
                lociMAF=lociMAF
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
    testCasesMats <- getTestCases(
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
                    testCase <- getTestCase(
                        testCasesMats=testCasesMats, 
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
                    test_seeds <- sample(
                        1:(10^6), 
                        length(unique(test_aseIDs))
                    )
                    ## get joint results:
                    jointResults <- runMBASED1sTest(
                        lociAllele1Counts=testCase$x, 
                        lociAllele2Counts=testCase$n-testCase$x, 
                        lociAllele1NoASEProbs=testCase$mu, 
                        lociRhos=testCase$rho, 
                        aseIDs=test_aseIDs,
                        numSim=test_numSim, 
                        BPPARAM=SerialParam(),
                        isPhased=test_isPhased,
                        tieBreakRandom=test_tieBreakRandom,
                        testSeeds=test_seeds
                    )
                    ## get individual results:
                    runMBASED1s1aseIDOutputs <- lapply(
                        unique(test_aseIDs), 
                        function(test_aseID) {
                            test_aseID_subv <- (test_aseIDs==test_aseID)
                            test_aseID_ind <- which(
                                unique(test_aseIDs)==test_aseID
                            )
                            runMBASED1s1aseIDOutput <- 
                                runMBASED1s1aseIDTest(
                                    lociAllele1Counts=testCase$x[
                                        test_aseID_subv
                                    ], 
                                    lociAllele2Counts=(testCase$n-testCase$x)[
                                        test_aseID_subv
                                    ], 
                                    lociAllele1NoASEProbs=testCase$mu[
                                        test_aseID_subv
                                    ], 
                                    lociRhos=testCase$rho[
                                        test_aseID_subv
                                    ], 
                                    numSim=test_numSim, 
                                    isPhased=test_isPhased,
                                    tieBreakRandom=test_tieBreakRandom,
                                    testSeed=test_seeds[
                                       test_aseID_ind
                                    ]
                                )
                            return(runMBASED1s1aseIDOutput)
                        }
                    )
                    individualResults <- list(
                        ASEResults=data.frame(
                            majorAlleleFrequency=sapply(
                                runMBASED1s1aseIDOutputs, 
                                function(el) {
                                    el$majorAlleleFrequency
                                }
                            ),
                            pValueASE=sapply(
                                runMBASED1s1aseIDOutputs, 
                                function(el) {
                                    el$pValueASE
                                }
                            ),
                            heterogeneityQ=sapply(
                                runMBASED1s1aseIDOutputs, 
                                function(el) {
                                    el$heterogeneityQ
                                }
                            ),
                            pValueHeterogeneity=sapply(
                                runMBASED1s1aseIDOutputs, 
                                function(el) {
                                    el$pValueHeterogeneity
                                }
                            ),
                            row.names=unique(test_aseIDs)
                        ),
                        allele1IsMajor=unlist(
                            lapply(
                                runMBASED1s1aseIDOutputs, 
                                function(el) {
                                    el$lociAllele1IsMajor
                                }
                            )
                        ),
                        lociMAF=unlist(
                            lapply(
                                runMBASED1s1aseIDOutputs, 
                                function(el) {
                                    el$lociMAF
                                }
                            )
                        )
                    )
                    checkEquals(
                        jointResults, 
                        individualResults, 
                        msg='test_runMBASED1s: failed test checking that runMBASED1s does the same job as individual calls to runMBASED1s1aseID'
                    )
                }, BPPARAM=myBPPARAM)
            }
        }
    }
    
    return(TRUE)
}    
    

    
    

