#' Function that runs between-sample (differential) ASE calling using data from loci (SNVs) within a single unit of ASE (gene). The i-th entry of each of vector arguments 'lociAllele1CountsSample1', 'lociAllele2CountsSample1', 'lociAllele1NoASEProbsSample1', 'lociRhosSample1', 'lociAllele1CountsSample2', 'lociAllele2CountsSample2', 'lociAllele1NoASEProbsSample2', and 'lociRhosSample2' should correspond to the i-th locus. If argument 'isPhased' (see below) is true, then entries corresponding to allele1 at each locus must represent the same haplotype. Note: for each locus in each sample, at least one allele should have >0 supporting reads.
#'
#' @param lociAllele1CountsSample1,lociAllele2CountsSample1,lociAllele1CountsSample2,lociAllele2CountsSample2 vectors of counts of allele1 (e.g. reference) and allele2 (e.g. alternative) at individiual loci in sample1 and sample2. Allele counts are not necessarily phased (see argument 'isPhased'),  so allele1 counts may not represent the same haplotype. However, the two alleles (allele1 and allele2) must be defined identically for both samples at each locus. All 4 arguments must be vectors of non-negative integers.
#' @param lociAllele1NoASEProbsSample1,lociAllele1NoASEProbsSample2 probabilities of observing allele1-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus) in sample1 and sample2, respectively. Note that these probabilities are allowed to be sample-specific. Each argument must be a vector with entries >0 and <1.
#' @param lociRhosSample1,lociRhosSample2 dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial). Note that the dispersions are allowed to be sample-specific. Each argument must be a vector with entries >=0 and <1.  
#' @param numSim number of simulations to perform. Must be a non-negative integer. If 0 (DEFAULT), no simulations are performed.
#' @param isPhased single boolean specifying whether the phasing has already been performed, in which case the lociAllele1CountsSample1 (and, therefore, lociAllele1CountsSample2) represent the same haplotype. DEFAULT: FALSE.
#' @param tieBreakRandom single boolean specifying how ties should be broken during pseudo-phasing in cases of unphased data (isPhased=FALSE). If TRUE, each of the two allele will be assigned to major haplotype with probability=0.5. If FALSE (DEFAULT), allele1 will be assigned to major haplotype and allele2 to minor haplotype. 
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return list with 7 elements
#' \item{majorAlleleFrequencyDifference}{Estimate of major allele frequency difference for this unit of ASE (gene). 'Major' here refers to the allelic imbalance within sample1, and the difference is defined as Frequency(major, sample1)-Frequency(major, sample2).}
#' \item{pValueASE}{Estimate of p-value for observed extent of ASE (nominal if no simulations are performed, simulations-based otherwise).}
#' \item{heterogeneityQ}{Statistic summarizing variability of locus-specific estimates of major allele frequency difference if >1 locus is present. Set to NA for single-locus cases.}
#' \item{pValueHeterogeneity}{Estimate of p-value for observed extent of variability of locus-specific estimates of major allele frequency difference if >1 locus is present. Set to NA for single-locus cases.}
#' \item{lociAllele1IsMajor}{Vector of booleans, specifying for each locus whether allele1 is assigned to major (TRUE) or minor (FALSE) haplotype (where 'major' and 'minor' refer to abundances in sample1). If the data is phased (isPhased=TRUE), then all elements of the vector are TRUE if haplotype 1 is found to be major in sample1, and are all FALSE if haplotype 1 is found to be minor. In cases of unphased data (isPhased=FALSE), the assignment is provided by the pseudo-phasing procedure within sample1.}
#' \item{nullHypothesisMAF}{Estimate of major allele frequency under the null hypothesis that allelic frequencies are the same in both samples. This estimate is obtained by maximum likelihood, and, in case of unphased data (isPhased=FALSE), the likelihood is further maximized over all possible assignments of alleles to haplotypes.}
#' \item{lociMAFDifference}{Estimate of the difference of major allele (haplotype) frequency at individual loci. Note that 'major' and 'minor' distinction is made at the level of gene haplotype in sample1.}
#'
#' @examples
#' \donttest{
#' SNVCoverageTumor=sample(10:100, 5) ## gene with 5 loci
#' SNVCoverageNormal=sample(10:100, 5)
#' SNVAllele1CountsTumor=rbinom(length(SNVCoverageTumor), SNVCoverageTumor, 0.5)
#' SNVAllele1CountsNormal=rbinom(length(SNVCoverageNormal), SNVCoverageNormal, 0.5)
#' MBASED:::runMBASED2s1aseID(lociAllele1CountsSample1=SNVAllele1CountsTumor, lociAllele2CountsSample1=SNVCoverageTumor-SNVAllele1CountsTumor, lociAllele1CountsSample2=SNVAllele1CountsNormal, lociAllele2CountsSample2=SNVCoverageNormal-SNVAllele1CountsNormal, lociAllele1NoASEProbsSample1=rep(0.5, length(SNVCoverageTumor)), lociAllele1NoASEProbsSample2=rep(0.5, length(SNVCoverageNormal)), lociRhosSample1=rep(0, length(SNVCoverageTumor)), lociRhosSample2=rep(0, length(SNVCoverageNormal)), numSim=10^6, isPhased=FALSE) 
#' }
#'
runMBASED2s1aseID <- function (
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
            stop('MBASED:runMBASED2s1aseID: argument lociAllele1CountsSample1 must be a vector of non-negative integers')
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
            stop('MBASED:runMBASED2s1aseID: argument lociAllele1CountsSample2 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociAllele2CountsSample1) || 
            any(is.na(lociAllele2CountsSample1)) || 
            !is.numeric(lociAllele2CountsSample1) || 
            !isTRUE(
                all.equal(
                    lociAllele2CountsSample1, 
                    round(lociAllele2CountsSample1)
                )
            ) || 
            any(lociAllele2CountsSample1<0) 
        ) {
            stop('MBASED:runMBASED2s1aseID: argument lociAllele2CountsSample1 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociAllele2CountsSample2) || 
            any(is.na(lociAllele2CountsSample2)) || 
            !is.numeric(lociAllele2CountsSample2) || 
            !isTRUE(
                all.equal(
                    lociAllele2CountsSample2, 
                    round(lociAllele2CountsSample2)
                )
            ) || 
            any(lociAllele2CountsSample2<0) 
        ) {
            stop('MBASED:runMBASED2s1aseID: argument lociAllele2CountsSample2 must be a vector of non-negative integers')
        }
        if (
            !is.vector(lociAllele1NoASEProbsSample1) || 
            any(is.na(lociAllele1NoASEProbsSample1)) || 
            !is.numeric(lociAllele1NoASEProbsSample1) || 
            any(lociAllele1NoASEProbsSample1<=0) || 
            any(lociAllele1NoASEProbsSample1>=1)
        ) {
            stop('MBASED:runMBASED2s1aseID: argument lociAllele1NoASEProbsSample1 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociAllele1NoASEProbsSample2) || 
            any(is.na(lociAllele1NoASEProbsSample2)) || 
            !is.numeric(lociAllele1NoASEProbsSample2) || 
            any(lociAllele1NoASEProbsSample2<=0) || 
            any(lociAllele1NoASEProbsSample2>=1)
        ) {
            stop('MBASED:runMBASED2s1aseID: argument lociAllele1NoASEProbsSample2 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhosSample1) || 
            any(is.na(lociRhosSample1)) || 
            !is.numeric(lociRhosSample1) || 
            any(lociRhosSample1<0) || 
            any(lociRhosSample1>=1)
        ) {
            stop('MBASED:runMBASED2s1aseID: argument lociRhosSample1 must be a vector with entries >=0 and <1')
        }
        if (
            !is.vector(lociRhosSample2) || 
            any(is.na(lociRhosSample2)) || 
            !is.numeric(lociRhosSample2) || 
            any(lociRhosSample2<0) || 
            any(lociRhosSample2>=1)
        ) {
            stop('MBASED:runMBASED2s1aseID: argument lociRhosSample2 must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociAllele1CountsSample1), 
                length(lociAllele2CountsSample1), 
                length(lociAllele1NoASEProbsSample1), 
                length(lociRhosSample1), 
                length(lociAllele1CountsSample2), 
                length(lociAllele2CountsSample2), 
                length(lociAllele1NoASEProbsSample2), 
                length(lociRhosSample2)
            )))>1
        ) {
            stop('MBASED:runMBASED2s1aseID: arguments lociAllele1CountsSample1, lociAllele2CountsSample1, lociAllele1NoASEProbsSample1, lociRhosSample1, lociAllele1CountsSample2, lociAllele2CountsSample2, lociAllele1NoASEProbsSample2, and lociRhosSample2 must be of same length')
        }
        if (
            any((lociAllele1CountsSample1+lociAllele2CountsSample1)<1) || 
            any((lociAllele1CountsSample2+lociAllele2CountsSample2)<1)
        ) {
            stop('MBASED:runMBASED2s1aseID: for each locus, at least 1 allele should have >0 supporting reads in each sample')
        }
        if ( 
            !(is.vector(numSim)) || 
            !is.numeric(numSim) || 
            length(numSim)!=1 || 
            is.na(numSim) || 
            !isTRUE(all.equal(numSim, round(numSim))) || 
            (numSim<0) 
        ) {
            stop('MBASED:runMBASED2s1aseID: argument numSim must be a single non-negative integer')
        }
        if ( 
            !(is.vector(isPhased)) || 
            !is.logical(isPhased) || 
            length(isPhased)!=1
        ) {
            stop('MBASED:runMBASED2s1aseID: argument isPhased must be a single TRUE or FALSE value')
        }
        if ( 
            !(is.vector(tieBreakRandom)) || 
            !is.logical(tieBreakRandom) || 
            length(tieBreakRandom)!=1 
        ) {
            stop('MBASED:runMBASED2s1aseID: argument tieBreakRandom must be a single TRUE or FALSE value')
        }
    }
    ## get the total coverage at each SNV
    lociTotalAlleleCountsSample1 <- (
        lociAllele1CountsSample1+lociAllele2CountsSample1
    )
    lociTotalAlleleCountsSample2 <- (
        lociAllele1CountsSample2+lociAllele2CountsSample2
    )
    ## which of the alleles is major: based on sample1.
    lociAllele1IsMajor <- runMBASED1s1aseID(
        lociAllele1Counts=lociAllele1CountsSample1, 
        lociAllele2Counts=lociAllele2CountsSample1, 
        lociAllele1NoASEProbs=lociAllele1NoASEProbsSample1, 
        lociRhos=lociRhosSample1, 
        ## no need to do simulations: we only require information 
        ## about which allele is major
        numSim=0, 
        isPhased=isPhased,
        tieBreakRandom=tieBreakRandom
    )$lociAllele1IsMajor
    lociMajorAlleleCountsSample1 <- ifelse(
        lociAllele1IsMajor, 
        lociAllele1CountsSample1, 
        lociAllele2CountsSample1
    )
    lociMajorAlleleCountsSample2 <- ifelse(
        lociAllele1IsMajor, 
        lociAllele1CountsSample2, 
        lociAllele2CountsSample2
    )
    lociMajorAlleleNoASEProbsSample1 <- ifelse(
        lociAllele1IsMajor, 
        lociAllele1NoASEProbsSample1, 
        1-lociAllele1NoASEProbsSample1
    )
    lociMajorAlleleNoASEProbsSample2 <- ifelse(
        lociAllele1IsMajor, 
        lociAllele1NoASEProbsSample2, 
        1-lociAllele1NoASEProbsSample2
    )
    ## get results without simulations
    noSimResults <- MBASEDVectorizedPropDiffTest( 
        countsMatSample1=matrix(
            lociMajorAlleleCountsSample1, 
            ncol=1
        ), 
        totalsMatSample1=matrix(
            lociTotalAlleleCountsSample1, 
            ncol=1
        ), 
        countsMatSample2=matrix(
            lociMajorAlleleCountsSample2, 
            ncol=1
        ), 
        totalsMatSample2=matrix(
            lociTotalAlleleCountsSample2, 
            ncol=1
        ),
        probsMatSample1=matrix(
            lociMajorAlleleNoASEProbsSample1, 
            ncol=1
        ), 
        probsMatSample2=matrix(
            lociMajorAlleleNoASEProbsSample2, 
            ncol=1
        ),
        rhosMatSample1=matrix(
            lociRhosSample1, 
            ncol=1
        ),
        rhosMatSample2=matrix(
            lociRhosSample2, 
            ncol=1
        ),
        alternative='two.sided'
    )
    ## calculate MAF under null hypothesis
    MAF <- estimateMAF2s(
            lociAllele1CountsSample1=lociAllele1CountsSample1, 
            lociTotalCountsSample1=lociTotalAlleleCountsSample1, 
            lociAllele1CountsSample2=lociAllele1CountsSample2,
            lociTotalCountsSample2=lociTotalAlleleCountsSample2, 
            lociAllele1NoASEProbsSample1=lociAllele1NoASEProbsSample1, 
            lociAllele1NoASEProbsSample2=lociAllele1NoASEProbsSample2, 
            lociRhosSample1=lociRhosSample1, 
            lociRhosSample2=lociRhosSample2, 
            isPhased=isPhased
        )$MAF 
    ## initialize results to be returned to the user
    returnResults <- list(
        majorAlleleFrequencyDifference=as.vector(
            noSimResults$propDifferenceFinal
        ), 
        pValueASE=as.vector(noSimResults$pValue),
        heterogeneityQ=as.vector(noSimResults$hetQ),
        pValueHeterogeneity=as.vector(noSimResults$hetPVal),
        lociAllele1IsMajor=lociAllele1IsMajor,
        nullHypothesisMAF=MAF,
        lociMAFDifference=as.vector(noSimResults$propDifferenceLoci)
    )
    ######## estimate distribution of MAF estimates and Q statistics ####
    ######## under null model and get empirical p-values ######
    if (numSim>0) {        
        observedMAFDifference <- returnResults$majorAlleleFrequencyDifference
        observedQ <- returnResults$heterogeneityQ
        ## get vector of total allele counts, one full set for each simulation     
         lociTotalAlleleNullCountsSample1 <- rep(
             lociTotalAlleleCountsSample1, 
             numSim
         )
         lociTotalAlleleNullCountsSample2 <- rep(
             lociTotalAlleleCountsSample2, 
             numSim
         )
         ## get vector of allele1 probabilities (under conditions of no ASE), 
         ## one full set for each simulation     
         lociAllele1NullNoASEProbsSample1 <- rep(
             lociAllele1NoASEProbsSample1, 
             numSim
         )
         lociAllele1NullNoASEProbsSample2 <- rep(
             lociAllele1NoASEProbsSample2, 
             numSim
         )
         ## get vector of overdispersion parameters, 
         ## one full set for each simulation     
         lociNullRhosSample1 <- rep(lociRhosSample1, numSim)
         lociNullRhosSample2 <- rep(lociRhosSample2, numSim) 
         ## at each locus, one of Allele1 or Allele2 will be 
         ## randomly chosen to belong to haplotype1 
         ## ('major' haplotype, but in in the phasing sense of the word).
         if (isPhased) {
            lociAllele1IsMajorGenerating <- rep(lociAllele1IsMajor, numSim)
        } else {
            lociAllele1IsMajorGenerating <- sample(
                c(TRUE,FALSE), 
                size=length(lociTotalAlleleCountsSample1)*numSim, 
                replace=TRUE
            )
        }
        ## adjust these vectors of probabilities to reflect the estimated MAF.
        lociAllele1NullGeneratingProbsSample1 <- ifelse(
             lociAllele1IsMajorGenerating, 
             getPFinal(
                 trueAF=MAF, 
                 noASEAF=lociAllele1NullNoASEProbsSample1
             ),
             getPFinal(
                 trueAF=1-MAF, 
                 noASEAF=lociAllele1NullNoASEProbsSample1
             )
         )
         lociAllele1NullGeneratingProbsSample2 <- ifelse(
             lociAllele1IsMajorGenerating, 
             getPFinal(
                 trueAF=MAF, 
                 noASEAF=lociAllele1NullNoASEProbsSample2
             ), 
             getPFinal(
                 trueAF=1-MAF, 
                 noASEAF=lociAllele1NullNoASEProbsSample2
             )
         )
        ## get allele1 counts based on null distribution
        lociAllele1NullCountsSample1 <- vectorizedRbetabinomMR(
            n=length(lociAllele1CountsSample1)*numSim, 
            size=lociTotalAlleleNullCountsSample1, 
            mu=lociAllele1NullGeneratingProbsSample1,
            rho=lociNullRhosSample1
        )
        lociAllele1NullCountsSample2 <- vectorizedRbetabinomMR(
            n=length(lociAllele1CountsSample2)*numSim, 
            size=lociTotalAlleleNullCountsSample2, 
            mu=lociAllele1NullGeneratingProbsSample2,
            rho=lociNullRhosSample2
        )         
        ## 'phase' the allele counts if needed; 
        ## otherwise initialize 'major' haplotype to be haplotype 1
        if (!isPhased) { 
            lociAllele1IsMajorNull <- runMBASED1s1aseID(
                 lociAllele1Counts=lociAllele1NullCountsSample1, 
                lociAllele2Counts=(
                    lociTotalAlleleNullCountsSample1-lociAllele1NullCountsSample1
                ), 
                lociAllele1NoASEProbs=lociAllele1NullNoASEProbsSample1, 
                lociRhos=lociNullRhosSample1, 
                ## no need to do simulations: we only require information 
                ## about which allele is major
                numSim=0, 
                isPhased=isPhased,
                tieBreakRandom=tieBreakRandom
            )$lociAllele1IsMajor
        } else {
            lociAllele1IsMajorNull <- lociAllele1IsMajorGenerating
        }
        lociMajorAlleleNullCountsSample1 <- ifelse(
            lociAllele1IsMajorNull, 
            lociAllele1NullCountsSample1, 
            lociTotalAlleleNullCountsSample1-lociAllele1NullCountsSample1
        )
        lociMajorAlleleNullCountsSample2 <- ifelse(
            lociAllele1IsMajorNull, 
            lociAllele1NullCountsSample2, 
            lociTotalAlleleNullCountsSample2-lociAllele1NullCountsSample2
        )
        lociMajorAlleleNullNoASEProbsSample1<- ifelse(
            lociAllele1IsMajorNull, 
            lociAllele1NullNoASEProbsSample1, 
            1-lociAllele1NullNoASEProbsSample1
        )
        lociMajorAlleleNullNoASEProbsSample2 <- ifelse(
            lociAllele1IsMajorNull, 
            lociAllele1NullNoASEProbsSample2, 
            1-lociAllele1NullNoASEProbsSample2
        )
        ##get the null-distribution of nominal p-values and obtain the ASE p-value 
        simNullResults <- MBASEDVectorizedPropDiffTest(
            countsMatSample1=matrix(
                lociMajorAlleleNullCountsSample1, 
                ncol=numSim
            ), 
            totalsMatSample1=matrix(
                lociTotalAlleleNullCountsSample1, 
                ncol=numSim
            ), 
            countsMatSample2=matrix(
                lociMajorAlleleNullCountsSample2, 
                ncol=numSim
            ), 
            totalsMatSample2=matrix(
                lociTotalAlleleNullCountsSample2, 
                ncol=numSim
            ),
            probsMatSample1=matrix(
                lociMajorAlleleNullNoASEProbsSample1, 
                ncol=numSim
            ), 
            probsMatSample2=matrix(
                lociMajorAlleleNullNoASEProbsSample2, 
                ncol=numSim
            ),
            rhosMatSample1=matrix(
                lociNullRhosSample1, 
                ncol=numSim
            ),
            rhosMatSample2=matrix(
                lociNullRhosSample2, 
                ncol=numSim
            ),
            alternative='two.sided'
        )
        simNullMAFDifferences <- as.vector(simNullResults$propDifferenceFinal)
        simNullQs <- as.vector(simNullResults$hetQ)
        returnResults$pValueASE <- getSimulationPvalue(
            observedVal=observedMAFDifference, 
            simulatedVals=simNullMAFDifferences, 
            direction='greater'
        )
        if (!is.na(observedQ)) {
            returnResults$pValueHeterogeneity <- getSimulationPvalue(
                observedVal=observedQ, 
                simulatedVals=simNullQs, 
                direction='greater'
            )
        }
       }
    return(returnResults)  
 }

#' Function that runs between-sample (differential) ASE calling using data from individual loci (SNVs) within units of ASE (genes). Vector arguments 'lociAllele1CountsSample1', 'lociAllele2CountsSample1', 'lociAllele1NoASEProbsSample1', 'lociRhosSample1', 'lociAllele1CountsSample2', 'lociAllele2CountsSample2', 'lociAllele1NoASEProbsSample2', 'lociRhosSample2', and 'aseIDs' should all be of the same length. Letting i1, i2, .., iN denote the indices corresponding to entries within aseIDs equal to a given aseID, the entries at those indices in the other vector arguments provide information for the loci within that aseID for the respective samples. This information is then used by runMBASED2s1aseID. It is assumed that for any i, the i-th entries of all vector arguments correspond to the same locus, and that the entries corresponding to allele1 in sample1 and sample2 provide information on the same allele. If argument 'isPhased' (see below) is true, then entries corresponding to allele1 at each locus must represent the same haplotype.
#'
#' @param lociAllele1CountsSample1,lociAllele2CountsSample1,lociAllele1CountsSample2,lociAllele2CountsSample2 vectors of counts of allele1 (e.g. reference) and allele2 (e.g. alternative) at individiual loci in sample1 and sample2. Allele counts are not necessarily phased (see argument 'isPhased'),  so allele1 counts may not represent the same haplotype. However, the two alleles (allele1 and allele2) must be defined identically for both samples at each locus. All 4 arguments must be vectors of non-negative integers.
#' @param lociAllele1NoASEProbsSample1,lociAllele1NoASEProbsSample2 probabilities of observing allele1-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus) in sample1 and sample2, respectively. Note that these probabilities are allowed to be sample-specific. Each argument must be a vector with entries >0 and <1.
#' @param lociRhosSample1,lociRhosSample2 dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial). Note that the dispersions are allowed to be sample-specific. Each argument must be a vector with entries >=0 and <1.
#' @param aseIDs the IDs of ASE units corresponding to the individual loci (e.g. gene names).  
#' @param numSim number of simulations to perform. Must be a non-negative integer. If 0 (DEFAULT), no simulations are performed.
#' @param BPPARAM argument to be passed to bplapply(), when parallel achitecture is used to speed up simulations (parallelization is done over aseIDs).  DEFAULT: SerialParam() (no parallelization).
#' @param isPhased single boolean specifying whether the phasing has already been performed, in which case the lociAllele1CountsSample1 (and, therefore, lociAllele1CountsSample2) represent the same haplotype. DEFAULT: FALSE.
#' @param tieBreakRandom single boolean specifying how ties should be broken during pseudo-phasing in cases of unphased data (isPhased=FALSE). If TRUE, each of the two allele will be assigned to major haplotype with probability=0.5. If FALSE (DEFAULT), allele1 will be assigned to major haplotype and allele2 to minor haplotype. 
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return list with 3 elements:
#' \item{ASEResults}{Data frame with each row reporting MBASED results for a given aseID (aseIDs are provided as row names of this data frame). The columns of the data frame are: majorAlleleFrequencyDifference, pValueASE, heterogeneityQ, and pValueHeterogeneity.}
#' \item{allele1IsMajor}{Vector of TRUE/FALSE of length equal to the number of supplied SNVs, reporting for each SNV whether allele1 represents major (TRUE) or minor (FALSE) haplotype of the corresponding aseID.}
#' \item{lociMAFDifference}{Vector of locus-specific estimates of the difference of major allele (haplotype) frequency between the two samples. Note that 'major' and 'minor' distinction is made at the level of gene haplotype in sample1.}
#'
#' @examples
#' \donttest{
#' SNVCoverageTumor=sample(10:100, 5)
#' SNVCoverageNormal=sample(10:100, 5)
#' SNVAllele1CountsTumor=rbinom(length(SNVCoverageTumor), SNVCoverageTumor, 0.5)
#' SNVAllele1CountsNormal=rbinom(length(SNVCoverageNormal), SNVCoverageNormal, 0.5)
#' MBASED:::runMBASED2s(lociAllele1CountsSample1=SNVAllele1CountsTumor, lociAllele2CountsSample1=SNVCoverageTumor-SNVAllele1CountsTumor, lociAllele1CountsSample2=SNVAllele1CountsNormal, lociAllele2CountsSample2=SNVCoverageNormal-SNVAllele1CountsNormal, lociAllele1NoASEProbsSample1=rep(0.5, length(SNVCoverageTumor)), lociAllele1NoASEProbsSample2=rep(0.5, length(SNVCoverageNormal)), lociRhosSample1=rep(0, length(SNVCoverageTumor)), lociRhosSample2=rep(0, length(SNVCoverageNormal)), aseIDs=c(rep('gene1',4), 'gene2'), numSim=10^6,  BPPARAM=SerialParam(), isPhased=FALSE) 
#' }
#'
runMBASED2s <- function (
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
            stop('MBASED:runMBASED2s: argument lociAllele1CountsSample1 must be a vector of non-negative integers')
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
            stop('MBASED:runMBASED2s: argument lociAllele1CountsSample2 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociAllele2CountsSample1) || 
            any(is.na(lociAllele2CountsSample1)) || 
            !is.numeric(lociAllele2CountsSample1) || 
            !isTRUE(
                all.equal(
                    lociAllele2CountsSample1, 
                    round(lociAllele2CountsSample1)
                )
            ) || 
            any(lociAllele2CountsSample1<0) 
        ) {
            stop('MBASED:runMBASED2s: argument lociAllele2CountsSample1 must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociAllele2CountsSample2) || 
            any(is.na(lociAllele2CountsSample2)) || 
            !is.numeric(lociAllele2CountsSample2) || 
            !isTRUE(
                all.equal(
                    lociAllele2CountsSample2, 
                    round(lociAllele2CountsSample2)
                )
            ) || 
            any(lociAllele2CountsSample2<0) 
        ) {
            stop('MBASED:runMBASED2s: argument lociAllele2CountsSample2 must be a vector of non-negative integers')
        }
        if (
            !is.vector(lociAllele1NoASEProbsSample1) || 
            any(is.na(lociAllele1NoASEProbsSample1)) || 
            !is.numeric(lociAllele1NoASEProbsSample1) || 
            any(lociAllele1NoASEProbsSample1<=0) || 
            any(lociAllele1NoASEProbsSample1>=1)
        ) {
            stop('MBASED:runMBASED2s: argument lociAllele1NoASEProbsSample1 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociAllele1NoASEProbsSample2) || 
            any(is.na(lociAllele1NoASEProbsSample2)) || 
            !is.numeric(lociAllele1NoASEProbsSample2) || 
            any(lociAllele1NoASEProbsSample2<=0) || 
            any(lociAllele1NoASEProbsSample2>=1)
        ) {
            stop('MBASED:runMBASED2s: argument lociAllele1NoASEProbsSample2 must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhosSample1) || 
            any(is.na(lociRhosSample1)) || 
            !is.numeric(lociRhosSample1) || 
            any(lociRhosSample1<0) || 
            any(lociRhosSample1>=1)
        ) {
            stop('MBASED:runMBASED2s: argument lociRhosSample1 must be a vector with entries >=0 and <1')
        }
        if (
            !is.vector(lociRhosSample2) || 
            any(is.na(lociRhosSample2)) || 
            !is.numeric(lociRhosSample2) || 
            any(lociRhosSample2<0) || 
            any(lociRhosSample2>=1)
        ) {
            stop('MBASED:runMBASED2s: argument lociRhosSample2 must be a vector with entries >=0 and <1')
        }
        if (
            !is.vector(aseIDs) || 
            any(is.na(aseIDs)) 
        ) {
            stop('MBASED:runMBASED2s: argument aseIDs must be a vector with no NAs')
        }
        if (
            length(unique(c(
                length(lociAllele1CountsSample1), 
                length(lociAllele2CountsSample1), 
                length(lociAllele1NoASEProbsSample1), 
                length(lociRhosSample1), 
                length(lociAllele1CountsSample2), 
                length(lociAllele2CountsSample2), 
                length(lociAllele1NoASEProbsSample2), 
                length(lociRhosSample2), 
                length(aseIDs)
            )))>1
        ) {
            stop('MBASED:runMBASED2s: arguments lociAllele1CountsSample1, lociAllele2CountsSample1, lociAllele1NoASEProbsSample1, lociRhosSample1, lociAllele1CountsSample2, lociAllele2CountsSample2, lociAllele1NoASEProbsSample2, lociRhosSample2, and aseIDs must be of same length')
        }
        if (
            any((lociAllele1CountsSample1+lociAllele2CountsSample1)<1) || 
            any((lociAllele1CountsSample2+lociAllele2CountsSample2)<1)
        ) {
            stop('MBASED:runMBASED2s: for each locus, at least 1 allele should have >0 supporting reads in each sample')
        }
        if ( 
            !(is.vector(numSim)) || 
            !is.numeric(numSim) || 
            length(numSim)!=1 || 
            is.na(numSim) || 
            !isTRUE(all.equal(numSim, round(numSim))) || 
            (numSim<0) 
        ) {
            stop('MBASED:runMBASED2s: argument numSim must be a single non-negative integer')
        }
        if (!(is(BPPARAM, 'BiocParallelParam'))) {
            stop('MBASED:runMBASED2s: argument BPPARAM must be an instance of BiocParallelParam')
        }
        if ( 
            !(is.vector(isPhased)) || 
            !is.logical(isPhased) || 
            length(isPhased)!=1
        ) {
            stop('MBASED:runMBASED2s: argument isPhased must be a single TRUE or FALSE value')
        }
        if ( 
            !(is.vector(tieBreakRandom)) || 
            !is.logical(tieBreakRandom) || 
            length(tieBreakRandom)!=1 
        ) {
            stop('MBASED:runMBASED2s: argument tieBreakRandom must be a single TRUE or FALSE value')
        }
    }
    ASEResults <- bplapply(unique(aseIDs), function(aseID) {
        aseIDsubv <- (aseIDs==aseID)
        res <- runMBASED2s1aseID(
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
                tieBreakRandom=tieBreakRandom
        ) 
        return(res)
    }, BPPARAM=BPPARAM)
    names(ASEResults) <- unique(aseIDs)
    allele1IsMajor <- unname(unlist(lapply(ASEResults, function(el) {
        el$lociAllele1IsMajor
    })))
    lociMAFDifference <- unname(unlist(lapply(ASEResults, function(el) {
        el$lociMAFDifference
    })))
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



