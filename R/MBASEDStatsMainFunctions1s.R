#' Function that runs single-sample ASE calling using data from loci (SNVs) within a single unit of ASE (gene). The i-th entry of each of vector arguments 'lociAllele1Counts', 'lociAllele2Counts', 'lociAllele1NoASEProbs', 'lociRhos' should correspond to the i-th locus. If argument 'isPhased' (see below) is true, then entries corresponding to allele1 at each locus must represent the same haplotype. Note: for each locus, at least one allele should have >0 supporting reads.
#'
#' @param lociAllele1Counts,lociAllele2Counts vectors of counts of allele1 (e.g. reference) and allele2 (e.g., alternative) at individual loci. Allele counts are not necessarily phased (see argument 'isPhased'), so allele1 counts may not represent the same haplotype. Both arguments must be vectors of non-negative integers.
#' @param lociAllele1NoASEProbs probabilities of observing allele1-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus). Must be a vector with entries >0 and <1.
#' @param lociRhos dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial). Must be a vector with entries >=0 and <1.
#' @param numSim number of simulations to perform. Must be a non-negative integer. If 0 (DEFAULT), no simulations are performed.
#' @param isPhased single boolean specifying whether the phasing has already been performed, in which case the lociAllele1Counts represent the same haplotype. DEFAULT: FALSE.
#' @param tieBreakRandom single boolean specifying how ties should be broken during pseudo-phasing in cases of unphased data (isPhased=FALSE). If TRUE, each of the two allele will be assigned to major haplotype with probability=0.5. If FALSE (DEFAULT), allele1 will be assigned to major haplotype and allele2 to minor haplotype. 
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return list with 6 elements
#' \item{majorAlleleFrequency}{Estimate of major allele frequency for this unit of ASE (gene).}
#' \item{pValueASE}{Estimate of p-value for observed extent of ASE (nominal if no simulations are performed, simulations-based otherwise).}
#' \item{heterogeneityQ}{Statistic summarizing variability of locus-specific estimates of major allele frequency if >1 locus is present. Set to NA for single-locus cases.}
#' \item{pValueHeterogeneity}{Estimate of p-value for observed extent of variability of locus-specific estimates of major allele frequency if >1 locus is present. Set to NA for single-locus cases.}
#' \item{lociAllele1IsMajor}{Vector of booleans, specifying for each locus whether allele1 is assigned to major (TRUE) or minor (FALSE) haplotype. If the data is phased (isPhased=TRUE), then all elements of the vector are TRUE if haplotype 1 is found to be major, and are all FALSE if haplotype 1 is found to be minor. In cases of unphased data (isPhased=FALSE), the assignment is provided by the pseudo-phasing procedure.}
#' \item{lociMAF}{Estimate of major allele (haplotype) frequency at individual loci. Note that since 'major' and 'minor' distinction is made at the level of gene haplotype, there may be some loci where the frequency of the 'major' haplotype is <0.5.}
#'
#' @examples
#' \donttest{
#' SNVCoverage <- sample(10:100,5) ## gene with 5 loci
#' SNVAllele1Counts <- rbinom(length(SNVCoverage), SNVCoverage, 0.5)
#' MBASED:::runMBASED1s1aseID(lociAllele1Counts=SNVAllele1Counts, lociAllele2Counts=SNVCoverage-SNVAllele1Counts, lociAllele1NoASEProbs=rep(0.5, length(SNVCoverage)), lociRhos=rep(0, length(SNVCoverage)), numSim=0, isPhased=FALSE, tieBreakRandom=FALSE) ## data is not phased, no simulations
#' MBASED:::runMBASED1s1aseID(lociAllele1Counts=SNVAllele1Counts, lociAllele2Counts=SNVCoverage-SNVAllele1Counts, lociAllele1NoASEProbs=rep(0.5, length(SNVCoverage)), lociRhos=rep(0, length(SNVCoverage)), numSim=10^6, isPhased=FALSE, tieBreakRandom=FALSE) ## data is not phased, simulations
#' MBASED:::runMBASED1s1aseID(lociAllele1Counts=SNVAllele1Counts, lociAllele2Counts=SNVCoverage-SNVAllele1Counts, lociAllele1NoASEProbs=rep(0.5, length(SNVCoverage)), lociRhos=rep(0, length(SNVCoverage)),  numSim=0, isPhased=TRUE, tieBreakRandom=FALSE) ## data is phased, no simulations
#' MBASED:::runMBASED1s1aseID(lociAllele1Counts=SNVAllele1Counts, lociAllele2Counts=SNVCoverage-SNVAllele1Counts, lociAllele1NoASEProbs=rep(0.5, length(SNVCoverage)), lociRhos=rep(0, length(SNVCoverage)),  numSim=10^6, isPhased=TRUE, tieBreakRandom=FALSE) ## data is phased, simulations
#' }
#'
runMBASED1s1aseID <- function (
    lociAllele1Counts, 
    lociAllele2Counts, 
    lociAllele1NoASEProbs, 
    lociRhos, 
    numSim=0, 
    isPhased=FALSE, 
    tieBreakRandom=FALSE, 
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
            stop('MBASED:runMBASED1s1aseID: argument lociAllele1Counts must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociAllele2Counts) || 
            any(is.na(lociAllele2Counts)) || 
            !is.numeric(lociAllele2Counts) || 
            !isTRUE(all.equal(lociAllele2Counts, round(lociAllele2Counts))) || 
            any(lociAllele2Counts<0) 
        ) {
            stop('MBASED:runMBASED1s1aseID: argument lociAllele2Counts must be a vector of non-negative integers')
        }
        if (
            !is.vector(lociAllele1NoASEProbs) || 
            any(is.na(lociAllele1NoASEProbs)) || 
            !is.numeric(lociAllele1NoASEProbs) || 
            any(lociAllele1NoASEProbs<=0) || 
            any(lociAllele1NoASEProbs>=1)
        ) {
            stop('MBASED:runMBASED1s1aseID: argument lociAllele1NoASEProbs must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhos) || 
            any(is.na(lociRhos)) || 
            !is.numeric(lociRhos) || 
            any(lociRhos<0) || 
            any(lociRhos>=1)
        ) {
            stop('MBASED:runMBASED1s1aseID: argument lociRhos must be a vector with entries >=0 and <1')
        }
        if (
            length(unique(c(
                length(lociAllele1Counts), 
                length(lociAllele2Counts), 
                length(lociAllele1NoASEProbs), 
                length(lociRhos)
            )))>1
        ) {
            stop('MBASED:runMBASED1s1aseID: arguments lociAllele1Counts, lociAllele2Counts, lociAllele1NoASEProbs, and lociRhos must be of same length')
        }
        if (any((lociAllele1Counts+lociAllele2Counts)<1)) {
            stop('MBASED:runMBASED1s1aseID: for each locus, at least 1 allele should have >0 supporting reads')
        }
        if ( 
            !(is.vector(numSim)) || 
            !is.numeric(numSim) || 
            length(numSim)!=1 || 
            is.na(numSim) || 
            !isTRUE(all.equal(numSim, round(numSim))) || 
            (numSim<0) 
        ) {
            stop('MBASED:runMBASED1s1aseID: argument numSim must be a single non-negative integer')
        }
        if ( 
            !(is.vector(isPhased)) || 
            !is.logical(isPhased) || 
            length(isPhased)!=1
        ) {
            stop('MBASED:runMBASED1s1aseID: argument isPhased must be a single TRUE or FALSE value')
        }
        if ( 
            !(is.vector(tieBreakRandom)) || 
            !is.logical(tieBreakRandom) || 
            length(tieBreakRandom)!=1 
        ) {
            stop('MBASED:runMBASED1s1aseID: argument tieBreakRandom must be a single TRUE or FALSE value')
        }
    }
    ## get the total coverage at each SNV
    lociTotalAlleleCounts <- (lociAllele1Counts+lociAllele2Counts)
    ## if alleles are not phased, assign alleles to major haplotype, 
    ## based on observed read counts (with appropriate adjustments) 
    ## if alelles ARE phased, initialize major haplotype to be based 
    ## on observed read counts at first SNV  and then 
    ## use estimate of MAF to re-adjust, if the other haplotype 
    ## turns out to correspond to higher allele frequency.
    if (!isPhased) { 
        lociAllele1IsMajor <- isCountMajorFT(
            x=lociAllele1Counts, 
            n=lociTotalAlleleCounts, 
            p=lociAllele1NoASEProbs, 
            tieBreakRandom=tieBreakRandom
        )
    } else {
        lociAllele1IsMajor <- rep(
            isCountMajorFT(
                x=lociAllele1Counts[1], 
                n=lociTotalAlleleCounts[1], 
                p=lociAllele1NoASEProbs[1], 
                tieBreakRandom=tieBreakRandom
            ), length(lociAllele1Counts)
        )
    }
    lociMajorAlleleCounts <- ifelse(
        lociAllele1IsMajor, 
        lociAllele1Counts, 
        lociAllele2Counts
    )
    lociMajorAlleleNoASEProbs <- ifelse(
        lociAllele1IsMajor, 
        lociAllele1NoASEProbs, 
        1-lociAllele1NoASEProbs
    )
    ## get the results without simulations
    noSimResults <- MBASEDVectorizedMetaprop(
        countsMat=matrix(
            lociMajorAlleleCounts,
            ncol=1
        ), 
        totalsMat=matrix(
            lociTotalAlleleCounts,
            ncol=1
        ), 
        probsMat=matrix(
            lociMajorAlleleNoASEProbs,
            ncol=1
        ),
        rhosMat=matrix(
            lociRhos, 
            ncol=1
        ),
        alternative='two.sided'
    )
    ## initialize results to be returned to the user
    returnResults <- list(
        majorAlleleFrequency=as.vector(noSimResults$propFinal), 
        pValueASE=as.vector(noSimResults$pValue),
        heterogeneityQ=as.vector(noSimResults$hetQ),
        pValueHeterogeneity=as.vector(noSimResults$hetPVal),
        lociAllele1IsMajor=lociAllele1IsMajor,
        lociMAF=as.vector(noSimResults$propLoci)
    )
    ## adjust assignment of major/minor haplotypes, 
    ## in case the data was already phased and 
    ## we used haplotype1 as 'major'
       if (isPhased) {
           if (returnResults$majorAlleleFrequency<0.5) {
               returnResults$majorAlleleFrequency <- 
                   1-returnResults$majorAlleleFrequency
               returnResults$lociAllele1IsMajor <- !(returnResults$lociAllele1IsMajor)
               returnResults$lociMAF <- 1-returnResults$lociMAF
              }
       }
######## estimate distribution of MAF estimates and Q statistics #######
######## under null model and get empirical p-values ############
    if (numSim>0) {   
        observedMAF <-  returnResults$majorAlleleFrequency
           observedQ <-  returnResults$heterogeneityQ
         ## get vector of total allele counts, one full set for each simulation     
         lociTotalAlleleNullCounts <- rep(lociTotalAlleleCounts, numSim)
         ## get vector of allele1 probabilities (under conditions of no ASE), 
         ## one full set for each simulation     
         lociAllele1NullNoASEProbs <- rep(lociAllele1NoASEProbs, numSim)
         ## get vector of overdispersion parameters, 
         ## one full set for each simulation 
         lociNullRhos <- rep(lociRhos, numSim)
           ## get allele1 counts based on null distribution
        lociAllele1NullCounts <- vectorizedRbetabinomMR(
            n=length(lociAllele1Counts)*numSim, 
            size=lociTotalAlleleNullCounts, 
            mu=lociAllele1NullNoASEProbs,
            rho=lociNullRhos
        )
        ## pseudo-phase the allele counts, if needed; 
        ## otherwise initialize 'major' haplotype to be haplotype 1 
        if (!isPhased) {
            lociAllele1IsMajorNull <- isCountMajorFT(
                x=lociAllele1NullCounts, 
                n=lociTotalAlleleNullCounts, 
                p=lociAllele1NullNoASEProbs, 
                tieBreakRandom=tieBreakRandom
            )
        } else {
            lociAllele1IsMajorNull <- rep(TRUE, length(lociAllele1NullCounts))
        }
        lociMajorAlleleNullCounts <- ifelse(
            lociAllele1IsMajorNull, 
            lociAllele1NullCounts, 
            lociTotalAlleleNullCounts-lociAllele1NullCounts
        )
        lociMajorAlleleNullNoASEProbs <- ifelse(
            lociAllele1IsMajorNull, 
            lociAllele1NullNoASEProbs, 
            1-lociAllele1NullNoASEProbs
        )
        ##get the null-distribution of MAF estimates and 
        ## the heterogeneity Q statistics and obtain the corresponding p-values 
        simNullResults <- MBASEDVectorizedMetaprop(
            countsMat=matrix(
                lociMajorAlleleNullCounts,
                ncol=numSim
            ), 
            totalsMat=matrix(
                lociTotalAlleleNullCounts,
                ncol=numSim
            ),
            probsMat=matrix(
                lociMajorAlleleNullNoASEProbs, 
                ncol=numSim
            ),
            rhosMat=matrix(
                lociNullRhos, 
                ncol=numSim
            ),
            alternative='two.sided'
        )
        simNullMAFs <- as.vector(simNullResults$propFinal)
        simNullQs <- as.vector(simNullResults$hetQ)
        if (isPhased) { 
        	## if phased, some of the null MAFs will be <0.5. 
        	## We're doing two-sided test, so need to make them all >0.5
            simNullMAFs <- pmax(simNullMAFs, 1-simNullMAFs)
        }
        returnResults$pValueASE <- getSimulationPvalue(
            observedVal=observedMAF, 
            simulatedVals=simNullMAFs, 
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

#' Function that runs single-sample ASE calling using data from individual loci (SNVs) within units of ASE (genes). Vector arguments 'lociAllele1Counts', 'lociAllele2Counts', 'lociAllele1NoASEProbs', 'lociRhos', and 'aseIDs' should all be of the same length. Letting i1, i2, .., iN denote the indices corresponding to entries within aseIDs equal to a given aseID, the entries at those indices in the other vector arguments provide information for the loci within that aseID. This information is then used by runMBASED1s1aseID. It is assumed that for any i, the i-th entries of all vector arguments correspond to the same locus. If argument 'isPhased' (see below) is true, then entries corresponding to allele1 at each locus must represent the same haplotype.
#'
#' @param lociAllele1Counts,lociAllele2Counts vectors of counts of allele1 (e.g. reference) and allele2 (e.g., alternative) at individual loci. Allele counts are not necessarily phased (see argument 'isPhased'), so allele1 counts may not represent the same haplotype. Both arguments must be vectors of non-negative integers.
#' @param lociAllele1NoASEProbs probabilities of observing allele1-supporting reads at individual loci under conditions of no ASE (e.g., vector with all entries set to 0.5, if there is no pre-existing allelic bias at any locus). Must be a vector with entries >0 and <1.
#' @param lociRhos dispersion parameters of beta distribution at individual loci (set to 0 if the read count-generating distribution at the locus is binomial).  Must be a vector with entries >=0 and <1.
#' @param aseIDs the IDs of ASE units corresponding to the individual loci (e.g. gene names).
#' @param numSim number of simulations to perform.  Must be a non-negative integer. If 0 (DEFAULT), no simulations are performed.
#' @param BPPARAM argument to be passed to bplapply(), when parallel achitecture is used to speed up simulations (parallelization is done over aseIDs).  DEFAULT: SerialParam() (no parallelization).
#' @param isPhased single boolean specifying whether the phasing has already been performed, in which case the lociAllele1Counts represent the same haplotype. DEFAULT: FALSE.
#' @param tieBreakRandom single boolean specifying how ties should be broken during pseudo-phasing in cases of unphased data (isPhased=FALSE). If TRUE, each of the two allele will be assigned to major haplotype with probability=0.5. If FALSE (DEFAULT), allele1 will be assigned to major haplotype and allele2 to minor haplotype. 
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return list with 3 elements:
#' \item{ASEResults}{Data frame with each row reporting MBASED results for a given aseID (aseIDs are provided as row names of this data frame). The columns of the data frame are: majorAlleleFrequency, pValueASE, heterogeneityQ, and pValueHeterogeneity.}
#' \item{allele1IsMajor}{Vector of TRUE/FALSE of length equal to the number of supplied SNVs, reporting for each SNV whether allele1 represents major (TRUE) or minor (FALSE) haplotype of the corresponding aseID.}
#' \item{lociMAF}{Vector of locus-specific estimates of the frequency of major allele, where 'major' refers to the haplotype of the gene found to be major by the ASE analysis. Note that since the determination of the major/minor status is done at the level of the gene, there may be loci with locus-specific MAF < 0.5.}
#' 
#' @examples
#' \donttest{
#' SNVCoverage1 <- sample(10:100,5) ## gene with 5 loci
#' SNVAllele1Counts1 <- rbinom(length(SNVCoverage1), SNVCoverage1, 0.5)
#' SNVCoverage2 <- sample(10:100,5) ## gene with 5 loci
#' SNVAllele1Counts2 <- rbinom(length(SNVCoverage2), SNVCoverage2, 0.5)
#' MBASED:::runMBASED1s(lociAllele1Counts=c(SNVAllele1Counts1, SNVAllele1Counts2), lociAllele2Counts=c(SNVCoverage1-SNVAllele1Counts1, SNVCoverage2-SNVAllele1Counts2), lociAllele1NoASEProbs=rep(0.5, 10), lociRhos=rep(0, 10), aseIDs=rep(c('gene1','gene2'), each=5), numSim=10^6, BPPARAM=SerialParam(), isPhased=FALSE, tieBreakRandom=FALSE) 
#' }
#'
runMBASED1s <- function (
    lociAllele1Counts, 
    lociAllele2Counts, 
    lociAllele1NoASEProbs, 
    lociRhos, 
    aseIDs, 
    numSim=0, 
    BPPARAM=SerialParam(), 
    isPhased=FALSE, 
    tieBreakRandom=FALSE, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(lociAllele1Counts) || 
            any(is.na(lociAllele1Counts)) || 
            !is.numeric(lociAllele1Counts) || 
            !isTRUE(all.equal(lociAllele1Counts, round(lociAllele1Counts))) 
            || any(lociAllele1Counts<0) 
        ) {
            stop('MBASED:runMBASED1s: argument lociAllele1Counts must be a vector of non-negative integers')
        }
        if ( 
            !is.vector(lociAllele2Counts) || 
            any(is.na(lociAllele2Counts)) || 
            !is.numeric(lociAllele2Counts) || 
            !isTRUE(all.equal(lociAllele2Counts, round(lociAllele2Counts))) || 
            any(lociAllele2Counts<0) 
        ) {
            stop('MBASED:runMBASED1s: argument lociAllele2Counts must be a vector of non-negative integers')
        }
        if (
            !is.vector(lociAllele1NoASEProbs) || 
            any(is.na(lociAllele1NoASEProbs)) || 
            !is.numeric(lociAllele1NoASEProbs) || 
            any(lociAllele1NoASEProbs<=0) || 
            any(lociAllele1NoASEProbs>=1)
        ) {
            stop('MBASED:runMBASED1s: argument lociAllele1NoASEProbs must be a vector with entries >0 and <1')
        }
        if (
            !is.vector(lociRhos) || 
            any(is.na(lociRhos)) || 
            !is.numeric(lociRhos) || 
            any(lociRhos<0) || 
            any(lociRhos>=1)
        ) {
            stop('MBASED:runMBASED1s: argument lociRhos must be a vector with entries >=0 and <1')
        }
        if (
            !is.vector(aseIDs) || 
            any(is.na(aseIDs)) 
        ) {
            stop('MBASED:runMBASED1s: argument aseIDs must be a vector with no NAs')
        }
        if (
            length(unique(c(
                length(lociAllele1Counts), 
                length(lociAllele2Counts), 
                length(lociAllele1NoASEProbs), 
                length(lociRhos), 
                length(aseIDs)
            )))>1
        ) {
            stop('MBASED:runMBASED1s: arguments lociAllele1Counts, lociAllele2Counts, lociAllele1NoASEProbs, lociRhos, and aseIDs must be of same length')
        }
        if (any((lociAllele1Counts+lociAllele2Counts)<1)) {
            stop('MBASED:runMBASED1s: for each locus, at least 1 allele should have >0 supporting reads')
        }
        if ( 
            !(is.vector(numSim)) || 
            !is.numeric(numSim) || 
            length(numSim)!=1 || 
            is.na(numSim) || 
            !isTRUE(all.equal(numSim, round(numSim))) || 
            (numSim<0) 
        ) {
            stop('MBASED:runMBASED1s: argument numSim must be a single non-negative integer')
        }
        if (!(is(BPPARAM, 'BiocParallelParam'))) {
            stop('MBASED:runMBASED1s: argument BPPARAM must be an instance of BiocParallelParam')
        }
        if ( 
            !(is.vector(isPhased)) || 
            !is.logical(isPhased) || 
            length(isPhased)!=1
        ) {
            stop('MBASED:runMBASED1s: argument isPhased must be a single TRUE or FALSE value')
        }
        if ( 
            !(is.vector(tieBreakRandom)) || 
            !is.logical(tieBreakRandom) || 
            length(tieBreakRandom)!=1 
        ) {
            stop('MBASED:runMBASED1s: argument tieBreakRandom must be a single TRUE or FALSE value')
        }
    }
    ASEResults <- bplapply(unique(aseIDs), function(aseID) {
        aseIDsubv <- (aseIDs==aseID)
        res <- runMBASED1s1aseID(
                lociAllele1Counts=lociAllele1Counts[aseIDsubv],
                lociAllele2Counts=lociAllele2Counts[aseIDsubv],
                lociAllele1NoASEProbs=lociAllele1NoASEProbs[aseIDsubv],
                lociRhos=lociRhos[aseIDsubv],
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
    lociMAF <- unname(unlist(lapply(ASEResults, function(el) {
        el$lociMAF
    })))
    ASEResults <- lapply(ASEResults, function(el) {
        el$lociAllele1IsMajor <- NULL
        el$lociMAF <- NULL
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
            lociMAF=lociMAF
        )
    )
}
