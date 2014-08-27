#' Main function that implements MBASED.  
#'
#' @param ASESummarizedExperiment SummarizedExperiment object containing information on read counts to be used for ASE detection. Rows represent individual heterozygous loci (SNVs), while columns represent individual samples. There should be either one or two columns, depending on whether one- or two-sample analysis is to be performed. Joint analysis of multiple samples or replicates is currently not supported, and one-sample analysis of multiple samples must be done through independent series of calls to runMBASED(). Note that for two-sample analysis, only loci which are heterozygous in both samples must be supplied (this excludes, e.g., tumor-specific mutations in cases of tumor/normal comparisons). For two-sample analysis, it is assumed that the first column corresponds to 'sample1' and the second column to 'sample2' in the sample1-vs-sample2 comparison. This is important, since differential ASE assessment is not symmetric and sample1-vs-sample2 comparison may yield different results from sample2-vs-sample1 comparison (the relationship is set up by assuming that only instances of ASE greater in sample1 than in sample2 are of interest). assays(ASESummarizedExperiment) must contain matrices lociAllele1Counts and lociAllele2Counts of non-negative integers, containing counts of allele1 (e.g. reference) and allele2 (e.g. alternative) at individual loci. All supplied loci must have total read count (across both alleles) greater than 0 (in each of the two samples, in the case of two-sample analysis). Allele counts are not necessarily phased (see 'isPhased' argument below), so allele1 counts may not represent the same haplotype. assays(ASESummarizedExperiment) may also contain matrix lociAllele1CountsNoASEProbs with entries >0 and <1, containing probabilities of observing allele1-supporting reads at individual loci under conditions of no ASE (which may differ for individual samples in the two-sample analysis). If this matrix is not provided, it is constructed such that every entry in the matrix is set to 0.5 (no pre-existing allelic bias at any locus in any sample). assays(ASESummarizedExperiment) may also contain matrix lociCountsDispersions with entries >=0 and <1, containing dispersion parameters of beta-binomial read count distribution at individual loci (which may differ for individual samples in the two-sample analysis). If this matrix is not provided, it is constructed such that every entry in the matrix is set to 0 (read count-generating distribution at each locus in each sample is binomial). Any other matrices in assays(ASESummarizedExperiment) are ignored by MBASED. rowData(ASESummarizedExperiment) must be supplied by the user, containing additional information about SNVs, including a required column 'aseID', specifying for each locus the unique unit of expression that it belongs to (e.g., gene; must be non-NA). MBASED uses names(rowData(ASESummarizedExperiment)), when specified, to give a unique identifier to each SNV; if no names are provided, the SNVs are labeled 'locus1', 'locus2', ..., in the row order. 
#' @param isPhased specifies whether the true haplotypes are known, in which case the lociAllele1Counts are assumed to represent allelic counts along the same haplotype (and the same is true of lociAllele2Counts). Must be either TRUE or FALSE (DEFAULT).
#' @param numSim number of simulations to perform to estimate statistical signficance of observed ASE. Must be a non-negative integer. If set to 0 (DEFAULT), no simulations are performed and nominal p-values are reported. 
#' @param BPPARAM argument to be passed to function bplapply(), when parallel achitecture is used to speed up simulations (parallelization is done over aseIDs).  DEFAULT: SerialParam() (no parallelization).
#'
#' @return SummarizedExperiment object with rows representing individual aseIDs (genes) and a single column. assays(returnObject) includes single-column matrices 'majorAlleleFrequency' (1-sample analysis only), 'majorAlleleFrequencyDifference' (2-sample analysis only), 'pValueASE' (unadjusted ASE p-value), 'pValueHeterogeneity' (unadjusted inter-loci variability p-value, set to NA for aseIDs with only 1 locus). Note that p-values are not adjusted for multiple hypothesis testing, and the users should carry out such an adjustment themselves, e.g. by employing the utilities in the multtest package. In addition, exptData(returnObject) is a list containing a SummarizedExperiment object names 'locusSpecificResults', with rows corresponding to individual loci (SNVs) and a single column, that provides information on locus-level MBASED analysis results. assays(exptData(returnObject)$locusSpecificResults) contains single-column matrices 'MAF' (estimate of allele frequency for gene-wide major allele at the locus, 1-sample analysis only), 'MAFDifference' (estimate of allele frequency difference for gene-wide major allele at the locus, 2-sample analysis only), and 'allele1IsMajor' (whether allele1 is assigned to major haplotype by MBASED).
#'
#' @export
#'
#' @examples
#' \donttest{
#' mySNVs <- GRanges(
#'     seqnames=c('chr1', 'chr2', 'chr2', 'chr2'),
#'      ranges=IRanges(start=c(1000, 20020, 20285, 21114), width=1),
#'      aseID=c('gene1', rep('gene2', 3)),
#'     allele1=c('G', 'A', 'C', 'A'),
#'     allele2=c('T', 'C', 'T', 'G')    
#' )
#' names(mySNVs) <- paste0('SNV', 1:4)
#' ## SummarizedExperiment object with data to run tumor vs. normal comparison
#' mySE_TumorVsNormal <- SummarizedExperiment(
#'      assays=list(
#'         lociAllele1Counts=matrix(
#'              c(
#'                 c(25,10,22,14),
#'                 c(18,17,14,28)
#'             ),
#'             ncol=2,
#'             dimnames=list(
#'                 names(mySNVs), 
#'                 c('tumor', 'normal')
#'             )
#'          ),
#'         lociAllele2Counts=matrix(
#'              c(
#'                 c(20,16,15,16),
#'                  c(23,9,24,17)
#'             ),
#'             ncol=2, 
#'             dimnames=list(
#'                 names(mySNVs), 
#'                  c('tumor', 'normal')
#'             )
#'          ),
#'          lociAllele1CountsNoASEProbs=matrix(
#'              c(
#'                  c(0.48, 0.51, 0.55, 0.45),
#'                  c(0.52, 0.43, 0.52, 0.43)
#'              ),
#'              ncol=2, 
#'              dimnames=list(
#'                  names(mySNVs), 
#'                  c('tumor', 'normal')
#'              )
#'          ),
#'         lociCountsDispersions=matrix(
#'             c(
#'                  c(0.005, 0.007, 0.003, 0.01),
#'                  c(0.001, 0.004, 0.02, 0.006)
#'              ),
#'             ncol=2,
#'             dimnames=list(
#'                  names(mySNVs), 
#'                  c('tumor', 'normal')
#'             )
#'          )
#'     ),
#'      rowData=mySNVs
#' )
#' twoSampleAnalysisTumorVsNormal <- runMBASED(
#'     ASESummarizedExperiment=mySE_TumorVsNormal,
#'      numSim=10^6, 
#'      BPPARAM=SerialParam(),
#'      isPhased=FALSE
#' ) 
#' rowData(twoSampleAnalysisTumorVsNormal)
#' assays(twoSampleAnalysisTumorVsNormal)$majorAlleleFrequencyDifference
#' assays(twoSampleAnalysisTumorVsNormal)$pValueASE
#' assays(twoSampleAnalysisTumorVsNormal)$pValueHeterogeneity
#' assays(exptData(twoSampleAnalysisTumorVsNormal)$locusSpecificResults)$MAFDifference
#' assays(exptData(twoSampleAnalysisTumorVsNormal)$locusSpecificResults)$allele1IsMajor
#' 
#' ## exchanging the order of the columns will allow us to run normal vs. tumor comparison
#' ## Note that while results are the same for single-locus gene1, they differ for multi-locus gene2
#' mySE_NormalVsTumor <- SummarizedExperiment(
#'     assays=lapply(names(assays(mySE_TumorVsNormal)), function(matName) {
#'         curMat <- assays(mySE_TumorVsNormal)[[matName]]
#'         modifiedMat <- curMat[,c('normal','tumor')]
#'         return(modifiedMat)
#'     }),
#'     colData=colData(mySE_TumorVsNormal)[2:1,],
#'     rowData=rowData(mySE_TumorVsNormal)
#' )
#' names(assays(mySE_NormalVsTumor )) <- names(assays(mySE_TumorVsNormal))
#' twoSampleAnalysisNormalVsTumor <- runMBASED(
#'     ASESummarizedExperiment=mySE_NormalVsTumor,
#'      numSim=10^6, 
#'      BPPARAM=SerialParam(),
#'      isPhased=FALSE
#' ) 
#' rowData(twoSampleAnalysisNormalVsTumor)
#' assays(twoSampleAnalysisNormalVsTumor)$majorAlleleFrequencyDifference
#' assays(twoSampleAnalysisNormalVsTumor)$pValueASE
#' assays(twoSampleAnalysisNormalVsTumor)$pValueHeterogeneity
#' assays(exptData(twoSampleAnalysisNormalVsTumor)$locusSpecificResults)$MAFDifference
#' assays(exptData(twoSampleAnalysisNormalVsTumor)$locusSpecificResults)$allele1IsMajor
#' 
#' ## we can also do separate one-sample analysis on tumor and normal samples
#' mySE_Tumor <- SummarizedExperiment(
#'     assays=lapply(names(assays(mySE_TumorVsNormal)), function(matName) {
#'         curMat <- assays(mySE_TumorVsNormal)[[matName]]
#'         modifiedMat <- curMat[,'tumor',drop=FALSE]
#'         return(modifiedMat)
#'     }),
#'     colData=colData(mySE_TumorVsNormal)[1,],
#'     rowData=rowData(mySE_TumorVsNormal)
#' )
#' names(assays(mySE_Tumor)) <- names(assays(mySE_TumorVsNormal))
#' oneSampleAnalysisTumor <- runMBASED(
#'     ASESummarizedExperiment=mySE_Tumor,
#'      numSim=10^6, 
#'      BPPARAM=SerialParam(),
#'      isPhased=FALSE
#' )
#' rowData(oneSampleAnalysisTumor)
#' assays(oneSampleAnalysisTumor)$majorAlleleFrequency
#' assays(oneSampleAnalysisTumor)$pValueASE
#' assays(oneSampleAnalysisTumor)$pValueHeterogeneity
#' assays(exptData(oneSampleAnalysisTumor)$locusSpecificResults)$MAF
#' assays(exptData(oneSampleAnalysisTumor)$locusSpecificResults)$allele1IsMajor
#' 
#' mySE_Normal <- SummarizedExperiment(
#'     assays=lapply(names(assays(mySE_TumorVsNormal)), function(matName) {
#'         curMat <- assays(mySE_TumorVsNormal)[[matName]]
#'         modifiedMat <- curMat[,'normal',drop=FALSE]
#'         return(modifiedMat)
#'     }),
#'     colData=colData(mySE_TumorVsNormal)[1,],
#'     rowData=rowData(mySE_TumorVsNormal)
#' )
#' names(assays(mySE_Normal)) <- names(assays(mySE_TumorVsNormal))
#' oneSampleAnalysisNormal <- runMBASED(
#'     ASESummarizedExperiment=mySE_Normal,
#'      numSim=10^6, 
#'      BPPARAM=SerialParam(),
#'      isPhased=FALSE
#' )
#' rowData(oneSampleAnalysisNormal)
#' assays(oneSampleAnalysisNormal)$majorAlleleFrequency
#' assays(oneSampleAnalysisNormal)$pValueASE
#' assays(oneSampleAnalysisNormal)$pValueHeterogeneity
#' assays(exptData(oneSampleAnalysisNormal)$locusSpecificResults)$MAF
#' assays(exptData(oneSampleAnalysisNormal)$locusSpecificResults)$allele1IsMajor
#' }
#' 
runMBASED <- function (
    ASESummarizedExperiment,
    isPhased=FALSE,
    numSim=0, 
    BPPARAM=SerialParam()
){
    if ( 
        !(is.vector(isPhased)) || 
        !is.logical(isPhased) || 
        length(isPhased)!=1 ||
        is.na(isPhased)
    ) {
        stop('MBASED:runMBASED: argument isPhased must be a single TRUE or FALSE value')
    }
    if ( 
        !(is.vector(numSim)) || 
        !is.numeric(numSim) || 
        length(numSim)!=1 || 
        is.na(numSim) || 
        !isTRUE(all.equal(numSim, round(numSim))) || 
        (numSim<0) 
    ) {
        stop('MBASED:runMBASED: argument numSim must be a single non-negative integer')
    }
    if (!(is(BPPARAM, 'BiocParallelParam'))) {
        stop('MBASED:runMBASED: argument BPPARAM must be an instance of BiocParallelParam')
    }
    lociAllele1Counts <- assays(ASESummarizedExperiment)$lociAllele1Counts
    lociAllele2Counts <- assays(ASESummarizedExperiment)$lociAllele2Counts
    if (
        is.null(lociAllele1Counts) || 
        is.null(lociAllele2Counts)
    ) {
        stop('MBASED:runMBASED: assay slot of ASESummarizedExperiment must contain matrices lociAllele1Counts  and lociAllele2Counts')
    }
    if (
        !is.matrix(lociAllele1Counts) ||
        !(ncol(lociAllele1Counts)%in%c(1,2)) ||
        !is.numeric(lociAllele1Counts) ||
        !isTRUE(all.equal(lociAllele1Counts, round(lociAllele1Counts))) ||
        any(lociAllele1Counts<0) 
    ) {
        stop('MBASED:runMBASED: lociAllele1Counts must be a matrix of non-negative integers with 1 or 2 columns')
    }
    if (
        !is.matrix(lociAllele2Counts) ||
        !(ncol(lociAllele2Counts)%in%c(1,2)) ||
        !is.numeric(lociAllele2Counts) ||
        !isTRUE(all.equal(lociAllele2Counts, round(lociAllele2Counts))) ||
        any(lociAllele2Counts<0)
    ) {
        stop('MBASED:runMBASED: lociAllele2Counts must be a matrix of non-negative integers with 1 or 2 columns')
    }
    ## determine if one-sample or two-sample analysis is to be performed
    if (ncol(lociAllele1Counts)==2) {
        pairedAnalysis <- TRUE
    } else {
        pairedAnalysis <- FALSE
    }
    ## set noASE allelic probabilities to default (0.5), if not supplied
    if (!(
        'lociAllele1CountsNoASEProbs'%in%
            names(
                assays(
                    ASESummarizedExperiment
                )
            )
        )
    ) {
    	## matrix with every entry equal to 0.5
        lociAllele1CountsNoASEProbs <- lociAllele1Counts-lociAllele1Counts+0.5 
    } else {
        lociAllele1CountsNoASEProbs <- assays(
            ASESummarizedExperiment
        )$lociAllele1CountsNoASEProbs
    }
    if (
        !is.matrix(lociAllele1CountsNoASEProbs) ||
        !(ncol(lociAllele1CountsNoASEProbs)%in%c(1,2)) ||
        !is.numeric(lociAllele1CountsNoASEProbs) ||
        any(lociAllele1CountsNoASEProbs<=0) ||
        any(lociAllele1CountsNoASEProbs>=1)
    ) {
        stop('MBASED:runMBASED: lociAllele1CountsNoASEProbs must be a matrix with 1 or 2 columns and with entries >0 and <1')
    }
    ## set dispersion values to default (0), if not supplied
    if (!(
        'lociCountsDispersions'%in%
            names(
                assays(
                    ASESummarizedExperiment
                )
            )
        )
    ) {
    	## matrix with every entry equal to 0
        lociCountsDispersions <- lociAllele1Counts-lociAllele1Counts 
    } else {
        lociCountsDispersions <- assays(
            ASESummarizedExperiment
        )$lociCountsDispersions
    }
    if (
        !is.matrix(lociCountsDispersions) ||
        !(ncol(lociCountsDispersions)%in%c(1,2)) ||
        !is.numeric(lociCountsDispersions) ||
        any(lociCountsDispersions<0) ||
        any(lociCountsDispersions>=1)
    ) {
        stop('MBASED:runMBASED: lociCountsDispersions must be a matrix with 1 or 2 columns and with entries >=0 and <1')
    }
    locusInfo <- rowData(ASESummarizedExperiment)
    if (
        is.null(locusInfo) ||
        !(is(locusInfo, 'GRanges')) ||
        is.null(locusInfo$aseID)
    ) {
        stop('MBASED:runMBASED: rowData(ASESummarizedExperiment) must be non-null GRanges object with a column aseID')
    }
    ## assign names to loci, if not already provided
    if (is.null(names(locusInfo))) {
        names(locusInfo) <- paste0('locus', 1:length(locusInfo))
    }
    aseIDs <- locusInfo$aseID
    if (
        !is.vector(aseIDs) ||
        any(is.na(aseIDs))
    ) {
        stop('MBASED:runMBASED: aseIDs must be a vector with no NAs')
    }
    aseIDs <- as.character(aseIDs)
    aseIDsUnique <- unique(aseIDs)
    ## run ASE on supplied samples
    if (pairedAnalysis) { ## two-sample analysis
        MBASEDResults <- runMBASED2s(
            lociAllele1CountsSample1=lociAllele1Counts[,1], 
            lociAllele2CountsSample1=lociAllele2Counts[,1], 
            lociAllele1CountsSample2=lociAllele1Counts[,2], 
            lociAllele2CountsSample2=lociAllele2Counts[,2],  
            lociAllele1NoASEProbsSample1=lociAllele1CountsNoASEProbs[,1], 
            lociAllele1NoASEProbsSample2=lociAllele1CountsNoASEProbs[,2], 
            lociRhosSample1=lociCountsDispersions[,1], 
            lociRhosSample2=lociCountsDispersions[,2], 
            aseIDs=aseIDs, 
            numSim=numSim, 
            BPPARAM=BPPARAM, 
            isPhased=isPhased, 
            tieBreakRandom=TRUE
        )
    } else { ## one-sample analysis
        MBASEDResults <- runMBASED1s(
            lociAllele1Counts=lociAllele1Counts[,1], 
            lociAllele2Counts=lociAllele2Counts[,1], 
            lociAllele1NoASEProbs=lociAllele1CountsNoASEProbs[,1],
            lociRhos=lociCountsDispersions[,1],   
            aseIDs=aseIDs, 
            numSim=numSim, 
            BPPARAM=BPPARAM, 
            isPhased=isPhased, 
            tieBreakRandom=TRUE
        ) 
    } 
    ## create output object
    names(MBASEDResults) <- gsub('loci', '', names(MBASEDResults))
    if (pairedAnalysis) {
        aseIDMatrixNames <- c(
            'majorAlleleFrequencyDifference',
            'pValueASE',
            'pValueHeterogeneity'
        )
        lociMatrixNames <- c(
            'allele1IsMajor',
            'MAFDifference'
        )
    } else {
        aseIDMatrixNames <- c(
            'majorAlleleFrequency',
            'pValueASE',
            'pValueHeterogeneity'
        )
        lociMatrixNames <- c(
            'allele1IsMajor',
            'MAF'
        )
    }
    if (pairedAnalysis) {
        outputName <- paste(colnames(lociAllele1Counts), collapse='/')
    } else {
        outputName <- colnames(lociAllele1Counts)
    }
    lociSpecificResults <- lapply(lociMatrixNames, function(lociMatrixName) {
        matrix(
            MBASEDResults[[lociMatrixName]],
            ncol=1,
            nrow=length(locusInfo),
            dimnames=list(
                names(locusInfo),
                outputName
            )
        )
    })
    names(lociSpecificResults) <- lociMatrixNames
    aseIDSpecificResults <- lapply(
        aseIDMatrixNames, 
        function(aseIDMatrixName) {
            matrix(
                MBASEDResults$ASEResults[,aseIDMatrixName],
                ncol=1,
                nrow=length(aseIDsUnique),
                dimnames=list(
                    aseIDsUnique,
                    outputName
                )
            )
        }
    )
    names(aseIDSpecificResults) <- aseIDMatrixNames
    outputSummarizedExperiment <- SummarizedExperiment(
        assays=aseIDSpecificResults,
        rowData <- split(locusInfo, factor(locusInfo$aseID, levels=aseIDsUnique)),
        exptData=SimpleList(
            locusSpecificResults=SummarizedExperiment(
                assays=lociSpecificResults,
                rowData <- locusInfo
            )
        )
    )
    return(outputSummarizedExperiment)
}
