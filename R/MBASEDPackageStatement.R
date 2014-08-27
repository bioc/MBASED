#' MBASED
#'
#' Package that contains functions to process sets of SNVs and determine the genes that show allele-specific expression (ASE)
#'
#' @docType package
#'
#' @name MBASED
#'
#' @aliases MBASED MBASED-package
#'
#' @import RUnit BiocGenerics BiocParallel GenomicRanges
#'
#' @details The package implements MBASED method for detecting allele-specific gene expression.  The main workhorse function is runMBASED which is used to run both 1-sample and 2-sample (allelic imbalance) analyses.  Please consult the accompanying vignette and the runMBASED help page for more details.
#'
#' @author Oleg Mayba <maybao@@gene.com> Houston Gilbert <gilbert.houston@@gene.com>
NULL
