#' Functions to convert between shape parameters a and b for beta distribution and parameters mu (mean) and rho (dispersion).
#'
#' @details getMuRho takes in shape parameters a and b and returns list with parameters mu (a/(a+b)) and rho (1/(a+b+1)). The function is vectorized (both a and b can be vectors (of the same length) or matrices (of the same dimension)).
#'
#' @param a,b shape parameters for beta distribution. Must be >0.
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return getMuRho returns a list with 2 elements: mu and rho (vectors, if the arguments a and b are vectors).
#' 
#' @rdname betaParametrizationConverters
#'
#' @family bbFunctions
#'
#' @aliases getMuRho getAB
#'
#' @examples
#' MBASED:::getMuRho(a=1, b=1)
#'
getMuRho <- function (
    a, 
    b, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            (!is.vector(a) && !is.matrix(a)) || 
            any(is.na(a)) || 
            !is.numeric(a) || 
            any(a<=0) 
        ) {
            stop('MBASED:getMuRho: argument a must be (vector/matrix of) positive number(s)')
        }
        if ( 
            (!is.vector(b) && !is.matrix(b)) || 
            any(is.na(b)) || 
            !is.numeric(b) || 
            any(b<=0) 
        ) {
            stop('MBASED:getMuRho: argument b must be (vector/matrix of) positive number(s)')
        }
        if (
            !isTRUE(all.equal(length(a),length(b))) || 
            !isTRUE(all.equal(dim(a), dim(b)))
        ) {
            stop('MBASED:getMuRho: arguments a and b must be of same length/dimension')
        }
    }
    return(
        list(
            mu=a/(a+b), 
            rho=1/(a+b+1)
        )
    )
}

#'
#' @details getAB takes in shape mean and dispersion parameters mu and rho and returns shape parameters a (mu*(1/rho-1)) and b ((1-mu)*(1/rho-1)). The function is vectorized (both mu and rho can be vectors (of the same length) or matrices (of the same dimension)).
#'
#' @param mu,rho mean and dispersion parameters for beta distribution, respectively. Must be in (0,1) interval, although rho is allowed to take on value of 0 (binomial distribution). 
#'
#' @return getAB returns a list with 2 elements: a and b (vectors, if arguments mu and rho are vectors). For values of rho=0, the resulting entries are NA.
#' 
#' @rdname betaParametrizationConverters
#'
#' @family bbFunctions
#'
#' @examples
#' MBASED:::getAB(mu=1/2, rho=1/3)
#' MBASED:::getMuRho(MBASED:::getAB(mu=0.7, rho=0.0045)$a, MBASED:::getAB(mu=0.7, rho=0.0045)$b) 
#' MBASED:::getAB(MBASED:::getMuRho(a=0.2, b=4)$mu, MBASED:::getMuRho(a=0.2, b=4)$rho)
#'
getAB <- function (
    mu, 
    rho, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            (!is.vector(mu) && !is.matrix(mu)) || 
            any(is.na(mu)) || 
            !is.numeric(mu) || 
            any(mu<=0) || 
            any(mu>=1)
        ) {
            stop('MBASED:getAB: argument mu must be (vector/matrix of) value(s) >0 and <1')
        }
        if ( 
            (!is.vector(rho) && !is.matrix(rho)) || 
            any(is.na(rho)) || 
            !is.numeric(rho) || 
            any(rho<0) || 
            any(rho>=1)
        ) {
            stop('MBASED:getAB: argument rho must be (vector/matrix of) value(s) >=0 and <1')
        }
        if (
            !isTRUE(all.equal(length(mu),length(rho))) || 
            !isTRUE(all.equal(dim(mu), dim(rho)))
        ) {
            stop('MBASED:getAB: arguments mu and rho must be of same length/dimension')
        }
    }
    return(
          list(
              a=ifelse(rho==0, NA, mu*(1/rho-1)), 
              b=ifelse(rho==0, NA, (1-mu)*(1/rho-1))
          )
      )
}

#' Functions to generate beta-binomial random variables.
#'
#' @details vectorizedRbetabinomAB is the same function as rbetabinom.ab from VGAM package but it avoids a lot of overhang and requires that arguments size, a (shape1), and b (shape2) be of length equal to argument n.
#'
#' @param n sample size, must be a single positive integer
#' @param size number of trials for each count to be generated in the sample, must be a vector of positive integers
#' @param a,b vectors of shape parameters for beta distributions used to generate probability of success for each count to be generated in the sample, must be >0
#' @param checkArgs single boolean specifying whether arguments should be checked for adherence to specifications. DEFAULT: FALSE
#'
#' @return a numeric vector of betabinomial random variables.
#' 
#' @rdname generatingBetaBinomialCounts
#'
#' @family bbFunctions
#'
#' @aliases vectorizedRbetabinomAB vectorizedRbetabinomMR
#'
#' @examples
#' set.seed(111)
#' MBASED:::vectorizedRbetabinomAB(n=10, size=rep(50,10), a=rep(1,10), b=rep(1,10))
#' 
vectorizedRbetabinomAB <- function (
    n, 
    size, 
    a, 
    b, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if (
            !is.vector(n) || 
            length(n)!=1 || 
            is.na(n) || 
            !is.numeric(n) || 
            !isTRUE(all.equal(n, round(n))) || 
            n<1
        ) {
            stop('MBASED:vectorizedRbetabinomAB: argument n must be a single positive integer')
        }
        if ( 
            !is.vector(size) || 
            length(size)!=n || 
            any(is.na(size)) || 
            !is.numeric(size) || 
            !isTRUE(all.equal(size, round(size))) || 
            any(size<1)
        ) {
            stop('MBASED:vectorizedRbetabinomAB: argument size must be a vector of positive integers of length equal to argument n')
        }
        if ( 
            !is.vector(a) || 
            length(a)!=n || 
            any(is.na(a)) || 
            !is.numeric(a) || 
            any(a<=0)
        ) {
            stop('MBASED:vectorizedRbetabinomAB: argument a must be a vector of positive numbers of length equal to argument n')
        }
        if ( 
            !is.vector(b) || 
            length(b)!=n || 
            any(is.na(b)) || 
            !is.numeric(b) || 
            any(b<=0)
        ) {
            stop('MBASED:vectorizedRbetabinomAB: argument b must be a vector of positive numbers of length equal to argument n')
        }
    }
    probs <- rbeta(
        n=n, 
        shape1=a, 
        shape2=b
    )
      return(
          rbinom(
              n=n, 
              size=size, 
              prob=probs
          )
      )
}    
  
#'
#' @details vectorizedRbetabinomMR is a wrapper around vectorizedRbetabinomAB using mu/rho parametrization.  Requires that arguments size, mu, and rho be of length equal to argument n.
#'
#' @param mu,rho mean (a/(a+b)) and dispersion (1/(a+b+1)) parameters for beta distribution, must be in (0,1). Value of 0 is allowed for rho and implies binomial distribution.
#'
#' @rdname generatingBetaBinomialCounts
#'
#' @family bbFunctions
#'
#' @examples
#' set.seed(111)
#' MBASED:::vectorizedRbetabinomMR(n=10, size=rep(50,10), mu=rep(1/2,10), rho=rep(1/3,10))
#'
vectorizedRbetabinomMR <- function (
    n, 
    size, 
    mu, 
    rho, 
    checkArgs=FALSE
) {
    if (checkArgs) {
        if ( 
            !is.vector(n) || 
            length(n)!=1 || 
            is.na(n) || 
            !is.numeric(n) || 
            !isTRUE(all.equal(n, round(n))) || 
            n<1
        ) {
            stop('MBASED:vectorizedRbetabinomMR: argument n must be a single positive integer')
        }
        if ( 
            !is.vector(size) || 
            length(size)!=n || 
            any(is.na(size)) || 
            !is.numeric(size) || 
            !isTRUE(all.equal(size, round(size))) || 
            any(size<1)
        ) {
            stop('MBASED:vectorizedRbetabinomMR: argument size must be a vector of positive integers of length equal to argument n')
        }
        if ( 
            !is.vector(mu) || 
            length(mu)!=n || 
            any(is.na(mu)) || 
            !is.numeric(mu) || 
            any(mu<=0) || 
            any(mu>=1)
        ) {
            stop('MBASED:vectorizedRbetabinomMR: argument mu must be a vector of positive numbers >0 and <1 of length equal to argument n')
        }
        if ( 
            !is.vector(rho) || 
            length(rho)!=n || 
            any(is.na(rho)) || 
            !is.numeric(rho) ||  
            any(rho<0) || 
            any(rho>=1)
        ) {
            stop('MBASED:vectorizedRbetabinomMR: argument rho must be a vector of positive numbers >=0 and <1 of length equal to argument n')
        }
    }
    probs <- mu
    betabinSubv <- (rho>0)
    ab <- getAB(
        mu=mu[betabinSubv], 
        rho=rho[betabinSubv]
    )
    probs[betabinSubv] <- rbeta(
        n=sum(betabinSubv), 
        shape1=ab$a, 
        shape2=ab$b
    )
    return(
        rbinom(
            n=n, 
            size=size, 
            prob=probs
        )
    )
}
 