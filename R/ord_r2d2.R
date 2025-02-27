#' @title Ordinal regression in Stan with R2D2 prior
#' @description This function carries out a Bayesian ordinal regression model
#' in Stan using the proposed psuedo-R2D2 prior
#' @param x covariate matrix
#' @param y response variables
#' @param K number of response categories
#' @param a hyper-parameter of prior for R2 ~ Beta(a,b)
#' @param b hyper-parameter of prior for R2 ~ Beta(a,b)
#' @param hyper hyper-parameters for W prior
#' @param alpha prior hyper-parameters for prior Dirichlet distribution on response probabilities
#' @param nsims number of times to simulate data
#' @param nreps number of times to run the algorithm (default = 5)
#' @param no_cores number of cores to parallelize data-generation process
#' @param progress logical. if TRUE, shows the progress bars from the posterior sampling.
#' @param ... optional hyper-parameters for Stan fitting
#' @return Stan model fit
#' @examples
#' \donttest{
#' # X are covariates, Y are responses, K is number of response categories
#' # This example will yield low R2 values as the response are independent of the covariates.
#' set.seed(1234)
#' n = 100
#' p = 5
#' X = matrix(rnorm(n*p), nrow = n, ncol=p)
#' K = 3
#' Y = sample(1:K, 100, replace=TRUE)
#' a = 1
#' b = 5
#' # Pre-computed hyperparameters
#' fit <- ord_r2d2(X, Y, K, hyper=c(0.002, 0.989, 1.013), no_cores=1)
#' out <- rstan::extract(fit)
#' # Plot histogram of posterior W
#' hist(out$W, xlab="W")
#' }
#' @export
ord_r2d2 <- function(x, y, K, a=1, b=10, hyper=NULL,
                     alpha=rep(1,K), nsims=1000, nreps=5, no_cores=10, progress=FALSE, ...){

  #### Check that the inputs are valid.
  if(length(y)!=nrow(x)){
    stop("The number of observations in y must equal the number of rows in x.")
  }

  n = length(y)
  p = ncol(x)

  # Only compute hyper-parameters if not pre-specified
  if(is.null(hyper)){
    hyper = R2D2ordinal::find_param(a, b, n, K, alpha, nsims, nreps, no_cores)
  }

  standata  = list(Y=y, X=x, K=K, N=n, p=p, hyper = hyper)

  init=list()
  for(i in 1:4){
    init[[i]] <-list("beta" = rep(0,p), "tau" = sort(stats::runif(K-1, -3, 3)), "W" = 1,
                     "xi"=1, "phi"=rep(1/p, p))
  }

  if(progress){
    out <- rstan::sampling(stanmodels$r2d2ord, data = standata, init=init, ...)
  }else{
    out <- rstan::sampling(stanmodels$r2d2ord, data = standata, init=init, refresh=0, ...)
  }

  return(out)
}




