#' @title PDF of cut-points
#' @description This function computes the value of the
#' probability density function for the cut-points. The
#' distribution is induced by a Dirichlet distribution
#' on the prior probabilities of the response.
#' @param tau cut-points
#' @param W global variance
#' @param alpha concentration parameters for prior probabilities of Y
#' @param log logical; if TRUE, returns log pdf
#' @return value of pdf at tau
#' @examples
#' tau = c(-Inf, -1, 1, Inf) # set cut points
#' W = 1 # set value of global variance
#' alpha = c(1,1,1) #concentration parameters
#' dcut(tau, W, alpha, log=FALSE)
#' @export
dcut <- function(tau, W, alpha, log=FALSE){

  #### Check that the inputs are valid.
  if((length(tau)-1)!=length(alpha)){
    stop("The length of tau must be one more than the length of alpha.")
  }

  if(W < 0){
    stop("W must be positive")
  }

  if(min(alpha) < 0){
    stop("All entries of alpha must be positive.")
  }

  if(is.unsorted(tau)){
    stop("tau must be sorted.")
  }
  
  # # Compute prior probabilities for PI based on cut-points tau
  # K = length(tau) + 1 # Number of categories
  # PI = numeric(K)  # Probability of each category: P(Y=k) = PI[k]
  # 
  # PI[1] = stats::pnorm(tau[1], 0, sqrt(1+W))
  # PI[K] = 1 - stats::pnorm(tau[K-1], 0, sqrt(1+W))
  # 
  # for(k in 2:(K-1)){
  #   PI[k] = stats::pnorm(tau[k], 0, sqrt(1+W)) - stats::pnorm(tau[k-1], 0, sqrt(1+W))
  # }

  # Compute prior probabilities for PI based on cut-points tau
  K = length(tau) - 1 # Number of categories
  PI = numeric(K)  # Probability of each category: P(Y=k) = PI[k]

  for(k in 1:K){
    PI[k] = stats::pnorm(tau[k+1], 0, sqrt(1+W)) - stats::pnorm(tau[k], 0, sqrt(1+W))
  }

  J = matrix(0, nrow = K-1, ncol = K-1)

  for(k in 1:(K-1)){
    rho = stats::dnorm(tau[k+1], 0, sqrt(1+W))
    J[k,k]   = rho
    J[k-1,k] = -rho
  }

  # Compute log-pdf
  lpdf = LaplacesDemon::ddirichlet(PI, alpha, log=TRUE) + as.numeric(determinant(J, logarithm = TRUE)$modulus)


  if(log){
    return(lpdf)
  }else{
    return(exp(lpdf))
  }

}

#' @title Log-Likelihood for ordinal regression
#' @description This function evaluates the log-likelihood
#' of the response for a given value of the parameters.
#' @param Y ordinal response
#' @param W global variance
#' @param tau cut-points
#' @return value of log-likelihood at Y, W and tau
#' @examples
#' set.seed(1234)
#' K = 3 # number of response categories
#' Y = sample(1:K, 10, replace=TRUE) # generate responses
#' W = 1 
#' tau = c(-Inf, -1, 1, Inf) # set parameter values
#' llike(Y, W, tau)
#' 
#' @export
llike <- function(Y, W, tau){
  #### Check that the inputs are valid.
  if(sum(Y-floor(Y)==0)!=length(Y)){
    stop("Y must have integer response.") 
  }
  
  if(min(Y) <= 0){
    stop("Y must be positive.") 
  }
  
  if(max(Y)> (length(tau)-1)){
    stop("The length of tau must at least one less than the maximum of Y") 
  }
  
  if(W < 0){
    stop("W must be positive") 
  }    
  
  if(tau[1]!=-Inf || tau[length(tau)]!=Inf){
    stop("The first and last entries in tau must be negative and positive infinity, respectively.")
  }
  
  if(is.unsorted(tau)){
    stop("tau must be sorted.") 
  }
  
  return(
    sum(log(stats::pnorm(tau[Y+1], 0, sqrt(1+W)) - stats::pnorm(tau[Y], 0, sqrt(1+W))))  
  )
}

#' @title Posterior distribution of McFadden's R2
#' @description This function finds the posterior distribution of
#' McFadden's R2 given the posterior samples from a Stan model fit
#' @param Y ordinal response
#' @param out posterior samples from R2D2 model fit in Stan
#' @return Posterior samples from McFadden's R2
#' @examples
#' \donttest{
#' # Obtain output from ord_r2d2() model fit
#' set.seed(1234)
#' # X are covariates, Y are responses, K is number of response categories
#' # This example will yield low R2 values as the response are independent of the covariates.
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
#' # Plot histogram of posterior R2
#' hist(r2_mc(Y, out), xlab="R2")
#' }
#' @export
r2_mc <- function(Y, out){
  #### Check that the inputs are valid.
  if(sum(Y-floor(Y)==0)!=length(Y)){
    stop("Y must have integer response.")
  }

  if(min(Y) <= 0){
    stop("Y must be positive.")
  }

  D <- cbind(out$W, out$tau)
  K = ncol(D)

  apply_fun <- function(d){
    1 - R2D2ordinal::llike(Y, d[1], c(-Inf, d[2:K], Inf)) / R2D2ordinal::llike(Y, 0, c(-Inf, d[2:K], Inf))
  }

  unlist(apply(D, 1, apply_fun))
}

#' @title Find optimal GIG parameters for W prior
#' @description This function finds the optimal GIG parameters for the prior on W
#' which induces a beta prior distribution on McFadden's R2.
#' @param a hyper-parameter of prior for R2 ~ Beta(a,b)
#' @param b hyper-parameter of prior for R2 ~ Beta(a,b)
#' @param n number of observations
#' @param K number of response categories
#' @param alpha prior hyper-parameters for prior Dirichlet distribution on response probabilities
#' @param nsims number of times to simulate data
#' @param nreps number of times to run the algorithm (default = 5)
#' @param no_cores number of cores to parallelize data-generation process
#' @return Optimal GIG parameters
#' @examples
#' \donttest{
#' a = 1
#' b = 5
#' n = 100
#' K = 3
#' find_param(a, b, n, K, no_cores=1)
#' }
#' @export
find_param <- function(a, b, n, K, alpha=rep(2,K), nsims=1000, nreps=5, no_cores=10){

  #### Check that the inputs are valid.
  if(a <= 0 || b <= 0){
    stop("The paramters a and b must both be positive.")
  }

  if(a/(a+b) > 0.5 ){
    warning("Approximation is most accurate for small prior R2 values.")
  }

  if(n <= 0){
    stop("The number of data points must be positive.")
  }

  if(K!=round(K)){
    stop("The number of categories, K, must be an integer.")
  }

  if(length(alpha)!=K){
    stop("The length of alpha must equal K.")
  }

  if(nsims <= 0 || nreps <= 0 || no_cores <= 0){
    stop("The number of simulations, repititions and cores must all be positive.")
  }



  # Objective function to minimize
  # Difference between desired and observed quantiles of R2 distribution
  obj_fun <- function(params){

    findR2 <- function(i){
      ### Simulate data-generating process

      # Sample W (global variance)
      W <- GIGrvg::rgig(1, params[1], exp(params[2]), exp(params[3]))

      # Sample z, latent variable
      z <- stats::rnorm(n, 0, sqrt(1+W))

      # Sample tau, cut-points
      # First sample pi's from Dirichlet.
      PI <- extraDistr::rdirichlet(1, alpha)
      tau <- numeric(K-1)

      # Transform pi's to tau's
      for(k in 1:(K-1)){
        tau[k] = stats::qnorm(sum(PI[1:k]), 0, sqrt(1+W))
      }
      tau <- c(-Inf, tau, Inf)

      # Find observations based on latent variable and cut-points.
      Y <- findInterval(z, tau)

      # Return McFadden's R2 with base model of W=0.
      return(1 - R2D2ordinal::llike(Y, W, tau) / R2D2ordinal::llike(Y, 0, tau)  )
    }

    R2 <- unlist(parallel::mclapply(1:nsims, findR2, mc.cores = no_cores))

    # Compute quantiles of R2
    qx = stats::quantile(R2, seq(0.005, 0.995, 0.005), na.rm = TRUE)

    # Compute quantiles of desired distribution
    qy = stats::qbeta(seq(0.005, 0.995, 0.005), a, b)

    # Find L2 distance between quantiles
    return(sum((qx-qy)^2))
  }

  # Find optimal parameters
  optp = c(0, 0, 0)
  optl = obj_fun(optp)
  # Repeat nreps times to ensure we don't just find a local max
  for(i in 1:nreps){
    # Run optimization routine
    out <- stats::optim(c(0,0,0), obj_fun, method="Nelder-Mead")$par
    loss = obj_fun(out)

    if(loss < optl){
      optp = out
      optl = loss
    }
    print(paste(i, "/", nreps, "iteration(s) complete of hyper-parameter computation."))
  }

  # Return optimal parameters on regular scale
  return(c(optp[1], exp(optp[2:3])))
}





