## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
    n = 100 # Number of data points
    p = 5   # Number of covariates
    K = 3   # Number of response categories
    
    # Regression coefficients
    beta <- c(1, -1, 2, -2, 0)

    # Generate covariates
    X <- matrix(rnorm(n*p), nrow=n, ncol=p)
    
    # Generate latent term
    eta = X%*%beta
    z <- rnorm(n, eta, 1)
    
    # Set cut-off values
    tau <- c(-Inf, 0, 2/3, Inf)

    # Find response values
    Y <- findInterval(z, tau)

## -----------------------------------------------------------------------------
a = 1
b = 10
# In practice, this function should be run.
# hyper <- R2D2ordinal::find_param(a, b, n, K, nreps = 5, no_cores=10)
hyper = c(0.045, 1.05, 1.05)

## -----------------------------------------------------------------------------
out <- R2D2ordinal::ord_r2d2(x=X, y=Y, K=K, hyper=hyper, warmup=500, iter=1500, verbose=FALSE)
fit <- rstan::extract(out)

## -----------------------------------------------------------------------------
df = dplyr::tibble(W = fit$W,
            beta1 = fit$beta[,1],
            beta2 = fit$beta[,2])

ggplot2::ggplot(df, ggplot2::aes(x=beta1))+
  ggplot2::geom_histogram(bins=25)+
  ggplot2::geom_vline(xintercept = beta[1], color="red")+
  ggplot2::xlab("beta1")+
  ggplot2::ylab("Count")+
  ggplot2::theme_bw()
ggplot2::ggplot(df, ggplot2::aes(x=beta2))+
  ggplot2::geom_histogram(bins=25)+
  ggplot2::geom_vline(xintercept = beta[2], color="red")+
  ggplot2::xlab("beta2")+
  ggplot2::ylab("Count")+
  ggplot2::theme_bw()
ggplot2::ggplot(df, ggplot2::aes(x=W))+
  ggplot2::geom_histogram(bins=25)+
  ggplot2::xlab("W")+
  ggplot2::ylab("Count")+
  ggplot2::theme_bw()

