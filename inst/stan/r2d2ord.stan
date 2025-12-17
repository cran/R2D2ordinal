functions {
  // Induced distribution on cutpoints, tau
  real induced_dirichlet_lpdf(vector tau, real W, vector alpha) {
    int K = num_elements(tau) + 1;
    vector[K] PI;
    matrix[K, K] J = rep_matrix(0, K, K);

    // Induced ordinal probabilities
    PI[1] = normal_cdf(tau[1]| 0, sqrt(1+W));
    for (k in 2:(K - 1))
      PI[k] = normal_cdf(tau[k]| 0, sqrt(1+W)) - normal_cdf(tau[k-1]| 0, sqrt(1+W));
    PI[K] = 1 - normal_cdf(tau[K-1]| 0, sqrt(1+W));

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = exp(normal_lpdf(tau[k-1]| 0, sqrt(1+W)));
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    return   dirichlet_lpdf(PI | alpha) + log_determinant(J);
  }

  // Inverse Gaussian distribution for global variance, W
  real inv_gauss_lpdf(real x, real mu, real lambda){
    real lpdf = 0.5*log(lambda) - 0.5*log(2*pi()) - 1.5*log(x) - lambda*(x-mu)^2/(2*mu^2*x);
    return lpdf;
  }

}

data {
  int<lower=2> K;
  int<lower=0> N;
  int<lower=0> p;
  array[N] row_vector[p] X;
  array[N] int<lower=1,upper=K> Y;
  //row_vector[p] X[N]; // deprecated syntax
  //int<lower=1,upper=K> Y[N];
  row_vector[3] hyper;
  vector[K] alpha;
  real<lower=0> xi0;
}
parameters {
  ordered[K-1] tau;    // Internal cutpoints
  vector[p] beta;      // Regression coefficients
  real<lower=0> W;     // Global variance
  real<lower=0> xi; // Latent term for global variance
  simplex[p] phi;      // Variance allocation parameters
}
model {
  tau ~ induced_dirichlet(W, alpha);
  for(j in 1:p){beta[j] ~ normal(0, sqrt(phi[j]*W));}
  phi ~ dirichlet(rep_vector(xi0,p));

  if(hyper[1] < -0.5){
    W  ~ inv_gauss(sqrt(hyper[3]/(hyper[2]+2*xi)), hyper[3]);
    xi ~ gamma(-(hyper[1]+0.5), W);
  }else{
    W  ~ inv_gauss(sqrt((hyper[3]+2*xi)/(hyper[2])), hyper[3]+2*xi);
    xi ~ gamma(hyper[1]+0.5, 1/W);
  }


  for (n in 1:N) {
    vector[K] theta;
    real gamma = X[n]*beta;

    theta[1] = Phi(tau[1]-gamma);
    for (k in 2:(K-1))
      theta[k] = Phi(tau[k]-gamma) - Phi(tau[k-1]-gamma);
    theta[K] = 1-Phi(tau[K-1]-gamma);
    Y[n] ~ categorical(theta);
  }
}
