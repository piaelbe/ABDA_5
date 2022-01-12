data {
  int<lower=0> J; //observations
  int<lower=0> P; //participants
  vector[J] y;                    // Vector of observations.
  int<lower = 1, upper = P> g[J]; // Vector of group assignments. 
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[P] beta;                 // Vector of group intercepts.
  real<lower=0> sigma;
}

model { 
  // Hyperpriors.
  mu ~ normal(0, 5);
  tau ~ normal(0, 5);
  
    // Prior.
  sigma ~ normal(0, 5);
 
  for (j in 1:P) {
    beta[j] ~ normal(mu, tau);
  }

   // Population model and likelihood.
  for (j in 1:J) {
    y[j] ~ normal(beta[g[j]], sigma); //different beta for each participant

  }
}

generated quantities {
  real beta_predicted = normal_rng(mu, tau);
  real logy_predicted = normal_rng(beta_predicted,sigma);
}
