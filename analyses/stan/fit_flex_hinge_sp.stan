data {
  int<lower=1> N;
  int<lower=0> Nspp;
  real x[N];
  real x0;
  real y[N];
  int species[N];
  
}

parameters {
  real alpha[Nspp];
  real beta1[Nspp];
  real beta2[Nspp];
  real<lower=0> sigma;
  
  real mu_alpha;
  real<lower = 0> sigma_a;
  real mu_beta1;
  real mu_beta2;
  real<lower = 0> sigma_b1;
  real<lower = 0> sigma_b2;
}

model {
  alpha ~ normal(mu_alpha, sigma_a);
  beta1 ~ normal(mu_beta1, sigma_b1);
  beta2 ~ normal(mu_beta2, sigma_b2);
  
  mu_alpha ~ normal(188, 50);
  sigma_a ~normal(0,50);
  mu_beta1 ~ normal(0, 10);
  mu_beta2 ~ normal(0, 10);
  sigma_b1 ~ normal(0, 10);
  sigma_b2 ~ normal(0, 10);
  
  sigma ~normal(0,10); 
 
  
  for (n in 1:N) {
    if (x[n] <= x0) {
      y[n] ~ normal(alpha[species[n]] + beta1[species[n]] * (x[n] - x0), sigma);
    } else {
      y[n] ~ normal(alpha[species[n]] + beta2[species[n]] * (x[n] - x0), sigma);
    }
  }
}

