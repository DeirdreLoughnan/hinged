// Started Aug 15, 2024 by Deirdre
// Aim of this code is to have the two hinges with species and study level intercepts
data {
  int<lower=1> N;
  int<lower=0> Nspp;
  int<lower=0> Nstudy;
  real x[N];
  real x0;
  real y[N];
  int species[N];
  int study[N];
  
}

parameters {
  
  real mu_grand;
  
  real alpha_sp[Nspp];
  real<lower=0> sigma_sp;
  
  real alpha_study[Nstudy];
  real<lower=0> sigma_study;
  
  real beta1[Nspp];
  real beta2[Nspp];
  real<lower=0> sigma;
  
  real mu_beta1;
  real mu_beta2;
  real<lower = 0> sigma_b1;
  real<lower = 0> sigma_b2;
}

model {
  
  alpha_sp ~ normal(0, sigma_sp);
  alpha_study ~ normal(0, sigma_study);
  beta1 ~ normal(mu_beta1, sigma_b1);
  beta2 ~ normal(mu_beta2, sigma_b2);
  
  mu_grand ~ normal(188, 50);
  sigma_sp ~normal(0,50);
  sigma_study ~normal(0,50);
  mu_beta1 ~ normal(0, 10);
  mu_beta2 ~ normal(0, 10);
  sigma_b1 ~ normal(0, 10);
  sigma_b2 ~ normal(0, 10);
  
  sigma ~normal(0,10); 
 
  
  for (n in 1:N) {
    if (x[n] <= x0) {
      y[n] ~ normal(mu_grand + alpha_sp[species[n]] + alpha_study[study[n]] + beta1[species[n]] * (x[n] - x0), sigma);
    } else {
      y[n] ~ normal(mu_grand + alpha_sp[species[n]] + alpha_study[study[n]] + beta2[species[n]] * (x[n] - x0), sigma);
    }
  }
}

