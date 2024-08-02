data {
  int<lower=1> N;
  real x[N];
  real x0;
  real y[N];
  
}

parameters {
  real alpha;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  for (n in 1:N) {
    if (x[n] <= x0) {
      y[n] ~ normal(alpha + beta1 * (x[n] - x0), sigma);
    } else {
      y[n] ~ normal(alpha + beta2 * (x[n] - x0), sigma);
    }
  }
}

