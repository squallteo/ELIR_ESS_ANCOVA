data {
  int<lower=0> N;
  vector[N] base;
  vector[N] trt;
  vector[N] y;
  real<lower=0> sigma;
  real delta_m;
  real delta_s;
}
parameters {
  real beta0;
  real beta1;
  real delta;
}
model {
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  delta ~ normal(delta_m, delta_s);
  y ~ normal(beta0 + beta1*base + delta*trt, sigma);
}
