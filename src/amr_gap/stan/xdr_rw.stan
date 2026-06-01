// Beta-binomial random-walk logit model for XDR prevalence projection.
//
// Model:
//   y_t ~ Binomial(n_t, inv_logit(theta_t))
//   theta_t = theta_{t-1} + epsilon_t,   epsilon_t ~ Normal(0, sigma)
//
// Forecast for future years uses the same random-walk transition with the
// posterior sigma; uncertainty grows with horizon as 2030 - 2024 = 6 years out.
// The mean projection matches the latest observed level (random-walk property);
// the credible interval widens predictably.
//
// Inputs:
//   T      : number of observed years
//   T_pred : number of years to forecast forward
//   y[T]   : observed XDR counts
//   n[T]   : observed total isolate counts
//
// Outputs (generated quantities):
//   p_hist[T]   : posterior prevalence at each historical year (smoothed fit)
//   p_pred[T_pred] : projected prevalence at each future year (with growing UI)

data {
  int<lower=1> T;
  int<lower=1> T_pred;
  array[T] int<lower=0> n;
  array[T] int<lower=0> y;
}

parameters {
  real theta0;             // initial logit-prevalence
  vector[T-1] z;           // standardized innovations (non-centered)
  real<lower=0> sigma;     // random-walk step SD on the logit scale
}

transformed parameters {
  vector[T] theta;
  theta[1] = theta0;
  for (t in 2:T)
    theta[t] = theta[t-1] + sigma * z[t-1];
}

model {
  // Weakly informative priors
  theta0 ~ normal(0, 2);     // logit-scale; allows prevalence ~ 1% to 99%
  z ~ std_normal();
  sigma ~ normal(0, 0.5);    // half-normal on RW step (years-scale)

  // Likelihood
  for (t in 1:T)
    y[t] ~ binomial_logit(n[t], theta[t]);
}

generated quantities {
  vector[T] p_hist;
  vector[T_pred] p_pred;
  vector[T_pred] theta_pred;

  for (t in 1:T) p_hist[t] = inv_logit(theta[t]);

  // Forecast: continue the random walk from theta[T] with the posterior sigma.
  real theta_t = theta[T];
  for (h in 1:T_pred) {
    theta_t = theta_t + normal_rng(0, sigma);
    theta_pred[h] = theta_t;
    p_pred[h] = inv_logit(theta_t);
  }
}
