// =====================================================================
// Type 3 moult model (Underhill & Zucchini 1988) for Question (d)
// Only birds that are actively moulting are used.
//
// Likelihood for bird i (caught on day t_i with moult score y_i):
//     tau * f_T(t_i - y_i * tau) / Q(t_i)
// where f_T is the N(mu, sigma^2) density and
//     Q(t_i) = Phi((t_i - mu)/sigma) - Phi((t_i - tau - mu)/sigma)
// =====================================================================

data {
  int<lower=1> m;                       // number of moulting birds
  vector[m] t;                          // capture day of each bird
  vector<lower=0, upper=1>[m] y;        // moult index (strictly between 0 and 1)
}

parameters {
  // The bounds do two jobs: they match the uniform priors' support,
  // and they tell Stan to sample on an unconstrained (logit / log) scale
  // and add the Jacobian automatically, which is what we did by hand in (c).
  real<lower=100, upper=160> mu;        // mean start date of moult
  real<lower=0> sigma2;                 // variance of start dates
  real<lower=50, upper=180> tau;        // duration of moult (days)
}

transformed parameters {
  // Stan's normal() is parameterised by the standard deviation, not the variance
  real<lower=0> sigma = sqrt(sigma2);
}

model {
  // ---------- Priors ----------
  mu     ~ uniform(100, 160);
  sigma2 ~ inv_gamma(1, 1);
  tau    ~ uniform(50, 180);

  // ---------- Likelihood ----------
  for (i in 1:m) {
    // log Q(t_i) = log( Phi(a) - Phi(b) ), computed stably on the log scale.
    // log_diff_exp(x, z) = log(exp(x) - exp(z)) and avoids log(0) underflow.
    real logQ = log_diff_exp(normal_lcdf(t[i]       | mu, sigma),
                             normal_lcdf(t[i] - tau | mu, sigma));

    // "target +=" adds a term to the log posterior directly.
    // We need it because the Type 3 likelihood is not a built-in distribution.
    target += log(tau) + normal_lpdf(t[i] - y[i] * tau | mu, sigma) - logQ;
  }
}
