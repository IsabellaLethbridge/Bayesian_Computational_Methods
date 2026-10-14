library(mvtnorm)

## ---------------------------------------------------------------
## 1. Type 3 log-likelihood (sigma2 is a variance)
## ---------------------------------------------------------------
log_diff_exp <- function(la, lb) la + log1p(-exp(lb - la))

Type_3_Log_Likelihood <- function(y, t, mu, sigma2, tau) {
  sd <- sqrt(sigma2)
  b  <- (t - mu) / sd
  a  <- (t - tau - mu) / sd
  
  logQ <- ifelse(
    a > 0,
    log_diff_exp(pnorm(a, lower.tail = FALSE, log.p = TRUE),
                 pnorm(b, lower.tail = FALSE, log.p = TRUE)),
    log_diff_exp(pnorm(b, log.p = TRUE),
                 pnorm(a, log.p = TRUE))
  )
  
  out <- sum(log(tau) + dnorm(t - y * tau, mu, sd, log = TRUE) - logQ)
  if (!is.finite(out)) -Inf else out
}

## ---------------------------------------------------------------
## 2. Parameter transforms: (mu, sigma2, tau) <-> real line
## ---------------------------------------------------------------
to_real <- function(theta) {
  c(qlogis((theta[1] - 100) / 60),
    log(theta[2]),
    qlogis((theta[3] - 50) / 130))
}

from_real <- function(phi) {
  c(100 + 60 * plogis(phi[1]),
    exp(phi[2]),
    50 + 130 * plogis(phi[3]))
}

## ---------------------------------------------------------------
## 3. Log posterior on the transformed scale (priors + Jacobians)
## ---------------------------------------------------------------
log_target <- function(phi, y, t) {
  theta <- from_real(phi)
  ll <- Type_3_Log_Likelihood(y, t, theta[1], theta[2], theta[3])
  if (!is.finite(ll)) return(-Inf)
  
  # Uniform priors x logit Jacobians (constants dropped)
  lp_mu  <- plogis(phi[1], log.p = TRUE) + plogis(-phi[1], log.p = TRUE)
  lp_tau <- plogis(phi[3], log.p = TRUE) + plogis(-phi[3], log.p = TRUE)
  # IG(1,1) prior on sigma2 times Jacobian sigma2
  lp_s2  <- -phi[2] - exp(-phi[2])
  
  ll + lp_mu + lp_s2 + lp_tau
}

## ---------------------------------------------------------------
## 4. Core MH loop: one step per row of `jitters` (symmetric proposals)
## ---------------------------------------------------------------
mh_run <- function(phi, jitters, y, t) {
  n   <- nrow(jitters)
  cur <- log_target(phi, y, t)
  if (!is.finite(cur)) stop("Starting value has zero density; redraw init.")
  
  phi_draws <- matrix(NA_real_, n, 3)
  accepted  <- logical(n)
  
  for (i in seq_len(n)) {
    prop      <- phi + jitters[i, ]
    prop_val  <- log_target(prop, y, t)
    log_alpha <- if (is.finite(prop_val)) min(0, prop_val - cur) else -Inf
    
    if (log(runif(1)) < log_alpha) {
      phi <- prop
      cur <- prop_val
      accepted[i] <- TRUE
    }
    phi_draws[i, ] <- phi
  }
  list(phi_draws = phi_draws, last = phi, accept_rate = mean(accepted))
}

## ---------------------------------------------------------------
## 5. Two-phase sampler
##    Phase 1: burn-in with isotropic jitter -> covariance of phi draws
##    Phase 2: all iterations retained, proposals ~ N(0, scale * Sigma)
## ---------------------------------------------------------------
MH_Type_3 <- function(y, t, init, iterations = 50000,
                      burn = 3000, burn_sd = 0.5,
                      cov_frac = 0.5,
                      scale = 2.38^2 / 3) {
  
  phi0 <- to_real(init)
  
  # Phase 1
  burn_jit <- rmvnorm(burn, mean = rep(0, 3), sigma = burn_sd^2 * diag(3))
  b <- mh_run(phi0, burn_jit, y, t)
  
  keep  <- seq.int(floor(burn * (1 - cov_frac)) + 1, burn)
  Sigma <- cov(b$phi_draws[keep, , drop = FALSE])
  if (any(!is.finite(Sigma)) || det(Sigma) < 1e-12)
    stop("Burn-in covariance is degenerate (chain barely moved?). Reduce burn_sd and retry.")
  
  # Phase 2
  main_jit <- rmvnorm(iterations, mean = rep(0, 3), sigma = scale * Sigma)
  m <- mh_run(b$last, main_jit, y, t)
  
  draws <- t(apply(m$phi_draws, 1, from_real))
  colnames(draws) <- c("mu", "sigma2", "tau")
  
  list(draws       = draws,          # original scale, retained
       phi_draws   = m$phi_draws,    # transformed scale
       Sigma       = Sigma,          # tuned covariance (transformed space)
       burn_accept = b$accept_rate,
       accept_rate = m$accept_rate)
}

## ---------------------------------------------------------------
## 6. Random starting values from the prior ranges (redraw if bad)
## ---------------------------------------------------------------
random_init <- function(y, t) {
  repeat {
    init <- c(runif(1, 100, 160), 1 / rgamma(1, 1, 1), runif(1, 50, 180))
    if (is.finite(log_target(to_real(init), y, t))) return(init)
  }
}

## ---------------------------------------------------------------
## 7. Simulated Type 3 data (replace with your own y and t)
## ---------------------------------------------------------------
set.seed(1)
mu_true <- 130; sigma2_true <- 100; tau_true <- 100
y <- numeric(0); t <- numeric(0)
while (length(y) < 200) {
  ti <- runif(1, 80, 260)                      # capture date
  Ti <- rnorm(1, mu_true, sqrt(sigma2_true))   # start of moult
  yi <- (ti - Ti) / tau_true                   # moult index
  if (yi > 0 && yi < 1) {                      # only moulting birds recorded
    y <- c(y, yi); t <- c(t, ti)
  }
}

## ---------------------------------------------------------------
## 8. Run three chains from random starts
## ---------------------------------------------------------------
n_chains <- 3
fits <- lapply(seq_len(n_chains), function(k) {
  MH_Type_3(y, t, init = random_init(y, t), iterations = 50000)
})

sapply(fits, function(f) c(burn_accept = f$burn_accept, accept = f$accept_rate))

## Quick look at the first chain
colMeans(fits[[1]]$draws)
plot(fits[[1]]$draws[, "tau"], type = "l", ylab = "tau", main = "Trace: tau (chain 1)")

library(moult)

data(sanderlings) 
dat <- sanderlings[sanderlings$MIndex!=0 & sanderlings$MIndex!=1,]
run <- MH_Type_3(dat$MIndex,dat$Day,init = c(110,150,100))

results <- data.frame(
  Parameter = c("mu", "sigma2", "tau"),
  Mean = apply(run$draws, 2, mean),
  Lower_95 = apply(run$draws, 2, quantile, probs = 0.025),
  Upper_95 = apply(run$draws, 2, quantile, probs = 0.975)
)

print(results, row.names = FALSE, digits = 4)
