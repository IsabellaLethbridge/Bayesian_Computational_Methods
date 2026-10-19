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

library(moult)
library(coda)

data(sanderlings) 
dat <- sanderlings[sanderlings$MIndex!=0 & sanderlings$MIndex!=1,]

## ---------------------------------------------------------------
## 1. Run three chains from random (prior-based) starting values
## ---------------------------------------------------------------
set.seed(2026)
n_chains <- 3
n_iter   <- 50000
y_obs <- dat$MIndex
t_obs <- dat$Day

# Retry with a fresh start if the burn-in covariance comes out degenerate
run_chain <- function() {
  repeat {
    init <- random_init(y_obs, t_obs)
    fit  <- tryCatch(MH_Type_3(y_obs, t_obs, init = init, iterations = n_iter),
                     error = function(e) NULL)
    if (!is.null(fit)) { fit$init <- init; return(fit) }
  }
}
fits <- lapply(seq_len(n_chains), function(i) run_chain())

pars <- c("mu", "sigma2", "tau")
cols <- c("firebrick", "steelblue", "darkgreen")

## ---------------------------------------------------------------
## 2. Acceptance rates and starting values
## ---------------------------------------------------------------
accept_tab <- data.frame(
  chain       = seq_len(n_chains),
  burn_accept = sapply(fits, `[[`, "burn_accept"),
  main_accept = sapply(fits, `[[`, "accept_rate"),
  init_mu     = sapply(fits, function(f) f$init[1]),
  init_sigma2 = sapply(fits, function(f) f$init[2]),
  init_tau    = sapply(fits, function(f) f$init[3])
)
print(accept_tab, digits = 3, row.names = FALSE)

## ---------------------------------------------------------------
## 3. Convergence: Gelman-Rubin (R-hat) and effective sample size
## ---------------------------------------------------------------
chains <- mcmc.list(lapply(fits, function(f) mcmc(f$draws)))

gelman.diag(chains, autoburnin = FALSE, multivariate = FALSE)  # R-hat per parameter
effectiveSize(chains)                                          # ESS summed over chains
gelman.plot(chains, autoburnin = FALSE)                        # R-hat as iterations accumulate

## ---------------------------------------------------------------
## 4. Trace plots (all chains overlaid; thinned for plotting only)
## ---------------------------------------------------------------
idx <- seq(1, n_iter, by = 10)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
for (p in pars) {
  plot(idx, fits[[1]]$draws[idx, p], type = "l", col = cols[1],
       ylim = range(sapply(fits, function(f) range(f$draws[idx, p]))),
       xlab = "Iteration", ylab = p, main = paste("Trace:", p))
  for (k in 2:n_chains) lines(idx, fits[[k]]$draws[idx, p], col = cols[k])
}

## ---------------------------------------------------------------
## 5. Marginal posteriors (per-chain densities + pooled histogram)
## ---------------------------------------------------------------
pooled <- do.call(rbind, lapply(fits, `[[`, "draws"))

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (p in pars) {
  hist(pooled[, p], breaks = 60, freq = FALSE, col = "grey90", border = "grey70",
       main = paste("Posterior:", p), xlab = p)
  for (k in seq_len(n_chains)) lines(density(fits[[k]]$draws[, p]), col = cols[k], lwd = 1.5)
  abline(v = quantile(pooled[, p], c(0.025, 0.975)), lty = 2)
}

## ---------------------------------------------------------------
## 6. Bivariate plots: how the parameters covary
## ---------------------------------------------------------------
thin_pool <- pooled[seq(1, nrow(pooled), by = 20), ]

panel_cor <- function(x, y, ...) {
  usr <- par("usr"); on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.5, sprintf("r = %.2f", cor(x, y)), cex = 1.6)
}
panel_pts <- function(x, y, ...) points(x, y, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.15))
panel_hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr = usr))
  h <- hist(x, plot = FALSE, breaks = 40)
  par(usr = c(usr[1:2], 0, max(h$counts) * 1.5))
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts, col = "grey80")
}
pairs(thin_pool, lower.panel = panel_pts, upper.panel = panel_cor,
      diag.panel = panel_hist, labels = pars,
      main = "Pairwise joint posterior (pooled chains, thinned)")

## ---------------------------------------------------------------
## 7. Autocorrelation plots (every chain overlaid)
## ---------------------------------------------------------------
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (p in pars) {
  plot(NA, xlim = c(0, 100), ylim = c(-0.1, 1), xlab = "Lag", ylab = "ACF",
       main = paste("ACF:", p))
  abline(h = 0, lty = 3)
  for (k in seq_len(n_chains)) {
    a <- acf(fits[[k]]$draws[, p], lag.max = 100, plot = FALSE)
    lines(a$lag[, 1, 1], a$acf[, 1, 1], col = cols[k], lwd = 1.5)
  }
}

## ---------------------------------------------------------------
## 8. Summary table (pooled) with ESS
## ---------------------------------------------------------------
ess <- effectiveSize(chains)
summary_tab <- data.frame(
  Parameter = pars,
  Mean      = colMeans(pooled),
  SD        = apply(pooled, 2, sd),
  Q2.5      = apply(pooled, 2, quantile, 0.025),
  Median    = apply(pooled, 2, quantile, 0.5),
  Q97.5     = apply(pooled, 2, quantile, 0.975),
  ESS       = as.numeric(ess[pars])
)
print(summary_tab, digits = 4, row.names = FALSE)
