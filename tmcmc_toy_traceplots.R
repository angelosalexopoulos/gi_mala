rm(list = ls())

# ============================================================
# tmcmc_toy_traceplots.R  --  TOY version of tmcmc.R
#
# Runs the GI-MALA sampler on a univariate Student-t target (same drift
# as mala_t in tmcmc.R) and plots the log-likelihood (log-target) trace
# for a few nu values, tuned to a target acceptance rate.
#
# NOTE on acceptance: the native GI-MALA proposal is near-exact on these
# (near-)Gaussian-ish targets, so its acceptance is very high (>0.89). To
# hit a requested acceptance we inflate the proposal standard deviation by
# a factor `scale` (bigger scale -> bigger steps -> lower acceptance). The
# GI drift and a consistent MH correction are kept, so it stays valid MCMC.
#
# `scale` is pre-calibrated per nu by bisection (acceptance is monotone in
# scale), which is robust -- unlike online adaptation it cannot collapse.
#
# OUTPUT
#   tmcmc_loglik_traceplots.pdf   (one panel per nu)
# ============================================================

set.seed(1234)

# --- TOY settings ---
nu_list    <- c(1, 5, 30, 100)  # degrees of freedom to show
n_samples  <- 3000              # plotted chain length
gam_gi     <- 0.85              # GI-MALA step (fixed; native drift)
x0         <- 5                 # start off in the tail
acc_target <- 0.85              # aim for ~0.85 acceptance

out_pdf <- file.path(
  "/Users/angelos/Library/CloudStorage/OneDrive-aueb.gr/gi-mala-new",
  "tmcmc_loglik_traceplots.pdf")

# ------------------------------------------------------------------
# One GI-MALA chain on Student-t(nu) at a fixed proposal scale.
# Returns the log-target trace and the acceptance rate.
# ------------------------------------------------------------------
run_chain <- function(n_samples, nu, gam_gi, scale, x0 = 0) {
  log_p      <- function(x) -((nu + 1) / 2) * log(1 + x^2 / nu)
  grad_log_p <- function(x) -((nu + 1) * x) / (nu + x^2)

  precond <- nu / (nu + 1)
  prop_sd <- scale * sqrt(precond * (2 * gam_gi - gam_gi^2))

  loglik <- numeric(n_samples)
  acc    <- 0
  x <- x0
  for (i in seq_len(n_samples)) {
    grad   <- grad_log_p(x)
    mu     <- x + gam_gi * grad * precond
    x_prop <- rnorm(1, mean = mu, sd = prop_sd)

    grad_prop <- grad_log_p(x_prop)
    mu_prop   <- x_prop + gam_gi * grad_prop * precond

    log_q_xp_x <- dnorm(x_prop, mean = mu,      sd = prop_sd, log = TRUE)
    log_q_x_xp <- dnorm(x,      mean = mu_prop, sd = prop_sd, log = TRUE)

    log_alpha <- (log_p(x_prop) + log_q_x_xp) - (log_p(x) + log_q_xp_x)
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      x <- x_prop; acc <- acc + 1
    }
    loglik[i] <- log_p(x)
  }
  list(loglik = loglik, acc = acc / n_samples)
}

# ------------------------------------------------------------------
# Calibrate `scale` so acceptance ~ acc_target, by bisection on log(scale).
# `scale` multiplies only the proposal NOISE (not the GI drift), so
# acceptance is monotone decreasing in scale ONLY for scale >= 1 (the
# native setting). Below 1 the proposal becomes a near-deterministic
# gradient jump and acceptance collapses, so we restrict the bracket to
# [1, 50] where bisection is valid.
# ------------------------------------------------------------------
calibrate_scale <- function(nu, gam_gi, acc_target, x0,
                            n_pilot = 1500, n_iter = 20) {
  lo <- 1; hi <- 50               # acc high at lo (native), low at hi
  # If even the native setting accepts below target, we can't lower it.
  if (run_chain(n_pilot, nu, gam_gi, scale = lo, x0 = x0)$acc <= acc_target)
    return(lo)
  for (k in seq_len(n_iter)) {
    mid <- sqrt(lo * hi)          # geometric midpoint
    a   <- run_chain(n_pilot, nu, gam_gi, scale = mid, x0 = x0)$acc
    if (a > acc_target) lo <- mid else hi <- mid   # too accepting -> bigger steps
  }
  sqrt(lo * hi)
}

# ------------------------------------------------------------------
# Calibrate per nu, run the chain, draw the log-likelihood traceplots
# ------------------------------------------------------------------
pdf(out_pdf, width = 9, height = 7)
op <- par(mfrow = c(2, 2), mar = c(4, 4.2, 3, 1))

for (nu in nu_list) {
  sc  <- calibrate_scale(nu, gam_gi, acc_target, x0)
  res <- run_chain(n_samples, nu = nu, gam_gi = gam_gi, scale = sc, x0 = x0)
  cat(sprintf("nu = %4d : calibrated scale = %.2f,  acc = %.3f\n",
              nu, sc, res$acc))

  plot(seq_len(n_samples), res$loglik, type = "l", col = "black",
       xlab = "Iteration", ylab = "log-target  log p(x)",
       main = sprintf("Student-t, nu = %d  (acc = %.2f)", nu, res$acc))
}

par(op)
dev.off()
cat("\nWritten:", out_pdf, "\n")
