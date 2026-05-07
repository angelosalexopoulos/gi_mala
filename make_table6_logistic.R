rm(list = ls())
library(MASS)

# ============================================================
# make_table6_logistic.R
#
# PURPOSE
#   Reproduces Table 6 (Section 5.2): variance reduction factors
#   for Bayesian logistic regression (Heart, Australian).
#   Reports the range [min - max] across all d dimensions of
#   Var(mu_n(F)) / Var(mu_{n,G}(F)) for the GI-MALA estimator
#   with two control variates H1, H2.
#
# VARIANCE REDUCTION METHOD
#   G(x_j) = x_j / gamma  (Poisson solution for F(x) = x_j)
#   H1 = alpha * (G(y_j) - G(x_j)) = alpha * (y_j - x_j) / gamma
#   H2 = G(y_j) - E_{q(.|x)}[G_j]
#      = (y_j - x_j) / gamma - (A_x * grad_log_pi)_j
#   Optimal coefficients (theta1, theta2) solve a 2x2 linear system
#   estimated online from the Markov chain.
#
# DATASETS
#   Heart (N=270) and Australian (N=690)  [Table 6 rows]
#
# OUTPUT
#   table6_logistic.tex  written to the working directory
#
# RUNTIME WARNING
#   Runs 100 repeats x 4 chain lengths x 2 datasets.
#   T=200000 is expensive; reduce Repeats or T_list for a quick check.
#   For a smoke-test: set Repeats <- 20 and T_list <- c(1000).
#
# DATA PATH
#   Set data_dir below to the folder containing Logistic_Heart/,
#   Logistic_Australian/ subdirectories.
# ============================================================

data_dir <- "/Users/angelos/Desktop/vr_pCN/data"   # <-- edit if needed
out_file <- "table6_logistic.tex"

set.seed(121088)

Burnin  <- 5000
Repeats <- 100
T_list  <- c(1000, 10000, 50000, 200000)

datasets <- list(
  list(name = "Heart",      sep = ","),
  list(name = "Australian", sep = "")
)

sigmoid     <- function(x) exp(x) / (1 + exp(x))
log_sigmoid <- function(x) {
  y <- numeric(length(x))
  pos <- x > 0
  y[ pos] <- -log(1 + exp(-x[ pos]))
  y[!pos] <-  x[!pos] - log(1 + exp(x[!pos]))
  y
}
ergMean <- function(x) cumsum(x) / seq_along(x)

run_VR_chain <- function(myX, Y, its, Burnin, seed_offset = 0) {
  # Returns list(means = numeric(d), means_cv = numeric(d))
  # for the GI-MALA estimator with two control variates.
  set.seed(121088 + seed_offset)
  d  <- ncol(myX)
  tY <- 2 * Y - 1

  log_target <- function(z) {
    Xz      <- myX %*% z
    log_lik <- sum(log_sigmoid(tY * Xz))
    grad    <- t(Y - sigmoid(Xz)) %*% myX
    list(log_lik = log_lik, grad = grad)
  }

  fit       <- suppressWarnings(glm(Y ~ myX - 1, family = binomial("logit")))
  Sprop     <- vcov(fit)
  R         <- t(chol(Sprop))
  xcurrent  <- fit$coefficients
  gamma_tune <- 1
  target_acc <- 0.774
  acceptance_rates <- numeric(Burnin)

  x_store     <- matrix(NA, its, d)
  yx_store    <- matrix(NA, its, d)   # y - x at each accepted/rejected step
  alpha_store <- numeric(its)
  gradS_store <- matrix(NA, its, d)  # A_x * grad at each step
  ii          <- 0

  lp <- log_target(xcurrent)

  for (i in seq_len(its + Burnin)) {
    lp       <- log_target(xcurrent)
    fx       <- lp$log_lik
    grad_x   <- lp$grad
    A_grad_x <- R %*% (t(R) %*% t(grad_x))

    noise <- rnorm(d)
    y     <- xcurrent + gamma_tune * A_grad_x +
             sqrt(2 * gamma_tune - gamma_tune^2) * (R %*% noise)

    v    <- y - xcurrent - 0.5 * gamma_tune * A_grad_x
    h_yx <- (1 / (2 - gamma_tune)) * sum(v * t(grad_x))

    lp_y     <- log_target(y)
    fy       <- lp_y$log_lik
    grad_y   <- lp_y$grad
    A_grad_y <- R %*% (t(R) %*% t(grad_y))

    v    <- xcurrent - y - 0.5 * gamma_tune * A_grad_y
    h_xy <- (1 / (2 - gamma_tune)) * sum(v * t(grad_y))

    log_a  <- (fy - fx) + (h_xy - h_yx)
    a_prob <- min(exp(log_a), 1)
    accept <- (runif(1) < a_prob)

    xcurrent  <- if (accept) y else xcurrent
    grad_curr <- if (accept) grad_y else grad_x

    if (i <= Burnin) {
      acceptance_rates[i] <- a_prob
      if (i %% 200 == 0 && i > 200) {
        avg_acc <- mean(acceptance_rates[(i - 199):i])
        if (!(avg_acc >= 0.7 && avg_acc <= 0.8)) {
          if (avg_acc < target_acc)
            gamma_tune <- gamma_tune / (1 + 0.1 * (target_acc - avg_acc) / target_acc)
          else
            gamma_tune <- gamma_tune * (1 + 0.1 * (avg_acc - target_acc) / target_acc)
        }
      }
    } else {
      ii <- ii + 1
      x_store[ii, ]     <- xcurrent
      yx_store[ii, ]    <- y - xcurrent          # y - x_current (after accept/reject)
      alpha_store[ii]   <- a_prob
      # A_x * grad = Sprop %*% grad (using fixed preconditioning Sprop)
      gradS_store[ii, ] <- as.vector(Sprop %*% t(grad_curr))
    }
  }

  # ------------------------------------------------------------------
  # Compute VR-corrected mean for each dimension using two CVs
  # H1 = alpha * (y - x) / gamma
  # H2 = (y - x) / gamma - A_x * grad
  # Optimal theta from online 2x2 linear system
  # ------------------------------------------------------------------
  lambda    <- 1e-8
  tol_rcond <- 1e-12
  means    <- numeric(d)
  means_cv <- numeric(d)

  for (dd in seq_len(d)) {
    g1 <- alpha_store * (yx_store[, dd] / gamma_tune)
    g2 <- yx_store[, dd] / gamma_tune - gradS_store[, dd]

    K11 <- ergMean(g1^2)
    K12 <- ergMean(g1 * g2)
    K22 <- ergMean(g2^2)
    th1 <- ergMean(x_store[, dd] * g1) - ergMean(x_store[, dd]) * ergMean(g1)
    th2 <- ergMean(x_store[, dd] * g2) - ergMean(x_store[, dd]) * ergMean(g2)

    theta <- c(NA_real_, NA_real_)
    for (i in 20:its) {
      Kmat <- matrix(c(K11[i], K12[i], K12[i], K22[i]), 2, 2)
      Kmat <- (Kmat + t(Kmat)) / 2
      if (rcond(Kmat) < tol_rcond) {
        next
      }
      theta <- tryCatch(
        solve(Kmat + lambda * diag(2)) %*% c(th1[i], th2[i]),
        error = function(e) theta)
    }

    means[dd]    <- mean(x_store[, dd])
    means_cv[dd] <- mean(x_store[, dd]) -
                    (if (!is.na(theta[1])) theta[1] * mean(g1) else 0) -
                    (if (!is.na(theta[2])) theta[2] * mean(g2) else 0)
  }

  list(means = means, means_cv = means_cv)
}

# ------------------------------------------------------------------
# Main loop: datasets x T values x Repeats
# ------------------------------------------------------------------
ratio_min <- matrix(NA, length(datasets), length(T_list))
ratio_max <- matrix(NA, length(datasets), length(T_list))

for (di in seq_along(datasets)) {
  ds <- datasets[[di]]
  cat("Dataset:", ds$name, "\n")

  sep_char <- if (ds$sep == "") " " else ds$sep
  myX <- as.matrix(read.table(
    file.path(data_dir, paste0("Logistic_", ds$name), "myX.txt"),
    sep = sep_char))
  Y   <- as.matrix(read.table(
    file.path(data_dir, paste0("Logistic_", ds$name), "Y.txt"),
    sep = sep_char))

  d <- ncol(myX)
  myX <- myX[, c(d, seq_len(d - 1))]

  for (ti in seq_along(T_list)) {
    T <- T_list[ti]
    cat("  T =", T, "\n")

    means_mat    <- matrix(NA, Repeats, d)
    means_cv_mat <- matrix(NA, Repeats, d)

    for (rep in seq_len(Repeats)) {
      res <- run_VR_chain(myX, Y, its = T, Burnin = Burnin,
                         seed_offset = rep + (di - 1) * 10000 + (ti - 1) * 1000)
      means_mat[rep, ]    <- res$means
      means_cv_mat[rep, ] <- res$means_cv
    }

    var_plain <- apply(means_mat,    2, var)
    var_vr    <- apply(means_cv_mat, 2, var)
    ratio     <- var_plain / pmax(var_vr, 1e-12)

    ratio_min[di, ti] <- min(ratio)
    ratio_max[di, ti] <- max(ratio)
    cat("    ratio range: [", round(ratio_min[di, ti], 2), ",",
                              round(ratio_max[di, ti], 2), "]\n")
  }
}

# ------------------------------------------------------------------
# Write LaTeX table (Table 6 format)
# ------------------------------------------------------------------
fid <- file(out_file, "w")

writeLines("% Table 6 (Section 5.2): variance reduction factors, Bayesian logistic regression", fid)
writeLines("% Rows: dataset; Columns: n values", fid)
writeLines("\\begin{tabular}{lrrrr}", fid)
writeLines("\\toprule", fid)
header <- paste(
  "Dataset",
  paste(sprintf("$n = %s$", formatC(T_list, format = "d", big.mark = ",")),
        collapse = " & "),
  sep = " & ")
writeLines(paste0(header, " \\\\"), fid)
writeLines("\\midrule", fid)

ds_labels <- c(
  Heart      = "Heart ($d=13$)",
  Australian = "Australian ($d=14$)"
)

for (di in seq_along(datasets)) {
  nm  <- datasets[[di]]$name
  lbl <- ds_labels[nm]
  vals <- paste(
    sapply(seq_along(T_list), function(ti) {
      if (!is.na(ratio_min[di, ti]))
        sprintf("%.2f-%.2f", ratio_min[di, ti], ratio_max[di, ti])
      else "--"
    }),
    collapse = " & ")
  writeLines(paste0(lbl, " & ", vals, " \\\\"), fid)
}

writeLines("\\bottomrule", fid)
writeLines("\\end{tabular}", fid)
close(fid)

cat("\nWritten:", out_file, "(Table 6, Section 5.2)\n")
