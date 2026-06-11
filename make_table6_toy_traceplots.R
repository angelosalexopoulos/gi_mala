rm(list = ls())
library(MASS)

# ============================================================
# make_table6_toy_traceplots.R  --  TOY companion to make_table6_logistic.R
#
# Runs the same GI-MALA sampler used in make_table6_logistic.R on the
# Bayesian logistic-regression posteriors (Heart, Australian) and plots
# the log-likelihood (log-target) trace. This is a quick smoke-test: it
# drops the variance-reduction control variates and the 100-rep loop, and
# just records log p(Y | x) along one chain per dataset.
#
# OUTPUT
#   make_table6_loglik_traceplots.pdf   (one panel per dataset)
# ============================================================

set.seed(1234)

data_dir <- "/Users/angelos/Desktop/vr_pCN/data"   # same as make_table6_logistic.R
out_pdf  <- file.path(
  "/Users/angelos/Library/CloudStorage/OneDrive-aueb.gr/gi-mala-new",
  "make_table6_loglik_traceplots.pdf")

n_burnin <- 2000
n_sample <- 3000

datasets <- list(
  list(name = "Heart",      sep = ","),
  list(name = "Australian", sep = "")
)

sigmoid     <- function(x) exp(x) / (1 + exp(x))
log_sigmoid <- function(x) {
  y <- numeric(length(x)); pos <- x > 0
  y[ pos] <- -log(1 + exp(-x[ pos]))
  y[!pos] <-  x[!pos] - log(1 + exp(x[!pos]))
  y
}

# ------------------------------------------------------------------
# GI-MALA on the logistic posterior; records the log-likelihood trace.
# (Same proposal as run_VR_chain() in make_table6_logistic.R, with
#  gamma adapted during burn-in to target ~0.774 acceptance; CVs removed.)
# ------------------------------------------------------------------
gimala_logistic_trace <- function(myX, Y, n_burnin, n_sample) {
  d  <- ncol(myX)
  tY <- 2 * Y - 1

  log_target <- function(z) {
    Xz   <- myX %*% z
    list(log_lik = sum(log_sigmoid(tY * Xz)),
         grad    = t(Y - sigmoid(Xz)) %*% myX)
  }

  fit      <- suppressWarnings(glm(Y ~ myX - 1, family = binomial("logit")))
  R        <- t(chol(vcov(fit)))      # preconditioner Cholesky
  xcurrent <- fit$coefficients
  gamma_tune <- 1
  target_acc <- 0.774

  n_total <- n_burnin + n_sample
  loglik  <- numeric(n_total)
  acc_hist <- numeric(n_total)

  for (i in seq_len(n_total)) {
    lp     <- log_target(xcurrent)
    fx     <- lp$log_lik
    grad_x <- lp$grad
    A_gx   <- R %*% (t(R) %*% t(grad_x))

    y    <- xcurrent + gamma_tune * A_gx +
            sqrt(2 * gamma_tune - gamma_tune^2) * (R %*% rnorm(d))
    h_yx <- (1 / (2 - gamma_tune)) *
            sum((y - xcurrent - 0.5 * gamma_tune * A_gx) * t(grad_x))

    lp_y   <- log_target(y)
    fy     <- lp_y$log_lik
    grad_y <- lp_y$grad
    A_gy   <- R %*% (t(R) %*% t(grad_y))
    h_xy   <- (1 / (2 - gamma_tune)) *
              sum((xcurrent - y - 0.5 * gamma_tune * A_gy) * t(grad_y))

    log_a  <- (fy - fx) + (h_xy - h_yx)
    accept <- (log(runif(1)) < log_a)
    if (accept) xcurrent <- y

    acc_hist[i] <- as.numeric(accept)
    loglik[i]   <- log_target(xcurrent)$log_lik

    # adapt gamma during burn-in toward target acceptance
    if (i <= n_burnin && i %% 200 == 0 && i > 200) {
      avg_acc <- mean(acc_hist[(i - 199):i])
      if (!(avg_acc >= 0.7 && avg_acc <= 0.8)) {
        if (avg_acc < target_acc)
          gamma_tune <- gamma_tune / (1 + 0.1 * (target_acc - avg_acc) / target_acc)
        else
          gamma_tune <- gamma_tune * (1 + 0.1 * (avg_acc - target_acc) / target_acc)
      }
    }
  }
  list(loglik = loglik,
       acc    = mean(acc_hist[(n_burnin + 1):n_total]),
       n_burnin = n_burnin)
}

# ------------------------------------------------------------------
# Load each dataset, run one chain, draw the log-likelihood traceplots
# ------------------------------------------------------------------
pdf(out_pdf, width = 10, height = 4.5)
op <- par(mfrow = c(1, 2), mar = c(4, 4.4, 3, 1))

for (ds in datasets) {
  sep_char <- if (ds$sep == "") "\t" else ds$sep
  raw_X <- read.table(file.path(data_dir, paste0("Logistic_", ds$name), "myX.txt"),
                      sep = sep_char, header = FALSE)
  myX <- matrix(as.numeric(unlist(raw_X)), nrow = nrow(raw_X), ncol = ncol(raw_X))
  raw_Y <- read.table(file.path(data_dir, paste0("Logistic_", ds$name), "Y.txt"),
                      sep = sep_char, header = FALSE)
  Y <- as.numeric(unlist(raw_Y))
  d <- ncol(myX); myX <- myX[, c(d, seq_len(d - 1))]

  res <- gimala_logistic_trace(myX, Y, n_burnin, n_sample)
  cat(sprintf("%-11s : N = %d, d = %d, post-burnin acc = %.3f\n",
              ds$name, nrow(myX), ncol(myX), res$acc))

  n_total <- length(res$loglik)
  plot(seq_len(n_total), res$loglik, type = "l", col = "black",
       xlab = "Iteration", ylab = "log-likelihood  log p(Y | x)",
       main = sprintf("%s  (N=%d, d=%d, acc=%.2f)",
                      ds$name, nrow(myX), ncol(myX), res$acc))
  abline(v = res$n_burnin, col = "red", lty = 2)   # burn-in end
}

par(op)
dev.off()
cat("\nWritten:", out_pdf, "\n")
