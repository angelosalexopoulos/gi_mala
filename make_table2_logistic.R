rm(list = ls())
library(MASS)
library(coda)

# ============================================================
# make_table2_logistic.R
#
# PURPOSE
#   Reproduces Table 2 (Section 5.1.1): ESS scores for Bayesian
#   logistic regression (GI-MALA vs. preconditioned MALA).
#
# DATASETS
#   Heart (N=270), Australian (N=690), German (N=1000)
#
# OUTPUT
#   table2_logistic.tex  written to the working directory
#   Rows: dataset; Cols: method x (Min, Median, Max ESS)
#   All numbers are averages across 10 random repeats.
#
# MCMC SETTINGS (matching Section 5.1.1)
#   its    = 10000  post-burn-in samples
#   Burnin = 5000   burn-in iterations
#   Target acceptance rate: pMALA -> 0.574, GI-MALA -> 0.774
#
# DATA PATH
#   Set data_dir below to the folder containing Logistic_Heart/,
#   Logistic_Australian/, Logistic_German/ subdirectories.
#   Each subfolder must contain myX.txt and Y.txt.
# ============================================================

data_dir <- "/Users/angelos/Desktop/vr_pCN/data"   # <-- edit if needed
out_file <- "table2_logistic.tex"

set.seed(1234)

its    <- 10000
Burnin <- 5000
Repeats <- 10

datasets <- list(
  list(name = "Heart",      sep = ",",  N = 270),
  list(name = "Australian", sep = "",   N = 690),
  list(name = "German",     sep = "",   N = 1000)
)

sigmoid     <- function(x) exp(x) / (1 + exp(x))
log_sigmoid <- function(x) {
  y <- numeric(length(x))
  pos <- x > 0
  y[ pos] <- -log(1 + exp(-x[ pos]))
  y[!pos] <-  x[!pos] - log(1 + exp(x[!pos]))
  y
}

run_chain <- function(myX, Y, its, Burnin, pMALA, seed_offset = 0) {
  set.seed(1234 + seed_offset)
  d  <- ncol(myX)
  tY <- 2 * Y - 1

  log_target <- function(z) {
    Xz       <- myX %*% z
    log_lik  <- sum(log_sigmoid(tY * Xz))
    grad     <- t(Y - sigmoid(Xz)) %*% myX
    list(log_lik = log_lik, grad = grad)
  }

  fit      <- suppressWarnings(glm(Y ~ myX - 1, family = binomial("logit")))
  Sprop    <- vcov(fit)
  R        <- t(chol(Sprop))

  xcurrent      <- fit$coefficients
  gamma_tune    <- 1
  target_acc    <- if (pMALA) 0.574 else 0.774
  lim_acc       <- if (pMALA) c(0.5, 0.6) else c(0.7, 0.8)
  acceptance_rates <- numeric(Burnin)

  lp      <- log_target(xcurrent)
  x_store <- matrix(NA, its, d)
  ii      <- 0

  for (i in seq_len(its + Burnin)) {
    lp      <- log_target(xcurrent)
    fx      <- lp$log_lik
    grad_x  <- lp$grad

    A_grad_x <- R %*% (t(R) %*% t(grad_x))
    noise    <- rnorm(d)

    if (pMALA) {
      y    <- xcurrent + gamma_tune * A_grad_x + sqrt(2 * gamma_tune) * (R %*% noise)
      v    <- y - xcurrent - 0.5 * gamma_tune * A_grad_x
      h_yx <- 0.5 * sum(v * t(grad_x))
    } else {
      y    <- xcurrent + gamma_tune * A_grad_x +
              sqrt(2 * gamma_tune - gamma_tune^2) * (R %*% noise)
      v    <- y - xcurrent - 0.5 * gamma_tune * A_grad_x
      h_yx <- (1 / (2 - gamma_tune)) * sum(v * t(grad_x))
    }

    lp_y     <- log_target(y)
    fy       <- lp_y$log_lik
    grad_y   <- lp_y$grad
    A_grad_y <- R %*% (t(R) %*% t(grad_y))

    v <- xcurrent - y - 0.5 * gamma_tune * A_grad_y
    h_xy <- if (pMALA) 0.5 * sum(v * t(grad_y)) else
            (1 / (2 - gamma_tune)) * sum(v * t(grad_y))

    log_a <- (fy - fx) + (h_xy - h_yx)
    a_prob <- min(exp(log_a), 1)
    accept <- (runif(1) < a_prob)

    xcurrent <- if (accept) y else xcurrent

    if (i <= Burnin) {
      acceptance_rates[i] <- a_prob
      if (i %% 200 == 0 && i > 200) {
        avg_acc <- mean(acceptance_rates[(i - 199):i])
        if (!(avg_acc >= lim_acc[1] && avg_acc <= lim_acc[2])) {
          if (avg_acc < target_acc)
            gamma_tune <- gamma_tune / (1 + 0.1 * (target_acc - avg_acc) / target_acc)
          else
            gamma_tune <- gamma_tune * (1 + 0.1 * (avg_acc - target_acc) / target_acc)
        }
      }
    } else {
      ii <- ii + 1
      x_store[ii, ] <- xcurrent
    }
  }
  x_store
}

# ------------------------------------------------------------------
# Main loop: datasets x methods
# ------------------------------------------------------------------
results <- list()

for (ds in datasets) {
  cat("Dataset:", ds$name, "\n")

  sep_char <- if (ds$sep == "") "\t" else ds$sep
  raw_X <- read.table(
    file.path(data_dir, paste0("Logistic_", ds$name), "myX.txt"),
    sep = sep_char, header = FALSE)
  myX <- matrix(as.numeric(unlist(raw_X)), nrow = nrow(raw_X), ncol = ncol(raw_X))
  raw_Y <- read.table(
    file.path(data_dir, paste0("Logistic_", ds$name), "Y.txt"),
    sep = sep_char, header = FALSE)
  Y <- as.numeric(unlist(raw_Y))

  d <- ncol(myX)
  # Reorder columns: put last column first (intercept convention)
  myX <- myX[, c(d, seq_len(d - 1))]

  ess_pmala   <- matrix(NA, 3, Repeats)   # rows: min/med/max
  ess_gimala  <- matrix(NA, 3, Repeats)

  for (rep in seq_len(Repeats)) {
    cat("  rep", rep, "\n")

    x_pmala  <- run_chain(myX, Y, its, Burnin, pMALA = TRUE,  seed_offset = rep)
    x_gimala <- run_chain(myX, Y, its, Burnin, pMALA = FALSE, seed_offset = rep + 100)

    ess_p <- effectiveSize(x_pmala)
    ess_g <- effectiveSize(x_gimala)

    ess_pmala[, rep]  <- c(min(ess_p),  median(ess_p),  max(ess_p))
    ess_gimala[, rep] <- c(min(ess_g),  median(ess_g),  max(ess_g))
  }

  results[[ds$name]] <- list(
    pmala  = rowMeans(ess_pmala),
    gimala = rowMeans(ess_gimala)
  )
}

# ------------------------------------------------------------------
# Write LaTeX table (Table 2 format)
# ------------------------------------------------------------------
fid <- file(out_file, "w")

writeLines("% Table 2 (Section 5.1.1): ESS in Bayesian logistic regression", fid)
writeLines("% Rows: dataset; each has two sub-rows (GI-MALA, pMALA)", fid)
writeLines("% Columns: Min ESS, Median ESS, Max ESS (averaged over 10 repeats)", fid)
writeLines("\\begin{tabular}{llrrr}", fid)
writeLines("\\toprule", fid)
writeLines("Dataset & Method & Min & Median & Max \\\\", fid)
writeLines("\\midrule", fid)

ds_labels <- list(
  Heart      = "Heart ($N=270, d=13$)",
  Australian = "Australian ($N=690, d=14$)",
  German     = "German ($N=1000, d=24$)"
)

for (ds in datasets) {
  nm  <- ds$name
  res <- results[[nm]]
  lbl <- ds_labels[[nm]]
  writeLines(sprintf("%s & GI-MALA & %.1f & %.1f & %.1f \\\\",
                     lbl, res$gimala[1], res$gimala[2], res$gimala[3]), fid)
  writeLines(sprintf("         & MALA   & %.1f & %.1f & %.1f \\\\",
                     res$pmala[1], res$pmala[2], res$pmala[3]), fid)
  writeLines("\\addlinespace", fid)
}

writeLines("\\bottomrule", fid)
writeLines("\\end{tabular}", fid)
close(fid)

cat("\nWritten:", out_file, "(Table 2, Section 5.1.1)\n")
