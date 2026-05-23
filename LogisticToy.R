rm(list = ls())

# ---- install missing packages ----------------------------------------
pkgs_needed <- c("MASS", "mvtnorm", "coda", "ggplot2", "gridExtra")
pkgs_new    <- pkgs_needed[!pkgs_needed %in% installed.packages()[, "Package"]]
if (length(pkgs_new)) install.packages(pkgs_new)
# ----------------------------------------------------------------------

library(MASS)
library(mvtnorm)
library(coda)
library(ggplot2)
library(gridExtra)

# ============================================================
# LogisticToy.R  --  Figure 3
#
# Runs pMALA and GI-MALA on the German logistic regression
# (d = 25), 10 independent repeats.
# Produces Figure 3: trace plots (panel 1) + boxplots of
# mean log-likelihood over 10 reps (panel 2).
#
# DATA PATH
#   Set data_dir to the folder containing Logistic_German/
#   with myX.txt and Y.txt inside.
#
# OUTPUT
#   Figure3.pdf  saved next to this script
# ============================================================

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

data_dir   <- "/Users/angelos/Desktop/vr_pCN/data"   # <-- edit if needed
fig3_path  <- file.path(script_dir, "Figure3.pdf")

set.seed(121088)

data_name    <- "German"
d            <- 25
initial_inds <- c(d, 1:(d - 1))

####MCMC settings
its    <- 10000
Burnin <- 5000

###auxiliary functions
sigmoid <- function(x) exp(x) / (1 + exp(x))

log_sigmoid <- function(x) {
  y        <- rep(0, length(x))
  idxPos   <- which(x >  0)
  idxNeg   <- which(x <= 0)
  y[idxPos] <- -log(1 + exp(-x[idxPos]))
  y[idxNeg] <-  x[idxNeg] - log(1 + exp(x[idxNeg]))
  y
}

log_target <- function(z) {
  n     <- length(z)
  Xz    <- myX %*% z
  log_lik <- sum(log_sigmoid(tY * Xz))
  grad    <- t(Y - sigmoid(Xz)) %*% myX
  list(log_lik = log_lik, grad = grad)
}

# ---- load data -------------------------------------------------------
myX <- matrix(
  as.numeric(unlist(read.table(
    file.path(data_dir, paste0("Logistic_", data_name), "myX.txt"),
    sep = ""
  ))),
  ncol = d
)[, initial_inds]

Y   <- as.numeric(unlist(read.table(
  file.path(data_dir, paste0("Logistic_", data_name), "Y.txt"),
  sep = ""
)))
tY  <- 2 * Y - 1

# ---- obtain MLE for proposal -----------------------------------------
fit   <- glm(Y ~ myX - 1, family = binomial("logit"))
Sprop <- vcov(fit)
R     <- t(chol(Sprop))

# ---- storage ----------------------------------------------------------
mysamples  <- array(NA, c(its, d, 10))
myESSstats <- matrix(NA, 3, 10)
samp       <- 0
results    <- list()
resultsESS <- array(NA, c(3, 10, 2))
likStore2  <- array(NA, c(its, 10, 2))

for (pMALA in c(TRUE, FALSE)) {
  samp       <- samp + 1
  acc_store  <- rep(NA, 10)

  for (rep in 1:10) {

    yx <- grad_store <- x <- matrix(NA, its, d)
    alpha <- logLik <- rep(NA, its)
    acceptance_rates <- rep(NA, Burnin)
    xcurrent         <- fit$coefficients
    gamma_tune       <- 1
    log_post_eval    <- log_target(xcurrent)
    log_post         <- log_post_eval[[1]]
    log_grad         <- log_post_eval[[2]]

    if (pMALA) {
      target_acc <- 0.574; lim_acc1 <- 0.5; lim_acc2 <- 0.6
    } else {
      target_acc <- 0.774; lim_acc1 <- 0.7; lim_acc2 <- 0.8
    }

    ii <- 0; acc <- 0
    for (i in 1:(its + Burnin)) {
      logp_grad  <- log_target(xcurrent)
      fx         <- logp_grad[[1]]
      grad_x     <- logp_grad[[2]]
      RT_grad_x  <- t(R) %*% t(grad_x)
      A_grad_x   <- R %*% RT_grad_x
      noise      <- rnorm(d, mean = 0, sd = 1)

      if (pMALA) {
        y <- xcurrent + gamma_tune * A_grad_x +
          sqrt(2 * gamma_tune) * (R %*% noise)
      } else {
        y <- xcurrent + gamma_tune * A_grad_x +
          sqrt(2 * gamma_tune - gamma_tune * gamma_tune) * (R %*% noise)
      }

      v <- y - xcurrent - 0.5 * gamma_tune * A_grad_x
      if (pMALA) {
        h_y_x <- 0.5 * sum(v * t(grad_x))
      } else {
        h_y_x <- (1 / (2 - gamma_tune)) * sum(v * t(grad_x))
      }

      logp_grad_y <- log_target(y)
      fy          <- logp_grad_y[[1]]
      grad_y      <- logp_grad_y[[2]]
      RT_grad_y   <- t(R) %*% t(grad_y)
      A_grad_y    <- R %*% RT_grad_y

      v <- xcurrent - y - 0.5 * gamma_tune * A_grad_y
      if (pMALA) {
        h_x_y <- 0.5 * sum(v * t(grad_y))
      } else {
        h_x_y <- (1 / (2 - gamma_tune)) * sum(v * t(grad_y))
      }

      delta_f  <- fy - fx
      delta_q  <- h_x_y - h_y_x
      logp_a   <- delta_f + delta_q
      a_prob   <- min(exp(logp_a), 1)
      a        <- as.integer(runif(1) < a_prob)

      xcurrent <- a * y + (1 - a) * xcurrent
      fcurrent <- a * fy + (1 - a) * fx
      grad_x_store <- a * grad_y + (1 - a) * grad_x
      if (i > Burnin) acc <- acc + a

      if (i <= Burnin) {
        acceptance_rates[i] <- a_prob
        if (i %% 200 == 0 && i > 200) {
          avg_acc <- mean(acceptance_rates[(i - 199):i])
          if (!(avg_acc >= lim_acc1 && avg_acc <= lim_acc2)) {
            if (avg_acc < target_acc) {
              gamma_tune <- gamma_tune / (1 + 0.1 * (target_acc - avg_acc) / target_acc)
            } else {
              gamma_tune <- gamma_tune * (1 + 0.1 * (avg_acc - target_acc) / target_acc)
            }
          }
        }
      }

      if (i > Burnin) {
        ii          <- ii + 1
        x[ii, ]     <- xcurrent
        yx[ii, ]    <- y - x[ii, ]
        logLik[ii]  <- fcurrent
        alpha[ii]   <- a_prob
        grad_store[ii, ] <- grad_x_store
      }
    }

    acc_store[rep]     <- acc / its
    mysamples[, , rep] <- x
    myESS              <- effectiveSize(x)
    myESSstats[, rep]  <- c(min(myESS), median(myESS), max(myESS))
    cat(c(gamma_tune, samp, rep, avg_acc), "\n")
    likStore2[, rep, samp] <- logLik
    cat(rep, "\n")
  }

  results[[samp]]       <- mysamples
  resultsESS[, , samp]  <- myESSstats
}

# ---- ESS summary table -----------------------------------------------
ESSmat2          <- matrix(NA, 2, 3)
rownames(ESSmat2) <- c("pMALA", "GI-MALA")
colnames(ESSmat2) <- c("min", "median", "max")
ESSmat2[1, ] <- apply(resultsESS[, , 1], 1, mean)
ESSmat2[2, ] <- apply(resultsESS[, , 2], 1, mean)

# =====================================================================
# Figure 3 -- two panels:
#   Left  : trace plots of log-likelihood (rep 1, post-burnin)
#   Right : boxplots of mean log-likelihood over 10 reps
# =====================================================================

method_labels <- c("pMALA", "GI-MALA")
method_lts    <- c("pMALA" = "dashed", "GI-MALA" = "solid")
method_factor <- factor(method_labels, levels = method_labels)

# --- Panel 1: trace plot (rep = 1) ------------------------------------
n_show <- its   # show all post-burnin iterations
trace_list <- lapply(seq_along(method_labels), function(mi) {
  data.frame(
    iter   = seq_len(n_show),
    logL   = likStore2[seq_len(n_show), 1, mi],
    method = method_labels[mi]
  )
})
trace_df <- do.call(rbind, trace_list)
trace_df$method <- factor(trace_df$method, levels = method_labels)

p1 <- ggplot(trace_df, aes(x = iter, y = logL,
                            group = method, linetype = method)) +
  geom_line(colour = "black", linewidth = 0.5) +
  scale_linetype_manual(values = method_lts, name = NULL) +
  labs(x = "Iteration", y = "Log-likelihood") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"))

# --- Panel 2: boxplots of mean log-likelihood over 10 reps -----------
box_list <- lapply(seq_along(method_labels), function(mi) {
  means <- apply(likStore2[, , mi], 2, mean)   # 10-vector
  data.frame(method = method_labels[mi], mean_logL = means)
})
box_df <- do.call(rbind, box_list)
box_df$method <- factor(box_df$method, levels = method_labels)

p2 <- ggplot(box_df, aes(x = method, y = mean_logL)) +
  geom_boxplot(fill = "white", colour = "black", width = 0.5) +
  labs(x = NULL, y = "Mean log-likelihood") +
  theme_bw(base_size = 11)

# --- combine and save -------------------------------------------------
fig3 <- grid.arrange(p1, p2, ncol = 2, widths = c(2, 1))
ggsave(fig3_path, plot = fig3, width = 9, height = 4)
cat("Written:", fig3_path, "\n")
