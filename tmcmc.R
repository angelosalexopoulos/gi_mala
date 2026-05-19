rm(list = ls())
library(ggplot2)

# ============================================================
# tmcmc.R -- Variance-reduction experiment for tail probability estimation
#            (Section 5.2.1 of the paper / Figure 1)
#
# PURPOSE
#   Compares the variance of a plain GI-MALA tail-probability estimator
#   E[1(X > c)] with a variance-reduced version using two control variates
#   (H1, H2 defined in Section 4 of the paper).  Computes the ratio
#   Var(plain) / Var(VR) across combinations of nu and c, for two values
#   of N (number of Poisson approximation terms): N=2 and N=5.
#
# TARGET DISTRIBUTION
#   Univariate Student-t with nu degrees of freedom.
#
# PARAMETERS (edit below to change the experiment)
#   mynu_list   degrees-of-freedom values  (default: 1,2,5,30,100,300,1000)
#   c_list      threshold values c for P(X > c)  (default: 0,1,2,3)
#   N_list      Poisson approximation terms       (default: 2, 5)
#
# OUTPUT
#   Figure1.pdf  Two-panel ggplot2 figure (left: N=2, right: N=5).
#                Black lines with different styles per threshold c.
#                Corresponds to Figure 1 (Section 5.2.1).
#   Table13.tex  Tabular version for all b values.
#                Corresponds to Table 13 (Appendix A.2).
#
# RUNTIME: approx. 5-10 minutes (100 repetitions x 10000 samples per cell).
#
# R REQUIREMENTS
#   ggplot2 (>= 3.4.0)
# ============================================================

set.seed(121088)

# ------------------------------------------------------------------
# Numerically stable log(Phi(z))
# ------------------------------------------------------------------
log_pnorm <- function(z) {
  result <- numeric(length(z))
  low <- z < -10
  result[ low] <- pnorm(z[ low], log.p = TRUE)
  result[!low] <- log(pnorm(z[!low]))
  result
}

# ------------------------------------------------------------------
# GI-MALA sampler on univariate Student-t with N-term Poisson approx.
#   N        number of approximation terms in G (beta^1 ... beta^N)
# ------------------------------------------------------------------
mala_t <- function(n_samples, nu, gam_gi, x0 = 0, N = 5, b) {

  log_p      <- function(x) -((nu + 1) / 2) * log(1 + x^2 / nu)
  grad_log_p <- function(x) -((nu + 1) * x) / (nu + x^2)

  samples      <- numeric(n_samples)
  H1           <- numeric(n_samples)
  H2           <- numeric(n_samples)
  accept_count <- 0
  Kn           <- matrix(0, 2, 2)

  beta    <- 1 - gam_gi
  idx_N   <- seq_len(N)               # 1 .. N
  idx_N1  <- c(0L, idx_N)             # 0 .. N  (for adjustments)

  x <- x0

  for (i in seq_len(n_samples)) {

    precond  <- nu / (nu + 1)
    grad     <- grad_log_p(x)
    mu       <- x + gam_gi * grad * precond
    prop_var <- 2 * gam_gi - gam_gi^2

    x_prop    <- rnorm(1, mean = mu, sd = sqrt(precond * prop_var))
    grad_prop <- grad_log_p(x_prop)
    mu_prop   <- x_prop + gam_gi * grad_prop * precond

    log_p_x      <- log_p(x)
    log_p_xp     <- log_p(x_prop)
    log_q_xp_x   <- dnorm(x_prop, mean = mu,      sd = sqrt(precond * prop_var), log = TRUE)
    log_q_x_xp   <- dnorm(x,      mean = mu_prop,  sd = sqrt(precond * prop_var), log = TRUE)

    log_alpha <- (log_p_xp + log_q_x_xp) - (log_p_x + log_q_xp_x)
    acc_ratio <- min(1, max(0, if (log_alpha >= 0) 1 else
                                if (log_alpha < -700) 0 else exp(log_alpha)))

    if (log(runif(1)) < log_alpha) {
      x    <- x_prop
      grad <- grad_prop
      accept_count <- accept_count + 1
    }
    samples[i] <- x

    # --- Control variates (N-term Poisson approximation) ---
    beta_pow  <- beta^idx_N
    beta_pow2 <- beta_pow^2

    vPhi2 <- precond * (1 - beta_pow2)
    vPhi  <- sqrt(pmax(vPhi2, .Machine$double.eps))

    z_Gx <- (beta_pow * x      - b) / vPhi
    z_Gy <- (beta_pow * x_prop - b) / vPhi
    Gx   <- as.numeric(x      > b) + sum(exp(log_pnorm(z_Gx)))
    Gy   <- as.numeric(x_prop > b) + sum(exp(log_pnorm(z_Gy)))

    muq     <- x + gam_gi * precond * grad
    s2q     <- pmax(precond * prop_var, .Machine$double.eps)
    sqrt_sq <- sqrt(s2q)

    Eqa   <- exp(log_pnorm((muq - b) / sqrt_sq))
    denom <- sqrt(pmax(vPhi2 + beta_pow2 * s2q, .Machine$double.eps))
    Eqb   <- sum(exp(log_pnorm((beta_pow * muq - b) / denom)))

    H1[i] <- acc_ratio * (Gy - Gx)
    H2[i] <- Gy - Eqa - Eqb

    Kn[1, 1] <- Kn[1, 1] + H1[i]^2
    Kn[1, 2] <- Kn[1, 2] + H1[i] * H2[i]
    Kn[2, 1] <- Kn[2, 1] + H2[i] * H1[i]
    Kn[2, 2] <- Kn[2, 2] + H2[i]^2
  }

  list(samples = samples, H1 = H1, H2 = H2,
       accept_rate = accept_count / n_samples,
       Kn = Kn / n_samples)
}

# ------------------------------------------------------------------
# Experiment parameters
# ------------------------------------------------------------------
mynu_list <- c(1, 2, 5, 30, 100, 300, 1000)
c_list    <- c(0, 1, 2, 3)
N_list    <- c(2, 5)
n_reps    <- 100
n_samples <- 10000

# ratio_array[nu, c, N_panel]
ratio_array <- array(NA,
  dim      = c(length(mynu_list), length(c_list), length(N_list)),
  dimnames = list(paste0("nu=", mynu_list), paste0("c=", c_list),
                  paste0("N=", N_list)))

for (ni in seq_along(N_list)) {
  N_val <- N_list[ni]
  cat(sprintf("\n=== N = %d ===\n", N_val))

  for (i in seq_along(mynu_list)) {
    for (j in seq_along(c_list)) {
      mynu    <- mynu_list[i]
      c_value <- c_list[j]

      tp_plain <- numeric(n_reps)
      tp_vr    <- numeric(n_reps)

      for (rep_iter in seq_len(n_reps)) {
        res <- mala_t(n_samples = n_samples, nu = mynu, gam_gi = 0.85,
                      N = N_val, b = c_value)

        tp_plain[rep_iter] <- mean(res$samples > c_value)

        vec1 <- c(mean((res$samples > c_value) * res$H1),
                  mean((res$samples > c_value) * res$H2))
        vec2 <- c(mean(res$samples > c_value) * mean(res$H1),
                  mean(res$samples > c_value) * mean(res$H2))

        betaVR    <- tryCatch(solve(res$Kn) %*% (vec1 - vec2),
                              error = function(e) matrix(c(0, 0), 2, 1))
        betaVR[2] <- 1
        betaVR[1] <- -1

        tp_vr[rep_iter] <- mean(res$samples > c_value) +
                           betaVR[2] * mean(res$H1) +
                           betaVR[1] * mean(res$H2)
      }

      ratio_array[i, j, ni] <- var(tp_plain) / max(var(tp_vr), 1e-12)
      cat(sprintf("  nu = %4d, c = %d  ->  ratio = %.2f\n",
                  mynu, c_value, ratio_array[i, j, ni]))
    }
  }
}

# ------------------------------------------------------------------
# Figure 1: two-panel ggplot2, black lines with different styles
# ------------------------------------------------------------------
line_types <- c("solid", "dashed", "dotted", "dotdash")

# Build long-format data frame
plot_df <- do.call(rbind, lapply(seq_along(N_list), function(ni) {
  do.call(rbind, lapply(seq_along(c_list), function(j) {
    data.frame(
      nu_idx  = seq_along(mynu_list),
      nu_lab  = factor(mynu_list, levels = mynu_list),
      ratio   = ratio_array[, j, ni],
      c_val   = paste0("c = ", c_list[j]),
      N_panel = paste0("N = ", N_list[ni])
    )
  }))
}))
plot_df$c_val   <- factor(plot_df$c_val,   levels = paste0("c = ", c_list))
plot_df$N_panel <- factor(plot_df$N_panel, levels = paste0("N = ", N_list))

p <- ggplot(plot_df, aes(x = nu_lab, y = ratio,
                         group = c_val, linetype = c_val)) +
  geom_line(colour = "black", linewidth = 0.8) +
  geom_point(colour = "black", size = 1.8) +
  scale_linetype_manual(values = line_types, name = NULL) +
  scale_x_discrete(name = expression(nu)) +
  scale_y_continuous(name = "Variance reduction factor") +
  facet_wrap(~ N_panel, ncol = 2) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "bottom",
    legend.key.width  = unit(1.8, "cm"),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 12, face = "bold"),
    panel.grid.minor  = element_blank()
  )

ggsave("Figure1.pdf", plot = p, width = 8, height = 4)
cat("Written: Figure1.pdf  (Figure 1, two panels N=2 and N=5)\n")

# ------------------------------------------------------------------
# Table 13: tabular version (all c values, both N panels)
# ------------------------------------------------------------------
out_file <- "Table13.tex"
fid <- file(out_file, "w")

writeLines("% Table 13 (Appendix A.2): variance-reduction ratios Var(plain)/Var(VR)", fid)
writeLines("% Rows: nu; Columns: threshold c; repeated for N=2 and N=5", fid)

for (ni in seq_along(N_list)) {
  writeLines(sprintf("\n%% --- N = %d ---", N_list[ni]), fid)
  writeLines(paste0("\\begin{tabular}{l",
                    paste(rep("c", length(c_list)), collapse = ""), "}"), fid)
  writeLines("\\toprule", fid)
  writeLines(paste0("$\\nu$ & ",
                    paste(sprintf("$c=%d$", c_list), collapse = " & "),
                    " \\\\"), fid)
  writeLines("\\midrule", fid)
  for (i in seq_along(mynu_list)) {
    writeLines(paste0(mynu_list[i], " & ",
                      paste(sprintf("%.2f", ratio_array[i, , ni]),
                            collapse = " & "),
                      " \\\\"), fid)
  }
  writeLines("\\bottomrule", fid)
  writeLines("\\end{tabular}", fid)
}

close(fid)
cat("Written: Table13.tex  (Table 13, Appendix A.2)\n")
