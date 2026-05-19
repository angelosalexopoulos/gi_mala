rm(list=ls())
library(ggplot2)

# Resolve output directory to the script's own folder so output always
# lands next to tmcmc_ess.R regardless of R's working directory.
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

set.seed(121088)

# ============================================================
# Multi-target GI-MALA tuning experiment for univariate Student t_nu target
#
# Methods:
#   1) GI-MALA  (target acceptance varied over gi_targets)
#   2) MALA     (target fixed at 0.60)
#   3) RWM      (target fixed at 0.275)
#
# For EACH GI-MALA target, outputs are saved in a dedicated folder:
#   base_out_dir/GI_target_0.25/
#   base_out_dir/GI_target_0.35/
#   ...
#
# Outputs (per GI target folder):
#   (A) Table14.tex                        : ESS, ESS/s, Acc(post), Gamma(final)
#   (B) traceplots_*.pdf                   : traces for nu = 1,2,5,100 (single run per nu)
#   (C) density_vs_true_*.pdf              : averaged KDE (10 reps) vs true dt()
#   (D) table_density_metrics.tex          : ISE and L1 density errors (mean(sd) over reps)
#   (E) table_tail_metrics.tex             : tail prob P(X>thr) accuracy (mean(sd) + errors)
#
# Notes:
# - Robbins–Monro adapts log(gamma) during burn-in only, then gamma is frozen.
# - ESS computed via initial-positive-sequence estimator on post-burnin samples.
# ============================================================

# ---------------------------
# Settings
# ---------------------------
mynu_list <- c(1,2,5,30,100,1000)
methods   <- c("GI-MALA","MALA","RWM")

nrep       <- 10
n_samples  <- 7000
burnin     <- 3000

nu_plot <- c(1,2,5,30,100,1000)

# GI-MALA target acceptance rates to iterate over
gi_targets <- c(0.70,0.75,0.80,0.85,0.90,0.95)

# Base output directory — subfolders per GI target are created here
base_out_dir <- file.path(script_dir, "tmcmc_ess_output")

# Initial gammas (rough; tuning will adjust)
gamma_init_map <- c("GI-MALA"=0.7, "MALA"=0.5, "RWM"=0.2)

# Tail event threshold
thr <- 3

# ---------------------------
# Utilities
# ---------------------------
ess_1d <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  x <- x - mean(x)
  v <- sum(x^2) / n
  if (v == 0) return(n)
  
  if (is.null(max_lag)) max_lag <- min(2000, n - 1)
  acf_vals <- acf(x, plot = FALSE, lag.max = max_lag)$acf[-1]  # drop lag 0
  
  s <- 0
  for (k in seq(1, length(acf_vals), by = 2)) {
    if (k + 1 > length(acf_vals)) break
    pair_sum <- acf_vals[k] + acf_vals[k + 1]
    if (pair_sum <= 0) break
    s <- s + pair_sum
  }
  tau <- 1 + 2 * s
  n / tau
}

log_p_t <- function(x, nu) -((nu + 1) / 2) * log(1 + (x^2 / nu))
grad_log_p_t <- function(x, nu) -((nu + 1) * x) / (nu + x^2)

density_errors_from_y <- function(f_hat, f_true, dx) {
  ISE <- sum((f_hat - f_true)^2) * dx
  L1  <- sum(abs(f_hat - f_true)) * dx
  c(ISE = ISE, L1 = L1)
}

# --- LaTeX helpers ---
latex_table_simple <- function(char_mat, row_labels, col_labels, caption, label, file) {
  # char_mat must be a CHARACTER matrix (already formatted for LaTeX)
  stopifnot(is.matrix(char_mat))
  stopifnot(length(row_labels) == nrow(char_mat))
  stopifnot(length(col_labels) == (ncol(char_mat) + 1))
  
  lines <- c()
  lines <- c(lines, "\\begin{table}[H]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, paste0("\\caption{", caption, "}"))
  lines <- c(lines, paste0("\\label{", label, "}"))
  lines <- c(lines, paste0("\\begin{tabular}{l", paste(rep("c", ncol(char_mat)), collapse=""), "}"))
  lines <- c(lines, "\\toprule")
  lines <- c(lines, paste(col_labels, collapse=" & "), " \\\\")
  lines <- c(lines, "\\midrule")
  for (i in seq_len(nrow(char_mat))) {
    lines <- c(lines, paste0(row_labels[i], " & ", paste(char_mat[i, ], collapse=" & "), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}")
  lines <- c(lines, "\\end{table}")
  writeLines(lines, file)
}

fmt_mean_sd_mat <- function(mean_mat, sd_mat, digits_mean = 6, digits_sd = 6) {
  stopifnot(all(dim(mean_mat) == dim(sd_mat)))
  out <- matrix("", nrow=nrow(mean_mat), ncol=ncol(mean_mat))
  for (i in 1:nrow(mean_mat)) {
    for (j in 1:ncol(mean_mat)) {
      out[i,j] <- paste0(
        formatC(mean_mat[i,j], format="f", digits=digits_mean),
        " (", formatC(sd_mat[i,j], format="f", digits=digits_sd), ")"
      )
    }
  }
  out
}

# ---------------------------
# Self-tuned samplers
# ---------------------------
run_method_tuned <- function(method_name, n_samples, burnin, nu,
                             gamma_init = 0.5, x0 = 0,
                             target_override = NULL) {
  
  target_default <- switch(method_name,
                           "GI-MALA" = 0.75,
                           "MALA"    = 0.60,
                           "RWM"     = 0.275,
                           stop("Unknown method"))
  target <- if (!is.null(target_override) && method_name == "GI-MALA") target_override else target_default
  
  gamma_min <- 1e-6
  gamma_max <- if (method_name == "GI-MALA") 1.99999 else 100
  
  log_gamma <- log(gamma_init)
  
  samples <- numeric(n_samples)
  x <- x0
  
  acc_total <- 0
  acc_post  <- 0
  
  t0 <- proc.time()[["elapsed"]]
  
  for (i in 1:n_samples) {
    gamma <- exp(log_gamma)
    gamma <- min(max(gamma, gamma_min), gamma_max)
    
    # Preconditioning: keep your factor (same for GI-MALA/MALA), none for RWM
    precond <- if (method_name == "RWM") 1 else nu / (nu + 1)
    
    if (method_name == "RWM") {
      x_prop <- x + rnorm(1, 0, sd = sqrt(2 * gamma))
      log_alpha <- log_p_t(x_prop, nu) - log_p_t(x, nu)
      
    } else {
      grad <- grad_log_p_t(x, nu)
      mu <- x + gamma * grad * precond
      
      prop_var <- if (method_name == "GI-MALA") (2 * gamma - gamma^2) else (2 * gamma)
      prop_var <- pmax(prop_var, .Machine$double.eps)
      
      x_prop <- rnorm(1, mean = mu, sd = sqrt(precond * prop_var))
      
      grad_p <- grad_log_p_t(x_prop, nu)
      mu_p <- x_prop + gamma * grad_p * precond
      
      log_q_xp_given_x <- dnorm(x_prop, mean = mu,  sd = sqrt(precond * prop_var), log = TRUE)
      log_q_x_given_xp <- dnorm(x,      mean = mu_p, sd = sqrt(precond * prop_var), log = TRUE)
      
      log_alpha <- (log_p_t(x_prop, nu) + log_q_x_given_xp) - (log_p_t(x, nu) + log_q_xp_given_x)
    }
    
    accepted <- (log(runif(1)) < log_alpha)
    if (accepted) x <- x_prop
    samples[i] <- x
    
    if (accepted) acc_total <- acc_total + 1
    if (i > burnin && accepted) acc_post <- acc_post + 1
    
    if (i <= burnin) {
      eta <- 0.9 / (50 + i)^0.6
      log_gamma <- log_gamma + eta * ((accepted * 1.0) - target)
      log_gamma <- log(min(max(exp(log_gamma), gamma_min), gamma_max))
    }
  }
  
  t1 <- proc.time()[["elapsed"]]
  
  list(
    samples = samples,
    gamma_final = exp(log_gamma),
    acc_total = acc_total / n_samples,
    acc_post  = acc_post / max(1, (n_samples - burnin)),
    elapsed_sec = (t1 - t0),
    target = target
  )
}


# ============================================================
# Main loop over GI-MALA targets
# ============================================================
for (gi_target in gi_targets) {
  
  target_tag <- sprintf("GI_target_%0.2f", gi_target)
  out_dir <- file.path(base_out_dir, target_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tex_out         <- file.path(out_dir, "Table14.tex")
  pdf_density_out <- file.path(out_dir, "Figure3.pdf")
  
  label_tag <- gsub("\\.", "p", sprintf("%0.2f", gi_target))  # for LaTeX labels
  
  cat("\n====================================================\n")
  cat("Running GI-MALA target =", gi_target, "\nOutput folder:", out_dir, "\n")
  
  # ============================================================
  # Run experiments and store post-burnin samples
  # ============================================================
  ESS_arr  <- array(NA_real_, dim=c(length(mynu_list), nrep, length(methods)),
                    dimnames=list(paste0("nu=",mynu_list), paste0("rep",1:nrep), methods))
  TIME_arr <- ESS_arr
  ACC_arr  <- ESS_arr
  GAM_arr  <- ESS_arr
  
  # Store post-burnin samples: list indexed by nu -> method -> rep
  samples_post <- setNames(vector("list", length(mynu_list)), paste0("nu=", mynu_list))
  for (ii in seq_along(mynu_list)) {
    samples_post[[ii]] <- setNames(vector("list", length(methods)), methods)
    for (mm in seq_along(methods)) {
      samples_post[[ii]][[mm]] <- vector("list", nrep)
    }
  }
  
  for (ii in seq_along(mynu_list)) {
    nu <- mynu_list[ii]
    for (r in 1:nrep) {
      for (mm in seq_along(methods)) {
        meth <- methods[mm]
        targ <- if (meth == "GI-MALA") gi_target else NULL
        
        out <- run_method_tuned(meth, n_samples=n_samples, burnin=burnin, nu=nu,
                                gamma_init=gamma_init_map[meth], x0=0,
                                target_override = targ)
        
        x_post <- out$samples[(burnin+1):n_samples]
        samples_post[[ii]][[meth]][[r]] <- x_post
        
        ESS_arr[ii, r, mm]  <- ess_1d(x_post)
        TIME_arr[ii, r, mm] <- out$elapsed_sec
        ACC_arr[ii, r, mm]  <- out$acc_post
        GAM_arr[ii, r, mm]  <- out$gamma_final
      }
    }
    cat("Done nu =", nu, "\n")
  }
  
  ESSps_arr <- ESS_arr / pmax(TIME_arr, .Machine$double.eps)
  
  ESS_mean   <- apply(ESS_arr,   c(1,3), mean, na.rm=TRUE)
  ESSps_mean <- apply(ESSps_arr, c(1,3), mean, na.rm=TRUE)
  ACC_mean   <- apply(ACC_arr,   c(1,3), mean, na.rm=TRUE)
  GAM_mean   <- apply(GAM_arr,   c(1,3), mean, na.rm=TRUE)
  
  # ============================================================
  # LaTeX table: ESS, ESS/s, Acc(post), Gamma(final)
  # ============================================================
  fmt2 <- function(x) sprintf("%.2f", x)
  fmt3 <- function(x) sprintf("%.3f", x)
  fmt4 <- function(x) sprintf("%.4f", x)
  
  lines <- c()
  lines <- c(lines, "\\begin{table}[H]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, paste0("\\caption{ESS, ESS/s, post-burnin acceptance, and final step size $\\gamma$ ",
                           "for self-tuned samplers targeting a univariate Student's $t_{\\nu}$ distribution. ",
                           "GI-MALA targets acceptance ", sprintf("%.2f", gi_target),
                           "; MALA targets 0.60; RWM targets 0.275. ",
                           "Each method adapts $\\gamma$ during burn-in (", burnin, " iterations) and then freezes it. ",
                           "Results are averages over ", nrep, " independent runs of length ", n_samples, ".}"))
  lines <- c(lines, paste0("\\label{tab:ess_t_selftuned_gi_", label_tag, "}"))
  lines <- c(lines, paste0("\\begin{tabular}{ll", paste(rep("c", length(methods)), collapse=""), "}"))
  lines <- c(lines, "\\toprule")
  lines <- c(lines, paste0(" &  & ", paste(methods, collapse=" & "), " \\\\"))
  lines <- c(lines, "\\midrule")
  
  for (ii in seq_along(mynu_list)) {
    nu <- mynu_list[ii]
    lines <- c(lines,
               paste0("$\\nu=", nu, "$ & ESS      & ", paste(fmt2(ESS_mean[ii, ]),   collapse=" & "), " \\\\"))
    lines <- c(lines,
               paste0(" & ESS/s    & ", paste(fmt2(ESSps_mean[ii, ]), collapse=" & "), " \\\\"))
    lines <- c(lines,
               paste0(" & Acc(post)& ", paste(fmt3(ACC_mean[ii, ]),   collapse=" & "), " \\\\"))
    lines <- c(lines,
               paste0(" & $\\gamma$ & ", paste(fmt4(GAM_mean[ii, ]),   collapse=" & "), " \\\\"))
  }
  
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}")
  lines <- c(lines, "\\end{table}")
  
  writeLines(lines, tex_out)
  cat("Saved LaTeX ESS table:", tex_out, "\n")
  
  cat("\nMean post-burnin acceptance rates:\n")
  print(round(ACC_mean, 3))
  cat("\nMean final gammas:\n")
  print(round(GAM_mean, 4))
  
  # ============================================================
  # Figure 3: averaged KDE over nrep reps vs true dt(), for nu in nu_plot
  # Two rows x three columns panel (one panel per nu value).
  # Black lines with different styles per method.
  # ============================================================
  all_line_types <- c("True" = "solid",
                      "GI-MALA" = "dashed",
                      "MALA"    = "dotted",
                      "RWM"     = "dotdash")

  fig3_df <- do.call(rbind, lapply(nu_plot, function(nu) {
    ii <- match(nu, mynu_list)
    if (is.na(ii)) stop("nu_plot contains a nu not in mynu_list")

    pooled <- unlist(lapply(methods,
                 function(m) unlist(samples_post[[ii]][[m]], use.names = FALSE)),
                     use.names = FALSE)
    lo   <- as.numeric(quantile(pooled, 0.001, na.rm = TRUE))
    hi   <- as.numeric(quantile(pooled, 0.999, na.rm = TRUE))
    grid <- seq(lo, hi, length.out = 600)

    # True density
    rows <- data.frame(x = grid, y = dt(grid, df = nu),
                       method = "True",
                       nu_lab = paste0("nu == ", nu),
                       stringsAsFactors = FALSE)

    # Averaged KDE per method
    for (meth in methods) {
      fhat_reps <- sapply(seq_len(nrep), function(r) {
        kd <- density(samples_post[[ii]][[meth]][[r]], from = lo, to = hi, n = 600)
        kd$y
      })
      rows <- rbind(rows, data.frame(
        x      = grid,
        y      = rowMeans(fhat_reps, na.rm = TRUE),
        method = meth,
        nu_lab = paste0("nu == ", nu),
        stringsAsFactors = FALSE
      ))
    }
    rows
  }))

  all_methods_ordered <- c("True", methods)
  fig3_df$method <- factor(fig3_df$method, levels = all_methods_ordered)
  fig3_df$nu_lab <- factor(fig3_df$nu_lab,
                           levels = paste0("nu == ", nu_plot))

  p3 <- ggplot(fig3_df, aes(x = x, y = y,
                             group = method, linetype = method)) +
    geom_line(colour = "black", linewidth = 0.75) +
    scale_linetype_manual(values = all_line_types, name = NULL) +
    scale_x_continuous(name = "x") +
    scale_y_continuous(name = "Density") +
    facet_wrap(~ nu_lab, nrow = 2, ncol = 3,
               labeller = label_parsed, scales = "free") +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.key.width = unit(1.5, "cm"),
      strip.background = element_blank(),
      strip.text       = element_text(size = 11, face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(pdf_density_out, plot = p3, width = 10, height = 6)
  cat("Saved Figure3.pdf:", pdf_density_out, "\n")
}

cat("\nAll GI target runs completed.\n")
