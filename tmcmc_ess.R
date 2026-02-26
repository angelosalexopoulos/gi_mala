rm(list=ls())

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
#   (A) table_ess_t_selftuned.tex         : ESS, ESS/s, Acc(post), Gamma(final)
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

# Base output directory (each GI target gets its subfolder)
base_out_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/pCN_revision/t-example_new"

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

# ---------------------------
# Plot formatting
# ---------------------------
cols <- grDevices::hcl.colors(length(methods), "Dark 3")
cex_lab  <- 2.0
cex_axis <- 1.6
cex_leg  <- 1.3
lwd_line <- 2.2
ltys <- rep(c(1, 2, 3, 4, 5, 6), length.out = length(methods))
names(ltys) <- methods

# ============================================================
# Main loop over GI-MALA targets
# ============================================================
for (gi_target in gi_targets) {
  
  target_tag <- sprintf("GI_target_%0.2f", gi_target)
  out_dir <- file.path(base_out_dir, target_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tex_out         <- file.path(out_dir, "table_ess_t_selftuned.tex")
  pdf_trace_out   <- file.path(out_dir, sprintf("traceplots_t_selftuned_%s_nu1_2_5_100.pdf", target_tag))
  pdf_density_out <- file.path(out_dir, sprintf("density_vs_true_t_target_avg10reps_%s.pdf", target_tag))
  
  dens_table_out  <- file.path(out_dir, "table_density_metrics.tex")
  tail_table_out  <- file.path(out_dir, "table_tail_metrics.tex")
  
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
  # Trace plots PDF for nu=1,2,5,100 (single fresh run per nu)
  # ============================================================
  pdf(pdf_trace_out, width = 12, height = 5.6)
  par(mfrow=c(2,2))
  
  for (nu in nu_plot) {
    chains <- vector("list", length(methods))
    names(chains) <- methods
    
    for (mm in seq_along(methods)) {
      meth <- methods[mm]
      targ <- if (meth == "GI-MALA") gi_target else NULL
      
      out <- run_method_tuned(meth, n_samples=n_samples, burnin=burnin, nu=nu,
                              gamma_init=gamma_init_map[meth], x0=0,
                              target_override = targ)
      chains[[mm]] <- out$samples
    }
    
    ylim_all <- range(unlist(chains), finite = TRUE)
    
    par(mar = c(4.8, 5.8, 2.2, 9.0) + 0.1, xpd = NA,
        cex.lab = cex_lab, cex.axis = cex_axis)
    
    plot(chains[[1]], type="l", col=cols[1], lwd=lwd_line, lty=ltys[1],
         xlab="Iteration", ylab="Sample value",
         main=bquote(nu == .(nu)),
         ylim=ylim_all, xlim=c(1, n_samples))
    
    for (mm in 2:length(methods)) {
      lines(chains[[mm]], col=cols[mm], lwd=lwd_line, lty=ltys[mm])
    }
    
    legend("right", inset=c(-0.25,0),
           legend=methods, col=cols, lwd=lwd_line, lty=ltys,
           cex=cex_leg, bty="n")
  }
  dev.off()
  cat("Saved trace plots:", pdf_trace_out, "\n")
  
  # ============================================================
  # Convergence-to-target diagnostics USING ALL 10 REPLICATES:
  #  (C1) PDF: averaged KDE over 10 reps vs true dt(), for nu in nu_plot
  #  (C2) LaTeX: density errors (ISE, L1) mean(sd) across reps
  #  (C3) LaTeX: tail prob estimation mean(sd) + abs/rel errors across reps
  # ============================================================
  
  # containers for LaTeX tables (rows nu_plot, cols methods)
  ISE_mean_mat <- matrix(NA_real_, nrow=length(nu_plot), ncol=length(methods), dimnames=list(nu_plot, methods))
  ISE_sd_mat   <- ISE_mean_mat
  L1_mean_mat  <- ISE_mean_mat
  L1_sd_mat    <- ISE_mean_mat
  
  p_true_vec      <- rep(NA_real_, length(nu_plot))
  p_hat_mean_mat  <- ISE_mean_mat
  p_hat_sd_mat    <- ISE_mean_mat
  abs_err_mean_mat <- ISE_mean_mat
  rel_err_mean_mat <- ISE_mean_mat
  
  pdf(pdf_density_out, width = 12, height = 5.6)
  par(mfrow = c(2,2))
  
  for (nu in nu_plot) {
    
    ii <- match(nu, mynu_list)
    if (is.na(ii)) stop("nu_plot contains a nu not in mynu_list")
    
    # pool all post-burnin samples across methods and reps to define grid bounds
    pooled <- unlist(lapply(methods, function(m) unlist(samples_post[[ii]][[m]], use.names = FALSE)),
                     use.names = FALSE)
    
    lo <- as.numeric(quantile(pooled, 0.001, na.rm = TRUE))
    hi <- as.numeric(quantile(pooled, 0.999, na.rm = TRUE))
    grid <- seq(lo, hi, length.out = 1200)
    dx <- grid[2] - grid[1]
    f_true <- dt(grid, df = nu)
    
    fhat_rep <- array(NA_real_, dim = c(length(grid), nrep, length(methods)),
                      dimnames = list(NULL, paste0("rep",1:nrep), methods))
    
    err_ISE <- matrix(NA_real_, nrow = nrep, ncol = length(methods), dimnames = list(paste0("rep",1:nrep), methods))
    err_L1  <- err_ISE
    p_hat   <- err_ISE
    
    for (mm in seq_along(methods)) {
      meth <- methods[mm]
      for (r in 1:nrep) {
        x_post <- samples_post[[ii]][[meth]][[r]]
        
        kd <- density(x_post, from = lo, to = hi, n = length(grid))
        fhat <- kd$y
        fhat_rep[, r, meth] <- fhat
        
        ee <- density_errors_from_y(fhat, f_true, dx)
        err_ISE[r, meth] <- ee["ISE"]
        err_L1[r, meth]  <- ee["L1"]
        
        p_hat[r, meth] <- mean(x_post > thr)
      }
    }
    
    # average KDE for plotting
    fhat_mean <- sapply(methods, function(m) rowMeans(fhat_rep[,,m, drop=FALSE][,,1], na.rm=TRUE))
    
    par(mar = c(4.5, 5.2, 1.8, 9.0) + 0.1, xpd = NA,
        cex.lab = cex_lab, cex.axis = cex_axis)
    
    plot(grid, f_true, type="l", lwd=2.5, col="black",
         xlab="x", ylab="Density", main=bquote(nu == .(nu)))
    
    for (mm in seq_along(methods)) {
      lines(grid, fhat_mean[, mm], col = cols[mm], lwd = lwd_line, lty = ltys[mm])
    }
    
    legend("right", inset=c(-0.25,0),
           legend = c("True", methods),
           col = c("black", cols),
           lwd = c(2.5, rep(lwd_line, length(methods))),
           lty = c(1, ltys),
           cex = cex_leg, bty="n")
    
    # --- store summaries for LaTeX ---
    nu_idx <- match(nu, nu_plot)
    
    ISE_mean_mat[nu_idx, ] <- apply(err_ISE, 2, mean)
    ISE_sd_mat[nu_idx, ]   <- apply(err_ISE, 2, sd)
    
    L1_mean_mat[nu_idx, ] <- apply(err_L1, 2, mean)
    L1_sd_mat[nu_idx, ]   <- apply(err_L1, 2, sd)
    
    p_true <- 1 - pt(thr, df = nu)
    p_true_vec[nu_idx] <- p_true
    
    p_hat_mean_mat[nu_idx, ] <- apply(p_hat, 2, mean)
    p_hat_sd_mat[nu_idx, ]   <- apply(p_hat, 2, sd)
    
    abs_err_mean_mat[nu_idx, ] <- apply(abs(p_hat - p_true), 2, mean)
    rel_err_mean_mat[nu_idx, ] <- apply(abs(p_hat - p_true) / p_true, 2, mean)
    
    # --- console summaries (optional) ---
    cat("\n=============================\n")
    cat("Convergence-to-target summary (post-burnin), nu =", nu, "\n")
    cat("GI-MALA target =", gi_target, "\n")
    
    err_tab <- data.frame(
      Method = methods,
      ISE_mean = apply(err_ISE, 2, mean),
      ISE_sd   = apply(err_ISE, 2, sd),
      L1_mean  = apply(err_L1,  2, mean),
      L1_sd    = apply(err_L1,  2, sd)
    )
    print(round(err_tab[,-1], 6), row.names = TRUE)
    
    tail_tab <- data.frame(
      Method = methods,
      p_true = rep(p_true, length(methods)),
      p_hat_mean = apply(p_hat, 2, mean),
      p_hat_sd   = apply(p_hat, 2, sd),
      abs_err_mean = apply(abs(p_hat - p_true), 2, mean),
      rel_err_mean = apply(abs(p_hat - p_true) / p_true, 2, mean)
    )
    cat("\nTail event accuracy for P(X >", thr, "):\n")
    print(round(tail_tab[,-1], 6), row.names =TRUE)
  }
  
  dev.off()
  cat("Saved averaged density-vs-true plots:", pdf_density_out, "\n")
  
  # ============================================================
  # Write LaTeX tables: density metrics + tail metrics
  # ============================================================
  
  # --- Density metrics table (ISE and L1): mean(sd) ---
  ISE_ms <- fmt_mean_sd_mat(ISE_mean_mat, ISE_sd_mat, digits_mean=6, digits_sd=6)
  L1_ms  <- fmt_mean_sd_mat(L1_mean_mat,  L1_sd_mat,  digits_mean=6, digits_sd=6)
  
  row_labels_dens <- as.vector(t(sapply(nu_plot, function(nu) {
    c(paste0("$\\nu=", nu, "$ ISE"), paste0("$\\nu=", nu, "$ L1"))
  })))
  
  dens_char <- matrix("", nrow=2*length(nu_plot), ncol=length(methods))
  for (k in seq_along(nu_plot)) {
    dens_char[2*k - 1, ] <- ISE_ms[k, ]
    dens_char[2*k,     ] <- L1_ms[k, ]
  }
  
  latex_table_simple(
    char_mat   = dens_char,
    row_labels = row_labels_dens,
    col_labels = c("", methods),
    caption    = paste0("Convergence-to-target density error metrics for $t_{\\nu}$ targets (GI-MALA target ",
                        sprintf("%.2f", gi_target),
                        "). Entries are mean(sd) across ", nrep, " post-burnin replicates. Lower is better."),
    label      = paste0("tab:density_metrics_gi_", label_tag),
    file       = dens_table_out
  )
  cat("Saved density metrics LaTeX table:", dens_table_out, "\n")
  
  # --- Tail probability table ---
  p_hat_ms <- fmt_mean_sd_mat(p_hat_mean_mat, p_hat_sd_mat, digits_mean=6, digits_sd=6)
  
  abs_err_char <- apply(abs_err_mean_mat, 2, function(x) formatC(x, format="f", digits=6))
  rel_err_char <- apply(rel_err_mean_mat, 2, function(x) formatC(x, format="f", digits=6))
  p_true_char  <- formatC(p_true_vec, format="f", digits=6)
  
  tail_char <- cbind(
    p_true_char,
    p_hat_ms,
    abs_err_char,
    rel_err_char
  )
  
  col_labels_tail <- c(
    "",
    "$p_{true}$",
    paste0("$\\hat p$ ", methods, " (mean(sd))"),
    paste0("abs err ", methods),
    paste0("rel err ", methods)
  )
  
  row_labels_tail <- paste0("$\\nu=", nu_plot, "$")
  
  latex_table_simple(
    char_mat   = tail_char,
    row_labels = row_labels_tail,
    col_labels = col_labels_tail,
    caption    = paste0("Tail probability estimation for $p=\\Pr(X>", thr, ")$ under $t_{\\nu}$ targets (GI-MALA target ",
                        sprintf("%.2f", gi_target),
                        "). $\\hat p$ is mean(sd) across ", nrep,
                        " post-burnin replicates. Smaller absolute/relative error is better."),
    label      = paste0("tab:tail_metrics_gi_", label_tag),
    file       = tail_table_out
  )
  cat("Saved tail metrics LaTeX table:", tail_table_out, "\n")
}

cat("\nAll GI target runs completed.\n")
