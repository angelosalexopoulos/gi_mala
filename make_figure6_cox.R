rm(list = ls())

# ---- install missing packages ----------------------------------------
pkgs_needed <- c("R.matlab", "ggplot2", "gridExtra")
pkgs_new    <- pkgs_needed[!pkgs_needed %in% installed.packages()[, "Package"]]
if (length(pkgs_new)) install.packages(pkgs_new)
# ----------------------------------------------------------------------

library(R.matlab)
library(ggplot2)
library(gridExtra)

# ============================================================
# make_figure6_cox.R  --  Figure 6
#
# PURPOSE
#   Produces Figure 6 for the GI-MALA paper:
#     Panel 1 (left)  : trace plots of log-likelihood from rep 1
#                       for all competing methods on the
#                       Log-Gaussian Cox process (64×64 grid).
#     Panel 2 (right) : boxplots of the mean log-likelihood
#                       over 10 independent repetitions.
#
# METHODS
#   GI-MALA, pMALA(M), mGrad, Ellipt, pCNL, pCN
#
# INPUT
#   results/Cox_regression/
#     logGaussianCoxGirolami_Marg_repeat<r>_gimala.mat
#     logGaussianCoxGirolami_Marg_repeat<r>_Marg_pMALA.mat
#     logGaussianCoxGirolami_Marg_repeat<r>.mat        (mGrad,  cell array)
#     logGaussianCoxGirolami_Ellipt_repeat<r>.mat      (Ellipt, cell array)
#     logGaussianCoxGirolami_pCNL_repeat<r>.mat        (pCNL,  cell array)
#     logGaussianCoxGirolami_pCN_repeat<r>.mat         (pCN,   cell array)
#
# NOTE
#   The mGrad/Ellipt/pCN/pCNL files save results inside a MATLAB
#   cell array ({1} indexed).  The find_field() helper handles
#   both scalar-struct and cell-of-struct layouts automatically.
#
# OUTPUT
#   Figure6.pdf  saved next to this script
# ============================================================

script_dir  <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
results_dir <- file.path(script_dir, "results", "Cox_regression")
fig6_path   <- file.path(script_dir, "Figure6.pdf")

Repeats <- 10

# ---- helper: deep-search a named field in nested R.matlab output -----
find_field <- function(obj, field) {
  if (is.null(obj)) return(NULL)
  if (is.list(obj)) {
    if (!is.null(names(obj)) && field %in% names(obj)) {
      v <- obj[[field]]
      while (is.list(v) && length(v) == 1) v <- v[[1]]
      return(as.vector(v))
    }
    for (item in obj) {
      res <- find_field(item, field)
      if (!is.null(res)) return(res)
    }
  }
  NULL
}

# ---- method specifications -------------------------------------------
#   flat = TRUE  : LogL stored as a top-level variable in the .mat file
#   flat = FALSE : LogL inside struct_name (scalar struct or cell{1})
methods <- list(
  list(name = "GI-MALA",  flat = TRUE,
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Marg_repeat%d_gimala.mat", r))),
  list(name = "pMALA(M)", flat = TRUE,
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Marg_repeat%d_Marg_pMALA.mat", r))),
  list(name = "mGrad",    flat = FALSE, struct_name = "summaryMarg",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Marg_repeat%d.mat", r))),
  list(name = "Ellipt",   flat = FALSE, struct_name = "summaryEllipt",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Ellipt_repeat%d.mat", r))),
  list(name = "pCNL",     flat = FALSE, struct_name = "summarypCNL",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_pCNL_repeat%d.mat", r))),
  list(name = "pCN",      flat = FALSE, struct_name = "summarypCN",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_pCN_repeat%d.mat", r)))
)

level_order <- c("pCN", "pCNL", "Ellipt", "mGrad", "pMALA(M)", "GI-MALA")
linetypes   <- c("GI-MALA"  = "solid",   "pMALA(M)" = "dashed",
                 "mGrad"    = "dotted",  "Ellipt"   = "dotdash",
                 "pCNL"     = "longdash","pCN"      = "twodash")

# ---- read LogL from a .mat file --------------------------------------
read_logL <- function(m, rep) {
  path <- m$file_fn(rep)
  if (!file.exists(path)) { warning("File not found: ", path); return(NULL) }
  mat <- readMat(path)
  if (m$flat) {
    v <- mat[["LogL"]]
    return(as.vector(Re(v)))   # take real part (Cox samples can be complex)
  }
  v <- find_field(mat[[m$struct_name]], "LogL")
  if (!is.null(v)) return(as.vector(Re(v)))
  NULL
}

# ---- trace-plot data frame (rep = 1) ---------------------------------
trace_list <- lapply(methods, function(m) {
  logL <- read_logL(m, rep = 1)
  if (is.null(logL)) return(NULL)
  data.frame(iter   = seq_along(logL),
             logL   = logL,
             method = m$name,
             stringsAsFactors = FALSE)
})
trace_df <- do.call(rbind, Filter(Negate(is.null), trace_list))
trace_df$method <- factor(trace_df$method, levels = level_order)

# ---- boxplot data frame (all reps) -----------------------------------
box_list <- lapply(methods, function(m) {
  means <- vapply(seq_len(Repeats), function(r) {
    logL <- read_logL(m, rep = r)
    if (is.null(logL)) return(NA_real_)
    mean(logL, na.rm = TRUE)
  }, numeric(1))
  data.frame(method = m$name, mean_logL = means, stringsAsFactors = FALSE)
})
box_df <- do.call(rbind, box_list)
box_df$method <- factor(box_df$method, levels = level_order)

# ---- Panel 1: trace plot --------------------------------------------
p1 <- ggplot(trace_df,
             aes(x = iter, y = logL, group = method, linetype = method)) +
  geom_line(colour = "black", linewidth = 0.45) +
  scale_linetype_manual(values = linetypes, name = NULL,
                        breaks = level_order) +
  labs(title = "Log-Gaussian Cox Process",
       x     = "Iteration",
       y     = "Log-likelihood") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.4, "cm"),
        plot.title       = element_text(size = 11))

# ---- Panel 2: boxplot -----------------------------------------------
p2 <- ggplot(box_df, aes(x = method, y = mean_logL)) +
  geom_boxplot(fill = "white", colour = "black", width = 0.5,
               outlier.size = 1) +
  labs(x = NULL, y = "Mean log-likelihood") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# ---- combine and save -----------------------------------------------
fig6 <- grid.arrange(p1, p2, ncol = 2, widths = c(2, 1))
ggsave(fig6_path, plot = fig6, width = 9, height = 4)
cat("Written:", fig6_path, "\n")
