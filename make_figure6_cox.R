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
#                       Log-Gaussian Cox process (64x64 grid).
#     Panel 2 (right) : boxplots of the mean log-likelihood
#                       over 10 independent repetitions.
#
# METHODS
#   GI-MALA, pMALA(M), mGrad, Ellipt, pCNL, pCN
#
# PRE-REQUISITE
#   Run patch_logL.py once before running this script.
#
# INPUT
#   results/Cox_regression/
#     logGaussianCoxGirolami_Marg_repeat<r>_gimala.mat
#     logGaussianCoxGirolami_Marg_repeat<r>_Marg_pMALA.mat
#     logGaussianCoxGirolami_Marg_repeat<r>.mat
#     logGaussianCoxGirolami_Ellipt_repeat<r>.mat
#     logGaussianCoxGirolami_pCNL_repeat<r>.mat
#     logGaussianCoxGirolami_pCN_repeat<r>.mat
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

# ---- method specifications -------------------------------------------
methods <- list(
  list(name = "GI-MALA",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Marg_repeat%d_gimala.mat", r))),
  list(name = "pMALA(M)",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Marg_repeat%d_Marg_pMALA.mat", r))),
  list(name = "mGrad",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Marg_repeat%d.mat", r))),
  list(name = "Ellipt",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_Ellipt_repeat%d.mat", r))),
  list(name = "pCNL",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_pCNL_repeat%d.mat", r))),
  list(name = "pCN",
       file_fn = function(r)
         file.path(results_dir,
           sprintf("logGaussianCoxGirolami_pCN_repeat%d.mat", r)))
)

facet_order <- c("pCN", "pCNL", "Ellipt", "mGrad", "GI-MALA", "pMALA(M)")

method_colours <- c(
  "pCN"       = "#E41A1C",
  "pCNL"      = "#FF7F00",
  "Ellipt"    = "#4DAF4A",
  "mGrad"     = "#00B3B3",
  "GI-MALA"   = "#377EB8",
  "pMALA(M)"  = "#984EA3"
)

# ---- read LogL -----------------------------------------------------------
# Re() strips imaginary numerical artefacts from Cox process results.
read_logL <- function(m, rep) {
  path <- m$file_fn(rep)
  if (!file.exists(path)) { warning("File not found: ", path); return(NULL) }
  L <- tryCatch(as.numeric(readMat(path)[["LogL"]]),
                error = function(e) { warning(basename(path), ": ", e$message); NULL })
  if (is.null(L) || length(L) == 0) {
    warning("LogL not found in ", basename(path)); return(NULL)
  }
  as.numeric(Re(L))
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
trace_df$method <- factor(trace_df$method, levels = facet_order)

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
box_df$method <- factor(box_df$method, levels = facet_order)

# ---- Panel 1: faceted trace plots (2 rows x 3 cols) -----------------
p1 <- ggplot(trace_df, aes(x = iter, y = logL, colour = method)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~ method, nrow = 2, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = method_colours, guide = "none") +
  labs(x = "Iteration", y = "Log-likelihood") +
  theme_bw(base_size = 10) +
  theme(strip.text       = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        axis.text        = element_text(size = 7),
        panel.spacing    = unit(0.4, "lines"))

# ---- Panel 2: boxplots of mean log-likelihood over 10 reps ----------
p2 <- ggplot(box_df, aes(x = method, y = mean_logL)) +
  geom_boxplot(fill = "grey90", colour = "black", width = 0.55,
               outlier.size = 1) +
  labs(x = NULL, y = "Means of log-likelihood") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# ---- combine and save -----------------------------------------------
fig6 <- grid.arrange(p1, p2, ncol = 2, widths = c(3, 2))
ggsave(fig6_path, plot = fig6, width = 11, height = 5)
cat("Written:", fig6_path, "\n")
