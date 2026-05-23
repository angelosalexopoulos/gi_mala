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
# make_figure4_logistic.R  --  Figure 4
#
# PURPOSE
#   Produces Figure 4 for the GI-MALA paper:
#     Panel 1 (left)  : trace plots of log-likelihood from rep 1
#                       for all competing methods on the Heart
#                       GP binary classification benchmark.
#     Panel 2 (right) : boxplots of the mean log-likelihood
#                       over 10 independent repetitions.
#
# METHODS
#   GI-MALA, pMALA(M), mGrad, Ellipt, pCNL, pCN
#
# DATASET
#   Heart (N = 270, d = 270)  -- change dataset_name below.
#   Available: 'Australian', 'German', 'Heart', 'Pima', 'Ripley'
#
# INPUT
#   results/LogisticRegression_GP/
#     <Dataset>_repeat<r>_Marg_fixedhypers_gimala.mat
#     <Dataset>_repeat<r>_Marg_fixedhypers_pMALA.mat
#     <Dataset>_repeat<r>_Marg_fixedhypers.mat
#     <Dataset>_repeat<r>_Ellipt_fixedhypers.mat
#     <Dataset>_repeat<r>_pCNL_fixedhypers.mat
#     <Dataset>_repeat<r>_pCN_fixedhypers.mat
#
# OUTPUT
#   Figure4.pdf  saved next to this script
# ============================================================

script_dir   <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
results_dir  <- file.path(script_dir, "results", "LogisticRegression_GP")
fig4_path    <- file.path(script_dir, "Figure4.pdf")

dataset_name <- "Heart"   # <-- change dataset here
Repeats      <- 10

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
#   name      : display label
#   file_fn   : function(dataset, rep) -> file name
#   flat      : TRUE  -> LogL stored as flat variable in the .mat file
#               FALSE -> LogL inside struct named struct_name
#   struct_name: name of the struct variable (if flat = FALSE)

methods <- list(
  list(name = "GI-MALA",  flat = TRUE,
       file_fn = function(ds, r)
         file.path(results_dir,
                   sprintf("%s_repeat%d_Marg_fixedhypers_gimala.mat", ds, r))),
  list(name = "pMALA(M)", flat = TRUE,
       file_fn = function(ds, r)
         file.path(results_dir,
                   sprintf("%s_repeat%d_Marg_fixedhypers_pMALA.mat", ds, r))),
  list(name = "mGrad",    flat = FALSE, struct_name = "summaryMarg",
       file_fn = function(ds, r)
         file.path(results_dir,
                   sprintf("%s_repeat%d_Marg_fixedhypers.mat", ds, r))),
  list(name = "Ellipt",   flat = FALSE, struct_name = "summaryEllipt",
       file_fn = function(ds, r)
         file.path(results_dir,
                   sprintf("%s_repeat%d_Ellipt_fixedhypers.mat", ds, r))),
  list(name = "pCNL",     flat = FALSE, struct_name = "summarypCNL",
       file_fn = function(ds, r)
         file.path(results_dir,
                   sprintf("%s_repeat%d_pCNL_fixedhypers.mat", ds, r))),
  list(name = "pCN",      flat = FALSE, struct_name = "summarypCN",
       file_fn = function(ds, r)
         file.path(results_dir,
                   sprintf("%s_repeat%d_pCN_fixedhypers.mat", ds, r)))
)

method_names <- sapply(methods, `[[`, "name")

# Order for display (worst to best, roughly)
level_order  <- c("pCN", "pCNL", "Ellipt", "mGrad", "pMALA(M)", "GI-MALA")
linetypes    <- c("GI-MALA" = "solid", "pMALA(M)" = "dashed",
                  "mGrad"   = "dotted","Ellipt"   = "dotdash",
                  "pCNL"    = "longdash", "pCN"   = "twodash")

# ---- read LogL from a .mat file (rep 1) for trace plot ---------------
read_logL <- function(m, ds, rep) {
  path <- m$file_fn(ds, rep)
  if (!file.exists(path)) {
    warning("File not found: ", path); return(NULL)
  }
  mat <- readMat(path)
  if (m$flat) {
    return(as.vector(mat[["LogL"]]))
  } else {
    return(find_field(mat[[m$struct_name]], "LogL"))
  }
}

# ---- build trace-plot data frame (rep = 1) ---------------------------
trace_list <- lapply(methods, function(m) {
  logL <- read_logL(m, dataset_name, rep = 1)
  if (is.null(logL)) return(NULL)
  data.frame(iter   = seq_along(logL),
             logL   = logL,
             method = m$name,
             stringsAsFactors = FALSE)
})
trace_df <- do.call(rbind, Filter(Negate(is.null), trace_list))
trace_df$method <- factor(trace_df$method, levels = level_order)

# ---- build boxplot data frame (all reps) -----------------------------
box_list <- lapply(methods, function(m) {
  means <- vapply(seq_len(Repeats), function(r) {
    logL <- read_logL(m, dataset_name, rep = r)
    if (is.null(logL)) return(NA_real_)
    mean(logL, na.rm = TRUE)
  }, numeric(1))
  data.frame(method    = m$name,
             mean_logL = means,
             stringsAsFactors = FALSE)
})
box_df <- do.call(rbind, box_list)
box_df$method <- factor(box_df$method, levels = level_order)

# ---- Panel 1: trace plot --------------------------------------------
p1 <- ggplot(trace_df,
             aes(x = iter, y = logL, group = method, linetype = method)) +
  geom_line(colour = "black", linewidth = 0.45) +
  scale_linetype_manual(values = linetypes, name = NULL,
                        breaks = level_order) +
  labs(title   = dataset_name,
       x       = "Iteration",
       y       = "Log-likelihood") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.4, "cm"),
        plot.title = element_text(size = 11))

# ---- Panel 2: boxplot -----------------------------------------------
p2 <- ggplot(box_df, aes(x = method, y = mean_logL)) +
  geom_boxplot(fill = "white", colour = "black", width = 0.5,
               outlier.size = 1) +
  labs(x = NULL, y = "Mean log-likelihood") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# ---- combine and save -----------------------------------------------
fig4 <- grid.arrange(p1, p2, ncol = 2, widths = c(2, 1))
ggsave(fig4_path, plot = fig4, width = 9, height = 4)
cat("Written:", fig4_path, "\n")
