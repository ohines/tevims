source("R/algorithms.R")
require(tidyr)
require(tidyverse)
require(ggplot2)
require(glue)
require(SuperLearner)

RESULTS_DIR <- "Output/"
INVESTIGATION_NAME <- "actg175"

# use a parallel version of lapply with the same signature
# note does not parallelise algorithm 2 yet
lapply_func <- function(X, FUN, ...) {
  parallel::mclapply(X, FUN, ..., mc.cores = 10)
}

# Set up the data
data(ACTG175, package = "speff2trial")
covs <- c(
  "age", "wtkg", "karnof", "cd40", "cd80", "gender",
  "homo", "race", "symptom", "drugs", "str2", "hemo"
)
label_list <- list(
  cd40 = "CD4", homo = "Homosexual activity", age = "Age", wtkg = "Weight",
  str2 = "Antiretroviral history", race = "Race", hemo = "Hemophilia",
  symptom = "Symptomatic", cd80 = "CD8", karnof = "Karnofsky score",
  drugs = "IV drug use", gender = "Gender"
)
df <- as_tibble(ACTG175) %>%
  filter(arms %in% c(1, 3)) %>%
  mutate(
    A = as.numeric(arms == 1),
    Y = cd420
  ) %>%
  select(all_of(c(covs, "A", "Y")))
x <- select(df, all_of(covs))
a <- df$A
y <- df$Y
n <- nrow(df)

k_folds <- 10
set.seed(1234567)
folds <- sample(rep(1:k_folds, length.out = n), size = n)

# define covariate subsets for TE-VIM estimators
# Leave-one-out
covariate_groups <- as.list(covs)
names(covariate_groups) <- paste0(covs, "_LO")

# Keep-one-in
for (cov in covs) {
  name <- paste0(cov, "_KI")
  covariate_groups[[name]] <- setdiff(covs, cov)
}

# Fitting functions
propensity_score_constant <- mean(a)


fn_constant <- function(y_train, ...) {
  list(pred = rep_len(propensity_score_constant, length(y_train)))
}


ranger_learners <- create.Learner("SL.ranger",
  params = list(num.trees = 2000),
  tune = list(mtry = c(3, 4)),
  name_prefix = "RANGER"
)

# mtry arg must be less than covariate dimension
# let's make a version that we can use for low dimensional training
small_ranger_learners <- create.Learner("SL.ranger",
  params = list(num.trees = 2000, mtry = 1),
  name_prefix = "RANGER"
)

gam_learners <- create.Learner("SL.gam",
  tune = list(deg.gam = c(2, 3, 4)),
  name_prefix = "GAM"
)

xgb_learners <- create.Learner("SL.xgboost",
  params = list(minobspernode = 10, ntrees = 2000, shrinkage = 0.01),
  tune = list(max_depth = c(2, 3)),
  name_prefix = "XGB"
)

sl_library <- c(
  "SL.glmnet", "SL.glm",
  ranger_learners$names,
  gam_learners$names,
  xgb_learners$names
)

small_sl_library <- c(
  "SL.glmnet", "SL.glm",
  small_ranger_learners$names,
  gam_learners$names,
  xgb_learners$names
)


fit_sl <- function(y_train, x_train, x_new) {
  sl_lib <- sl_library

  low_dimensional <- ncol(x_train) < 4
  if (low_dimensional) {
    sl_lib <- small_sl_library
  }

  sl <- SuperLearner(
    Y = y_train,
    X = x_train,
    newX = x_new,
    family = gaussian(),
    cvControl = list(V = 10),
    SL.library = sl_lib
  )

  winning_algorithm <- which.min(sl$cvRisk)
  list(pred = sl$library.predict[, winning_algorithm])
}


fitfunc_outcome <- function(y_train, x_train, x_new) {
  t_learner(y_train, x_train, x_new, "A", fit_sl)
}
fitfunc_cate <- fit_sl
fit_func_ps <- fn_constant

results <- all_algorithms(
  y, a, x, folds,
  fitfunc_outcome,
  fit_func_ps,
  fitfunc_cate,
  covariate_groups
)

z_width <- 1.959964
atevte <- results$ate_vte %>%
  as_tibble() %>%
  mutate(
    ate_lb = ate - z_width * ate.std_err,
    ate_ub = ate + z_width * ate.std_err,
    ate_pval = 1 - pchisq((ate / ate.std_err)^2, df = 1),
    vte_lb = vte - z_width * vte.std_err,
    vte_ub = vte + z_width * vte.std_err,
    vte_pval = 1 - pchisq((vte / vte.std_err)^2, df = 1),
    rootvte = sqrt(vte),
  )
atevte %>%
  arrow::write_parquet(
    glue("{RESULTS_DIR}{INVESTIGATION_NAME}_atevte.parquet")
  )

res <- results$estimates %>%
  c() %>%
  unlist() %>%
  matrix(ncol = ncol(results$estimates)) %>%
  t()
colnames(res) <- rownames(results$estimates)

res <- as_tibble(res) %>%
  mutate(estimand = colnames(results$estimates)) %>%
  pivot_longer(
    cols = c(starts_with("tevim"), starts_with("std_err")),
    names_to = c(".value", "algorithm"),
    names_pattern = "(.{1,20})_(.{2})$"
  ) %>%
  rename(tevim_u = tevim, std_err_u = std_err) %>%
  pivot_longer(
    cols = c(starts_with("tevim"), starts_with("std_err")),
    names_to = c(".value", "scale"),
    names_pattern = "(.{1,20})_(.{1})$"
  )

res %>%
  arrow::write_parquet(
    glue("{RESULTS_DIR}{INVESTIGATION_NAME}_tevims.parquet")
  )

# make and save some pretty forest plots
te_vim_forest <- function(df, algo, scale = "u", x_intercept = 0) {
  df_plot <- filter(df, algorithm == !!algo, scale == !!scale) %>%
    arrange(tevim) %>%
    mutate(est = estimand) %>%
    separate_wider_delim(
      est,
      delim = "_",
      names = c("est", "strategy"),
    )
  df_plot$est <- label_list[df_plot$est]
  df_plot$strategy <- c(KI = "(KOI)", LO = "(LOO)", sep = " ")[df_plot$strategy]
  df_plot$Mode <- df_plot$strategy
  df_plot <- unite(df_plot, "estimand", c(est, strategy), sep = " ")
  df_plot$estimand <- factor(df_plot$estimand, levels = df_plot$estimand)

  ggplot(data = df_plot, aes(
    x = estimand,
    y = tevim,
    color = Mode,
    ymin = tevim - 1.96 * std_err,
    ymax = tevim + 1.96 * std_err
  )) +
    geom_pointrange(shape = 4, size = 0.4) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_hline(yintercept = x_intercept, lty = 2) +
    coord_flip() +
    labs(x = "Variable", y = "Importance") +
    theme_bw() +
    scale_color_manual(
      values = c(
        "(KOI)" = "red",
        "(LOO)" = "blue"
      )
    )
}

nice_display <- function(p1, p2, p3, p4, labels = "AUTO") {
  # https://github.com/wilkelab/cowplot/blob/master/vignettes/shared_legends.Rmd
  legend_b <- cowplot::get_legend(
    p2 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  prow <- cowplot::plot_grid(
    p1 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"),
    p3 + theme(legend.position = "none"),
    p4 + theme(legend.position = "none"),
    ncol = 2,
    align = "hv",
    axis = "b",
    labels = labels
  )
  cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
}

p_0T <- te_vim_forest(res, "0T")
p_0D <- te_vim_forest(res, "0D")
p_01 <- te_vim_forest(res, "01")
p_02 <- te_vim_forest(res, "02")
p_0T_s <- te_vim_forest(res, "0T", scale = "s", x_intercept = 1)
p_0D_s <- te_vim_forest(res, "0D", scale = "s", x_intercept = 1)
p_01_s <- te_vim_forest(res, "01", scale = "s", x_intercept = 1)
p_02_s <- te_vim_forest(res, "02", scale = "s", x_intercept = 1)

pdf(file = glue("{RESULTS_DIR}{INVESTIGATION_NAME}_applied_plot_unscaled.pdf"), width = 8, height = 8)
plt <- nice_display(p_0T, p_0D, p_01, p_02, labels = c("1A", "1B", "2A", "2B"))
print(plt)
dev.off()

pdf(file = glue("{RESULTS_DIR}{INVESTIGATION_NAME}_applied_plot_scaled.pdf"), width = 8, height = 8)
plt <- nice_display(p_0T_s, p_0D_s, p_01_s, p_02_s, labels = c("1A", "1B", "2A", "2B"))
print(plt)
dev.off()
