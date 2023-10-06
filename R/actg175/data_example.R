source("R/modelling_utilities.R")
source("R/light_gbm.R")
require(tidyr)
require(tidyverse)
require(ggplot2)
require(glue)
require(gridExtra)

RESULTS_DIR <- "Output/"
SIM_NAME <- "actg175"

# Set up the data
data(ACTG175, package = "speff2trial")
covs <- c(
  "age", "wtkg", "karnof", "cd40", "cd80", "gender",
  "homo", "race", "symptom", "drugs", "str2", "hemo"
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
set.seed(123456)
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

fitfunc_outcome <- lgbm_gaussian
fitfunc_cate <- lgbm_gaussian
fit_func_ps <- wrapper_constant

results <- all_algorithms(
  y, a, x, folds,
  fitfunc_outcome,
  fit_func_ps,
  fitfunc_cate,
  covariate_groups
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

# make and save some pretty forest plots
te_vim_forest <- function(df, algo, scale = "s", x_intercept = 1) {
  df_plot <- filter(df, algorithm == !!algo, scale == !!scale) %>%
    arrange(tevim)
  df_plot$estimand <- factor(df_plot$estimand, levels = df_plot$estimand)

  ggplot(data = df_plot, aes(
    x = estimand,
    y = tevim,
    ymin = tevim - 1.96 * std_err,
    ymax = tevim + 1.96 * std_err
  )) +
    geom_pointrange() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_hline(yintercept = x_intercept, lty = 2) +
    coord_flip() +
    labs(x = "Variable", y = "Importance") +
    theme_bw()
}

p_0T <- te_vim_forest(res, "0T")
p_0D <- te_vim_forest(res, "0D")
p_01 <- te_vim_forest(res, "01")
p_02 <- te_vim_forest(res, "02")

pdf(file = glue("{RESULTS_DIR}applied_plot_nocv.pdf"), width = 8, height = 4)
plt <- cowplot::plot_grid(
  p_0T, p_0D, ncol = 2, align = "hv", axis = "b", labels = "AUTO"
)
print(plt)
dev.off()

pdf(file = glue("{RESULTS_DIR}applied_plot_cv.pdf"), width = 8, height = 4)
plt <- cowplot::plot_grid(
  p_01, p_02, ncol = 2, align = "hv", axis = "b", labels = "AUTO"
)
print(plt)
dev.off()
