source("R/simulation/utils.R")
source("R/algorithms.R")
source("R/modelling_utilities.R")
source("R/simulation/true_values_experiment_1")
require(tidyr)
require(dplyr)
require(ggplot2)
require(glue)

RESULTS_DIR <- "Output/"
SIM_NAME <- "experiment0"


dgp <- function(n) {
  noise = rnorm(n)
  tibble(
    X1 = runif(n, -1, 1),
    X2 = runif(n, -1, 1),
    A = rbinom(n, 1, plogis(-0.4 * X1 + 0.1 * X1 * X2)),
    cate = (X1 * X1 * (X1 + 1.4) + 25 * X2 * X2 / 9),
    mu0 = X1 * X2 + 2 * X2 * X2 - X1,
    Y1 = mu0 + (A * cate) + noise,
    Y2 = mu0 + (A * cate / 10) + noise
  )
}


# true values for TEVIM estimands under `dgp1`
true_vals1 <- tibble(
  estimand = c("a", "a", "b", "b"),
  scale = c("u", "s", "u", "s"),  # unscaled vs scaled
  dgp = "1",
  true_value = as.numeric(true_values(1.4, 25 / 9, 1))[3:6],
)
true_vals2 <- tibble(
  estimand = c("a", "a", "b", "b"),
  scale = c("u", "s", "u", "s"),  # unscaled vs scaled
  dgp = "2",
  true_value = as.numeric(true_values(0.14, 25 / 90, 0.1))[3:6],
)
true_vals <- bind_rows(true_vals1, true_vals2)


get_estimates <- function(df, k_folds) {
  y1 <- df$Y1
  y2 <- df$Y2
  a <- df$A
  x <- df[c("X1", "X2")]
  n <- nrow(df)
  covariate_groups <- c(a = "X1", b = "X2")
  # no need to randomize folds, since data already randomized
  folds <- rep(1:k_folds, length.out = n)

  res_1 <- all_algorithms(
    y1, a, x, folds,
    wrapper_gam_gaussian,
    wrapper_gam_binomial,
    wrapper_gam_gaussian,
    covariate_groups
  )$estimates

  res_2 <- all_algorithms(
    y2, a, x, folds,
    wrapper_gam_gaussian,
    wrapper_gam_binomial,
    wrapper_gam_gaussian,
    covariate_groups
  )$estimates

  list_1 <- c(res_1)
  names(list_1) <- outer(
    rownames(res_1),
    colnames(res_1),
    paste,
    sep = "_dgp1_"
  ) %>% c()

  list_2 <- c(res_2)
  names(list_2) <- outer(
    rownames(res_2),
    colnames(res_2),
    paste,
    sep = "_dgp2_"
  ) %>% c()

  c(list(n = n), list_1, list_2)
}

# The meat of this script
res <- run_simulation(
  n_datasets = 2,
  sample_sizes = c(500, 1000, 2000, 3000, 4000, 5000),
  generate_data = dgp,
  get_estimates = get_estimates,
  sim_name = SIM_NAME,
  results_directory = RESULTS_DIR,
  k_folds = 5,
  append = TRUE
)

res_a <- res %>%
  select(matches("^tevim.*_a$")) %>%
  rename_all(~ gsub("^tevim", "order", gsub("_a$", "", .)))
res_b <- res %>%
  select(matches("^tevim.*_b$")) %>%
  rename_all(~ gsub("^tevim", "order", gsub("_b$", "", .)))
ordering <- purrr::map2_dfc(res_b, res_a, ~ .x > .y) %>%
  bind_cols(select(res, n)) %>%
  pivot_longer(
    cols = c(starts_with("order")),
    names_to = c(".value", "algorithm", "dgp"),
    names_pattern = "(.{1,20})_(.{2})_dgp(.{1})$"
  )
ordering_summary <- ordering %>%
  group_by(n, algorithm, dgp) %>%
  summarise(ordering_u = mean(order), ordering_s = mean(order_s))

# process and tabulate results
df_r <- res %>%
  pivot_longer(
    cols = c(starts_with("tevim"), starts_with("std_err")),
    names_to = c(".value", "algorithm", "dgp", "estimand"),
    names_pattern = "(.{1,20})_(.{2})_dgp(.{1})_(.{1})$"
  ) %>%
  rename(tevim_u = tevim, std_err_u = std_err) %>%
  pivot_longer(
    cols = c(starts_with("tevim"), starts_with("std_err")),
    names_to = c(".value", "scale"),
    names_pattern = "(.{1,20})_(.{1})$"
  ) %>%
  mutate(
    algorithm = factor(algorithm, levels = c("0T", "0D", "01", "02")),
  ) %>%
  left_join(true_vals, by = join_by(estimand, scale, dgp)) %>%
  mutate(
    bias = tevim - true_value,
    coverage_yn = abs(bias) <= 1.959964 * std_err,
  )

summary_stats <- df_r %>%
  group_by(n, algorithm, dgp, estimand, scale) %>%
  summarise(
    n_samples = length(bias),
    coverage = mean(coverage_yn),
    root_n_bias = mean(bias * sqrt(n)),
    n_var  = mean(bias^2 * n) - root_n_bias^2,
    root_n_bias_std_err = sqrt(n_var / n_samples),
    bias_min = root_n_bias - 1.96 * root_n_bias_std_err,
    bias_max = root_n_bias + 1.96 * root_n_bias_std_err
  )


# helper function for producing plots
sim_plots <- function(df, var_name, algs, scaled) {
  scale_code <- ifelse(scaled, "s", "u")
  df_plot <- filter(df,
    estimand == var_name,
    algorithm %in% algs,
    scale == scale_code,
  ) %>%
    mutate(sv = sqrt(n_var))

  # base plots
  bias_plt <- ggplot(data = df_plot, aes(
    x = n, y = root_n_bias, ymin = bias_min, ymax = bias_max
  )) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_pointrange(aes(color = algorithm, shape = learner))

  var_plt <- ggplot(data = df_plot, aes(x = n, y = sv)) +
    geom_point(aes(color = algorithm, shape = learner))

  cov_plt <- ggplot(data = df_plot, aes(x = n, y = coverage)) +
    geom_hline(yintercept = 0.95, lty = 2) +
    geom_point(aes(color = algorithm, shape = learner))

  out <- list(Bias = bias_plt, Variance = var_plt, Coverage = cov_plt)

  # style settings
  plot_limits <- list(
    Bias = c(min(df_plot$bias_min, 0), max(df_plot$bias_max, 0)),
    Variance = c(min(df_plot$sv), max(df_plot$sv)),
    Coverage = c(min(df_plot$coverage, 0.95), max(df_plot$coverage, 0.95))
  )

  y_labels <- list(
    Bias = "$\\sqrt{n}$(Bias)",
    Variance = "$\\sqrt{n}$(Standard Deviation)",
    Coverage = "Coverage"
  )

  # adjust plot style
  for (plt in names(out)) {
    lims <- plot_limits[[plt]]
    out[[plt]] <- out[[plt]] +
      ylim(lims[1], lims[2]) +
      xlab("n") +
      ylab(latex2exp::TeX(y_labels[[plt]])) +
      theme_bw() +
      scale_shape_manual(values = c(4, 5)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(color = "Algorithm", shape = "Learner")
  }

  return(out)
}


nice_display <- function(sim_plts, labels = "AUTO") {
  # https://github.com/wilkelab/cowplot/blob/master/vignettes/shared_legends.Rmd
  p1 <- sim_plts$Bias
  p2 <- sim_plts$Variance
  p3 <- sim_plts$Coverage

  legend_b <- cowplot::get_legend(
    p2 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  prow <- cowplot::plot_grid(
    p1 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"),
    p3 + theme(legend.position = "none"),
    ncol = 3,
    align = "hv",
    axis = "b",
    labels = labels
  )
  cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
}

# Produce and save plots
plots <- list(
  X1_cv = sim_plots(summary_stats, "a", c("01", "02"), scaled = FALSE),
  X2_cv = sim_plots(summary_stats, "b", c("01", "02"), scaled = FALSE),
  X1_nocv = sim_plots(summary_stats, "a", c("0D", "0T"), scaled = FALSE),
  X2_nocv = sim_plots(summary_stats, "b", c("0D", "0T"), scaled = FALSE),
  X1_cv_scaled = sim_plots(summary_stats, "a", c("01", "02"), scaled = TRUE),
  X2_cv_scaled = sim_plots(summary_stats, "b", c("01", "02"), scaled = TRUE),
  X1_nocv_scaled = sim_plots(summary_stats, "a", c("0D", "0T"), scaled = TRUE),
  X2_nocv_scaled = sim_plots(summary_stats, "b", c("0D", "0T"), scaled = TRUE)

)

for (name in names(plots)) {
  file_name <- glue("{RESULTS_DIR}sim_plot_{SIM_NAME}_{name}.pdf")
  plot <- nice_display(plots[[name]])
  pdf(file = file_name, width = 8, height = 3.5)
  print(plot)
  dev.off()
}
