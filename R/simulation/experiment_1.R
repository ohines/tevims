source("R/simulation/utils.R")
source("R/algorithms.R")
source("R/simulation/true_values_experiment_1.R")
require(tidyr)
require(dplyr)
require(ggplot2)
require(glue)

RESULTS_DIR <- "Output/"
SIM_NAME <- "experiment1"


dgp <- function(n) {
  noise <- rnorm(n)
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
  scale = c("u", "s", "u", "s"), # unscaled vs scaled
  dgp = "1",
  true_value = as.numeric(true_values(1.4, 25 / 9, 1))[3:6],
)
true_vals2 <- tibble(
  estimand = c("a", "a", "b", "b"),
  scale = c("u", "s", "u", "s"), # unscaled vs scaled
  dgp = "2",
  true_value = as.numeric(true_values(0.14, 25 / 90, 0.1))[3:6],
)
true_vals <- bind_rows(true_vals1, true_vals2)

# fitting function wrappers
wrapper_gam <- function(y_train, x_train, x_new, family) {
  p <- dim(x_train)[2]

  if (p == 2) {
    gam_model <- as.formula(
      "y_train~ s(X1) + s(X2) + ti(X1, X2)"
    )
  } else {
    gam_model <- as.formula(
      paste("y_train~", paste(
        paste("s(", colnames(x_train), ")", sep = ""),
        collapse = "+"
      ))
    )
  }

  fit <- mgcv::gam(gam_model, data = data.frame(x_train), family = family)
  pred <- mgcv::predict.gam(
    fit,
    newdata = data.frame(x_new),
    type = "response"
  )
  list(pred = pred)
}


gam_propensity_score <- function(y_train, x_train, x_new) {
  wrapper_gam(y_train, x_train, x_new, family = binomial())
}

gam_cate <- function(y_train, x_train, x_new) {
  wrapper_gam(y_train, x_train, x_new, family = gaussian())
}

gam_outcome <- function(y_train, x_train, x_new) {
  t_learner(y_train, x_train, x_new, "A", wrapper_gam, family = gaussian())
}


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
    gam_outcome,
    gam_propensity_score,
    gam_cate,
    covariate_groups
  )$estimates

  res_2 <- all_algorithms(
    y2, a, x, folds,
    gam_outcome,
    gam_propensity_score,
    gam_cate,
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
  n_datasets = 1000,
  sample_sizes = c(500, 1000, 2000, 3000, 4000, 5000),
  generate_data = dgp,
  get_estimates = get_estimates,
  sim_name = SIM_NAME,
  results_directory = RESULTS_DIR,
  k_folds = 8,
  append = FALSE
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
  summarise(ordering_u = mean(order), ordering_s = mean(order_s)) %>%
  pivot_longer(
    cols = starts_with("ordering"),
    names_to = c(".value", "scale"),
    names_pattern = "(.{1,20})_(.{1})$"
  )

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
    n_var = mean(bias^2 * n) - root_n_bias^2,
    root_n_bias_std_err = sqrt(n_var / n_samples),
    bias_min = root_n_bias - 1.96 * root_n_bias_std_err,
    bias_max = root_n_bias + 1.96 * root_n_bias_std_err
  )


# helper function for producing plots
sim_plots <- function(df, dgp_n, algs) {
  j <- 80
  df_plot <- filter(
    df,
    (dgp == dgp_n) & (algorithm %in% algs)
  ) %>%
    mutate(
      sv = sqrt(n_var),
      algorithm = case_match(
        algorithm,
        "0T" ~ "1A",
        "0D" ~ "1B",
        "01" ~ "2A",
        "02" ~ "2B",
      ),
    ) %>%
    unite("est", estimand:scale) %>%
    mutate(
      # add a small amount of plot jitter
      n = n + case_match(
        est,
        "a_s" ~ -2 * j,
        "a_u" ~ -j,
        "b_s" ~ +j,
        "b_u" ~ 2 * j,
      ) + case_match(
        algorithm,
        "1A" ~ -j / 2,
        "1B" ~ +j / 2,
        "2A" ~ -j / 2,
        "2B" ~ +j / 2,
      )
    )

  # base plots
  bias_plt <- ggplot(data = df_plot, aes(
    x = n, y = root_n_bias, ymin = bias_min, ymax = bias_max
  )) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(aes(color = algorithm, shape = est))
  # geom_pointrange(aes(color = algorithm, shape = est))

  var_plt <- ggplot(data = df_plot, aes(x = n, y = sv)) +
    geom_point(aes(color = algorithm, shape = est))

  cov_plt <- ggplot(data = df_plot, aes(x = n, y = coverage)) +
    geom_hline(yintercept = 0.95, lty = 2) +
    geom_point(aes(color = algorithm, shape = est))

  out <- list(Bias = bias_plt, Variance = var_plt, Coverage = cov_plt)

  # style settings
  df_cov <- filter(df, dgp == dgp_n)
  plot_limits <- list(
    # Bias = c(min(df_plot$bias_min, 0), max(df_plot$bias_max, 0)),
    Bias = c(min(df_plot$root_n_bias, 0), max(df_plot$root_n_bias, 0)),
    Variance = c(min(df_plot$sv), max(df_plot$sv)),
    Coverage = c(min(df_cov$coverage, 0.95), max(df_cov$coverage, 0.95))
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
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(color = "Algorithm", shape = "Estimand") +
      scale_color_manual(
        values = c(
          "1A" = "red",
          "1B" = "blue",
          "2A" = "darkgreen",
          "2B" = "purple"
        ),
        limits = c("1A", "1B", "2A", "2B"),
        drop = FALSE
      ) +
      scale_shape_manual(
        values = c(
          "a_s" = 3,
          "a_u" = 22,
          "b_s" = 4,
          "b_u" = 21
        ),
        labels = c(
          "a_s" = latex2exp::TeX("\\Psi_1"),
          "a_u" = latex2exp::TeX("\\Theta_1"),
          "b_s" = latex2exp::TeX("\\Psi_2"),
          "b_u" = latex2exp::TeX("\\Theta_2")
        )
      )
  }

  return(out)
}


ordering_plots <- function(df_order) {
  j <- 80
  df_plot <- filter(df_order) %>%
    mutate(
      dgp = paste0("DGP ", dgp),
      algorithm = case_match(
        algorithm,
        "0T" ~ "1A",
        "0D" ~ "1B",
        "01" ~ "2A",
        "02" ~ "2B",
      ),
      # add a small amount of plot jitter
      n = n + case_match(
        scale,
        "s" ~ +1.5 * j,
        "u" ~ -1.5 * j,
      ) + case_match(
        algorithm,
        "1A" ~ -j,
        "1B" ~ -j / 2,
        "2A" ~ +j / 2,
        "2B" ~ +j,
      )
    )

  order_plt <- ggplot(data = df_plot, aes(x = n, y = ordering)) +
    geom_point(aes(color = algorithm, shape = scale)) +
    facet_wrap(~dgp) +
    xlab("n") +
    ylab(latex2exp::TeX("Proportion rank correct")) +
    theme_bw() +
    scale_color_manual(
      values = c(
        "1A" = "red",
        "1B" = "blue",
        "2A" = "darkgreen",
        "2B" = "purple"
      ),
      limits = c("1A", "1B", "2A", "2B"),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(
        "s" = 3,
        "u" = 22
      ),
      labels = c(
        "s" = latex2exp::TeX("\\Psi"),
        "u" = latex2exp::TeX("\\Theta")
      )
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(color = "Algorithm", shape = "Scale type")
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
  dgp_1_cv = sim_plots(summary_stats, "1", c("01", "02")),
  dgp_2_cv = sim_plots(summary_stats, "2", c("01", "02")),
  dgp_1_nocv = sim_plots(summary_stats, "1", c("0T", "0D")),
  dgp_2_nocv = sim_plots(summary_stats, "2", c("0T", "0D"))
)
for (name in names(plots)) {
  file_name <- glue("{RESULTS_DIR}sim_plot_{SIM_NAME}_{name}.pdf")
  plot <- nice_display(plots[[name]])
  pdf(file = file_name, width = 8, height = 3)
  print(plot)
  dev.off()
}

order_plots <- list(
  ordering = ordering_plots(ordering_summary)
)
for (name in names(order_plots)) {
  file_name <- glue("{RESULTS_DIR}sim_plot_{SIM_NAME}_{name}.pdf")
  plot <- order_plots[[name]]
  pdf(file = file_name, width = 8, height = 3)
  print(plot)
  dev.off()
}

nice_display2 <- function(sim_plts1, sim_plts2, labels = "AUTO") {
  # https://github.com/wilkelab/cowplot/blob/master/vignettes/shared_legends.Rmd
  p1 <- sim_plts1$Bias
  p2 <- sim_plts1$Variance
  p3 <- sim_plts1$Coverage
  p4 <- sim_plts2$Bias
  p5 <- sim_plts2$Variance
  p6 <- sim_plts2$Coverage

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
    p5 + theme(legend.position = "none"),
    p6 + theme(legend.position = "none"),
    ncol = 3,
    align = "hv",
    axis = "b",
    labels = labels
  )
  cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
}

nds <- list(
  dgp1 = nice_display2(plots$dgp_1_nocv, plots$dgp_1_cv),
  dgp2 = nice_display2(plots$dgp_2_nocv, plots$dgp_2_cv)
)

for (name in names(nds)) {
  file_name <- glue("{RESULTS_DIR}sim_plot_{SIM_NAME}_{name}.pdf")
  pdf(file = file_name, width = 8, height = 7)
  print(nds[[name]])
  dev.off()
}
