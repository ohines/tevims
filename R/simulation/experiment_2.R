source("R/simulation/utils.R")
source("R/algorithms.R")
source("R/simulation/true_values_experiment_2.R")
require(tidyr)
require(dplyr)
require(ggplot2)
require(glue)

RESULTS_DIR <- "Output/"
SIM_NAME <- "experiment2"


dgp3 <- function(n) {
    rho <- 0.5
    k <- sqrt(1 - rho^2)
    tibble(
        X1 = rnorm(n),
        X2 = rho * X1 + k * rnorm(n),
        X3 = rnorm(n),
        X4 = rho * X3 + k * rnorm(n),
        X5 = rnorm(n),
        X6 = rho * X5 + k * rnorm(n),
        A = rbinom(n, 1, plogis(-0.4 * X1 + 0.1 * X1 * X2 + 0.2 * X2)),
        Y = rnorm(
            n,
            X3 - X6 +
                A * (X1 + 2 * X2 + X3), # the cate
            s = 3
        )
    )
}


# specify tevim exclusion plan
covariate_groups <- list(
    # leave one out
    l1 = c("X1"),
    l2 = c("X2"),
    l3 = c("X3"),
    l4 = c("X4"),
    l5 = c("X5"),
    l6 = c("X6"),
    # keep one in
    k1 = c("X2", "X3", "X4", "X5", "X6"),
    k2 = c("X1", "X3", "X4", "X5", "X6"),
    k3 = c("X1", "X2", "X4", "X5", "X6"),
    k4 = c("X1", "X2", "X3", "X5", "X6"),
    k5 = c("X1", "X2", "X3", "X4", "X6"),
    k6 = c("X1", "X2", "X3", "X4", "X5")
)

# true values for TEVIM estimands under `dgp2`
true_vals <- tibble(
    estimand = names(covariate_groups),
    scale = "u",
    true_value = c(true_values(
        a = 1, b = 2, c = 1, rho = 0.5, scaled = FALSE
    )[1:12])
)
true_vals_scaled <- tibble(
    estimand = names(covariate_groups),
    scale = "s",
    true_value = c(true_values(
        a = 1, b = 2, c = 1, rho = 0.5, scaled = TRUE
    )[1:12])
)
true_vals <- bind_rows(true_vals, true_vals_scaled)

# fitting function wrappers
wrapper_gam <- function(y_train, x_train, x_new, family) {
    gam_model <- as.formula(
        paste("y_train~", paste(
            paste("s(", colnames(x_train), ")", sep = ""),
            collapse = "+"
        ))
    )
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


get_estimates <- function(df, k_folds, covariate_groups) {
    y <- df$Y
    a <- df$A
    x <- df[c("X1", "X2", "X3", "X4", "X5", "X6")]
    n <- nrow(df)
    # no need to randomize folds, since data already randomized
    folds <- rep(1:k_folds, length.out = n)

    res <- all_algorithms(
        y, a, x, folds,
        gam_outcome,
        gam_propensity_score,
        gam_cate,
        covariate_groups
    )$estimates

    res_list <- c(res)
    names(res_list) <- outer(
        rownames(res),
        colnames(res),
        paste,
        sep = "_xxx_"
    ) %>% c()

    c(list(n = n), res_list)
}

# The meat of this script
res <- run_simulation(
    n_datasets = 500,
    sample_sizes = c(5000),
    generate_data = dgp3,
    get_estimates = get_estimates,
    sim_name = SIM_NAME,
    results_directory = RESULTS_DIR,
    k_folds = 8,
    covariate_groups = covariate_groups,
    append = FALSE
)

# process and tabulate results
df_r <- res %>%
    pivot_longer(
        cols = c(starts_with("tevim"), starts_with("std_err")),
        names_to = c(".value", "algorithm", "learner", "estimand"),
        names_pattern = "(.{1,20})_(.{2})_(.{3})_(.{2})$"
    ) %>%
    rename(tevim_u = tevim, std_err_u = std_err) %>%
    pivot_longer(
        cols = c(starts_with("tevim"), starts_with("std_err")),
        names_to = c(".value", "scale"),
        names_pattern = "(.{1,20})_(.{1})$"
    ) %>%
    left_join(true_vals, by = join_by(estimand, scale)) %>%
    mutate(
        bias = tevim - true_value,
        coverage_yn = abs(bias) <= 1.959964 * std_err,
    )

summary_stats <- df_r %>%
    group_by(n, algorithm, learner, estimand, scale) %>%
    summarise(
        n_samples = length(bias),
        coverage = mean(coverage_yn),
        root_n_bias = mean(bias * sqrt(n)),
        n_var = mean(bias^2 * n) - root_n_bias^2,
        root_n_bias_std_err = sqrt(n_var / n_samples),
        bias_min = root_n_bias - 1.96 * root_n_bias_std_err,
        bias_max = root_n_bias + 1.96 * root_n_bias_std_err
    )


sim_plots <- function(df, n_plot) {
    j <- 0.075
    df_plot <- df %>%
        filter(n == n_plot) %>%
        mutate(
            sv = sqrt(n_var),
            # rescale bias and std deviation to match
            root_n_bias = root_n_bias * case_match(
                scale,
                "u" ~ 1,
                "s" ~ 8,
            ),
            sv = sv * case_match(
                scale,
                "u" ~ 1,
                "s" ~ 8,
            ),
            algorithm = case_match(
                algorithm,
                "0T" ~ "1A",
                "0D" ~ "1B",
                "01" ~ "2A",
                "02" ~ "2B",
            ),
        ) %>%
        separate_wider_position(estimand, c(Strategy = 1, Covariate = 1)) %>%
        unite("est", c(Strategy, scale)) %>%
        mutate(
            est = case_match(
                est,
                "k_u" ~ "Keep-one-in, unscaled",
                "l_u" ~ "Leave-one-out, unscaled",
                "k_s" ~ "Keep-one-in, scaled",
                "l_s" ~ "Leave-one-out, scaled",
            ),
            x = as.numeric(Covariate) + case_match(
                algorithm,
                "1A" ~ -3 * j / 2,
                "1B" ~ -j / 2,
                "2A" ~ j / 2,
                "2B" ~ 3 * j / 2,
            ),
        )
    pch <- 0
    bias_plt <- ggplot(data = df_plot, aes(
        x = x, y = root_n_bias, ymin = bias_min, ymax = bias_max
    )) +
        geom_hline(yintercept = 0, lty = 2) +
        geom_point(aes(color = algorithm), shape = pch)
    # geom_pointrange(aes(color = algorithm, shape = est))

    var_plt <- ggplot(data = df_plot, aes(x = x, y = sv)) +
        geom_point(aes(color = algorithm), shape = pch)

    cov_plt <- ggplot(data = df_plot, aes(x = x, y = coverage)) +
        geom_hline(yintercept = 0.95, lty = 2) +
        geom_point(aes(color = algorithm), shape = pch)

    out <- list(Bias = bias_plt, Variance = var_plt, Coverage = cov_plt)

    plot_limits <- list(
        # Bias = c(min(df_plot$bias_min, 0), max(df_plot$bias_max, 0)),
        Bias = c(min(df_plot$root_n_bias, 0), max(df_plot$root_n_bias, 0)),
        Variance = c(min(df_plot$sv), max(df_plot$sv)),
        Coverage = c(min(df$coverage, 0.95), max(df$coverage, 0.95))
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
            xlab("Covariate") +
            ylab(latex2exp::TeX(y_labels[[plt]])) +
            theme_bw() +
            facet_wrap(vars(est)) +
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
            theme(axis.text.x = element_text(
                angle = 0, vjust = 0.5, hjust = 1
            )) +
            scale_x_continuous(breaks = 1:6)
        labs(color = "Scale", shape = "Estimand")
    }

    return(out)
}

p5000 <- sim_plots(summary_stats, 5000)

# Produce and save plots
plots <- list(
    dgp_3_n5000_bias = p5000$Bias,
    dgp_3_n5000_variance = p5000$Variance,
    dgp_3_n5000_coverage = p5000$Coverage
)

for (name in names(plots)) {
    file_name <- glue("{RESULTS_DIR}sim_plot_{SIM_NAME}_{name}.pdf")
    pdf(file = file_name, width = 8, height = 4)
    print(plots[[name]])
    dev.off()
}
