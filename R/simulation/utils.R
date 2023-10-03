require(dplyr)
require(glue)


# modified log writing function from github.com/sellorm/rlog
# writes to log_file as well as stdout / stderr
log_to_file <- function(message, log_file, msg_level = "INFO") {
    level <- Sys.getenv("LOG_LEVEL", "INFO")
    levels <- list(
        "all" = 0,
        "trace" = 1,
        "debug" = 2,
        "info" = 3,
        "warn" = 4,
        "error" = 5,
        "fatal" = 6,
        "off" = 7
    )
    if (is.null(levels[tolower(level)])) {
        stop("The LOG_LEVEL environment variable is incorrectly set")
    }

    log_level_int <- as.integer(levels[tolower(level)])
    msg_level_int <- as.integer(levels[tolower(msg_level)])

    if (log_level_int <= msg_level_int) {
        msg <- paste0(Sys.time(), " [", msg_level, "] ", message, "\n")
        cat(msg, file = log_file)
        if (msg_level_int > 3) {
            cat(msg, file = stderr())
        } else {
            cat(msg, file = stdout())
        }
        invisible(TRUE)
    } else {
        invisible(FALSE)
    }
}


# convert list of lists to tibble, similar to base::simplify2array
simplify2tibble <- function(x, names) {
    out <- list()
    y <- unlist(x, recursive = FALSE, use.names = FALSE)
    l <- length(y) / length(x)
    s <- seq(by = l, length.out = length(x))

    for (i in seq_along(names)) {
        out[[i]] <- unlist(
            y[(s + (i - 1))],
            recursive = FALSE,
            use.names = FALSE
        )
    }

    names(out) <- names
    as_tibble(out)
}


run_simulation <- function(
    n_datasets,
    sample_sizes,
    generate_data,
    get_estimates,
    sim_name,
    results_directory = "",
    mc_cores = 10,
    append = TRUE,
    ...) {
    log_file_path <- glue("{results_directory}{sim_name}.log")
    data_file_path <- glue("{results_directory}{sim_name}.parquet")
    log_file <- file(log_file_path, "a")
    on.exit(close(log_file))

    if (file.exists(data_file_path) && !append) {
        log_to_file(
            sprintf(
                "Existing simulation file '%s' will be overwritten.",
                data_file_path
            ),
            msg_level = "WARN",
            log_file
        )
    }

    start_time <- Sys.time()
    log_to_file(sprintf(
        "Beginning simulation '%s' with %s cores.", sim_name, mc_cores
    ), log_file)
    log_to_file(sprintf(
        "Simulation size: %s datasets for %s sample size(s).",
        n_datasets, length(sample_sizes)
    ), log_file)

    sample_set <- rep(sample_sizes, times = n_datasets)
    results <- parallel::mclapply(
        # shuffle the sample set so sample sizes
        # are distributed evenly over cores
        sample(sample_set),
        function(n) {
            df <- generate_data(n = n)
            get_estimates(df, ...)
        },
        mc.cores = mc_cores
    )

    estimate_names <- names(results[[1]])
    results <- simplify2tibble(results, estimate_names)

    end_time <- Sys.time()
    sim_time <- round(difftime(end_time, start_time, units = "secs"), 2)

    log_to_file(sprintf(
        "Finished simulation in %s seconds.", sim_time
    ), log_file)

    log_to_file(sprintf(
        "Writing simulation results to '%s'", data_file_path
    ), log_file)
    if (file.exists(data_file_path) && append) {
        existing_results <- arrow::read_parquet(data_file_path)
        results <- bind_rows(existing_results, results)
    }
    arrow::write_parquet(results, data_file_path)
    (results)
}


## Example ##
#
# test <- run_simulation(
#   n_datasets = 10,
#   sample_sizes = c(250, 500, 750, 1000, 1250, 1500),
#   generate_data = function(n){tibble(X = rnorm(n))},
#   get_estimates = function(df) {list(nrow(df), mean(df$X))},
#   sim_name = "test",
#   results_directory = "simulation-studies/"
# )
#
