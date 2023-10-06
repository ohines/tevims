int_to_set <- function(set_number) {
    # each number 0 to 63 denotes a permutation of included/ excluded covariates
    # This function converts the number to a set of 6 True/ False values
    # which correspond to inclusion/ exclusion of each covariate
    as.logical(intToBits(set_number)[1:6])
}

value_cond <- function(set_number, a, b, c, rho) {
    # The Unscaled TE-VIM (Called Theta_s in the paper)
    # Theta_s has a relatively simple form for the DGP2
    # We abstracting the DGP slightly we consider
    # cate = a * X1 + b * X2 + c * X3
    # where (X1,X2), (X3,X4), (X5,X6)
    # are each joint unit normal with correlation coeficient rho
    conditioning_set <- int_to_set(set_number)

    if (conditioning_set[1] && conditioning_set[2]) {
        s1 <- 0
        s2 <- 0
        s4 <- 0
    } else if (conditioning_set[1]) {
        s1 <- 0
        s2 <- 1 - rho^2
        s4 <- 0
    } else if (conditioning_set[2]) {
        s1 <- 1 - rho^2
        s2 <- 0
        s4 <- 0
    } else {
        s1 <- 1
        s2 <- 1
        s4 <- rho
    }

    if (conditioning_set[3]) {
        s3 <- 0
    } else if (conditioning_set[4]) {
        s3 <- 1 - rho^2
    } else {
        s3 <- 1
    }

    a^2 * s1 + b^2 * s2 + c^2 * s3 + 2 * a * b * s4
}

single_shapley <- function(covariate_index, a, b, c, rho) {
    number_of_covariates <- 6
    step <- 2^(covariate_index - 1) # used to switch the index bit to True

    shaps <- sapply(0:63, function(set) {
        index_bit_is_false <- set %% (2 * step) < step

        if (index_bit_is_false) {

            v_without_i <- -value_cond(set, a, b, c, rho)
            v_with_i <- -value_cond(set + step, a, b, c, rho)

            conditioning_set_size <- sum(int_to_set(set))
            weight <- 1 / choose(
                number_of_covariates - 1, conditioning_set_size
            ) / number_of_covariates

            weight * (v_with_i - v_without_i)
        } else {
            NA
        }
    })
    sum(shaps, na.rm = TRUE)
}


true_values <- function(a, b, c, rho, scaled = TRUE) {
    set_with_no_covariates <- 0
    vte <- value_cond(set_with_no_covariates, a, b, c, rho)

    sets_with_single_false <- 63 - 2^(0:5)
    leave_one_out <- sapply(
        sets_with_single_false,
        value_cond,
        a = a, b = b, c = c, rho = rho
    )

    sets_with_single_true <- 2^(0:5)
    keep_one_in <- sapply(
        sets_with_single_true,
        value_cond,
        a = a, b = b, c = c, rho = rho
    )
    # keep_one_in <- vte - keep_one_in

    shapley_values <- sapply(
        1:6,
        single_shapley,
        a = a, b = b, c = c, rho = rho
    )

    out <- matrix(c(
        leave_one_out,
        keep_one_in,
        shapley_values
    ), ncol = 3)

    if (scaled) {
        out <- out / vte
    }

    colnames(out) <- c("Leave-one-out", "Keep-one-in", "Shapley")
    out
}
