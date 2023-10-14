require(dplyr)

# for convenience we isolate lapply
# this means we can that we can replace it with a parallel version if we want
lapply_func <- lapply

te_vim <- function(po, cate, sub_cate, scaled = FALSE) {
  # average treatment effect
  n <- length(po)
  ate <- sum(po) / n

  # three residual like terms
  r_ate <- (po - ate)^2
  r_cate <- (po - cate)^2
  r_subcate <- (po - sub_cate)^2

  # evaluate TE-VIM (Theta_s in the paper)
  tevim <- sum(r_subcate - r_cate) / n
  infl <- r_subcate - r_cate - tevim
  std_err <- sqrt(sum(infl^2)) / n

  # evaluate scaled TE-VIM (Psi_s in the paper)
  vte <- sum(r_ate - r_cate) / n
  tevim_s <- tevim / vte
  infl <- (r_subcate - tevim_s * r_ate + (tevim_s - 1) * r_cate) / vte
  std_err_s <- sqrt(sum(infl^2)) / n

  list(
    tevim = tevim,
    std_err = std_err,
    tevim_s = tevim_s,
    std_err_s = std_err_s
  )
}


ate_vte <- function(po, cate) {
  # average treatment effect
  n <- length(po)
  ate <- sum(po) / n

  r_ate <- (po - ate)^2
  r_cate <- (po - cate)^2
  vte <- sum(r_ate - r_cate) / n

  ate_std_err <- sqrt(sum(r_ate)) / n
  vte_std_err <- sqrt(sum((r_ate - r_cate - vte)^2)) / n
  list(
    ate = ate,
    ate.std_err = ate_std_err,
    vte = vte,
    vte.std_err = vte_std_err
  )
}


# Algorithm: no cross fitting
algorithm_0 <- function(
    y,
    a,
    x,
    y_model,
    ps_model,
    cate_model,
    covariate_groups) {
  # get pseudo-outcomes and T-learner CATE using outcome and ps model
  po_fit <- pseudo_outcome_model(y, a, x, y, a, x, y_model, ps_model)
  po <- po_fit$pred
  cate_t <- po_fit$cate

  # fit the DR-learner by regressing on pseudo-outcomes
  cate_dr <- cate_model(po, x, x)$pred

  # Fit reduced variable CATE and evaluate
  te_vims <- lapply_func(covariate_groups, function(covs) {
    x_s <- select(x, -all_of(covs))
    t_subcate <- cate_model(cate_t, x_s, x_s)$pred
    d_subcate <- cate_model(cate_dr, x_s, x_s)$pred
    list(
      te_vim_t = te_vim(po, cate_t, t_subcate),
      te_vim_dr = te_vim(po, cate_dr, d_subcate)
    )
  }) %>% purrr::list_transpose()

  out_t <- list(
    estiamtes = simplify2array(te_vims$te_vim_t),
    ate_vte = simplify2array(ate_vte(po, cate_t))
  )

  out_dr <- list(
    estiamtes = simplify2array(te_vims$te_vim_dr),
    ate_vte = simplify2array(ate_vte(po, cate_dr))
  )

  list(
    algorithm_0T = out_t,
    algorithm_0D = out_dr
  )
}


# Algorithm: T-learner with cross fitting
algorithm_1 <- function(
    y,
    a,
    x,
    folds,
    y_model,
    ps_model,
    cate_model,
    covariate_groups,
    return_ates = FALSE) {
  # cross fit to obtain pseudo-outcome and cate predictions
  fold_list <- unique(folds)
  n <- length(folds)
  cross_fits <- lapply_func(fold_list, function(fold) {
    in_train <- folds != fold
    in_test <- !in_train

    model_algorithm1(
      y[in_train],
      a[in_train],
      x[in_train, ],
      y[in_test],
      a[in_test],
      x[in_test, ],
      y_model,
      ps_model,
      cate_model,
      covariate_groups
    )
  })

  # merge predictions into output vectors
  po <- rep(NA, length.out = n)
  cate <- rep(NA, length.out = n)
  subcates <- matrix(nrow = n, ncol = length(covariate_groups))

  for (i in seq_along(fold_list)) {
    in_fold <- folds == fold_list[i]
    po[in_fold] <- cross_fits[[i]]$po
    cate[in_fold] <- cross_fits[[i]]$cate
    for (j in seq_along(covariate_groups)) {
      subcates[in_fold, j] <- cross_fits[[i]]$subcates[[j]]
    }
  }

  tevims <- apply(subcates, 2, function(subcate) te_vim(po, cate, subcate))
  names(tevims) <- names(covariate_groups)

  list(
    estiamtes = simplify2array(tevims),
    ate_vte = simplify2array(ate_vte(po, cate))
  )
}


# Algorithm: DR-learner with cross-fitting
algorithm_2 <- function(
    y,
    a,
    x,
    folds,
    y_model,
    ps_model,
    cate_model,
    covariate_groups) {
  n <- length(y)
  fold_list <- unique(folds)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

  # (double) cross-fit pseudo outcomes
  # double in the sense of leaving out two folds
  cross_fits <- lapply_func(fold_pairs, function(fold_pair) {
    in_train <- !folds %in% fold_pair
    pseudo_outcome_model(
      y[in_train],
      a[in_train],
      x[in_train, ],
      y[!in_train],
      a[!in_train],
      x[!in_train, ],
      y_model,
      ps_model
    )
  })

  # merge cross fitted pseudo-outcomes into a matrix
  # po is an n by k matrix, with missing values for each fold
  # i.e. each row has one missing value.
  po_matrix <- lapply(fold_list, function(fold) {
    po_hat <- rep(NA, length.out = n)
    for (i in seq_along(fold_pairs)) {
      fold_pair <- fold_pairs[[i]]
      if (fold %in% fold_pair) {
        next
      }
      po_hat[folds %in% fold_pair] <- cross_fits[[i]]$pred
    }
    po_hat[folds == fold] <- NA
    po_hat
  }) %>% simplify2array()

  # average predictions to get out of sample pseudo-outcome estimate
  po <- rowMeans(po_matrix, na.rm = TRUE)

  # fit CATE and subcate models for each fold
  cate <- rep(NA, length.out = n)
  subcates <- matrix(nrow = n, ncol = length(covariate_groups))

  for (i in seq_along(fold_list)) {
    fold <- fold_list[[i]]
    in_train <- folds != fold
    in_fold <- !in_train
    x_train <- x[in_train, ]
    x_new <- x[in_fold, ]

    # fit the cate
    cate_fit <- cate_model(
      po_matrix[in_train, i],
      x_train,
      x
    )$pred
    cate[in_fold] <- cate_fit[in_fold]
    cate_train <- cate_fit[in_train]

    # train subcate models
    for (j in seq_along(covariate_groups)) {
      covs <- covariate_groups[[j]]
      subcates[in_fold, j] <- cate_model(
        cate_train,
        select(x_train, -any_of(covs)),
        select(x_new, -any_of(covs))
      )$pred
    }
  }

  tevims <- apply(subcates, 2, function(subcate) te_vim(po, cate, subcate))
  names(tevims) <- names(covariate_groups)

  list(
    estiamtes = simplify2array(tevims),
    ate_vte = simplify2array(ate_vte(po, cate))
  )
}


all_algorithms <- function(
    y,
    a,
    x,
    folds,
    y_model,
    ps_model,
    cate_model,
    covariate_groups) {
  res_0 <- algorithm_0(
    y, a, x,
    y_model,
    ps_model,
    cate_model,
    covariate_groups
  )
  res_1 <- algorithm_1(
    y, a, x, folds,
    y_model,
    ps_model,
    cate_model,
    covariate_groups
  )
  res_2 <- algorithm_2(
    y, a, x, folds,
    y_model,
    ps_model,
    cate_model,
    covariate_groups
  )

  estimates <- do.call(
    rbind,
    list(
      res_0$algorithm_0T$estiamtes,
      res_0$algorithm_0D$estiamtes,
      res_1$estiamtes,
      res_2$estiamtes
    )
  )

  algorithm_suffixes <- c("_0T", "_0D", "_01", "_02")
  rownames(estimates) <- purrr::map2_chr(
    rownames(estimates),
    rep(algorithm_suffixes, each = 4),
    paste0
  )

  ate_vte <- do.call(
    rbind,
    list(
      res_0$algorithm_0T$ate_vte,
      res_0$algorithm_0D$ate_vte,
      res_1$ate_vte,
      res_2$ate_vte
    )
  )
  rownames(ate_vte) <- paste0("algorithm", algorithm_suffixes)

  list(
    estimates = estimates,
    ate_vte = ate_vte
  )
}


model_algorithm1 <- function(
    y_train,
    a_train,
    x_train,
    y_new,
    a_new,
    x_new,
    y_model,
    ps_model,
    cate_model,
    covariate_groups) {
  m <- length(y_new)
  n <- m + length(y_train)

  # outcome predictions for the all data
  x <- rbind(x_new, x_train, deparse.level = 0)
  y_pred <- y_model(
    y_train,
    cbind(A = a_train, x_train, deparse.level = 0),
    rbind(
      cbind(A = 1, x, deparse.level = 0),
      cbind(A = 0, x, deparse.level = 0),
      deparse.level = 0
    )
  )$pred
  y1_hat <- y_pred[1:n]
  y0_hat <- y_pred[(n + 1):(2 * n)]
  cate <- y1_hat - y0_hat

  # submodel cate predictions for new data only
  cate_train <- cate[(m + 1):n]
  subcates <- lapply(covariate_groups, function(covs) {
    cate_model(
      cate_train,
      select(x_train, -any_of(covs)),
      select(x_new, -any_of(covs))
    )$pred
  })

  # get pseudo outcomes for new data using propensity score model
  cate_new <- cate[1:m]
  ps_hat <- ps_model(a_train, x_train, x_new)$pred
  y_hat <- if_else(a_new > 0.5, y1_hat[1:m], y0_hat[1:m])
  po_new <- cate_new + (y_new - y_hat) *
    (a_new - ps_hat) / ps_hat / (1 - ps_hat)

  # return predictions for 'new' data only
  list(po = po_new, cate = cate_new, subcates = subcates)
}


pseudo_outcome_model <- function(
    y,
    a,
    x,
    y_new,
    a_new,
    x_new,
    y_model,
    ps_model) {
  n <- length(y_new)

  # outcome and propensity score predictions for the new data
  ps_hat <- ps_model(a, x, x_new)$pred
  y_pred <- y_model(
    y,
    cbind(A = a, x, deparse.level = 0),
    rbind(
      cbind(A = 1, x_new, deparse.level = 0),
      cbind(A = 0, x_new, deparse.level = 0),
      deparse.level = 0
    )
  )$pred
  y1_hat <- y_pred[1:n]
  y0_hat <- y_pred[(n + 1):(2 * n)]

  y_hat <- if_else(a_new > 0.5, y1_hat, y0_hat)
  cate <- y1_hat - y0_hat
  po <- cate + (y_new - y_hat) * (a_new - ps_hat) / ps_hat / (1 - ps_hat)

  list(pred = po, cate = cate)
}


t_learner <- function(
    y_train,
    x_train,
    x_new,
    treatment_var_name,
    fit_fn,
    ...) {
  a_train <- x_train[, treatment_var_name] == 1
  a_new <- x_new[, treatment_var_name] == 1

  non_treat <- colnames(x_train) != treatment_var_name

  fit_1 <- fit_fn(
    y_train[a_train],
    x_train[a_train, non_treat],
    x_new[a_new, non_treat],
    ...
  )
  fit_0 <- fit_fn(
    y_train[!a_train],
    x_train[!a_train, non_treat],
    x_new[!a_new, non_treat],
    ...
  )

  pred <- rep(NA, length(a_new))
  pred[a_new] <- fit_1$pred
  pred[!a_new] <- fit_0$pred
  list(pred = pred)
}
