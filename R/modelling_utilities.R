wrapper_gam <- function(y_train, x_train, x_new, family, interactions) {
  p <- dim(x_train)[2]

  if ((p == 3) && interactions) {
    gam_model <- as.formula(
      "y_train~ s(X1) + s(X2) + ti(X1, X2) + s(X1, by=A) + s(X2, by=A) + ti(X1, X2, by=A)"
    )
  } else if ((p == 2) && interactions) {
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
  pred <- mgcv::predict.gam(fit, newdata = data.frame(x_new), type = "response")
  list(pred = pred)
}


wrapper_gam_binomial <- function(y_train, x_train, x_new) {
  wrapper_gam(y_train, x_train, x_new, family = binomial(), TRUE)
}


wrapper_gam_gaussian <- function(y_train, x_train, x_new) {
  wrapper_gam(y_train, x_train, x_new, family = gaussian(), TRUE)
}
