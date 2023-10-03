true_values <- function(a, b, c) {
    # consider cate of the form: a X1^2 + b X2^2 + c X1^3
    ate <- (a + b) / 3
    te_sq <- c^2 / 7 + (a^2 + b^2) / 5 + 2 * a * b / 9
    vte <- te_sq - ate^2

    tevim_2 <- b^2 * (1 / 5 - 1 / 9)

    tevim_1 <- c^2 / 7 + a^2 * (1 / 5 - 1 / 9)

    list(
        ate = ate,
        vte = vte,
        tevim_1 = tevim_1,
        tevim_1_scaled = tevim_1 / vte,
        tevim_2 = tevim_2,
        tevim_2_scaled = tevim_2 / vte
    )
}
