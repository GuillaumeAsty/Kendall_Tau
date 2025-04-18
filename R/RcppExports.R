# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Calcul du tau de Kendall avec une approche naive
#'
#' @param x NumericVector, premier vecteur de donnees
#' @param y NumericVector, second vecteur de donnees
#' @return Une valeur numerique representant le tau de Kendall
#' @export
naive_tau_kendall_cpp <- function(x, y) {
    .Call(`_Kendall_naive_tau_kendall_cpp`, x, y)
}

#' Calcul du tau de Kendall optimise
#'
#' @param X NumericVector, premier vecteur de donnees
#' @param Y NumericVector, second vecteur de donnees
#' @return Une valeur numerique representant le tau de Kendall
#' @export
NDtau_kendall_cpp <- function(X, Y) {
    .Call(`_Kendall_NDtau_kendall_cpp`, X, Y)
}

