#' Calcul du tau de Kendall (version naive)
#'
#' Cette fonction calcule le coefficient de Kendall en O(n**2) en comparant chaque paire.
#'
#' @param X Un vecteur numerique.
#' @param Y Un vecteur numerique.
#' @return Le coefficient de Kendall tau.
#' @examples
#' X <- c(3, 1, 4, 2, 5)
#' Y <- c(2, 3, 1, 5, 4)
#' naive_tau_kendall_R(X, Y)
#' @export
naive_tau_kendall_R <- function(X, Y) {
  n <- length(X)
  if (n != length(Y)) stop("X et Y doivent avoir la meme longueur.")
  if (n < 2) return(0)  # Tau n'est pas defini pour n < 2
  
  concordant <- 0
  discordant <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dx <- X[i] - X[j]
      dy <- Y[i] - Y[j]
      if (dx * dy > 0) {
        concordant <- concordant + 1
      } else if (dx * dy < 0) {
        discordant <- discordant + 1
      }
    }
  }
  
  total_pairs <- n * (n - 1) / 2
  if (total_pairs == 0) return(0)  # Cas ou n < 2
  
  tau <- (concordant - discordant) / total_pairs
  
  return(tau)
}


#' Calcul du tau de Kendall optimise (O(n log n))
#'
#' Cette fonction utilise un tri et une structure de donnees optimisee pour calculer tau en O(n log n).
#'
#' @param X Un vecteur numerique.
#' @param Y Un vecteur numerique.
#' @return Le coefficient de Kendall tau.
#' @examples
#' X <- c(3, 1, 4, 2, 5)
#' Y <- c(2, 3, 1, 5, 4)
#' NDtau_kendall_R(X, Y)
#' @export
NDtau_kendall_R <- function(X, Y) {
  count_inv_merge <- function(vec) {
    len <- length(vec)
    if (len == 1) return(list(sorted=vec, inv=0))
    mid <- len %/% 2
    left <- count_inv_merge(vec[1:mid])
    right <- count_inv_merge(vec[(mid+1):len])
    inv_count <- left$inv + right$inv
    merged <- numeric(len)
    i <- j <- k <- 1
    left_vec <- left$sorted
    right_vec <- right$sorted
    while (i <= length(left_vec) && j <= length(right_vec)) {
      if (left_vec[i] <= right_vec[j]) {
        merged[k] <- left_vec[i]; i <- i + 1
      } else {
        merged[k] <- right_vec[j]; inv_count <- inv_count + (length(left_vec)-i+1); j <- j + 1
      }
      k <- k + 1
    }
    while (i <= length(left_vec)) { merged[k] <- left_vec[i]; i <- i + 1; k <- k + 1 }
    while (j <= length(right_vec)) { merged[k] <- right_vec[j]; j <- j + 1; k <- k + 1 }
    return(list(sorted=merged, inv=inv_count))
  }
  
  if (length(X) != length(Y)) stop("X et Y doivent avoir la meme longueur.")
  n <- length(X)
  if (n < 2) return(0)
  ord <- order(X)
  Y <- Y[ord]
  result <- count_inv_merge(Y)
  tau <- 1 - ((4 * result$inv) / (n * (n - 1)))
  return(tau)
  
}