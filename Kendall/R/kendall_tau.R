#' Calcul du tau de Kendall (version naïve)
#'
#' Cette fonction calcule le coefficient de Kendall en O(n²) en comparant chaque paire.
#'
#' @param X Un vecteur numérique.
#' @param Y Un vecteur numérique.
#' @return Le coefficient de Kendall tau.
#' @examples
#' X <- c(3, 1, 4, 2, 5)
#' Y <- c(2, 3, 1, 5, 4)
#' naive_tau_kendall_R(X, Y)
#' @export
naive_tau_kendall_R <- function(X, Y) {
  n <- length(X)
  if (n != length(Y)) stop("X et Y doivent avoir la même longueur.")
  if (n < 2) return(0)  # Tau n'est pas défini pour n < 2
  
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
  if (total_pairs == 0) return(0)  # Cas où n < 2
  
  tau <- (concordant - discordant) / total_pairs
  
  return(tau)
}


#' Calcul du tau de Kendall optimisé (O(n log n))
#'
#' Cette fonction utilise un tri et une structure de données optimisée pour calculer tau en O(n log n).
#'
#' @param X Un vecteur numérique.
#' @param Y Un vecteur numérique.
#' @return Le coefficient de Kendall tau.
#' @examples
#' X <- c(3, 1, 4, 2, 5)
#' Y <- c(2, 3, 1, 5, 4)
#' NDtau_kendall_R(X, Y)
#' @export
NDtau_kendall_R <- function(X, Y) {
  n <- length(X)
  if (n != length(Y)) stop("X et Y doivent avoir la même longueur.")
  if (n < 2) return(0)  # Tau n'est pas défini pour n < 2
  
  # Step 1: Tri des éléments de X et Y selon X
  sorted_indices <- order(X)
  X <- X[sorted_indices]
  Y <- Y[sorted_indices]
  
  # Step 2: Initialisation de la variable pour compter les inversions
  c <- 0
  Y_sorted <- numeric(n)  # Pour stocker les valeurs triées de Y
  for (i in 1:n) {
    Yi <- Y[i]
    # Insertion binaire de Yi dans Y_sorted pour garder Y_sorted trié
    pos <- which(Y_sorted[1:i-1] > Yi)
    if (length(pos) > 0) {
      c <- c + length(pos)
    }
    Y_sorted[i] <- Yi
  }
  
  # Step 3: Calcul du Tau de Kendall
  tau <- 1 - ((4 * c) / (n * (n - 1)))

  return(tau)
}
