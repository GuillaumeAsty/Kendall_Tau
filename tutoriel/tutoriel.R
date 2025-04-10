devtools::install_github("GuillaumeAsty/Kendall_Tau")
library(Kendall)
library(ggplot2)

# Fonction pour générer des données synthétiques suivant un tau fixé
generate_correlated_data <- function(n = 100) {
  tau <- runif(1, -1, 1)
  rho <- sin(tau * pi / 2)
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- rho * X + sqrt(1 - rho^2) * Z
  return(list(X = X, Y = Y, tau = tau, rho = rho))
}

# Génération des données
data <- generate_correlated_data(100)
x <- data$X
y <- data$Y

# Calculs de tau
tau_R_naif <- naive_tau_kendall_R(x, y)
tau_R_ND <- NDtau_kendall_R(x, y)
tau_cpp_naif <- naive_tau_kendall_cpp(x, y)
tau_cpp_ND <- NDtau_kendall_cpp(x, y)
tau_standard <- cor(x, y, method = "kendall")

# Affichage des résultats
cat("Tau Kendall R naif     :", round(tau_R_naif, 3), "\n")
cat("Tau Kendall R ND       :", round(tau_R_ND, 3), "\n")
cat("Tau Kendall C++ naif   :", round(tau_cpp_naif, 3), "\n")
cat("Tau Kendall C++ ND     :", round(tau_cpp_ND, 3), "\n")
cat("Tau Kendall STANDARD R :", round(tau_standard, 3), "\n")
