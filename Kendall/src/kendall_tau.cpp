#include <Rcpp.h>
#include <set>
#include <cmath>
#include <map>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
double naive_tau_kendall_cpp(NumericVector x, NumericVector y) {
    int n = x.size();
    if (n != y.size()) Rcpp::stop("x et y doivent avoir la même longueur.");
    if (n < 2) return 0.0;  // Cas trivial

    int concordant = 0, discordant = 0;

    // Parcours naïf O(n²)
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            if (dx * dy > 0) {
                concordant++;
            } else if (dx * dy < 0) {
                discordant++;
            }
        }
    }

    int total_pairs = n * (n - 1) / 2;
    if (total_pairs == 0) return 0.0; // Évite la division par zéro

    // Calcul final de tau
    double tau = (double)(concordant - discordant) / total_pairs;
    return tau;
}

// [[Rcpp::export]]
double NDtau_kendall_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  
  // Vérifier que X et Y ont la même longueur
  if (n != Y.size()) {
    stop("X et Y doivent avoir la même longueur.");
  }
  
  // Tau n'est pas défini pour n < 2
  if (n < 2) {
    return 0.0;
  }
  
  // Step 1: Tri des éléments de X et Y selon X
  NumericVector Xcopy = clone(X);
  NumericVector Ycopy = clone(Y);
  
  // Créer indices de tri
  IntegerVector indices = seq_len(n) - 1; // 0-based indices
  std::sort(indices.begin(), indices.end(), 
            [&](int i, int j) {return Xcopy[i] < Xcopy[j];});
  
  // Réorganiser Y selon les indices de tri de X
  NumericVector Ysorted(n);
  for (int i = 0; i < n; i++) {
    Ysorted[i] = Ycopy[indices[i]];
  }
  
  // Step 2: Initialisation du compteur d'inversions
  int c = 0;
  NumericVector Y_sorted(n);
  
  for (int i = 0; i < n; i++) {
    double Yi = Ysorted[i];
    
    // Compter les inversions pour l'élément courant
    for (int j = 0; j < i; j++) {
      if (Y_sorted[j] > Yi) {
        c++;
      }
    }
    
    // Insérer Yi dans Y_sorted pour garder Y_sorted trié
    // Trouver la position d'insertion
    int pos = i;
    while (pos > 0 && Y_sorted[pos-1] > Yi) {
      Y_sorted[pos] = Y_sorted[pos-1];
      pos--;
    }
    Y_sorted[pos] = Yi;
  }
  
  // Step 3: Calcul du Tau de Kendall
  double tau = 1.0 - ((4.0 * c) / (n * (n - 1)));
  
  return tau;
}