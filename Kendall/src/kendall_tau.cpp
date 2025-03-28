#include <Rcpp.h>
#include <set>
#include <cmath>
#include <map>
#include <vector>

using namespace Rcpp;
//' Calcul du tau de Kendall avec une approche naive
//'
//' @param x NumericVector, premier vecteur de donnees
//' @param y NumericVector, second vecteur de donnees
//' @return Une valeur numerique representant le tau de Kendall
//' @export
// [[Rcpp::export]]
double naive_tau_kendall_cpp(NumericVector x, NumericVector y) {
    int n = x.size();
    if (n != y.size()) Rcpp::stop("x et y doivent avoir la meme longueur.");
    if (n < 2) return 0.0;  // Cas trivial

    int concordant = 0, discordant = 0;

    // Parcours naif O(n**2)
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
    if (total_pairs == 0) return 0.0; // Evite la division par zero

    // Calcul final de tau
    double tau = (double)(concordant - discordant) / total_pairs;
    return tau;
}

// Fusion de deux sous-tableaux triés et comptage des inversions
long long merge_and_count(NumericVector &vec, NumericVector &temp, int left, int mid, int right) {
  int i = left, j = mid, k = left;
  long long inv_count = 0;
  
  while (i <= mid - 1 && j <= right) {
    if (vec[i] <= vec[j]) {
      temp[k++] = vec[i++];
    } else {
      temp[k++] = vec[j++];
      inv_count += (mid - i);
    }
  }
  
  while (i <= mid - 1) temp[k++] = vec[i++];
  while (j <= right) temp[k++] = vec[j++];
  
  for (i = left; i <= right; i++) vec[i] = temp[i];
  
  return inv_count;
}


// Tri-fusion récursif avec comptage d'inversions
long long merge_sort(NumericVector &vec, NumericVector &temp, int left, int right) {
  long long inv_count = 0;
  if (right > left) {
    int mid = (right + left) / 2;
    
    inv_count += merge_sort(vec, temp, left, mid);
    inv_count += merge_sort(vec, temp, mid + 1, right);
    
    inv_count += merge_and_count(vec, temp, left, mid + 1, right);
  }
  return inv_count;
}


//' Calcul du tau de Kendall optimise
//'
//' @param X NumericVector, premier vecteur de donnees
//' @param Y NumericVector, second vecteur de donnees
//' @return Une valeur numerique representant le tau de Kendall
//' @export
// [[Rcpp::export]]
double NDtau_kendall_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  if (n != Y.size()) stop("X et Y doivent avoir la même longueur.");
  if (n < 2) return 0.0;
  
  // Étape 1: Trier Y selon X
  NumericVector Xcopy = clone(X);
  NumericVector Ycopy = clone(Y);
  IntegerVector indices = seq_len(n) - 1;
  std::sort(indices.begin(), indices.end(), [&](int i, int j){ return Xcopy[i] < Xcopy[j]; });
  
  NumericVector Ysorted(n);
  for (int i = 0; i < n; i++) Ysorted[i] = Ycopy[indices[i]];
  
  // Étape 2: Comptage efficace des inversions par tri-fusion
  NumericVector temp(n);
  long long count = merge_sort(Ysorted, temp, 0, n - 1);
  
  // Étape 3: Calculer Tau de Kendall
  double tau = 1.0 - (4.0 * count) / (n * (n - 1));
  
  return tau;
}