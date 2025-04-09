# Kendall_Tau
Implementation of Kendall's Tau calculation algorithms in R and Rcpp

Un fichier .Rmd de tutoriel et un fichier.Rmd contenant des tests des fonctions du package sont disponibles dans le dossier Kendall/tutoriel/.

Les 4 fonctions disponibles sont présentées ci-dessous, avec x et y les 2 vecteurs dont on veut calculer la corrélation au sens de Kendall :

- naive_tau_kendall_R(x, y) : implémentation en R avec une complexité en N^2

- NDtau_kendall_R(x, y) : implémentation en R avec une complexité en N*logN

- naive_tau_kendall_cpp(x, y) : implémentation en Rcpp avec une complexité en N^2

- NDtau_kendall_cpp(x, y) : implémentation en Rcpp avec une complexité en N*logN

