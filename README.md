# Kendall_Tau
Implementation of Kendall's Tau calculation algorithms in Rcpp

Un .R et un .Rmd contenant des tests et une analyse des fonctions du package se trouve dans le dossier Kendall/tutoriel/.

Les 4 fonctions disponibles sont présentées ci-dessous, avec x et y les 2 vecteurs dont on veut calculer la corrélation au sens de Kendall :

- naive_tau_kendall_R(x, y) : implémentation en R avec une complexité en N^2

- NDtau_kendall_R(x, y) : implémentation en R avec une complexité en N*logN

- naive_tau_kendall_cpp(x, y) : implémentation en Rcpp avec une complexité en N^2

- NDtau_kendall_cpp(x, y) : implémentation en Rcpp avec une complexité en N*logN

