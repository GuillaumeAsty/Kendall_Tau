# Kendall_Tau
Implementation of Kendall's Tau calculation algorithms in Rcpp

Un .Rmd contenant des tests des fonctions du package se trouve dans le dossier Kendall/tests/.

Les 4 fonctions disponibles sont présentées ci-dessous, avec x et y les 2 vecteurs dont on veut calculer la corrélation au sens de Kendall :

naive_tau_kendall_R(x, y)
NDtau_kendall_R(x, y)
naive_tau_kendall_cpp(x, y)
NDtau_kendall_cpp(x, y)
