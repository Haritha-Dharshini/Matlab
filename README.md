Matlab
======


Ce projet effectué sous Matlab permet de calculer la solution par approximation numérique du problème -DELTA(U(x,y)) = f(x,y), où 'DELTA' désigne l'opérateur Laplacien et f(x,y)=-2cos(x)cos(y), dans deux cas : 

- le cas de Neuman homogène : dérivee de la vitesse U nulle à la frontière (étape 4).
- le cas de Dirichlet : condition aux limites sur la vitesse U sur les frontières du domaine (étape 5).

Les étapes 1, 2 et 3 sont sont des fonctions intermédiaires détaillées et mises en application dans les étapes 4 (fichier "Neumann.m") et 5 (fichier "Dirichlet.m") et suivent celles du rapport du projet.


Etape 1 :
--------

Création de tableaux contenant les valeurs des pas de discrétisation, les coordonnées des points et leurs indices.


Etape 2 :
--------

Calcul des matrices de masse et de rigidité.


Etape 3 :
--------

Calcul de la solution U du problème -Delta(u(x,y)) = f(x,y).


Etape 4 :
--------

Calcul de la solution U du problème -Delta(U(x,y)) = -2cos(x)cos(y) dans le cas de Neumann homogène et de l'erreur entre la solution exacte et approximée.


Etape 5 :
--------

Calcul de la solution U du problème -Delta(U(x,y)) = -2cos(x)cos(y) dans le cas de Dirichlet et de l'erreur entre la solution exacte et approximée
