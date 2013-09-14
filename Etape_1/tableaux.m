function [x, y] = tableaux(n, m, Lx, Ly)

% Cette fonction permet de construire les tableaux x et y � partir de la 
% donn�e de Lx, Ly, n et m.

% Entr�es : n : le nombre de points de discr�tisation selon l'axe x.
%           m : le nombre de points de discr�tisation selon l'axe y.
%           Lx : coordonn�e en x du dernier point sur l'axe des abscisses.
%           Ly : coordonn�e en y du dernier point sur l'axe des ordonn�es.

% Sorties : x : tableau des points xi
%           y : tableau des points yi
%==========================================================================

% CONSTANTES DE SIMULATION x1, y1, Dx et Dy

x(1) = 0;
y(1) = 0;
Dx = Lx/(n-1);
Dy = Ly/(m-1);

% Calcul du vecteur x qui a pour coordonn�es les diff�rents points xi

for i = 1:n-1
    x(i+1) = x(i) + Dx;
end

% Calcul du vecteur y qui a pour coordonn�es les diff�rents points yi

for j = 1:m-1
    y(j+1) = y(j) + Dy;
end