function [x, y] = tableaux(n, m, Lx, Ly)

% Cette fonction permet de construire les tableaux x et y à partir de la 
% donnée de Lx, Ly, n et m.

% Entrées : n : le nombre de points de discrétisation selon l'axe x.
%           m : le nombre de points de discrétisation selon l'axe y.
%           Lx : coordonnée en x du dernier point sur l'axe des abscisses.
%           Ly : coordonnée en y du dernier point sur l'axe des ordonnées.

% Sorties : x : tableau des points xi
%           y : tableau des points yi
%==========================================================================

% CONSTANTES DE SIMULATION x1, y1, Dx et Dy

x(1) = 0;
y(1) = 0;
Dx = Lx/(n-1);
Dy = Ly/(m-1);

% Calcul du vecteur x qui a pour coordonnées les différents points xi

for i = 1:n-1
    x(i+1) = x(i) + Dx;
end

% Calcul du vecteur y qui a pour coordonnées les différents points yi

for j = 1:m-1
    y(j+1) = y(j) + Dy;
end