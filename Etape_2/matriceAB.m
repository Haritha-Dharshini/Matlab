function [Ax, Ay, Bx, By] = matriceAB(X, Y)

% Cette fonction permet de calculer les matrices A et B.
% Entrées : X : vecteur contenant les coordonnées respectives des points xi
%           Y : vecteur contenant les coordonnées respectives des points yi
% Sorties : Ax : matrice de rigidité selon x
%           Ay : matrice de rigidité selon y
%           Bx : matrice de masse selon x
%           By : matrice de masse selon y.
%==========================================================================

% CONSTANTES DE SIMULATION a, b, n, m, Dx, Dy

[a] = matrice_a();
[b] = matrice_b();
[deltax, deltay] = delta(X, Y);
n = length(X);
m = length(Y);

% INITIALISATION DES MATRICES Ax, Ay, Bx et By

Ax = zeros(n);
Ay = zeros(m);
Bx = zeros(n);
By = zeros(m);

% CALCUL DES MATRICES A ET B PAR RECURRENCE

for i = 1:n-1
    Ax(i:i+1,i:i+1) = Ax(i:i+1,i:i+1) + (1/deltax(i))*a;
   
    Bx(i:i+1,i:i+1) = Bx(i:i+1,i:i+1) + deltax(i)*b;
    
end

for j = 1:m-1
    
    Ay(j:j+1,j:j+1) = Ay(j:j+1,j:j+1) + (1/deltay(j))*a;
    
    By(j:j+1,j:j+1) = By(j:j+1,j:j+1) + deltay(j)*b;
    
end
