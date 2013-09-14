function [U] = solutionU(X, Y, F)
 
% Cette fonction calcule la solution U du probleme -Delta(u(x,y)) = f(x,y).
% Entrées : X : vecteur contenant les coordonnées respectives des points xi
%           Y : vecteur contenant les coordonnées respectives des points yi
%           F : matrice nxm contenant les valeurs de la fonction du second
%               membre f
%
% Sorties : U : matrice solution du problème -Delta(u) = f.
%========================================================================== 
 
% CONSTANTES ET MATRICES DE SIMULATION 
 
[~, Adx, P, ~, Ady, Q] = matricesAxAy(X, Y);
n = length(X);
m = length(Y);
U1 = zeros(n,m); 
Pi = inv(P);
Qi = inv(Q);
VPx = diag(Adx);
VPy = diag(Ady);

% CALCUL DE LA MATRICE F1 = P^-1*F*(Q^-1)^T
 
F1 = Pi*(F*(Qi'));
 
% CALCUL DE U1 = P^-1*U*(Q^-1)^T
 
for i = 1:n
    for j = 1:m
        if (VPx(i)==0) && (VPy(j)==0)
            U1(i,j) = 0;
        else
            U1(i,j) = F1(i,j)/(VPx(i)+VPy(j));
        end
    end
end

% CALCUL DE LA MATRICE U A PARTIR DE U1
 
U = (P*U1)*(Q');
