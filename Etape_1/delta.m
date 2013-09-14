function [deltax, deltay] = delta(X,Y)

% Cette fonction calcule les tableaux contenant les valeurs Dx et Dy.
% 
% Entr�es : X : vecteur contenant les coordonn�es respectives des points xi
%           Y : vecteur contenant les coordonn�es respectives des points yi
% 
% Sorties : deltax : tableau des pas de discr�tisation Dxi
%           deltay : tableau des pas de discr�tisation Dyi. 
%==========================================================================

n = length(X);
m = length(Y);

for i = 1:n-1
    deltax(i) = X(i+1) - X(i);
end

for j = 1:m-1
    deltay(j) = Y(j+1) - Y(j);
end