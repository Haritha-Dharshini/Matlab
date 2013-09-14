function [Ax1, Adx, P, Ay1, Ady, Q] = matricesAxAy(X, Y)

% Cette fonction calcule les matrices Ax1 = Bx^-1*Ax et Ay1 = By^-1*Ay et 
% leur diagonales successives.
% Entr�es : X : vecteur contenant les coordonn�es respectives des points xi
%           Y : vecteur contenant les coordonn�es respectives des points yi
% 
% Sorties : Ax1 : matrice d�finie par (Bx)^-1*Ax
%           Adx : matrice obtenue apr�s la diagonalisation de Ax1
%           P : matrice de passage de Ax1
%           Ay1 : matrice d�finie par (By)^-1*Ay
%           Ady : matrice obtenue apr�s la diagonalisation de Ay1
%           Q : matrice de passage de Ay1
%==========================================================================

[Ax, Ay, Bx, By] = matriceAB(X, Y);

% Calcul des matrices Ax1 = (Bx)^-1*Ax et Ay1 = (By)^-1*Ay

Ax1 = Bx\Ax;
Ay1 = By\Ay;

% Diagonalisation de Ax1 et Ay1

[P, Adx] = eig(Ax1);
[Q, Ady] = eig(Ay1);