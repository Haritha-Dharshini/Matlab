function [U, V] = Neumann(n, m, Lx, Ly, plot)

 
% Cette fonction permet de calculer la solution par approximation numerique
% du probleme -DELTA(U(x,y)) = f(x,y), ou 'DELTA' designe l'operateur Laplacien
% et f(x,y)=-2cos(x)cos(y), dans le cas de Neumann homogene (derivee de U nulle a la
% frontiere).
 
% Entrees : n  : nombre de points de discretisation selon x.
%           Lx : longueur de l'intervalle selon x.
%           m  : nombre de points de discretisation selon y.
%           Ly : longueur de l'intervalle selon y.
%           plot : representation graphique des resultats (utile pour la
%                  fonction erreur)
%                  1 : Oui
%                  0 : Non 
 
% Sorties : V : solution exacte du probleme, donnee par V = cos(x)cos(y).
%           U : approximation numerique de V.

% Developpe par Farah Yasmina Houdroge.
% Login : 20805577


%==========================================================================
%                               Debut du programme
%==========================================================================


% Calcul des vecteurs x, y, dx, dy

dx = Lx/(n-1);
dy = Ly/(m-1);

x = zeros(1,n);
y = zeros(1,m);

for i = 1:n-1
    x(i+1) = x(i) + dx;
end

for j = 1:m-1
    y(j+1) = y(j) + dy;
end

[dx, dy] = delta(x,y);

% Calcul des matrices A et B

a = [1, -1; -1, 1];
b = [1/3, 1/6; 1/6, 1/3];
Ax = zeros(n);
Ay = zeros(m);
Bx = zeros(n);
By = zeros(m);

for i = 1:n-1
    Ax(i:i+1,i:i+1) = Ax(i:i+1,i:i+1) + (1/dx(i))*a;
    Bx(i:i+1,i:i+1) = Bx(i:i+1,i:i+1) + dx(i)*b;
end

for i = 1:m-1  
    Ay(i:i+1,i:i+1) = Ay(i:i+1,i:i+1) + (1/dy(i))*a;
    By(i:i+1,i:i+1) = By(i:i+1,i:i+1) + dy(i)*b;
end

% Calcul des matrices Ax1 = (Bx)^-1*Ax et Ay1 = (By)^-1*Ay

Bxi = inv(Bx);
Byi = inv(By);
Ax1 = Bxi*Ax;
Ay1 = Byi*Ay;

% Diagonalisation de Ax1 et Ay1. P et Q sont les matrices de passages de
% Ax1 et Ay1 respectivement. Adx et Ady sont les matrices Ax1 et Ay1
% diagonalisees, contenant les valeurs propres sur la diagonale, que l'on
% stocke dans deux vecteurs VPx et VPy.

[P, Adx] = eig(Ax1);
[Q, Ady] = eig(Ay1);
VPx = diag(Adx);
VPy = diag(Ady);

% Calcul du second membre F

F = zeros(n,m);

for i = 1:n
    for j = 1:m
        F(i,j) = 2*cos(x(i))*cos(y(j));
    end
end

% Calcul de la matrice F1 = P^-1*F*(Q^-1)^T

Pi = inv(P);
Qi = inv(Q);
F2 = F*(Qi');
F1 = Pi*F2;
 
% Calcul de U1 = P^-1*U*(Q^-1)^T

U1 = zeros(n,m);

for i = 1:n
    for j = 1:m
        if (VPx(i)<10^(-6)) && (VPy(j)<10^(-6))
            U1(i,j) = 0;
        else
            U1(i,j) = F1(i,j)/(VPx(i) + VPy(j));
        end
    end
end

% Calcul de la solution U � partir de U1
 
U = (P*U1)*(Q');

% Solution exacte

V = zeros(n,m);

for i = 1:n
    for j = 1:m
        V(i,j) = cos(x(i))*cos(y(j));
    end
end

% Representation graphique de U

if plot==1
    
    [X,Y] = meshgrid(x,y);
    
    figure
    surf(X,Y,V) ;
    xlabel('x')
    ylabel('y')
    zlabel('V = f(x,y)')
    title('Solution exacte V du probleme de Neumann pour n = m = 100','fontsize',12,'FontWeight','bold')
    colorbar
    
    figure
    surf(X,Y,U) ;
    xlabel('x')
    ylabel('y')
    zlabel('U = f(x,y)')
    title('Solution approximee U du probleme de Neumann pour n = m = 100','fontsize',12,'FontWeight','bold')
    colorbar
    
end

%==========================================================================
%                               Fin du programme
%==========================================================================
