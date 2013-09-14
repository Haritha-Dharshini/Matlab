function [U,V] = Dirichlet(n,m,Lx,Ly,plot) 

% Cette fonction permet de calculer la solution par approximation numerique
% du probleme -DELTA(U(x,y)) = f(x,y), ou 'DELTA' designe l'operateur Laplacien
% et f(x,y)=-2cos(x)cos(y), dans le cas de Dirichlet.

% Entrees : n : nombre de points de discretisation selon x.
%           Lx : longueur de l'intervalle selon x.
%           m : nombre de points de discretisation selon y.
%           Ly : longueur de l'intervalle selon y.
%           plot : representation graphique des resultats (utile pour la
%                  fonction erreur)
%                  1 : Oui
%                  0 : Non
 
% Sorties : V : solution exacte du probleme, donnee par V = cos(x)cos(y).
%           U : approximation numerique de V.

% Développé par Farah Yasmina Houdroge.
% Login : 20805577

%==========================================================================
%                               Début du programme
%==========================================================================

% Calcul des vecteurs x, y, dx et dy

dx = Lx/(n-1);
dy = Ly/(m-1);

x = zeros(1,n);
y = zeros(1,m);

for i = 1:(n-1)
    x(i+1) = x(i) + dx;
end

for j = 1:(m-1)
    y(j+1) = y(j) + dy;
end

[dx, dy] = delta(x,y);

% Calcul du second membre F

F = zeros(n,m);

for i = 1:n
    for j = 1:m
        F(i,j) = 2*cos(x(i))*cos(y(j));
    end
end


% Solution exacte

V = zeros(n,m);

for i = 1:n
    for j = 1:m        
        V(i,j) = cos(x(i))*cos(y(j));        
    end
end

% Construction des matrices A et B

a = [1, -1; -1, 1];
b = [1/3, 1/6; 1/6, 1/3];
Ax = zeros(n);
Bx = zeros(n);
Ay = zeros(m);
By = zeros(m);

for i = 1:n-1
    
    Ax(i:i+1,i:i+1) = Ax(i:i+1,i:i+1) + ((1/dx(i))*a);
    Bx(i:i+1,i:i+1) = Bx(i:i+1,i:i+1) + ((dx(i))*b);
    
end

for j = 1:m-1
    
    Ay(j:j+1,j:j+1) = Ay(j:j+1,j:j+1) + ((1/dy(j))*a);
    By(j:j+1,j:j+1) = By(j:j+1,j:j+1) + ((dy(j))*b);
    
end


% Résolution du système linéaire

Bxi = inv(Bx);
Byi = inv(By);
Ax1 = Bxi*Ax;
Ay1 = Byi*Ay;

% Condition aux limites

Ub = V;

for i = 1:m
    Ub(1,i) = 0;
    Ub(n,i) = 0;
end

for j = 1:n
    Ub(j,1) = 0;
    Ub(j,m) = 0;
end

Uk = V - Ub;

% Application de la Condition aux limites sur F, Axt et Ayt

F = F - Ax1(:,1)*Uk(1,:) - Ax1(:,n)*Uk(n,:) - Uk(:,1)*Ay1(:,1)'- Uk(:,m)*Ay1(:,m)';

Ax1(:,1) = 0; 
Ax1(:,n) = 0; 
Ay1(:,1) = 0; 
Ay1(:,m) = 0;

[P, Adx] = eig(Ax1);
[Q, Ady] = eig(Ay1);
VPx = diag(Adx);
VPy = diag(Ady);

Pi = inv(P);
Qi = inv(Q);

F2 = F*Qi';
F1 = Pi*F2;

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

U = (P*U1)*(Q');

% Application de la condition aux limites sur U

U(1,:) = 0; 
U(n,:) = 0; 
U(:,1) = 0; 
U(:,m) = 0;

U = U + Uk;

% Representation graphique de la solution

if plot==1
    
    [X,Y] = meshgrid(x,y);
    
    figure
    surf(X,Y,U);
    xlabel('x')
    ylabel('y')
    zlabel('U = f(x,y)')
    title('Solution approximée U du problème dans le cas de Dirichlet pour n = m = 100','fontsize',12,'FontWeight','bold')
    colorbar
    
    figure
    surf(X,Y,V);
    xlabel('x')
    ylabel('y')
    zlabel('V = f(x,y)')
    title('Solution exacte V du problème dans le cas de Dirichlet pour n = m = 100','fontsize',12,'FontWeight','bold')
    colorbar
    
end

%==========================================================================
%                               Fin du programme
%==========================================================================