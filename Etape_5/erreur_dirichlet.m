function [E] = erreur_dirichlet(L)

% Cette fonction calcule l'erreur E entre la solution exacte du probleme de Neuman
% homogene et la solution trouvee par approximation numérique

% Initialisation

E = [];
ddx = [];

% Calcul de l'erreur et de son logarithme

for n = 20:100

[U, V] = Dirichlet(n,n,L,L,0);
E = [E,log(max(max(abs(U-V))))];
ddx =[ddx,log(L/(n-1))];

end 

P = polyfit(ddx,E,1)

% La fonction polyfit utilisee ci-dessus nous donnera les coefficients A et B de
% la droite log(E) = f(log(ddx)) qui est de la forme y = Ax + B.


% Affichage de log(E) en fonction de log(ddx)

plot(ddx,E,'r')
xlabel('log(\Deltax)')
ylabel('log(E)')
title(['Dirichlet : graphe de log(E) = f(log(\Deltax))',10,'Pente de la droite : P = ',num2str(P(1))])
