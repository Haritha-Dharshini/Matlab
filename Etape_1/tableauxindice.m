function [IP, JP, IR, JR] = tableauxindice(n, m)

% Cette fonction d�finit les tableaux � un indice IP, JP, IR et JR.
 
% Entr�es : n : nombre de points de discr�tisation selon x.
%           m : nombre de points de discr�tisation selon y.
 
% Sorties : IP et JP : fonctions qui donnent les indices i et j du point de
%                      r�f�rence du rectangle numero k.
%           IR et JR : fonctions qui retrouvent les coordonn�es du point de
%                      r�f�rence du rectangle numero k.
%==========================================================================

Nr = (n-1)*(m-1);
Nt = n*m;

for j=1:Nr

IR(j) = j - (n-1)*floor((j-1)/(n-1));
JR(j) = 1 + floor((j-1)/(n-1));

end


for i=1:Nt

JP(i) = 1 + floor((i-1)/n);
IP(i) = i - n*floor((i-1)/n);

end

