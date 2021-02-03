function [ somme, A ] = Addition_RSB(S, B, RSB)
% [ somme, A ] = Addition_RSB(S, B, RSB) 
% additionne un son S � un bruit B avec un RSB donn� (en dB), et renvoie 
% le facteur A tel que somme = A*S + B;

if length(S)~=length(B)
    error('les deux sons � additioner doivent avoir la m�me longueur');
end

Ps = mean(S.^2);
Pb = mean(B.^2);

% On fixe A pour qu'il y ait le bon RSB dans l'�quation du stim
% 10*log10(PS/PN)=RSB;PS/PN=10^(RSB/10);A^2*Pson/Pbackground=10^(RSB/10);

A=sqrt((Pb/Ps)*10^(RSB/10));
somme = A*S + B;

end

