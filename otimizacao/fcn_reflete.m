%Fun��o para refletir o vetor de posi��es dos cabos antes de fazer os calculos:

function [x_refl] = fcn_reflete(x) %,y_refl ,y_s

x_o = x(1,1:size(x,2)/2);
y_o = x(1,(size(x,2)/2)+1:size(x,2));

if size(x,2) == 4 || size(x,2) == 10
    x_r = [ -flip(x_o(2:end)) x_o ]; %gera uma matriz com as posi��es x na 1� linha e o reflexo negativo na 2� linha
    y_r = [  flip(y_o(2:end)) y_o ]; %gera uma matriz com as posi��es y na 1� linha e o reflexo na 2� linha
    
    x_refl = [ x_r y_r ];
elseif size(x,2) == 6 || size(x,2) == 12 %pega o tamanho do vetor de posi��es (caso 2 )  
    x_r = [ -flip(x_o(1:end)) x_o ]; %gera uma matriz com as posi��es x na 1� linha e o reflexo negativo na 2� linha
    y_r = [  flip(y_o(1:end)) y_o ]; %gera uma matriz com as posi��es y na 1� linha e o reflexo na 2� linha
    
    x_refl = [ x_r y_r ];
end

%nota: tirei o ; entre os vetores