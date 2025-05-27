%Função para refletir o vetor de posições dos cabos antes de fazer os calculos:

function [x_refl] = fcn_reflete(x) %,y_refl ,y_s

x_o = x(1,1:size(x,2)/2);
y_o = x(1,(size(x,2)/2)+1:size(x,2));

if size(x,2) == 4 || size(x,2) == 10
    x_r = [ -flip(x_o(2:end)) x_o ]; %gera uma matriz com as posições x na 1ª linha e o reflexo negativo na 2ª linha
    y_r = [  flip(y_o(2:end)) y_o ]; %gera uma matriz com as posições y na 1ª linha e o reflexo na 2ª linha
    
    x_refl = [ x_r y_r ];
elseif size(x,2) == 6 || size(x,2) == 12 %pega o tamanho do vetor de posições (caso 2 )  
    x_r = [ -flip(x_o(1:end)) x_o ]; %gera uma matriz com as posições x na 1ª linha e o reflexo negativo na 2ª linha
    y_r = [  flip(y_o(1:end)) y_o ]; %gera uma matriz com as posições y na 1ª linha e o reflexo na 2ª linha
    
    x_refl = [ x_r y_r ];
end

%nota: tirei o ; entre os vetores