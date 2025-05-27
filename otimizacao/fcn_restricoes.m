function [Erms,armazena] = fcn_restricoes(x,V,r,dmin,Dmin)

if size(x,2) == 4 %CASO 1

%x = [ -3 0 3 14.01 14.01 14.01 ] posi��o de entrada original esperada

%gera 3 vari�veis de condi��es
cond1 = 0; cond2 = 0; cond3 = 0;

%% RESTRI��O 1: Dist�ncia m�nima entre condutores

%calcula a dist�ncia � direita entre condutores de fase diferentes
right1 = sqrt(((x(4)-x(3)).^2)+(((x(2)-x(1))).^2));
left1 = Dmin;

armazena.distancia1 = [ right1 left1 ]; %armazena as dist�ncias geradas em um vetor

if right1 < left1 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
	cond1 = cond1+1; %vari�vel condi��o 1 recebe um contador
end

%% RESTRI��O 2: Campo el�trico cr�tico

%armazena em um vetor o resultado da fun��o que calcula o campo el�trico
%superficial para este caso
x_reflet = fcn_reflete(x);
[armazenaEsup] = fcn_sup1(x_reflet,1.0662*V,r);

%armazena na struct 'armazenaAdj' o primeiro field 'Esuperficial' contendo o
%resultado adquirido na linha anterior
armazena.Esup = armazenaEsup;

armazenaEcrit = fcn_supcrit(r*10^(2));

armazena.Ecrit = armazenaEcrit;

cond2cabo = zeros(1,size(x_reflet,2)/2); %cria um vetor de zeros com a metade do tamanho do vetor x original

for t = 1:(size(x_reflet,2)/2) %contador que sai de 1 at� o n� de cabos
    if max(armazenaEsup(t,:)) > armazenaEcrit %caso o m�ximo valor de Esup do determinado cabo for superior ao cr�tico
        %a condi��o de campo eletrico superficial tenta jogar a solu��o para um valor elevado   
        cond2cabo(t) = 1; %cria um vetor que armazena 1 na posi��o t que violou esta restri��o
        cond2 = cond2+1;
        armazena.Esupcond(t) = cond2cabo(t);
    end
end

%% RESTRI��O 3: Campo el�trico ao n�vel do solo

limiteE = 8.33e3; % 8.33 kV/m resolu��o 398 da Aneel

%calcula o campo el�trico ao n�vel do solo para esta config
[armazenaEsolo] = fcn_solo(x_reflet,V,r);

maxE = max(armazenaEsolo); %armazena o maior n�vel de campo el�trico em uma vari�vel
if maxE > limiteE %se este m�ximo valor ultrapassar o limite
    cond3 = cond3+1;
end

%armazena na struct 'armazena' o quinto field 'limite_E' contendo o
%resultado adquirido dos m�ximos valores de E e o limite E da norma
armazena.Esolo = [ maxE limiteE ];

%% Adicionando os pesos

%dividir a fun��o objetivo por um valor para a penalidade ser mais efetiva

Esolonorm = armazenaEsolo/100;

% Segunda Maneira tentativa:
% fx2 = abs(sum(armazenaEsolo));
% coef = round(abs(floor(log10(fx2))));
% 
% %for i=1:length(coef)
% K = 10^coef;
% %end
% 
% Esolonorm = (fx2)/K;

armazena.Esolonorm = Esolonorm;

%cria um vetor :,1 contendo valores em escala de 10 chamados pesos
weights = [ 10e5; 10e8; 10000 ];
%a escala faz diferen�a

%cria um vetor :,1 que armazena os resultados das condi��es
restriction =  [ cond1; cond2; cond3 ];

%cria um vetor 1,: que armazena os resultados das penalidades
%(se a condi��o receber uma contagem, logo ela ter� um resultado > 0 que
%contar� como penalidade e sair� do resultado final)
penalidades =  [ cond1, cond2, cond3 ]*weights;

%armazena na struct 'armazenaAdj' o segundo field 'restriction' contendo o
%resultado adquirido no vetor de condi��es
armazena.restriction = restriction;

%armazena na struct 'armazenaAdj' o segundo field 'penalidades' contendo o
%resultado adquirido no vetor de condi��es*pesos
armazena.penalizacao = penalidades;

%soma a entrada de Esolo ao vetor de penalidades e envia ao main
% Erms = armazenaEsolo + penalidades;
Erms = Esolonorm + penalidades;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x,2) == 6 %CASO 2
    
%x = [ 0.1140 10.136 10.364 14.29 14.29 14.29 ] posi��o de entrada original esperada

%gera 8 vari�veis de condi��es
cond1 = 0; cond2 = 0; cond3 = 0; cond4 = 0; cond5 = 0; cond6 = 0; cond7 = 0; cond8 = 0;

%% RESTRI��O 1: Dist�ncia m�nima entre condutores

%calcula a dist�ncia � direita entre condutores de fase diferentes
right1 = sqrt(((x(4)-x(5)).^2)+(((x(1)-x(2))).^2));
left1 = Dmin;

armazena.distancia1 = [ right1 left1 ]; %armazena as dist�ncias geradas em um vetor

if right1 < left1 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
	cond1 = cond1+1; %vari�vel condi��o 1 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcula a dist�ncia � direita entre o �ltimo condutor e o primeiro condutor
right2 = sqrt(((x(6)-x(4)).^2)+((x(3)-x(1)).^2));
left2 = Dmin+dmin;

armazena.distancia2 = [ right2  left2 ]; %armazena as dist�ncias geradas em um vetor

if right2 < left2 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond2 = cond2+1; %vari�vel condi��o 2 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcula a dist�ncia � direita entre o primeiro condutor e o centro
%do eixo x
right3 = sqrt(((x(4)-x(4)).^2)+(((x(1)-0)).^2));
left3 = dmin/2;

armazena.distancia3 = [ right3  left3 ]; %tirei o Adj por enquanto

if right3 < left3 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond3 = cond3+1; %vari�vel condi��o 3 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right4 = sqrt(((x(6)-x(5)).^2)+(((x(3)-x(2))).^2));
left4 = dmin;

armazena.distancia4 = [ right4  left4 ]; %tirei o Adj por enquanto

if right4 < left4 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond4 = cond4+1; %vari�vel condi��o 4 recebe um contador
end

%% RESTRI��O 2: Campo el�trico cr�tico

%armazena em um vetor o resultado da fun��o que calcula o campo el�trico
%superficial para este caso
x_reflet = fcn_reflete(x);
[armazenaEsup] = fcn_sup2(x_reflet,1.0493*V,r);

%armazena na struct 'armazenaAdj' o primeiro field 'Esuperficial' contendo o
%resultado adquirido na linha anterior
armazena.Esup = armazenaEsup;

armazenaEcrit = fcn_supcrit(r*10^(2));

armazena.Ecrit = armazenaEcrit;

cond5cabo = zeros(1,size(x_reflet,2)/2); %cria um vetor de zeros com a metade do tamanho do vetor x original
 
for t = 1:(size(x_reflet,2)/2) %contador que sai de 1 at� o n� de cabos
    if max(armazenaEsup(t,:)) > armazenaEcrit %caso o m�ximo valor de Esup do determinado cabo for superior ao cr�tico
    %a condi��o de campo eletrico superficial tenta jogar a solu��o para um valor elevado   
        cond5cabo(t) = 1; %cria um vetor que armazena 1 na posi��o t que violou esta restri��o
        cond5 = cond5+1; %cond5 recebe um contador
        armazena.Esupcond(t) = cond5cabo(t);
    end
end

%% RESTRI��O 3: Campo el�trico ao n�vel do solo

limiteE = 8.33e3; % 8.33 kV/m resolu��o 398 da Aneel

%calcula o campo el�trico ao n�vel do solo para esta config
[armazenaEsolo] = fcn_solo(x_reflet,V,r);

maxE = max(armazenaEsolo); %armazena o maior n�vel de campo el�trico em uma vari�vel
if maxE > limiteE %se este m�ximo valor ultrapassar o limite
    cond6 = cond6+1; %condi��o 8 recebe um contador
end

%armazena na struct 'armazena' o quinto field 'limite_E' contendo o
%resultado adquirido dos m�ximos valores de E e o limite E da norma
armazena.Esolo = [ maxE limiteE ];

%% RESTRI��O 4: Dist�ncia m�xima entre subcondutores
dmax2 = 1.50;

right5 = dmax2/2; %armazena numa vari�vel o valor da dist�ncia m�xima entre subcondutores
left5  = sqrt(((x(4)-x(4)).^2)+(((x(1)-0)).^2)); %armazena o c�lculo da dist�ncia entre condutores de fase diferentes

%armazena na struct 'armazena' o segundo field 'distancia12max' contendo o
%resultado adquirido da dist�ncias
armazena.distancia5 = [ left5 right5 ];
 
if right5 < left5 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond7 = cond7+1; %vari�vel condi��o 7 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right6 = dmax2; %armazena a metade da dist�ncia m�xima entre subcondutores
left6  = sqrt(((x(6)-x(5)).^2)+(((x(3)-x(2))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia6 = [ left6 right6 ];
 
if right6 < left6 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond8 = cond8+1; %vari�vel condi��o 8 recebe um contador
end

%% Adicionando os pesos

Esolonorm = armazenaEsolo/1000;

%cria um vetor :,1 contendo valores em escala de 10 chamados pesos
weights = [ 10e2; 10e2; 10e2; 10e2; 10e7; 10e4; 10e2; 10e2 ];
%a escala faz diferen�a

%cria um vetor :,1 que armazena os resultados das condi��es
restriction =  [ cond1; cond2; cond3; cond4; cond5; cond6; cond7; cond8 ];

%cria um vetor 1,: que armazena os resultados das penalidades
%(se a condi��o receber uma contagem, logo ela ter� um resultado > 0 que
%contar� como penalidade e sair� do resultado final)
penalidades =  [ cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8 ]*weights;

%armazena na struct 'armazenaAdj' o segundo field 'restriction' contendo o
%resultado adquirido no vetor de condi��es
armazena.restriction = restriction;

%armazena na struct 'armazenaAdj' o segundo field 'penalidades' contendo o
%resultado adquirido no vetor de condi��es*pesos
armazena.penalizacao = penalidades;

%soma a entrada de Esolo ao vetor de penalidades e envia ao main
% Erms = armazenaEsolo + penalidades;
Erms = Esolonorm + penalidades;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x,2) == 10 %CASO 3
    
%x = [ 0 0.228 10.021 10.25 10.478 16.758 16.530 16.530 16.758 16.530 ]
%posi��o de entrada original esperada

%gera 16 vari�veis de condi��es
cond1 = 0; cond2 = 0; cond3 = 0; cond4 = 0; cond5 = 0; cond6 = 0; cond7 = 0; cond8 = 0; cond9 = 0; cond10 = 0; cond11 = 0; cond12 = 0; cond13 = 0; cond14 = 0; cond15 = 0; cond16 = 0;

%% RESTRI��O 1: Dist�ncia m�nima entre condutores

%calcula a dist�ncia da direita para esquerda entre condutores
right1 = sqrt(((x(7)-x(6)).^2)+(((x(2)-x(1))).^2));
left1 = dmin;

armazena.distancia1 = [ right1 left1 ]; %armazena as dist�ncias geradas em um vetor

if right1 < left1 %se a dist�ncia m�nima for maior que a dist�ncia calculada
	cond1 = cond1+1; %vari�vel condi��o 1 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcula a dist�ncia � direita entre o primeiro condutor e o terceiro da outra fase
right2 = sqrt(((x(8)-x(6)).^2)+((x(3)-x(1)).^2));
left2 = Dmin; %+dmin

armazena.distancia2 = [ right2  left2 ]; %armazena as dist�ncias geradas em um vetor

if right2 < left2 %se a dist�ncia m�nima entre fases for maior que a dist�ncia calculada
    cond2 = cond2+1; %vari�vel condi��o 2 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcula a dist�ncia � direita entre o primeiro condutor e o centro
%do eixo x
right3 = sqrt(((x(9)-x(6)).^2)+(((x(4)-x(1))).^2));
left3 = Dmin; %+2*dmin

armazena.distancia3 = [ right3  left3 ]; %tirei o Adj por enquanto

if right3 < left3 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond3 = cond3+1; %vari�vel condi��o 3 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right4 = sqrt(((x(10)-x(6)).^2)+(((x(5)-x(1))).^2));
left4 = Dmin; %+3*dmin

armazena.distancia4 = [ right4  left4 ]; %tirei o Adj por enquanto

if right4 < left4 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond4 = cond4+1; %vari�vel condi��o 4 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right5 = sqrt(((x(8)-x(7)).^2)+(((x(3)-x(2))).^2));
left5 = Dmin;

armazena.distancia5 = [ right5  left5 ]; %tirei o Adj por enquanto

if right5 < left5 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond5 = cond5+1; %vari�vel condi��o 5 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right6 = sqrt(((x(9)-x(7)).^2)+(((x(4)-x(2))).^2));
left6 = Dmin; %+dmin

armazena.distancia6 = [ right6  left6 ]; %tirei o Adj por enquanto

if right6 < left6 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond6 = cond6+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right7 = sqrt(((x(10)-x(7)).^2)+(((x(5)-x(2))).^2));
left7 = Dmin; %+2*dmin

armazena.distancia7 = [ right7  left7 ]; %tirei o Adj por enquanto

if right7 < left7 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond7 = cond7+1; %vari�vel condi��o 7 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right8 = sqrt(((x(9)-x(8)).^2)+(((x(4)-x(3))).^2));
left8 = dmin;

armazena.distancia8 = [ right8  left8 ]; %tirei o Adj por enquanto

if right8 < left8 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond8 = cond8+1; %vari�vel condi��o 8 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right9 = sqrt(((x(10)-x(8)).^2)+(((x(5)-x(3))).^2));
left9 = dmin;

armazena.distancia9 = [ right9  left9 ]; %tirei o Adj por enquanto

if right9 < left9 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond9 = cond9+1; %vari�vel condi��o 9 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right10 = sqrt(((x(10)-x(9)).^2)+(((x(5)-x(4))).^2));
left10 = dmin;

armazena.distancia10 = [ right10  left10 ]; %tirei o Adj por enquanto

if right10 < left10 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond10 = cond10+1; %vari�vel condi��o 10 recebe um contador
end

%% RESTRI��O 2: Campo el�trico cr�tico

%armazena em um vetor o resultado da fun��o que calcula o campo el�trico
%superficial para este caso
x_reflet = fcn_reflete(x);
[armazenaEsup] = fcn_sup3(x_reflet,1.05*V,r);

%armazena na struct 'armazenaAdj' o primeiro field 'Esuperficial' contendo o
%resultado adquirido na linha anterior
armazena.Esup = armazenaEsup;

armazenaEcrit = fcn_supcrit(r*10^(2));

armazena.Ecrit = armazenaEcrit;

cond5cabo = zeros(1,size(x_reflet,2)/2); %cria um vetor de zeros com a metade do tamanho do vetor x original
 
for t = 1:(size(x_reflet,2)/2) %contador que sai de 1 at� o n� de cabos
    if max(armazenaEsup(t,:)) > armazenaEcrit %caso o m�ximo valor de Esup do determinado cabo for superior ao cr�tico
    %a condi��o de campo eletrico superficial tenta jogar a solu��o para um valor elevado   
        cond5cabo(t) = 1; %cria um vetor que armazena 1 na posi��o t que violou esta restri��o
        cond11 = cond11+1; %cond9 recebe um contador
        armazena.Esupcond(t) = cond5cabo(t);
    end
end

%% RESTRI��O 3: Campo el�trico ao n�vel do solo
 
limiteE = 8.33e3; % 8.33 kV/m resolu��o 398 da Aneel

%calcula o campo el�trico ao n�vel do solo para esta config
[armazenaEsolo] = fcn_solo(x_reflet,V,r);

maxE = max(armazenaEsolo); %armazena o maior n�vel de campo el�trico em uma vari�vel
if maxE > limiteE %se este m�ximo valor ultrapassar o limite
    cond12 = cond12+1; %condi��o 8 recebe um contador
end

%armazena na struct 'armazena' o quinto field 'limite_E' contendo o
%resultado adquirido dos m�ximos valores de E e o limite E da norma
armazena.Esolo = [ maxE limiteE ];

%% RESTRI��O 4: Dist�ncia m�xima entre subcondutores
dmax3 = 2.20;

right11 = dmax3; %armazena numa vari�vel o valor da dist�ncia m�xima entre subcondutores
left11  = sqrt(((x(7)-x(6)).^2)+(((x(2)-x(1))).^2)); %armazena o c�lculo da dist�ncia entre condutores de fase diferentes

%armazena na struct 'armazena' o segundo field 'distancia12max' contendo o
%resultado adquirido da dist�ncias
armazena.distancia11 = [ left11 right11 ];
 
if right11 < left11 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond13 = cond13+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right12 = dmax3; %armazena a metade da dist�ncia m�xima entre subcondutores
left12  = sqrt(((x(9)-x(8)).^2)+(((x(4)-x(3))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia12 = [ left12 right12 ];
 
if right12 < left12 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond14 = cond14+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right13 = 2*dmax3; %armazena a metade da dist�ncia m�xima entre subcondutores
left13  = sqrt(((x(10)-x(8)).^2)+(((x(5)-x(3))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia13 = [ left13 right13 ];
 
if right13 < left13 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond15 = cond15+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right14 = dmax3; %armazena a metade da dist�ncia m�xima entre subcondutores
left14  = sqrt(((x(10)-x(9)).^2)+(((x(5)-x(4))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia14 = [ left14 right14 ];
 
if right14 < left14 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond16 = cond16+1; %vari�vel condi��o 6 recebe um contador
end

%% Adicionando os pesos

Esolonorm = armazenaEsolo/10000;

% Segunda Maneira tentativa:
% fx2 = abs(sum(armazenaEsolo));
% coef = round(abs(floor(log10(fx2))));
% 
% %for i=1:length(coef)
% K = 10^coef;
% %end
% 
% Esolonorm = (fx2)/K;

armazena.Esolonorm = Esolonorm;

%cria um vetor :,1 contendo valores em escala de 10 chamados pesos
weights = [ 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e8; 10e4; 10e4; 10e4; 10e4; 10e4 ];
%a escala faz diferen�a

%cria um vetor :,1 que armazena os resultados das condi��es
restriction =  [ cond1; cond2; cond3; cond4; cond5; cond6; cond7; cond8; cond9; cond10; cond11; cond12; cond13; cond15 ; cond16 ];

%cria um vetor 1,: que armazena os resultados das penalidades
%(se a condi��o receber uma contagem, logo ela ter� um resultado > 0 que
%contar� como penalidade e sair� do resultado final)
penalidades =  [ cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8, cond9, cond10, cond11, cond12, cond13, cond14, cond15, cond16 ]*weights;

%armazena na struct 'armazenaAdj' o segundo field 'restriction' contendo o
%resultado adquirido no vetor de condi��es
armazena.restriction = restriction;

%armazena na struct 'armazenaAdj' o segundo field 'penalidades' contendo o
%resultado adquirido no vetor de condi��es*pesos
armazena.penalizacao = penalidades;

%soma a entrada de Esolo ao vetor de penalidades e envia ao main
% Erms = armazenaEsolo + penalidades;
Erms = Esolonorm + penalidades;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x,2) == 12 %CASO 4

%x = [ 0.475 0.475 7.025 7.025 7.975 7.975 25.95 25 18.45 17.5 18.45 17.5 ]   
%posi��o de entrada original esperada (diferente da disposi��o usada nos
%c�lculos dos Esolo e Esup)

%gera 26 vari�veis de condi��es
cond1 = 0; cond2 = 0; cond3 = 0; cond4 = 0; cond5 = 0; cond6 = 0; cond7 = 0; cond8 = 0; cond9 = 0; cond10 = 0; cond11 = 0; cond12 = 0; cond13 = 0; cond14 = 0; cond15 = 0; cond16 = 0; cond17 = 0; cond18 = 0; cond19 = 0; cond20 = 0; cond21 = 0; cond22 = 0; cond23 = 0; cond24 = 0; cond25 = 0; cond26 = 0;

%% RESTRI��O 1: Dist�ncia m�nima entre condutores

%calcula a dist�ncia da direita para esquerda entre condutores
right1 = sqrt(((x(8)-x(7)).^2)+(((x(2)-x(1))).^2));
left1 = dmin;

armazena.distancia1 = [ right1 left1 ]; %armazena as dist�ncias geradas em um vetor

if right1 < left1 %se a dist�ncia m�nima for maior que a dist�ncia calculada
	cond1 = cond1+1; %vari�vel condi��o 1 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcula a dist�ncia � direita entre o primeiro condutor e o terceiro da outra fase
right2 = sqrt(((x(9)-x(7)).^2)+((x(3)-x(1)).^2));
left2 = Dmin;

armazena.distancia2 = [ right2  left2 ]; %armazena as dist�ncias geradas em um vetor

if right2 < left2 %se a dist�ncia m�nima entre fases for maior que a dist�ncia calculada
    cond2 = cond2+1; %vari�vel condi��o 2 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcula a dist�ncia � direita entre o primeiro condutor e o centro
%do eixo x
right3 = sqrt(((x(10)-x(7)).^2)+(((x(4)-x(1))).^2));
left3 = Dmin;

armazena.distancia3 = [ right3  left3 ]; %tirei o Adj por enquanto

if right3 < left3 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond3 = cond3+1; %vari�vel condi��o 3 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right4 = sqrt(((x(11)-x(7)).^2)+(((x(5)-x(1))).^2));
left4 = Dmin+dmin;

armazena.distancia4 = [ right4  left4 ]; %tirei o Adj por enquanto

if right4 < left4 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond4 = cond4+1; %vari�vel condi��o 4 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right5 = sqrt(((x(12)-x(7)).^2)+(((x(6)-x(1))).^2));
left5 = Dmin+dmin;

armazena.distancia5 = [ right5  left5 ]; %tirei o Adj por enquanto

if right5 < left5 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond5 = cond5+1; %vari�vel condi��o 5 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right6 = sqrt(((x(9)-x(8)).^2)+(((x(3)-x(2))).^2));
left6 = Dmin;

armazena.distancia6 = [ right6  left6 ]; %tirei o Adj por enquanto

if right6 < left6 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond6 = cond6+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right7 = sqrt(((x(10)-x(8)).^2)+(((x(4)-x(2))).^2));
left7 = Dmin;

armazena.distancia7 = [ right7  left7 ]; %tirei o Adj por enquanto

if right7 < left7 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond7 = cond7+1; %vari�vel condi��o 7 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right8 = sqrt(((x(11)-x(8)).^2)+(((x(5)-x(2))).^2));
left8 = Dmin+dmin;

armazena.distancia8 = [ right8  left8 ]; %tirei o Adj por enquanto

if right8 < left8 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond8 = cond8+1; %vari�vel condi��o 8 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right9 = sqrt(((x(12)-x(8)).^2)+(((x(6)-x(2))).^2));
left9 = Dmin+dmin;

armazena.distancia9 = [ right9  left9 ]; %tirei o Adj por enquanto

if right9 < left9 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond9 = cond9+1; %vari�vel condi��o 9 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right10 = sqrt(((x(10)-x(9)).^2)+(((x(4)-x(3))).^2));
left10 = dmin;

armazena.distancia10 = [ right10  left10 ]; %tirei o Adj por enquanto

if right10 < left10 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond10 = cond10+1; %vari�vel condi��o 10 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right11 = sqrt(((x(11)-x(9)).^2)+(((x(5)-x(3))).^2));
left11 = dmin;

armazena.distancia11 = [ right11  left11 ]; %tirei o Adj por enquanto

if right11 < left11 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond11 = cond11+1; %vari�vel condi��o 10 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right12 = sqrt(((x(12)-x(9)).^2)+(((x(6)-x(3))).^2));
left12 = dmin;

armazena.distancia12 = [ right12  left12 ]; %tirei o Adj por enquanto

if right12 < left12 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond12 = cond12+1; %vari�vel condi��o 10 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right13 = sqrt(((x(11)-x(10)).^2)+(((x(5)-x(4))).^2));
left13 = dmin;

armazena.distancia13 = [ right13  left13 ]; %tirei o Adj por enquanto

if right13 < left13 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond13 = cond13+1; %vari�vel condi��o 10 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right14 = sqrt(((x(12)-x(10)).^2)+(((x(6)-x(4))).^2));
left14 = dmin;

armazena.distancia14 = [ right14  left14 ]; %tirei o Adj por enquanto

if right14 < left14 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond14 = cond14+1; %vari�vel condi��o 10 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right15 = sqrt(((x(12)-x(11)).^2)+(((x(6)-x(5))).^2));
left15 = dmin;

armazena.distancia15 = [ right15  left15 ]; %tirei o Adj por enquanto

if right15 < left15 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond15 = cond15+1; %vari�vel condi��o 10 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right23 = sqrt(((x(7)-x(7)).^2)+(((x(1)-0)).^2));
left23 = dmin/2;

armazena.distancia23 = [ right23  left23 ]; %tirei o Adj por enquanto

if right23 < left23 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond25 = cond25+1; %vari�vel condi��o 10 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right24 = sqrt(((x(8)-x(8)).^2)+(((x(2)-0)).^2));
left24 = dmin/2;

armazena.distancia24 = [ right24  left24 ]; %tirei o Adj por enquanto

if right24 < left24 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond26 = cond26+1; %vari�vel condi��o 10 recebe um contador
end

%% RESTRI��O 2: Campo el�trico cr�tico

%armazena em um vetor o resultado da fun��o que calcula o campo el�trico
%superficial para este caso
x_reflet = fcn_reflete(x);
[armazenaEsup] = fcn_sup4(x_reflet,1.05*V,r);

%armazena na struct 'armazenaAdj' o primeiro field 'Esuperficial' contendo o
%resultado adquirido na linha anterior
armazena.Esup = armazenaEsup;

armazenaEcrit = fcn_supcrit(r*10^(2));

armazena.Ecrit = armazenaEcrit;

cond5cabo = zeros(1,size(x_reflet,2)/2); %cria um vetor de zeros com a metade do tamanho do vetor x original

for t = 1:(size(x_reflet,2)/2) %contador que sai de 1 at� o n� de cabos
    if max(armazenaEsup(t,:)) > armazenaEcrit %caso o m�ximo valor de Esup do determinado cabo for superior ao cr�tico
    %a condi��o de campo eletrico superficial tenta jogar a solu��o para um valor elevado   
        cond5cabo(t) = 1; %cria um vetor que armazena 1 na posi��o t que violou esta restri��o
        cond16 = cond16+1; %cond9 recebe um contador
        armazena.Esupcond(t) = cond5cabo(t);
    end
end


%% RESTRI��O 3: Campo el�trico ao n�vel do solo
 
limiteE = 8.33e3; % 8.33 kV/m resolu��o 398 da Aneel

%calcula o campo el�trico ao n�vel do solo para esta config
[armazenaEsolo] = fcn_solo(x_reflet,V,r);

maxE = max(armazenaEsolo); %armazena o maior n�vel de campo el�trico em uma vari�vel
if maxE > limiteE %se este m�ximo valor ultrapassar o limite
    cond17 = cond17+1; %condi��o 8 recebe um contador
end

%armazena na struct 'armazena' o quinto field 'limite_E' contendo o
%resultado adquirido dos m�ximos valores de E e o limite E da norma
armazena.Esolo = [ maxE limiteE ];

%% RESTRI��O 4: Dist�ncia m�xima entre subcondutores
dmax4 = 2.20;

right16 = dmax4; %armazena numa vari�vel o valor da dist�ncia m�xima entre subcondutores
left16 = sqrt(((x(8)-x(7)).^2)+(((x(2)-x(1))).^2)); %armazena o c�lculo da dist�ncia entre condutores de fase diferentes

%armazena na struct 'armazena' o segundo field 'distancia12max' contendo o
%resultado adquirido da dist�ncias
armazena.distancia16 = [ left16 right16 ];
 
if right16 < left16 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond18 = cond18+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right17 = dmax4; %armazena a metade da dist�ncia m�xima entre subcondutores
left17  = sqrt(((x(10)-x(9)).^2)+(((x(4)-x(3))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia17 = [ left17 right17 ];
 
if right17 < left17 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond19 = cond19+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right18 = dmax4; %armazena a metade da dist�ncia m�xima entre subcondutores
left18 = sqrt(((x(11)-x(9)).^2)+(((x(5)-x(3))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia18 = [ left18 right18 ];
 
if right18 < left18 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond20 = cond20+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right19 = dmax4; %armazena a metade da dist�ncia m�xima entre subcondutores
left19 = sqrt(((x(12)-x(9)).^2)+(((x(6)-x(3))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia19 = [ left19 right19 ];
 
if right19 < left19 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond21 = cond21+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right20 = dmax4; %armazena a metade da dist�ncia m�xima entre subcondutores
left20 = sqrt(((x(11)-x(10)).^2)+(((x(5)-x(4))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia20 = [ left20 right20 ];
 
if right20 < left20 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond22 = cond22+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right21 = dmax4; %armazena a metade da dist�ncia m�xima entre subcondutores
left21 = sqrt(((x(12)-x(10)).^2)+(((x(6)-x(4))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia21 = [ left21 right21 ];
 
if right21 < left21 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond23 = cond23+1; %vari�vel condi��o 6 recebe um contador
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
right22 = dmax4; %armazena a metade da dist�ncia m�xima entre subcondutores
left22 = sqrt(((x(12)-x(11)).^2)+(((x(6)-x(5))).^2)); %armazena o c�lculo da dist�ncia entre o condutor mais � direito e o centro

%armazena na struct 'armazena' o terceiro field 'distancia30max' contendo o
%resultado adquirido das dist�ncias
armazena.distancia22 = [ left22 right22 ];
 
if right22 < left22 %se a dist�ncia � esquerda for maior que a dist�ncia � direita
    cond24 = cond24+1; %vari�vel condi��o 6 recebe um contador
end

%% Adicionando os pesos

%penalizar um pouco mais e ver o que acontece (feito)
%3. abaixar a penaliza��o ou 2. aumentar o n�mero de itera��es/indiv�duos
%ou 1. manter os condutores ao centro parados (em y)
Esolonorm = armazenaEsolo/10000;

%cria um vetor :,1 contendo valores em escala de 10 chamados pesos
weights = [ 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e8; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4; 10e4 ];
%a escala faz diferen�a

%cria um vetor :,1 que armazena os resultados das condi��es
restriction =  [ cond1; cond2; cond3; cond4; cond5; cond6; cond7; cond8; cond9; cond10; cond11; cond12; cond13; cond15; cond16; cond17; cond18; cond19; cond20; cond21; cond22; cond23; cond24; cond25 ; cond26 ];

%cria um vetor 1,: que armazena os resultados das penalidades
%(se a condi��o receber uma contagem, logo ela ter� um resultado > 0 que
%contar� como penalidade e sair� do resultado final)
penalidades =  [ cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8, cond9, cond10, cond11, cond12, cond13, cond14, cond15, cond16, cond17, cond18, cond19, cond20, cond21, cond22, cond23, cond24, cond25, cond26 ]*weights;

%armazena na struct 'armazenaAdj' o segundo field 'restriction' contendo o
%resultado adquirido no vetor de condi��es
armazena.restriction = restriction;

%armazena na struct 'armazenaAdj' o segundo field 'penalidades' contendo o
%resultado adquirido no vetor de condi��es*pesos
armazena.penalizacao = penalidades;

%soma a entrada de Esolo ao vetor de penalidades e envia ao main
% Erms = armazenaEsolo + penalidades;
Erms = Esolonorm + penalidades;

end