% Raio menor

clear all; close all; clc

% Cálculo do campo elétrico superficial do caso 2 do sistema trifásico

e_0 = 8.854*(10^(-12));
r = 0.00733; % raio do condutor
n = 6; % número de condutores do sistema
nc = n*2; % número de condutores total (com imagens)
ci = (2*n-1); % número de cargas imagens

cond = 6;
%% matriz P

% posição dos condutores reais do sistema
xr = [-10.478 -10.25 -0.228 0 10.022 10.25];
yr = 14.290;
yi = -yr;

% cálculo matriz P
P = zeros(n,n);
for i = 1:n
     for j = 1:n
        if(i==j)
           P(i,j) = (1/(2*pi*e_0))*log((4*yr)/(2*r));
        else
           P(i,j) = (1/(2*pi*e_0))*log(sqrt((xr(j)-xr(i))^2+(yi-yr)^2)/(sqrt((xr(j)-xr(i))^2)));
        end
     end
end

%% Cálculo tensão por fase

V = 362*10^3;

%tensão condutores 1 e 2 fase a
V_ra = V/sqrt(3);
V_ia = 0;

%tensão condutores 3 e 4 fase b
V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

%tensão condutores 5 e 6 fase c
V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ];

%% Cálculo densidade de carga

rho = P\Vf;

%busco cada uma das posicoes de rho real e imaginario (ver indice)

%cabo 1 fase a
rho_r1 = real(rho(1)); 
rho_i1 = imag(rho(1));

%cabo 2 fase a
rho_r2 = real(rho(2));
rho_i2 = imag(rho(2));

%cabo 3 fase b
rho_r3 = real(rho(3));
rho_i3 = imag(rho(3));

%cabo 4 fase b
rho_r4 = real(rho(4));
rho_i4 = imag(rho(4));

%cabo 5 fase c
rho_r5 = real(rho(5));
rho_i5 = imag(rho(5));

%cabo 6 fase c
rho_r6 = real(rho(6));
rho_i6 = imag(rho(6));


%% Distância entre condutores

x = [-10.478 -10.25 -0.228 0 10.022 10.25 -10.478 -10.25 -0.228 0 10.022 10.25]; % posição x dos condutores (com imagens)
y = [yr yr yr yr yr yr -yr -yr -yr -yr -yr -yr]; % posição y dos condutores (com imagens)

%sendo: i -> o número do primeiro condutor e j -> o número do segundo condutor
%exemplo: 1,2 é a distância entre condutor 1 e o condutor 2, os próximos
%cálculos seguem esta mesma lógica
%sendo 7 a imagem do condutor 1, 9 a imagem do condutor 3 e 11 a imagem do
%condutor 5

distance = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j
            distance(i,j) = 0; %como não existe distância no centro do próprio condutor então é zero
        else
            distance(i,j) = sqrt(((x(j)-x(i))^2)+((y(j)-y(i))^2)); %cálculo distância entre dois pontos (entre centro de dois condutores)
        end
    end
end

%% Cálculo de delta

delta = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j
            delta(i,j) = 0; %delta para o centro do próprio condutor é zero
        else
            delta(i,j) = (r^2)./(distance(i,j)); %cálculo de delta utilizando a fórmula e a distância calculada acima
        end
    end
end

%% Cálculo das posições em x das cargas imagens

phix = zeros(nc,nc);
posx = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j || x(j)==x(i) %cálculo de posição para cargas posicionadas no centro dos condutores
            posx(i,j) = x(i);
        elseif y(j)==y(i) %cálculo para cargas de mesma altura
            if x(i)~=x(j) && i > j
                posx(i,j) = x(i) - delta(i,j);
            elseif x(i)~=x(j) && j > i
                posx(i,j) = x(i) + delta(i,j);
            end
        elseif y(i) > 0 %cálculo de posição para cargas de altura diferentes mas positiva
            if x(i)~=x(j) && x(i) > 0  || i==4 && j==7 || i==4 && j==8 || i==4 && j==9 %cargas do lado direito do plano e carga 2-4
                phix(i,j) = acos((2*y(i))/(distance(i,j))); %cálculo do ângulo phi entre o triângulo gerado
                posx(i,j) = x(i) - delta(i,j)*sin(phix(i,j)); %cálculo da posição de x para este caso, onde a projeção de delta em x é feita pelo seno do ângulo phi
            elseif x(i)~=x(j) && x(i) < 0 || i==4 && j==11 || i==4 && j==12 %condutores do lado esquerdo do plano e carga 2-6
                phix(i,j) = acos((2*y(i))/(distance(i,j)));
                posx(i,j) = x(i) + delta(i,j)*sin(phix(i,j));
            end
        elseif y(i) < 0 %cálculo de posição para cargas de altura diferentes mas negativa
            if x(i)~=x(j) && x(i) > 0 || i==10 && j==1 || i==10 && j==2 || i==10 && j==3 %cargas do lado direito do plano e carga 5-1
                phix(i,j) = asin((2*y(i))/(distance(i,j)));
                posx(i,j) = x(i) - delta(i,j)*cos(phix(i,j));
            elseif x(i)~=x(j) && x(i) < 0 || i==10 && j==5 || i==10 && j==6 %cargas do lado esquerdo do plano e carga 5-3
                phix(i,j) = asin((2*y(i))/(distance(i,j)));
                posx(i,j) = x(i) + delta(i,j)*cos(phix(i,j));
            end
        end
    end
end

%% Cálculo das posições em y das cargas imagens

phiy = zeros(nc,nc);
posy = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j || y(j)==y(i) %cálculo de posição para cargas posicionadas no centro dos condutores
            posy(i,j) = y(i);
        elseif x(i)==x(j) %cálculo de posição para cargas posicionadas na mesma altura dos condutores
            posy(i,j) = y(i) - delta(i,j);
        elseif x(i)~=x(j) %cálculo de posição para cargas posicionadas em altura diferentes dos condutores
            if y(i) > 0 %altura positiva (condutores reais)
                phiy(i,j) = acos((2*y(i))/(distance(i,j))); %cálculo de phi (segue mesma lógica da posição de x
                posy(i,j) = y(i) - delta(i,j)*cos(phiy(i,j));
            elseif y(i) < 0 %altura negativa (condutores imagens)
                phiy(i,j) = asin((2*y(i))/(distance(i,j)));
                posy(i,j) = y(i) + delta(i,j)*sin(phiy(i,j));
            end
        end
    end
end

% plot(posx,posy,'o');
% grid
% xlim([-2.8e-05,2.8e-05])
% ylim([-14.010003,-14.009]) 

%% Pontos de avaliação (toda superfície do condutor)

theta = linspace(0,2*pi,360); %gera a superfície do condutor em 360 pontos

xcj = x(cond); % posição x do centro condutor 1 fase a (para trocar o condutor avaliado deve-se mudar este valor)
ycj = y(cond); % posição y do centro condutor 1 fase a

xf = r.*cos(theta) + xcj; % eixo x do ponto de avaliação
yf = r.*sin(theta) + ycj; % eixo y do ponto de avaliação

% plot(xf,yf)

%% E_xr componente x real campo elétrico condutor 1 fase a assim segue:

%busca o valor da carga imagem no condutor 1 real, valor do eixo x dos
%pontos de avaliação e as posições geradas pelas matrizes posx e posy das
%cargas imagens
E_xr12 = (-rho_r2/(2*pi*e_0)).*((xf - posx(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_xr13 = (-rho_r3/(2*pi*e_0)).*((xf - posx(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_xr14 = (-rho_r4/(2*pi*e_0)).*((xf - posx(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_xr15 = (-rho_r5/(2*pi*e_0)).*((xf - posx(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_xr16 = (-rho_r6/(2*pi*e_0)).*((xf - posx(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_xr17 = (rho_r1/(2*pi*e_0)).*((xf - posx(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_xr18 = (rho_r2/(2*pi*e_0)).*((xf - posx(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_xr19 = (rho_r3/(2*pi*e_0)).*((xf - posx(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_xr110 = (rho_r4/(2*pi*e_0)).*((xf - posx(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_xr1_11 = (rho_r5/(2*pi*e_0)).*((xf - posx(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_xr1_12 = (rho_r6/(2*pi*e_0)).*((xf - posx(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));

E_xr21 = (-rho_r1/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xr23 = (-rho_r3/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xr24 = (-rho_r4/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xr25 = (-rho_r5/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xr26 = (-rho_r6/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_xr27 = (rho_r1/(2*pi*e_0)).*((xf - posx(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_xr28 = (rho_r2/(2*pi*e_0)).*((xf - posx(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_xr29 = (rho_r3/(2*pi*e_0)).*((xf - posx(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_xr210 = (rho_r4/(2*pi*e_0)).*((xf - posx(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_xr2_11 = (rho_r5/(2*pi*e_0)).*((xf - posx(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_xr2_12 = (rho_r6/(2*pi*e_0)).*((xf - posx(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));

E_xr31 = (-rho_r1/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xr32 = (-rho_r2/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xr34 = (-rho_r4/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xr35 = (-rho_r5/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xr36 = (-rho_r6/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_xr37 = (rho_r1/(2*pi*e_0)).*((xf - posx(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_xr38 = (rho_r2/(2*pi*e_0)).*((xf - posx(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_xr39 = (rho_r3/(2*pi*e_0)).*((xf - posx(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_xr310 = (rho_r4/(2*pi*e_0)).*((xf - posx(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_xr311 = (rho_r5/(2*pi*e_0)).*((xf - posx(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_xr312 = (rho_r6/(2*pi*e_0)).*((xf - posx(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));

E_xr41 = (-rho_r1/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xr42 = (-rho_r2/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xr43 = (-rho_r3/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xr45 = (-rho_r5/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xr46 = (-rho_r6/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_xr47 = (rho_r1/(2*pi*e_0)).*((xf - posx(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_xr48 = (rho_r2/(2*pi*e_0)).*((xf - posx(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_xr49 = (rho_r3/(2*pi*e_0)).*((xf - posx(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_xr410 = (rho_r4/(2*pi*e_0)).*((xf - posx(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_xr411 = (rho_r5/(2*pi*e_0)).*((xf - posx(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_xr412 = (rho_r6/(2*pi*e_0)).*((xf - posx(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));

E_xr51 = (-rho_r1/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xr52 = (-rho_r2/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xr53 = (-rho_r3/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xr54 = (-rho_r4/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xr56 = (-rho_r6/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_xr57 = (rho_r1/(2*pi*e_0)).*((xf - posx(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_xr58 = (rho_r2/(2*pi*e_0)).*((xf - posx(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_xr59 = (rho_r3/(2*pi*e_0)).*((xf - posx(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_xr510 = (rho_r4/(2*pi*e_0)).*((xf - posx(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_xr511 = (rho_r5/(2*pi*e_0)).*((xf - posx(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_xr512 = (rho_r6/(2*pi*e_0)).*((xf - posx(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));

E_xr61 = (-rho_r1/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xr62 = (-rho_r2/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xr63 = (-rho_r3/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xr64 = (-rho_r4/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xr65 = (-rho_r5/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_xr67 = (rho_r1/(2*pi*e_0)).*((xf - posx(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_xr68 = (rho_r2/(2*pi*e_0)).*((xf - posx(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_xr69 = (rho_r3/(2*pi*e_0)).*((xf - posx(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_xr610 = (rho_r4/(2*pi*e_0)).*((xf - posx(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_xr611 = (rho_r5/(2*pi*e_0)).*((xf - posx(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_xr612 = (rho_r6/(2*pi*e_0)).*((xf - posx(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));

E_xr71 = (-rho_r1/(2*pi*e_0)).*((xf - posx(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_xr72 = (-rho_r2/(2*pi*e_0)).*((xf - posx(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_xr73 = (-rho_r3/(2*pi*e_0)).*((xf - posx(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_xr74 = (-rho_r4/(2*pi*e_0)).*((xf - posx(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_xr75 = (-rho_r5/(2*pi*e_0)).*((xf - posx(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_xr76 = (-rho_r6/(2*pi*e_0)).*((xf - posx(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_xr78 = (rho_r2/(2*pi*e_0)).*((xf - posx(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_xr79 = (rho_r3/(2*pi*e_0)).*((xf - posx(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_xr710 = (rho_r4/(2*pi*e_0)).*((xf - posx(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_xr711 = (rho_r5/(2*pi*e_0)).*((xf - posx(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_xr712 = (rho_r6/(2*pi*e_0)).*((xf - posx(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));

E_xr81 = (-rho_r1/(2*pi*e_0)).*((xf - posx(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_xr82 = (-rho_r2/(2*pi*e_0)).*((xf - posx(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_xr83 = (-rho_r3/(2*pi*e_0)).*((xf - posx(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_xr84 = (-rho_r4/(2*pi*e_0)).*((xf - posx(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_xr85 = (-rho_r5/(2*pi*e_0)).*((xf - posx(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_xr86 = (-rho_r6/(2*pi*e_0)).*((xf - posx(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_xr87 = (rho_r1/(2*pi*e_0)).*((xf - posx(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_xr89 = (rho_r3/(2*pi*e_0)).*((xf - posx(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_xr810 = (rho_r4/(2*pi*e_0)).*((xf - posx(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_xr811 = (rho_r5/(2*pi*e_0)).*((xf - posx(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_xr812 = (rho_r6/(2*pi*e_0)).*((xf - posx(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));

E_xr91 = (-rho_r1/(2*pi*e_0)).*((xf - posx(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_xr92 = (-rho_r2/(2*pi*e_0)).*((xf - posx(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_xr93 = (-rho_r3/(2*pi*e_0)).*((xf - posx(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_xr94 = (-rho_r4/(2*pi*e_0)).*((xf - posx(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_xr95 = (-rho_r5/(2*pi*e_0)).*((xf - posx(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_xr96 = (-rho_r6/(2*pi*e_0)).*((xf - posx(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_xr97 = (rho_r1/(2*pi*e_0)).*((xf - posx(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_xr98 = (rho_r2/(2*pi*e_0)).*((xf - posx(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_xr910 = (rho_r4/(2*pi*e_0)).*((xf - posx(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_xr911 = (rho_r5/(2*pi*e_0)).*((xf - posx(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_xr912 = (rho_r6/(2*pi*e_0)).*((xf - posx(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));

E_xr101 = (-rho_r1/(2*pi*e_0)).*((xf - posx(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_xr102 = (-rho_r2/(2*pi*e_0)).*((xf - posx(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_xr103 = (-rho_r3/(2*pi*e_0)).*((xf - posx(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_xr104 = (-rho_r4/(2*pi*e_0)).*((xf - posx(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_xr105 = (-rho_r5/(2*pi*e_0)).*((xf - posx(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_xr106 = (-rho_r6/(2*pi*e_0)).*((xf - posx(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_xr107 = (rho_r1/(2*pi*e_0)).*((xf - posx(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_xr108 = (rho_r2/(2*pi*e_0)).*((xf - posx(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_xr109 = (rho_r3/(2*pi*e_0)).*((xf - posx(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_xr1011 = (rho_r5/(2*pi*e_0)).*((xf - posx(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_xr1012 = (rho_r6/(2*pi*e_0)).*((xf - posx(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));

E_xr11_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_xr11_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_xr113 = (-rho_r3/(2*pi*e_0)).*((xf - posx(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_xr114 = (-rho_r4/(2*pi*e_0)).*((xf - posx(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_xr115 = (-rho_r5/(2*pi*e_0)).*((xf - posx(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_xr116 = (-rho_r6/(2*pi*e_0)).*((xf - posx(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_xr117 = (rho_r1/(2*pi*e_0)).*((xf - posx(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_xr118 = (rho_r2/(2*pi*e_0)).*((xf - posx(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_xr119 = (rho_r3/(2*pi*e_0)).*((xf - posx(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_xr1110 = (rho_r4/(2*pi*e_0)).*((xf - posx(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_xr1112 = (rho_r6/(2*pi*e_0)).*((xf - posx(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));

E_xr12_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_xr12_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_xr123 = (-rho_r3/(2*pi*e_0)).*((xf - posx(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_xr124 = (-rho_r4/(2*pi*e_0)).*((xf - posx(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_xr125 = (-rho_r5/(2*pi*e_0)).*((xf - posx(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_xr126 = (-rho_r6/(2*pi*e_0)).*((xf - posx(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_xr127 = (rho_r1/(2*pi*e_0)).*((xf - posx(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_xr128 = (rho_r2/(2*pi*e_0)).*((xf - posx(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_xr129 = (rho_r3/(2*pi*e_0)).*((xf - posx(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_xr1210 = (rho_r4/(2*pi*e_0)).*((xf - posx(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_xr1211 = (rho_r5/(2*pi*e_0)).*((xf - posx(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));

%% E_xi componente x imaginario campo elétrico condutor 2 fase b assim segue:

E_xi12 = (-rho_i2/(2*pi*e_0)).*((xf - posx(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_xi13 = (-rho_i3/(2*pi*e_0)).*((xf - posx(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_xi14 = (-rho_i4/(2*pi*e_0)).*((xf - posx(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_xi15 = (-rho_i5/(2*pi*e_0)).*((xf - posx(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_xi16 = (-rho_i6/(2*pi*e_0)).*((xf - posx(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_xi17 = (rho_i1/(2*pi*e_0)).*((xf - posx(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_xi18 = (rho_i2/(2*pi*e_0)).*((xf - posx(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_xi19 = (rho_i3/(2*pi*e_0)).*((xf - posx(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_xi110 = (rho_i4/(2*pi*e_0)).*((xf - posx(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_xi1_11 = (rho_i5/(2*pi*e_0)).*((xf - posx(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_xi1_12 = (rho_i6/(2*pi*e_0)).*((xf - posx(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));

E_xi21 = (-rho_i1/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xi23 = (-rho_i3/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xi24 = (-rho_i4/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xi25 = (-rho_i5/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xi26 = (-rho_i6/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_xi27 = (rho_i1/(2*pi*e_0)).*((xf - posx(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_xi28 = (rho_i2/(2*pi*e_0)).*((xf - posx(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_xi29 = (rho_i3/(2*pi*e_0)).*((xf - posx(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_xi210 = (rho_i4/(2*pi*e_0)).*((xf - posx(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_xi2_11 = (rho_i5/(2*pi*e_0)).*((xf - posx(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_xi2_12 = (rho_i6/(2*pi*e_0)).*((xf - posx(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));

E_xi31 = (-rho_i1/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xi32 = (-rho_i2/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xi34 = (-rho_i4/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xi35 = (-rho_i5/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xi36 = (-rho_i6/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_xi37 = (rho_i1/(2*pi*e_0)).*((xf - posx(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_xi38 = (rho_i2/(2*pi*e_0)).*((xf - posx(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_xi39 = (rho_i3/(2*pi*e_0)).*((xf - posx(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_xi310 = (rho_i4/(2*pi*e_0)).*((xf - posx(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_xi311 = (rho_i5/(2*pi*e_0)).*((xf - posx(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_xi312 = (rho_i6/(2*pi*e_0)).*((xf - posx(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));

E_xi41 = (-rho_i1/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xi42 = (-rho_i2/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xi43 = (-rho_i3/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xi45 = (-rho_i5/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xi46 = (-rho_i6/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_xi47 = (rho_i1/(2*pi*e_0)).*((xf - posx(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_xi48 = (rho_i2/(2*pi*e_0)).*((xf - posx(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_xi49 = (rho_i3/(2*pi*e_0)).*((xf - posx(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_xi410 = (rho_i4/(2*pi*e_0)).*((xf - posx(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_xi411 = (rho_i5/(2*pi*e_0)).*((xf - posx(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_xi412 = (rho_i6/(2*pi*e_0)).*((xf - posx(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));

E_xi51 = (-rho_i1/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xi52 = (-rho_i2/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xi53 = (-rho_i3/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xi54 = (-rho_i4/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xi56 = (-rho_i6/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_xi57 = (rho_i1/(2*pi*e_0)).*((xf - posx(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_xi58 = (rho_i2/(2*pi*e_0)).*((xf - posx(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_xi59 = (rho_i3/(2*pi*e_0)).*((xf - posx(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_xi510 = (rho_i4/(2*pi*e_0)).*((xf - posx(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_xi511 = (rho_i5/(2*pi*e_0)).*((xf - posx(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_xi512 = (rho_i6/(2*pi*e_0)).*((xf - posx(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));

E_xi61 = (-rho_i1/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xi62 = (-rho_i2/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xi63 = (-rho_i3/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xi64 = (-rho_i4/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xi65 = (-rho_i5/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_xi67 = (rho_i1/(2*pi*e_0)).*((xf - posx(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_xi68 = (rho_i2/(2*pi*e_0)).*((xf - posx(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_xi69 = (rho_i3/(2*pi*e_0)).*((xf - posx(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_xi610 = (rho_i4/(2*pi*e_0)).*((xf - posx(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_xi611 = (rho_i5/(2*pi*e_0)).*((xf - posx(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_xi612 = (rho_i6/(2*pi*e_0)).*((xf - posx(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));

E_xi71 = (-rho_i1/(2*pi*e_0)).*((xf - posx(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_xi72 = (-rho_i2/(2*pi*e_0)).*((xf - posx(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_xi73 = (-rho_i3/(2*pi*e_0)).*((xf - posx(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_xi74 = (-rho_i4/(2*pi*e_0)).*((xf - posx(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_xi75 = (-rho_i5/(2*pi*e_0)).*((xf - posx(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_xi76 = (-rho_i6/(2*pi*e_0)).*((xf - posx(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_xi78 = (rho_i2/(2*pi*e_0)).*((xf - posx(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_xi79 = (rho_i3/(2*pi*e_0)).*((xf - posx(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_xi710 = (rho_i4/(2*pi*e_0)).*((xf - posx(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_xi711 = (rho_i5/(2*pi*e_0)).*((xf - posx(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_xi712 = (rho_i6/(2*pi*e_0)).*((xf - posx(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));

E_xi81 = (-rho_i1/(2*pi*e_0)).*((xf - posx(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_xi82 = (-rho_i2/(2*pi*e_0)).*((xf - posx(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_xi83 = (-rho_i3/(2*pi*e_0)).*((xf - posx(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_xi84 = (-rho_i4/(2*pi*e_0)).*((xf - posx(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_xi85 = (-rho_i5/(2*pi*e_0)).*((xf - posx(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_xi86 = (-rho_i6/(2*pi*e_0)).*((xf - posx(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_xi87 = (rho_i1/(2*pi*e_0)).*((xf - posx(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_xi89 = (rho_i3/(2*pi*e_0)).*((xf - posx(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_xi810 = (rho_i4/(2*pi*e_0)).*((xf - posx(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_xi811 = (rho_i5/(2*pi*e_0)).*((xf - posx(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_xi812 = (rho_i6/(2*pi*e_0)).*((xf - posx(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));

E_xi91 = (-rho_i1/(2*pi*e_0)).*((xf - posx(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_xi92 = (-rho_i2/(2*pi*e_0)).*((xf - posx(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_xi93 = (-rho_i3/(2*pi*e_0)).*((xf - posx(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_xi94 = (-rho_i4/(2*pi*e_0)).*((xf - posx(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_xi95 = (-rho_i5/(2*pi*e_0)).*((xf - posx(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_xi96 = (-rho_i6/(2*pi*e_0)).*((xf - posx(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_xi97 = (rho_i1/(2*pi*e_0)).*((xf - posx(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_xi98 = (rho_i2/(2*pi*e_0)).*((xf - posx(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_xi910 = (rho_i4/(2*pi*e_0)).*((xf - posx(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_xi911 = (rho_i5/(2*pi*e_0)).*((xf - posx(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_xi912 = (rho_i6/(2*pi*e_0)).*((xf - posx(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));

E_xi101 = (-rho_i1/(2*pi*e_0)).*((xf - posx(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_xi102 = (-rho_i2/(2*pi*e_0)).*((xf - posx(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_xi103 = (-rho_i3/(2*pi*e_0)).*((xf - posx(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_xi104 = (-rho_i4/(2*pi*e_0)).*((xf - posx(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_xi105 = (-rho_i5/(2*pi*e_0)).*((xf - posx(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_xi106 = (-rho_i6/(2*pi*e_0)).*((xf - posx(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_xi107 = (rho_i1/(2*pi*e_0)).*((xf - posx(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_xi108 = (rho_i2/(2*pi*e_0)).*((xf - posx(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_xi109 = (rho_i3/(2*pi*e_0)).*((xf - posx(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_xi1011 = (rho_i5/(2*pi*e_0)).*((xf - posx(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_xi1012 = (rho_i6/(2*pi*e_0)).*((xf - posx(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));

E_xi11_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_xi11_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_xi113 = (-rho_i3/(2*pi*e_0)).*((xf - posx(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_xi114 = (-rho_i4/(2*pi*e_0)).*((xf - posx(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_xi115 = (-rho_i5/(2*pi*e_0)).*((xf - posx(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_xi116 = (-rho_i6/(2*pi*e_0)).*((xf - posx(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_xi117 = (rho_i1/(2*pi*e_0)).*((xf - posx(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_xi118 = (rho_i2/(2*pi*e_0)).*((xf - posx(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_xi119 = (rho_i3/(2*pi*e_0)).*((xf - posx(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_xi1110 = (rho_i4/(2*pi*e_0)).*((xf - posx(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_xi1112 = (rho_i6/(2*pi*e_0)).*((xf - posx(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));

E_xi12_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_xi12_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_xi123 = (-rho_i3/(2*pi*e_0)).*((xf - posx(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_xi124 = (-rho_i4/(2*pi*e_0)).*((xf - posx(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_xi125 = (-rho_i5/(2*pi*e_0)).*((xf - posx(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_xi126 = (-rho_i6/(2*pi*e_0)).*((xf - posx(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_xi127 = (rho_i1/(2*pi*e_0)).*((xf - posx(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_xi128 = (rho_i2/(2*pi*e_0)).*((xf - posx(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_xi129 = (rho_i3/(2*pi*e_0)).*((xf - posx(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_xi1210 = (rho_i4/(2*pi*e_0)).*((xf - posx(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_xi1211 = (rho_i5/(2*pi*e_0)).*((xf - posx(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));

%% E_yr componente y real campo elétrico condutor 2 fase b

E_yr12 = (-rho_r2/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yr13 = (-rho_r3/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yr14 = (-rho_r4/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yr15 = (-rho_r5/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yr16 = (-rho_r6/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_yr17 = (rho_r1/(2*pi*e_0)).*((yf - posy(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_yr18 = (rho_r2/(2*pi*e_0)).*((yf - posy(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_yr19 = (rho_r3/(2*pi*e_0)).*((yf - posy(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_yr110 = (rho_r4/(2*pi*e_0)).*((yf - posy(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_yr1_11 = (rho_r5/(2*pi*e_0)).*((yf - posy(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_yr1_12 = (rho_r6/(2*pi*e_0)).*((yf - posy(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));

E_yr21 = (-rho_r1/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yr23 = (-rho_r3/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yr24 = (-rho_r4/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yr25 = (-rho_r5/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yr26 = (-rho_r6/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_yr27 = (rho_r1/(2*pi*e_0)).*((yf - posy(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_yr28 = (rho_r2/(2*pi*e_0)).*((yf - posy(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_yr29 = (rho_r3/(2*pi*e_0)).*((yf - posy(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_yr210 = (rho_r4/(2*pi*e_0)).*((yf - posy(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_yr2_11 = (rho_r5/(2*pi*e_0)).*((yf - posy(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_yr2_12 = (rho_r6/(2*pi*e_0)).*((yf - posy(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));

E_yr31 = (-rho_r1/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yr32 = (-rho_r2/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yr34 = (-rho_r4/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yr35 = (-rho_r5/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yr36 = (-rho_r6/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_yr37 = (rho_r1/(2*pi*e_0)).*((yf - posy(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_yr38 = (rho_r2/(2*pi*e_0)).*((yf - posy(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_yr39 = (rho_r3/(2*pi*e_0)).*((yf - posy(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_yr310 = (rho_r4/(2*pi*e_0)).*((yf - posy(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_yr311 = (rho_r5/(2*pi*e_0)).*((yf - posy(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_yr312 = (rho_r6/(2*pi*e_0)).*((yf - posy(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));

E_yr41 = (-rho_r1/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yr42 = (-rho_r2/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yr43 = (-rho_r3/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yr45 = (-rho_r5/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yr46 = (-rho_r6/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_yr47 = (rho_r1/(2*pi*e_0)).*((yf - posy(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_yr48 = (rho_r2/(2*pi*e_0)).*((yf - posy(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_yr49 = (rho_r3/(2*pi*e_0)).*((yf - posy(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_yr410 = (rho_r4/(2*pi*e_0)).*((yf - posy(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_yr411 = (rho_r5/(2*pi*e_0)).*((yf - posy(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_yr412 = (rho_r6/(2*pi*e_0)).*((yf - posy(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));

E_yr51 = (-rho_r1/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yr52 = (-rho_r2/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yr53 = (-rho_r3/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yr54 = (-rho_r4/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yr56 = (-rho_r6/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_yr57 = (rho_r1/(2*pi*e_0)).*((yf - posy(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_yr58 = (rho_r2/(2*pi*e_0)).*((yf - posy(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_yr59 = (rho_r3/(2*pi*e_0)).*((yf - posy(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_yr510 = (rho_r4/(2*pi*e_0)).*((yf - posy(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_yr511 = (rho_r5/(2*pi*e_0)).*((yf - posy(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_yr512 = (rho_r6/(2*pi*e_0)).*((yf - posy(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));

E_yr61 = (-rho_r1/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yr62 = (-rho_r2/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yr63 = (-rho_r3/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yr64 = (-rho_r4/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yr65 = (-rho_r5/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_yr67 = (rho_r1/(2*pi*e_0)).*((yf - posy(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_yr68 = (rho_r2/(2*pi*e_0)).*((yf - posy(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_yr69 = (rho_r3/(2*pi*e_0)).*((yf - posy(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_yr610 = (rho_r4/(2*pi*e_0)).*((yf - posy(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_yr611 = (rho_r5/(2*pi*e_0)).*((yf - posy(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_yr612 = (rho_r6/(2*pi*e_0)).*((yf - posy(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));

E_yr71 = (-rho_r1/(2*pi*e_0)).*((yf - posy(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_yr72 = (-rho_r2/(2*pi*e_0)).*((yf - posy(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_yr73 = (-rho_r3/(2*pi*e_0)).*((yf - posy(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_yr74 = (-rho_r4/(2*pi*e_0)).*((yf - posy(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_yr75 = (-rho_r5/(2*pi*e_0)).*((yf - posy(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_yr76 = (-rho_r6/(2*pi*e_0)).*((yf - posy(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_yr78 = (rho_r2/(2*pi*e_0)).*((yf - posy(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_yr79 = (rho_r3/(2*pi*e_0)).*((yf - posy(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_yr710 = (rho_r4/(2*pi*e_0)).*((yf - posy(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_yr711 = (rho_r5/(2*pi*e_0)).*((yf - posy(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_yr712 = (rho_r6/(2*pi*e_0)).*((yf - posy(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));

E_yr81 = (-rho_r1/(2*pi*e_0)).*((yf - posy(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_yr82 = (-rho_r2/(2*pi*e_0)).*((yf - posy(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_yr83 = (-rho_r3/(2*pi*e_0)).*((yf - posy(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_yr84 = (-rho_r4/(2*pi*e_0)).*((yf - posy(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_yr85 = (-rho_r5/(2*pi*e_0)).*((yf - posy(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_yr86 = (-rho_r6/(2*pi*e_0)).*((yf - posy(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_yr87 = (rho_r1/(2*pi*e_0)).*((yf - posy(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_yr89 = (rho_r3/(2*pi*e_0)).*((yf - posy(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_yr810 = (rho_r4/(2*pi*e_0)).*((yf - posy(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_yr811 = (rho_r5/(2*pi*e_0)).*((yf - posy(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_yr812 = (rho_r6/(2*pi*e_0)).*((yf - posy(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));

E_yr91 = (-rho_r1/(2*pi*e_0)).*((yf - posy(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_yr92 = (-rho_r2/(2*pi*e_0)).*((yf - posy(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_yr93 = (-rho_r3/(2*pi*e_0)).*((yf - posy(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_yr94 = (-rho_r4/(2*pi*e_0)).*((yf - posy(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_yr95 = (-rho_r5/(2*pi*e_0)).*((yf - posy(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_yr96 = (-rho_r6/(2*pi*e_0)).*((yf - posy(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_yr97 = (rho_r1/(2*pi*e_0)).*((yf - posy(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_yr98 = (rho_r2/(2*pi*e_0)).*((yf - posy(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_yr910 = (rho_r4/(2*pi*e_0)).*((yf - posy(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_yr911 = (rho_r5/(2*pi*e_0)).*((yf - posy(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_yr912 = (rho_r6/(2*pi*e_0)).*((yf - posy(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));

E_yr101 = (-rho_r1/(2*pi*e_0)).*((yf - posy(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_yr102 = (-rho_r2/(2*pi*e_0)).*((yf - posy(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_yr103 = (-rho_r3/(2*pi*e_0)).*((yf - posy(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_yr104 = (-rho_r4/(2*pi*e_0)).*((yf - posy(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_yr105 = (-rho_r5/(2*pi*e_0)).*((yf - posy(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_yr106 = (-rho_r6/(2*pi*e_0)).*((yf - posy(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_yr107 = (rho_r1/(2*pi*e_0)).*((yf - posy(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_yr108 = (rho_r2/(2*pi*e_0)).*((yf - posy(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_yr109 = (rho_r3/(2*pi*e_0)).*((yf - posy(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_yr1011 = (rho_r5/(2*pi*e_0)).*((yf - posy(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_yr1012 = (rho_r6/(2*pi*e_0)).*((yf - posy(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));

E_yr11_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_yr11_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_yr113 = (-rho_r3/(2*pi*e_0)).*((yf - posy(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_yr114 = (-rho_r4/(2*pi*e_0)).*((yf - posy(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_yr115 = (-rho_r5/(2*pi*e_0)).*((yf - posy(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_yr116 = (-rho_r6/(2*pi*e_0)).*((yf - posy(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_yr117 = (rho_r1/(2*pi*e_0)).*((yf - posy(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_yr118 = (rho_r2/(2*pi*e_0)).*((yf - posy(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_yr119 = (rho_r3/(2*pi*e_0)).*((yf - posy(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_yr1110 = (rho_r4/(2*pi*e_0)).*((yf - posy(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_yr1112 = (rho_r6/(2*pi*e_0)).*((yf - posy(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));

E_yr12_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_yr12_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_yr123 = (-rho_r3/(2*pi*e_0)).*((yf - posy(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_yr124 = (-rho_r4/(2*pi*e_0)).*((yf - posy(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_yr125 = (-rho_r5/(2*pi*e_0)).*((yf - posy(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_yr126 = (-rho_r6/(2*pi*e_0)).*((yf - posy(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_yr127 = (rho_r1/(2*pi*e_0)).*((yf - posy(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_yr128 = (rho_r2/(2*pi*e_0)).*((yf - posy(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_yr129 = (rho_r3/(2*pi*e_0)).*((yf - posy(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_yr1210 = (rho_r4/(2*pi*e_0)).*((yf - posy(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_yr1211 = (rho_r5/(2*pi*e_0)).*((yf - posy(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));


%% E_yi21 componente y imaginario campo elétrico condutor 2 fase b

E_yi12 = (-rho_i2/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yi13 = (-rho_i3/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yi14 = (-rho_i4/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yi15 = (-rho_i5/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yi16 = (-rho_i6/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_yi17 = (rho_i1/(2*pi*e_0)).*((yf - posy(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_yi18 = (rho_i2/(2*pi*e_0)).*((yf - posy(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_yi19 = (rho_i3/(2*pi*e_0)).*((yf - posy(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_yi110 = (rho_i4/(2*pi*e_0)).*((yf - posy(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_yi1_11 = (rho_i5/(2*pi*e_0)).*((yf - posy(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_yi1_12 = (rho_i6/(2*pi*e_0)).*((yf - posy(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));

E_yi21 = (-rho_i1/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yi23 = (-rho_i3/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yi24 = (-rho_i4/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yi25 = (-rho_i5/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yi26 = (-rho_i6/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_yi27 = (rho_i1/(2*pi*e_0)).*((yf - posy(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_yi28 = (rho_i2/(2*pi*e_0)).*((yf - posy(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_yi29 = (rho_i3/(2*pi*e_0)).*((yf - posy(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_yi210 = (rho_i4/(2*pi*e_0)).*((yf - posy(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_yi2_11 = (rho_i5/(2*pi*e_0)).*((yf - posy(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_yi2_12 = (rho_i6/(2*pi*e_0)).*((yf - posy(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));

E_yi31 = (-rho_i1/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yi32 = (-rho_i2/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yi34 = (-rho_i4/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yi35 = (-rho_i5/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yi36 = (-rho_i6/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_yi37 = (rho_i1/(2*pi*e_0)).*((yf - posy(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_yi38 = (rho_i2/(2*pi*e_0)).*((yf - posy(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_yi39 = (rho_i3/(2*pi*e_0)).*((yf - posy(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_yi310 = (rho_i4/(2*pi*e_0)).*((yf - posy(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_yi311 = (rho_i5/(2*pi*e_0)).*((yf - posy(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_yi312 = (rho_i6/(2*pi*e_0)).*((yf - posy(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));

E_yi41 = (-rho_i1/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yi42 = (-rho_i2/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yi43 = (-rho_i3/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yi45 = (-rho_i5/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yi46 = (-rho_i6/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_yi47 = (rho_i1/(2*pi*e_0)).*((yf - posy(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_yi48 = (rho_i2/(2*pi*e_0)).*((yf - posy(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_yi49 = (rho_i3/(2*pi*e_0)).*((yf - posy(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_yi410 = (rho_i4/(2*pi*e_0)).*((yf - posy(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_yi411 = (rho_i5/(2*pi*e_0)).*((yf - posy(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_yi412 = (rho_i6/(2*pi*e_0)).*((yf - posy(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));

E_yi51 = (-rho_i1/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yi52 = (-rho_i2/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yi53 = (-rho_i3/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yi54 = (-rho_i4/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yi56 = (-rho_i6/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_yi57 = (rho_i1/(2*pi*e_0)).*((yf - posy(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_yi58 = (rho_i2/(2*pi*e_0)).*((yf - posy(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_yi59 = (rho_i3/(2*pi*e_0)).*((yf - posy(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_yi510 = (rho_i4/(2*pi*e_0)).*((yf - posy(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_yi511 = (rho_i5/(2*pi*e_0)).*((yf - posy(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_yi512 = (rho_i6/(2*pi*e_0)).*((yf - posy(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));

E_yi61 = (-rho_i1/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yi62 = (-rho_i2/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yi63 = (-rho_i3/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yi64 = (-rho_i4/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yi65 = (-rho_i5/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_yi67 = (rho_i1/(2*pi*e_0)).*((yf - posy(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_yi68 = (rho_i2/(2*pi*e_0)).*((yf - posy(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_yi69 = (rho_i3/(2*pi*e_0)).*((yf - posy(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_yi610 = (rho_i4/(2*pi*e_0)).*((yf - posy(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_yi611 = (rho_i5/(2*pi*e_0)).*((yf - posy(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_yi612 = (rho_i6/(2*pi*e_0)).*((yf - posy(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));

E_yi71 = (-rho_i1/(2*pi*e_0)).*((yf - posy(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_yi72 = (-rho_i2/(2*pi*e_0)).*((yf - posy(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_yi73 = (-rho_i3/(2*pi*e_0)).*((yf - posy(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_yi74 = (-rho_i4/(2*pi*e_0)).*((yf - posy(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_yi75 = (-rho_i5/(2*pi*e_0)).*((yf - posy(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_yi76 = (-rho_i6/(2*pi*e_0)).*((yf - posy(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_yi78 = (rho_i2/(2*pi*e_0)).*((yf - posy(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_yi79 = (rho_i3/(2*pi*e_0)).*((yf - posy(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_yi710 = (rho_i4/(2*pi*e_0)).*((yf - posy(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_yi711 = (rho_i5/(2*pi*e_0)).*((yf - posy(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_yi712 = (rho_i6/(2*pi*e_0)).*((yf - posy(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));

E_yi81 = (-rho_i1/(2*pi*e_0)).*((yf - posy(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_yi82 = (-rho_i2/(2*pi*e_0)).*((yf - posy(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_yi83 = (-rho_i3/(2*pi*e_0)).*((yf - posy(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_yi84 = (-rho_i4/(2*pi*e_0)).*((yf - posy(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_yi85 = (-rho_i5/(2*pi*e_0)).*((yf - posy(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_yi86 = (-rho_i6/(2*pi*e_0)).*((yf - posy(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_yi87 = (rho_i1/(2*pi*e_0)).*((yf - posy(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_yi89 = (rho_i3/(2*pi*e_0)).*((yf - posy(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_yi810 = (rho_i4/(2*pi*e_0)).*((yf - posy(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_yi811 = (rho_i5/(2*pi*e_0)).*((yf - posy(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_yi812 = (rho_i6/(2*pi*e_0)).*((yf - posy(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));

E_yi91 = (-rho_i1/(2*pi*e_0)).*((yf - posy(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_yi92 = (-rho_i2/(2*pi*e_0)).*((yf - posy(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_yi93 = (-rho_i3/(2*pi*e_0)).*((yf - posy(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_yi94 = (-rho_i4/(2*pi*e_0)).*((yf - posy(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_yi95 = (-rho_i5/(2*pi*e_0)).*((yf - posy(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_yi96 = (-rho_i6/(2*pi*e_0)).*((yf - posy(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_yi97 = (rho_i1/(2*pi*e_0)).*((yf - posy(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_yi98 = (rho_i2/(2*pi*e_0)).*((yf - posy(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_yi910 = (rho_i4/(2*pi*e_0)).*((yf - posy(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_yi911 = (rho_i5/(2*pi*e_0)).*((yf - posy(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_yi912 = (rho_i6/(2*pi*e_0)).*((yf - posy(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));

E_yi101 = (-rho_i1/(2*pi*e_0)).*((yf - posy(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_yi102 = (-rho_i2/(2*pi*e_0)).*((yf - posy(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_yi103 = (-rho_i3/(2*pi*e_0)).*((yf - posy(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_yi104 = (-rho_i4/(2*pi*e_0)).*((yf - posy(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_yi105 = (-rho_i5/(2*pi*e_0)).*((yf - posy(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_yi106 = (-rho_i6/(2*pi*e_0)).*((yf - posy(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_yi107 = (rho_i1/(2*pi*e_0)).*((yf - posy(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_yi108 = (rho_i2/(2*pi*e_0)).*((yf - posy(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_yi109 = (rho_i3/(2*pi*e_0)).*((yf - posy(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_yi1011 = (rho_i5/(2*pi*e_0)).*((yf - posy(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_yi1012 = (rho_i6/(2*pi*e_0)).*((yf - posy(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));

E_yi11_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_yi11_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_yi113 = (-rho_i3/(2*pi*e_0)).*((yf - posy(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_yi114 = (-rho_i4/(2*pi*e_0)).*((yf - posy(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_yi115 = (-rho_i5/(2*pi*e_0)).*((yf - posy(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_yi116 = (-rho_i6/(2*pi*e_0)).*((yf - posy(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_yi117 = (rho_i1/(2*pi*e_0)).*((yf - posy(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_yi118 = (rho_i2/(2*pi*e_0)).*((yf - posy(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_yi119 = (rho_i3/(2*pi*e_0)).*((yf - posy(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_yi1110 = (rho_i4/(2*pi*e_0)).*((yf - posy(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_yi1112 = (rho_i6/(2*pi*e_0)).*((yf - posy(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));

E_yi12_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_yi12_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_yi123 = (-rho_i3/(2*pi*e_0)).*((yf - posy(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_yi124 = (-rho_i4/(2*pi*e_0)).*((yf - posy(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_yi125 = (-rho_i5/(2*pi*e_0)).*((yf - posy(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_yi126 = (-rho_i6/(2*pi*e_0)).*((yf - posy(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_yi127 = (rho_i1/(2*pi*e_0)).*((yf - posy(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_yi128 = (rho_i2/(2*pi*e_0)).*((yf - posy(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_yi129 = (rho_i3/(2*pi*e_0)).*((yf - posy(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_yi1210 = (rho_i4/(2*pi*e_0)).*((yf - posy(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_yi1211 = (rho_i5/(2*pi*e_0)).*((yf - posy(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));

%% aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr12 ; E_xr13 ; E_xr14 ; E_xr15 ; E_xr16 ; E_xr17 ; E_xr18 ; E_xr19 ; E_xr110 ; E_xr1_11 ; E_xr1_12 ; E_xr21 ; E_xr23 ; E_xr24 ; E_xr25 ; E_xr26 ; E_xr27 ; E_xr28 ; E_xr29 ; E_xr210 ; E_xr2_11 ; E_xr2_12 ; E_xr31 ; E_xr32 ; E_xr34 ; E_xr35 ; E_xr36 ; E_xr37 ; E_xr38 ; E_xr39 ; E_xr310 ; E_xr311 ; E_xr312 ; E_xr41 ; E_xr42 ; E_xr43 ; E_xr45 ; E_xr46 ; E_xr47 ; E_xr48 ; E_xr49 ; E_xr410 ; E_xr411 ; E_xr412 ; E_xr51 ; E_xr52 ; E_xr53 ; E_xr54 ; E_xr56 ; E_xr57 ; E_xr58 ; E_xr59 ; E_xr510 ; E_xr511 ; E_xr512 ; E_xr61 ; E_xr62 ; E_xr63 ; E_xr64 ; E_xr65 ; E_xr67 ; E_xr68 ; E_xr69 ; E_xr610 ; E_xr611 ; E_xr612 ; E_xr71 ; E_xr72 ; E_xr73 ; E_xr74 ; E_xr75 ; E_xr76 ; E_xr78 ; E_xr79 ; E_xr710 ; E_xr711 ; E_xr712 ; E_xr81 ; E_xr82 ; E_xr83 ; E_xr84 ; E_xr85 ; E_xr86 ; E_xr87 ; E_xr89 ; E_xr810 ; E_xr811 ; E_xr812 ; E_xr91 ; E_xr92 ; E_xr93 ; E_xr94 ; E_xr95 ; E_xr96 ; E_xr97 ; E_xr98 ; E_xr910 ; E_xr911 ; E_xr912 ; E_xr101 ; E_xr102 ; E_xr103 ; E_xr104 ; E_xr105 ; E_xr106 ; E_xr107 ; E_xr108 ; E_xr109 ; E_xr1011 ; E_xr1012 ; E_xr11_1 ; E_xr11_2 ; E_xr113 ; E_xr114 ; E_xr115 ; E_xr116 ; E_xr117 ; E_xr118 ; E_xr119 ; E_xr1110 ; E_xr1112 ; E_xr12_1 ; E_xr12_2 ; E_xr123 ; E_xr124 ; E_xr125 ; E_xr126 ; E_xr127 ; E_xr128 ; E_xr129 ; E_xr1210 ; E_xr1211 ];

E_xi = [ E_xi12 ; E_xi13 ; E_xi14 ; E_xi15 ; E_xi16 ; E_xi17 ; E_xi18 ; E_xi19 ; E_xi110 ; E_xi1_11 ; E_xi1_12 ; E_xi21 ; E_xi23 ; E_xi24 ; E_xi25 ; E_xi26 ; E_xi27 ; E_xi28 ; E_xi29 ; E_xi210 ; E_xi2_11 ; E_xi2_12 ; E_xi31 ; E_xi32 ; E_xi34 ; E_xi35 ; E_xi36 ; E_xi37 ; E_xi38 ; E_xi39 ; E_xi310 ; E_xi311 ; E_xi312 ; E_xi41 ; E_xi42 ; E_xi43 ; E_xi45 ; E_xi46 ; E_xi47 ; E_xi48 ; E_xi49 ; E_xi410 ; E_xi411 ; E_xi412 ; E_xi51 ; E_xi52 ; E_xi53 ; E_xi54 ; E_xi56 ; E_xi57 ; E_xi58 ; E_xi59 ; E_xi510 ; E_xi511 ; E_xi512 ; E_xi61 ; E_xi62 ; E_xi63 ; E_xi64 ; E_xi65 ; E_xi67 ; E_xi68 ; E_xi69 ; E_xi610 ; E_xi611 ; E_xi612 ; E_xi71 ; E_xi72 ; E_xi73 ; E_xi74 ; E_xi75 ; E_xi76 ; E_xi78 ; E_xi79 ; E_xi710 ; E_xi711 ; E_xi712 ; E_xi81 ; E_xi82 ; E_xi83 ; E_xi84 ; E_xi85 ; E_xi86 ; E_xi87 ; E_xi89 ; E_xi810 ; E_xi811 ; E_xi812 ; E_xi91 ; E_xi92 ; E_xi93 ; E_xi94 ; E_xi95 ; E_xi96 ; E_xi97 ; E_xi98 ; E_xi910 ; E_xi911 ; E_xi912 ; E_xi101 ; E_xi102 ; E_xi103 ; E_xi104 ; E_xi105 ; E_xi106 ; E_xi107 ; E_xi108 ; E_xi109 ; E_xi1011 ; E_xi1012 ; E_xi11_1 ; E_xi11_2 ; E_xi113 ; E_xi114 ; E_xi115 ; E_xi116 ; E_xi117 ; E_xi118 ; E_xi119 ; E_xi1110 ; E_xi1112 ; E_xi12_1 ; E_xi12_2 ; E_xi123 ; E_xi124 ; E_xi125 ; E_xi126 ; E_xi127 ; E_xi128 ; E_xi129 ; E_xi1210 ; E_xi1211 ];

E_yr = [ E_yr12 ; E_yr13 ; E_yr14 ; E_yr15 ; E_yr16 ; E_yr17 ; E_yr18 ; E_yr19 ; E_yr110 ; E_yr1_11 ; E_yr1_12 ; E_yr21 ; E_yr23 ; E_yr24 ; E_yr25 ; E_yr26 ; E_yr27 ; E_yr28 ; E_yr29 ; E_yr210 ; E_yr2_11 ; E_yr2_12 ; E_yr31 ; E_yr32 ; E_yr34 ; E_yr35 ; E_yr36 ; E_yr37 ; E_yr38 ; E_yr39 ; E_yr310 ; E_yr311 ; E_yr312 ; E_yr41 ; E_yr42 ; E_yr43 ; E_yr45 ; E_yr46 ; E_yr47 ; E_yr48 ; E_yr49 ; E_yr410 ; E_yr411 ; E_yr412 ; E_yr51 ; E_yr52 ; E_yr53 ; E_yr54 ; E_yr56 ; E_yr57 ; E_yr58 ; E_yr59 ; E_yr510 ; E_yr511 ; E_yr512 ; E_yr61 ; E_yr62 ; E_yr63 ; E_yr64 ; E_yr65 ; E_yr67 ; E_yr68 ; E_yr69 ; E_yr610 ; E_yr611 ; E_yr612 ; E_yr71 ; E_yr72 ; E_yr73 ; E_yr74 ; E_yr75 ; E_yr76 ; E_yr78 ; E_yr79 ; E_yr710 ; E_yr711 ; E_yr712 ; E_yr81 ; E_yr82 ; E_yr83 ; E_yr84 ; E_yr85 ; E_yr86 ; E_yr87 ; E_yr89 ; E_yr810 ; E_yr811 ; E_yr812 ; E_yr91 ; E_yr92 ; E_yr93 ; E_yr94 ; E_yr95 ; E_yr96 ; E_yr97 ; E_yr98 ; E_yr910 ; E_yr911 ; E_yr912 ; E_yr101 ; E_yr102 ; E_yr103 ; E_yr104 ; E_yr105 ; E_yr106 ; E_yr107 ; E_yr108 ; E_yr109 ; E_yr1011 ; E_yr1012 ; E_yr11_1 ; E_yr11_2 ; E_yr113 ; E_yr114 ; E_yr115 ; E_yr116 ; E_yr117 ; E_yr118 ; E_yr119 ; E_yr1110 ; E_yr1112 ; E_yr12_1 ; E_yr12_2 ; E_yr123 ; E_yr124 ; E_yr125 ; E_yr126 ; E_yr127 ; E_yr128 ; E_yr129 ; E_yr1210 ; E_yr1211 ];

E_yi = [ E_yi12 ; E_yi13 ; E_yi14 ; E_yi15 ; E_yi16 ; E_yi17 ; E_yi18 ; E_yi19 ; E_yi110 ; E_yi1_11 ; E_yi1_12 ; E_yi21 ; E_yi23 ; E_yi24 ; E_yi25 ; E_yi26 ; E_yi27 ; E_yi28 ; E_yi29 ; E_yi210 ; E_yi2_11 ; E_yi2_12 ; E_yi31 ; E_yi32 ; E_yi34 ; E_yi35 ; E_yi36 ; E_yi37 ; E_yi38 ; E_yi39 ; E_yi310 ; E_yi311 ; E_yi312 ; E_yi41 ; E_yi42 ; E_yi43 ; E_yi45 ; E_yi46 ; E_yi47 ; E_yi48 ; E_yi49 ; E_yi410 ; E_yi411 ; E_yi412 ; E_yi51 ; E_yi52 ; E_yi53 ; E_yi54 ; E_yi56 ; E_yi57 ; E_yi58 ; E_yi59 ; E_yi510 ; E_yi511 ; E_yi512 ; E_yi61 ; E_yi62 ; E_yi63 ; E_yi64 ; E_yi65 ; E_yi67 ; E_yi68 ; E_yi69 ; E_yi610 ; E_yi611 ; E_yi612 ; E_yi71 ; E_yi72 ; E_yi73 ; E_yi74 ; E_yi75 ; E_yi76 ; E_yi78 ; E_yi79 ; E_yi710 ; E_yi711 ; E_yi712 ; E_yi81 ; E_yi82 ; E_yi83 ; E_yi84 ; E_yi85 ; E_yi86 ; E_yi87 ; E_yi89 ; E_yi810 ; E_yi811 ; E_yi812 ; E_yi91 ; E_yi92 ; E_yi93 ; E_yi94 ; E_yi95 ; E_yi96 ; E_yi97 ; E_yi98 ; E_yi910 ; E_yi911 ; E_yi912 ; E_yi101 ; E_yi102 ; E_yi103 ; E_yi104 ; E_yi105 ; E_yi106 ; E_yi107 ; E_yi108 ; E_yi109 ; E_yi1011 ; E_yi1012 ; E_yi11_1 ; E_yi11_2 ; E_yi113 ; E_yi114 ; E_yi115 ; E_yi116 ; E_yi117 ; E_yi118 ; E_yi119 ; E_yi1110 ; E_yi1112 ; E_yi12_1 ; E_yi12_2 ; E_yi123 ; E_yi124 ; E_yi125 ; E_yi126 ; E_yi127 ; E_yi128 ; E_yi129 ; E_yi1210 ; E_yi1211 ];


%% aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

Exr = sum(E_xr).^2;

Exi = sum(E_xi).^2;

Eyr = sum(E_yr).^2;

Eyi = sum(E_yi).^2;


%% aqui retira-se o módulo dos quadrados:

% V/m -> kV/cm = 10^-3/10^2
Erms = ((Exr + Exi + Eyr + Eyi).^(1/2))*(10^-5);
thetad = rad2deg(theta); %converte o ângulo de radianos para graus

fs = 0.82; %fator de superfície (constante)
pa = 760; %pressão atmosférica [mmHg]
Ecrit = fcn_supcrit(r*10^(2));
E_crit = linspace(Ecrit,Ecrit,360);

Emed = mean(Erms);
Emax = max(Erms);
Kirreg = Emax/Emed;

%%

E2rms1 = [11.6647044468313 11.6649245793722 11.6658053669753 11.6673462728508 11.6695463808756 11.6724043967054 11.6759186493382 11.6800870931197 11.6849073101901 11.6903765133595 11.6964915494145 11.7032489028422 11.7106446999682 11.7186747134860 11.7273343673999 11.7366187423226 11.7465225811653 11.7570402951728 11.7681659703053 11.7798933739560 11.7922159619765 11.8051268860145 11.8186190011336 11.8326848736994 11.8473167895302 11.8625067622768 11.8782465420364 11.8945276241543 11.9113412582322 11.9286784572872 11.9465300070850 11.9648864755925 11.9837382225507 12.0030754091582 12.0228880078259 12.0431658120092 12.0638984460736 12.0850753752132 12.1066859153731 12.1287192431671 12.1511644057962 12.1740103309162 12.1972458364799 12.2208596404972 12.2448403707389 12.2691765743401 12.2938567273081 12.3188692439162 12.3442024859636 12.3698447719155 12.3957843858747 12.4220095864101 12.4485086152082 12.4752697055570 12.5022810906303 12.5295310115976 12.5570077255146 12.5846995130278 12.6125946858434 12.6406815939974 12.6689486328963 12.6973842501236 12.7259769520339 12.7547153100921 12.7835879669933 12.8125836425395 12.8416911392742 12.8708993478783 12.9001972523354 12.9295739348405 12.9590185804767 12.9885204816629 13.0180690423415 13.0476537819541 13.0772643391694 13.1068904753835 13.1365220779935 13.1661491634412 13.1957618800469 13.2253505106118 13.2549054748078 13.2844173313757 13.3138767800870 13.3432746635336 13.3726019687089 13.4018498283925 13.4310095223582 13.4600724784029 13.4890302731860 13.5178746329178 13.5465974338652 13.5751907027058 13.6036466167349 13.6319575039141 13.6601158427823 13.6881142622260 13.7159455411234 13.7436026078650 13.7710785397404 13.7983665622269 13.8254600481595 13.8523525167971 13.8790376328002 13.9055092050954 13.9317611856767 13.9577876683049 13.9835828871399 14.0091412152915 14.0344571633108 14.0595253776182 14.0843406388626 14.1088978602332 14.1331920857305 14.1572184883722 14.1809723683791 14.2044491513035 14.2276443861455 14.2505537434194 14.2731730132104 14.2954981032015 14.3175250366835 14.3392499505473 14.3606690932745 14.3817788229011 14.4025756049952 14.4230560106194 14.4432167142906 14.4630544919547 14.4825662189643 14.5017488680448 14.5205995072965 14.5391152981959 14.5572934936113 14.5751314358369 14.5926265546432 14.6097763653530 14.6265784669366 14.6430305401273 14.6591303455576 14.6748757219500 14.6902645842969 14.7052949220939 14.7199647976086 14.7342723441600 14.7482157644502 14.7617933289236 14.7750033741638 14.7878443013193 14.8003145745837 14.8124127196917 14.8241373224711 14.8354870274321 14.8464605363896 14.8570566071355 14.8672740521431 14.8771117373287 14.8865685808374 14.8956435518932 14.9043356696676 14.9126440022226 14.9205676654705 14.9281058221987 14.9352576811268 14.9420224960204 14.9483995648424 14.9543882289543 14.9599878723603 14.9651979210108 14.9700178421348 14.9744471436299 14.9784853734997 14.9821321193372 14.9853870078522 14.9882497044544 14.9907199128778 14.9927973748518 14.9944818698271 14.9957732147444 14.9966712638540 14.9971759085819 14.9972870774397 14.9970047359940 14.9963288868776 14.9952595698418 14.9937968618685 14.9919408773262 14.9896917681719 14.9870497242044 14.9840149733610 14.9805877820667 14.9767684556331 14.9725573386957 14.9679548157088 14.9629613114833 14.9575772917713 14.9518032638963 14.9456397774431 14.9390874249747 14.9321468428099 14.9248187118418 14.9171037584087 14.9090027551993 14.9005165222175 14.8916459277767 14.8823918895611 14.8727553757019 14.8627374059291 14.8523390527407 14.8415614426364 14.8304057573750 14.8188732352899 14.8069651726344 14.7946829249788 14.7820279086395 14.7690016021447 14.7556055477571 14.7418413530106 14.7277106923023 14.7132153085133 14.6983570146660 14.6831376956090 14.6675593097544 14.6516238908253 14.6353335496379 14.6186904759376 14.6016969402235 14.5843552956270 14.5666679798112 14.5486375168900 14.5302665193668 14.5115576900968 14.4925138242677 14.4731378113962 14.4534326373402 14.4334013863091 14.4130472429100 14.3923734941814 14.3713835316309 14.3500808532858 14.3284690657420 14.3065518861969 14.2843331444954 14.2618167851495 14.2390068693556 14.2159075769931 14.1925232086068 14.1688581873582 14.1449170609671 14.1207045036036 14.0962253177662 14.0714844361030 14.0464869232130 14.0212379773804 13.9957429322692 13.9700072585689 13.9440365655726 13.9178366026889 13.8914132608966 13.8647725741156 13.8379207205112 13.8108640236963 13.7836089538694 13.7561621288393 13.7285303149611 13.7007204279677 13.6727395336826 13.6445948486342 13.6162937405351 13.5878437286388 13.5592524839688 13.5305278294072 13.5016777396450 13.4727103409738 13.4436339109300 13.4144568777871 13.3851878198644 13.3558354646854 13.3264086879454 13.2969165122999 13.2673681059760 13.2377727811870 13.2081399923385 13.1784793340611 13.1488005390069 13.1191134754543 13.0894281446994 13.0597546782159 13.0301033346047 13.0004844963129 12.9709086661270 12.9413864634361 12.9119286202491 12.8825459769982 12.8532494780828 12.8240501671835 12.7949591823442 12.7659877507973 12.7371471835577 12.7084488697781 12.6799042708616 12.6515249143315 12.6233223874803 12.5953083307647 12.5674944309874 12.5398924142440 12.5125140386416 12.4853710868167 12.4584753582180 12.4318386612032 12.4054728049156 12.3793895909865 12.3536008050377 12.3281182080156 12.3029535273492 12.2781184479635 12.2536246031253 12.2294835651681 12.2057068360835 12.1823058379982 12.1592919035547 12.1366762661924 12.1144700503716 12.0926842617149 12.0713297771242 12.0504173348455 12.0299575245343 12.0099607773178 11.9904373558580 11.9713973444705 11.9528506392746 11.9348069384216 11.9172757323940 11.9002662944148 11.8837876709778 11.8678486724965 11.8524578641272 11.8376235567429 11.8233537981138 11.8096563642748 11.7965387511278 11.7840081662731 11.7720715211042 11.7607354231702 11.7500061688191 11.7398897361542 11.7303917782977 11.7215176169904 11.7132722365374 11.7056602780998 11.6986860343776 11.6923534446487 11.6866660902198 11.6816271902697 11.6772395981068 11.6735057978425 11.6704279014978 11.6680076465366 11.6662463938417 11.6651451261320 11.6647044468313];
E2rms2 = [15.2218745436352 15.2215570352969 15.2208199747662 15.2196634636016 15.2180876493585 15.2160927256859 15.2136789324747 15.2108465560485 15.2075959294037 15.2039274324865 15.1998414925220 15.1953385843844 15.1904192310146 15.1850840038699 15.1793335234437 15.1731684597973 15.1665895331619 15.1595975145712 15.1521932265387 15.1443775437848 15.1361513939962 15.1275157586402 15.1184716738169 15.1090202311479 15.0991625787169 15.0888999220437 15.0782335251114 15.0671647114163 15.0556948650796 15.0438254319763 15.0315579209297 15.0188939049251 15.0058350223640 14.9923829783701 14.9785395461160 14.9643065682023 14.9496859580500 14.9346797013566 14.9192898575682 14.9035185613813 14.8873680242964 14.8708405361784 14.8539384668735 14.8366642678290 14.8190204737623 14.8010097043428 14.7826346659075 14.7638981531968 14.7448030511060 14.7253523364767 14.7055490798835 14.6853964474555 14.6648977027035 14.6440562083689 14.6228754282695 14.6013589291771 14.5795103826745 14.5573335670486 14.5348323691517 14.5120107862913 14.4888729281014 14.4654230184088 14.4416653971017 14.4176045219691 14.3932449705408 14.3685914419032 14.3436487584923 14.3184218678618 14.2929158444327 14.2671358911953 14.2410873413828 14.2147756601143 14.1882064459712 14.1613854325522 14.1343184899567 14.1070116262197 14.0794709886839 14.0517028653033 14.0237136858879 13.9955100232581 13.9670985943213 13.9384862610838 13.9096800315374 13.8806870604837 13.8515146502417 13.8221702512478 13.7926614625503 13.7629960321943 13.7331818574642 13.7032269850241 13.6731396109046 13.6429280803582 13.6126008875795 13.5821666752650 13.5516342340206 13.5210125016110 13.4903105620415 13.4595376444767 13.4287031219656 13.3978165100050 13.3668874649052 13.3359257819648 13.3049413934603 13.2739443664135 13.2429449001809 13.2119533238128 13.1809800932124 13.1500357880651 13.1191311085626 13.0882768718957 13.0574840085174 13.0267635581795 12.9961266657477 12.9655845767645 12.9351486327985 12.9048302665361 12.8746409966601 12.8445924224704 12.8146962182812 12.7849641275789 12.7554079569460 12.7260395697454 12.6968708795897 12.6679138435611 12.6391804552269 12.6106827374221 12.5824327348139 12.5544425062632 12.5267241169830 12.4992896304754 12.4721511003023 12.4453205616556 12.4188100227539 12.3926314560696 12.3667967893946 12.3413178967630 12.3162065892291 12.2914746055165 12.2671336025462 12.2431951458823 12.2196707000541 12.1965716188198 12.1739091353667 12.1516943524460 12.1299382324837 12.1086515876678 12.0878450700305 12.0675291615327 12.0477141641917 12.0284101902322 12.0096271523103 11.9913747538120 11.9736624792394 11.9564995847157 11.9398950886078 11.9238577623087 11.9083961211660 11.8935184156078 11.8792326224437 11.8655464364055 11.8524672618952 11.8400022049954 11.8281580657330 11.8169413306310 11.8063581655462 11.7964144088197 11.7871155647481 11.7784667974013 11.7704729247766 11.7631384133231 11.7564673728421 11.7504635517710 11.7451303328599 11.7404707292613 11.7364873810269 11.7331825520269 11.7305581273016 11.7286156108443 11.7273561238257 11.7267804032601 11.7268888011142 11.7276812838744 11.7291574325568 11.7313164431625 11.7341571275890 11.7376779149857 11.7418768535505 11.7467516127687 11.7522994860814 11.7585173939889 11.7654018875731 11.7729491524303 11.7811550130149 11.7900149373759 11.7995240422785 11.8096770986988 11.8204685376914 11.8318924565917 11.8439426255646 11.8566124944752 11.8698952000702 11.8837835734464 11.8982701478074 11.9133471664716 11.9290065911450 11.9452401103991 11.9620391483892 11.9793948737479 11.9972982086746 12.0157398381728 12.0347102194449 12.0541995914065 12.0741979843173 12.0946952295012 12.1156809691390 12.1371446661360 12.1590756140105 12.1814629468268 12.2042956491291 12.2275625658717 12.2512524123231 12.2753537839485 12.2998551662240 12.3247449443893 12.3500114131403 12.3756427861983 12.4016272057967 12.4279527520461 12.4546074521703 12.4815792895974 12.5088562129077 12.5364261446194 12.5642769898050 12.5923966445338 12.6207730041131 12.6493939711632 12.6782474634706 12.7073214216379 12.7366038165301 12.7660826565027 12.7957459943964 12.8255819343273 12.8555786382271 12.8857243321679 12.9160073124478 12.9464159514409 12.9769387032087 13.0075641088873 13.0382808018166 13.0690775124547 13.0999430730370 13.1308664220185 13.1618366082656 13.1928427950248 13.2238742636651 13.2549204171904 13.2859707835199 13.3170150185634 13.3480429090658 13.3790443752515 13.4100094732389 13.4409283972751 13.4717914817492 13.5025892030195 13.5333121810496 13.5639511808464 13.5944971137360 13.6249410384459 13.6552741620233 13.6854878405853 13.7155735799117 13.7455230358843 13.7753280147650 13.8049804733394 13.8344725189282 13.8637964092444 13.8929445521481 13.9219095052594 13.9506839754598 13.9792608182893 14.0076330372344 14.0357937828996 14.0637363521167 14.0914541869271 14.1189408735005 14.1461901409750 14.1731958602057 14.1999520424557 14.2264528380170 14.2526925347685 14.2786655566857 14.3043664622779 14.3297899430041 14.3549308216219 14.3797840505050 14.4043447099331 14.4286080063325 14.4525692704990 14.4762239557954 14.4995676363253 14.5225960050834 14.5453048721049 14.5676901625870 14.5897479150103 14.6114742792577 14.6328655147118 14.6539179883821 14.6746281730000 14.6949926451490 14.7150080833751 14.7346712663276 14.7539790708972 14.7729284703757 14.7915165326207 14.8097404182553 14.8275973788625 14.8450847552196 14.8621999755477 14.8789405537820 14.8953040878743 14.9112882581098 14.9268908254710 14.9421096300064 14.9569425892539 14.9713876966712 14.9854430201186 14.9991067003704 15.0123769496423 15.0252520501829 15.0377303528776 15.0498102759052 15.0614903034155 15.0727689842566 15.0836449307457 15.0941168174556 15.1041833800698 15.1138434142519 15.1230957745788 15.1319393734952 15.1403731803216 15.1483962202961 15.1560075736668 15.1632063748219 15.1699918114568 15.1763631237964 15.1823196038489 15.1878605947106 15.1929854899144 15.1976937328086 15.2019848160070 15.2058582808497 15.2093137169347 15.2123507616815 15.2149690999415 15.2171684636501 15.2189486315307 15.2203094288349 15.2212507271334 15.2217724441479 15.2218745436352];
E2rms3 = [12.4629884061254 12.4633204158175 12.4643690071130 12.4661335764979 12.4686131135283 12.4718062020249 12.4757110217415 12.4803253505053 12.4856465668233 12.4916716529490 12.4983971984031 12.5058194039398 12.5139340859511 12.5227366812965 12.5322222525531 12.5423854936678 12.5532207360036 12.5647219547682 12.5768827758069 12.5896964827482 12.6031560244869 12.6172540229870 12.6319827813937 12.6473342924302 12.6633002470662 12.6798720434444 12.6970407960389 12.7147973450352 12.7331322659088 12.7520358791842 12.7714982603630 12.7915092499898 12.8120584638480 12.8331353032640 12.8547289654994 12.8768284542174 12.8994225899962 12.9225000208860 12.9460492329812 12.9700585609895 12.9945161987932 13.0194102099714 13.0447285382800 13.0704590180646 13.0965893845963 13.1231072843151 13.1500002849691 13.1772558856307 13.2048615265781 13.2328045990425 13.2610724547870 13.2896524155275 13.3185317821722 13.3476978438778 13.3771378869026 13.4068392032725 13.4367890992116 13.4669749033830 13.4973839748746 13.5280037109837 13.5588215547442 13.5898250022277 13.6210016095967 13.6523389999078 13.6838248696722 13.7154469951612 13.7471932384598 13.7790515532646 13.8110099904343 13.8430567032776 13.8751799525906 13.9073681114490 13.9396096697324 13.9718932384212 14.0042075536274 14.0365414803923 14.0688840162377 14.1012242944789 14.1335515873122 14.1658553086589 14.1981250167838 14.2303504167125 14.2625213623963 14.2946278586975 14.3266600631505 14.3586082875192 14.3904629991632 14.4222148222187 14.4538545385735 14.4853730886874 14.5167615722162 14.5480112484749 14.5791135367466 14.6100600164229 14.6408424269978 14.6714526679142 14.7018827982742 14.7321250364210 14.7621717593773 14.7920155021798 14.8216489570884 14.8510649726847 14.8802565528761 14.9092168557818 14.9379391925487 14.9664170260609 14.9946439695769 15.0226137852816 15.0503203827726 15.0777578174808 15.1049202890172 15.1318021394687 15.1583978516499 15.1847020472872 15.2107094851774 15.2364150592909 15.2618137968542 15.2869008563836 15.3116715257092 15.3361212199617 15.3602454795415 15.3840399680768 15.4075004703635 15.4306228902924 15.4534032487784 15.4758376816750 15.4979224376928 15.5196538763246 15.5410284657701 15.5620427808621 15.5826935010112 15.6029774081638 15.6228913847650 15.6424324117414 15.6615975665088 15.6803840209895 15.6987890396644 15.7168099776358 15.7344442787200 15.7516894735813 15.7685431778696 15.7850030904061 15.8010669914013 15.8167327406896 15.8319982760190 15.8468616113598 15.8613208352554 15.8753741092137 15.8890196661329 15.9022558087634 15.9150809082153 15.9274934025079 15.9394917951523 15.9510746537829 15.9622406088298 15.9729883522370 15.9833166362177 15.9932242720623 16.0027101289861 16.0117731330240 16.0204122659694 16.0286265643638 16.0364151185268 16.0437770716364 16.0507116188555 16.0572180065080 16.0632955312977 16.0689435395799 16.0741614266775 16.0789486362488 16.0833046596998 16.0872290356470 16.0907213494300 16.0937812326717 16.0964083628854 16.0986024631352 16.1003633017416 16.1016906920389 16.1025844921791 16.1030446049880 16.1030709778685 16.1026636027536 16.1018225161091 16.1005477989851 16.0988395771174 16.0966980210773 16.0941233464712 16.0911158141895 16.0876757307045 16.0838034484159 16.0794993660478 16.0747639290941 16.0695976303100 16.0640010102545 16.0579746578826 16.0515192111825 16.0446353578641 16.0373238360921 16.0295854352721 16.0214209968788 16.0128314153343 16.0038176389329 15.9943806708150 15.9845215699825 15.9742414523642 15.9635414919261 15.9524229218278 15.9408870356189 15.9289351884886 15.9165687985511 15.9037893481784 15.8905983853727 15.8769975251880 15.8629884511852 15.8485729169281 15.8337527475266 15.8185298412140 15.8029061709605 15.7868837861259 15.7704648141552 15.7536514622922 15.7364460193496 15.7188508574940 15.7008684340649 15.6825012934395 15.6637520689036 15.6446234845618 15.6251183572784 15.6052395986291 15.5849902168866 15.5643733190164 15.5433921126976 15.5220499083627 15.5003501212433 15.4782962734278 15.4558919959390 15.4331410308115 15.4100472331704 15.3866145733147 15.3628471388063 15.3387491365429 15.3143248948337 15.2895788654620 15.2645156257293 15.2391398804888 15.2134564641610 15.1874703427115 15.1611866156197 15.1346105177986 15.1077474214930 15.0806028381231 15.0531824201026 15.0254919625876 14.9975374051836 14.9693248335996 14.9408604812292 14.9121507306656 14.8832021151508 14.8540213199400 14.8246151835931 14.7949906991595 14.7651550152934 14.7351154372508 14.7048794277912 14.6744546079699 14.6438487578032 14.6130698168310 14.5821258845358 14.5510252206320 14.5197762452198 14.4883875387920 14.4568678420930 14.4252260558104 14.3934712401117 14.3616126140195 14.3296595545918 14.2976215959457 14.2655084280768 14.2333298954920 14.2010959956517 14.1688168772078 14.1365028380177 14.1041643229787 14.0718119216094 14.0394563654302 14.0071085251184 13.9747794074192 13.9424801518314 13.9102220270549 13.8780164271913 13.8458748677120 13.8138089811588 13.7818305126188 13.7499513149296 13.7181833436404 13.6865386517220 13.6550293840141 13.6236677714254 13.5924661248778 13.5614368289983 13.5305923355564 13.4999451566583 13.4695078576911 13.4392930500192 13.4093133834541 13.3795815384676 13.3501102182026 13.3209121402297 13.2920000281131 13.2633866027426 13.2350845734761 13.2071066290812 13.1794654284918 13.1521735913843 13.1252436885950 13.0986882323697 13.0725196664786 13.0467503561964 13.0213925781595 12.9964585101191 12.9719602206015 12.9479096584915 12.9243186425485 12.9011988508799 12.8785618103772 12.8564188861413 12.8347812709070 12.8136599744823 12.7930658132299 12.7730093995974 12.7535011317205 12.7345511831146 12.7161694924746 12.6983657536048 12.6811494054864 12.6645296225138 12.6485153049121 12.6331150693516 12.6183372397861 12.6041898385211 12.5906805775406 12.5778168501063 12.5656057226403 12.5540539269181 12.5431678525802 12.5329535399825 12.5234166733978 12.5145625745834 12.5063961967293 12.4989221188015 12.4921445402871 12.4860672763594 12.4806937534698 12.4760270053751 12.4720696696135 12.4688239844318 12.4662917861749 12.4644745071404 12.4633731739049 12.4629884061254];
E2rms4 = [16.1031120104670 16.1028687638313 16.1021917979976 16.1010811878429 16.0995370528628 16.0975595572933 16.0951489102793 16.0923053660953 16.0890292244124 16.0853208306171 16.0811805761763 16.0766088990536 16.0716062841739 16.0661732639343 16.0603104187687 16.0540183777559 16.0472978192769 16.0401494717241 16.0325741142541 16.0245725775900 16.0161457448711 16.0072945525482 15.9980199913309 15.9883231071757 15.9782050023208 15.9676668363737 15.9567098274348 15.9453352532739 15.9335444525469 15.9213388260562 15.9087198380614 15.8956890176246 15.8822479600012 15.8683983280758 15.8541418538346 15.8394803398829 15.8244156609937 15.8089497657055 15.7930846779531 15.7768224987304 15.7601654077991 15.7431156654213 15.7256756141354 15.7078476805565 15.6896343772113 15.6710383044012 15.6520621520971 15.6327087018563 15.6129808287605 15.5928815033918 15.5724137938118 15.5515808675736 15.5303859937435 15.5088325449436 15.4869239993962 15.4646639430022 15.4420560713943 15.4191041920346 15.3958122262752 15.3721842114560 15.3482243029730 15.3239367763563 15.2993260293361 15.2743965838905 15.2491530882877 15.2236003191044 15.1977431832233 15.1715867198017 15.1451361022190 15.1183966399799 15.0913737805837 15.0640731113600 15.0365003612382 15.0086614024915 14.9805622524042 14.9522090748905 14.9236081820435 14.8947660356138 14.8656892484236 14.8363845856851 14.8068589662359 14.7771194637068 14.7471733075530 14.7170278840216 14.6866907369932 14.6561695687075 14.6254722403792 14.5946067726931 14.5635813461494 14.5324043012992 14.5010841388149 14.4696295194210 14.4380492636795 14.4063523516035 14.3745479221097 14.3426452722971 14.3106538565517 14.2785832854717 14.2464433245893 14.2142438929158 14.1819950612765 14.1497070504418 14.1173902290562 14.0850551113361 14.0527123545677 14.0203727563653 13.9880472517100 13.9557469097484 13.9234829303621 13.8912666404975 13.8591094902445 13.8270230486736 13.7950189994342 13.7631091360858 13.7313053571921 13.6996196611514 13.6680641407836 13.6366509776494 13.6053924361314 13.5743008572500 13.5433886522298 13.5126682958208 13.4821523193674 13.4518533036304 13.4217838713748 13.3919566797113 13.3623844122091 13.3330797707854 13.3040554673691 13.2753242153459 13.2468987208022 13.2187916735711 13.1910157380800 13.1635835440187 13.1365076768408 13.1098006680962 13.0834749856214 13.0575430235824 13.0320170923949 13.0069094085396 12.9822320842625 12.9579971172010 12.9342163799396 12.9109016094992 12.8880643968008 12.8657161760949 12.8438682143862 12.8225316008715 12.8017172364021 12.7814358229868 12.7616978533638 12.7425136006515 12.7238931080945 12.7058461789291 12.6883823663855 12.6715109638433 12.6552409951574 12.6395812051788 12.6245400504818 12.6101256903200 12.5963459778263 12.5832084514784 12.5707203268418 12.5588884886108 12.5477194829629 12.5372195102435 12.5273944179936 12.5182496943382 12.5097904617466 12.5020214711827 12.4949470966508 12.4885713301550 12.4828977770810 12.4779296520105 12.4736697749747 12.4701205681613 12.4672840530765 12.4651618481713 12.4637551669363 12.4630648164699 12.4630911965230 12.4638342990212 12.4652937080684 12.4674686004274 12.4703577464816 12.4739595116710 12.4782718584024 12.4832923484277 12.4890181456875 12.4954460196083 12.5025723488540 12.5103931255165 12.5189039597386 12.5281000847596 12.5379763623739 12.5485272887859 12.5597470008545 12.5716292827086 12.5841675727255 12.5973549708497 12.6111842462441 12.6256478452527 12.6407378996632 12.6564462352453 12.6727643805560 12.6896835759873 12.7071947830446 12.7252886938290 12.7439557407175 12.7631861062101 12.7829697329338 12.8032963337800 12.8241554021634 12.8455362223754 12.8674278800166 12.8898192724977 12.9126991195809 12.9360559739473 12.9598782317747 12.9841541433132 13.0088718234249 13.0340192620993 13.0595843349002 13.0855548133442 13.1119183752040 13.1386626146969 13.1657750525675 13.1932431460492 13.2210542986770 13.2491958699581 13.2776551848748 13.3064195432185 13.3354762287463 13.3648125181376 13.3944156897534 13.4242730321947 13.4543718526407 13.4846994849634 13.5152432976168 13.5459907012985 13.5769291563643 13.6080461800091 13.6393293532000 13.6707663273580 13.7023448307940 13.7340526748942 13.7658777600426 13.7978080813032 13.8298317338344 13.8619369180623 13.8941119445867 13.9263452388525 13.9586253455533 13.9909409327954 14.0232807960212 14.0556338616801 14.0879891906613 14.1203359814957 14.1526635733184 14.1849614486122 14.2172192357097 14.2494267110975 14.2815738014869 14.3136505856849 14.3456472962586 14.3775543209918 14.4093622041628 14.4410616476204 14.4726435116790 14.5040988158401 14.5354187393404 14.5665946215376 14.5976179621256 14.6284804212031 14.6591738192031 14.6896901366573 14.7200215138501 14.7501602503231 14.7800988042588 14.8098297917542 14.8393459859782 14.8686403162040 14.8977058667760 14.9265358759416 14.9551237346145 14.9834629850557 15.0115473194569 15.0393705784616 15.0669267496144 15.0942099657402 15.1212145032761 15.1479347805250 15.1743653558849 15.2005009260108 15.2263363239438 15.2518665172059 15.2770866058501 15.3019918204901 15.3265775203001 15.3508391909944 15.3747724427828 15.3983730083195 15.4216367406362 15.4445596110596 15.4671377071421 15.4893672305597 15.5112444950499 15.5327659243104 15.5539280499427 15.5747275093715 15.5951610437966 15.6152254961475 15.6349178090577 15.6542350228502 15.6731742735548 15.6917327909286 15.7099078965141 15.7276970017197 15.7450976059206 15.7621072945920 15.7787237374708 15.7949446867510 15.8107679753048 15.8261915149451 15.8412132947135 15.8558313792134 15.8700439069750 15.8838490888537 15.8972452064769 15.9102306107232 15.9228037202448 15.9349630200279 15.9467070599973 15.9580344536676 15.9689438768267 15.9794340662722 15.9895038185897 15.9991519889722 16.0083774900914 16.0171792910060 16.0255564161223 16.0335079442015 16.0410330074074 16.0481307904070 16.0548005295157 16.0610415118906 16.0668530747727 16.0722346047730 16.0771855372107 16.0817053555000 16.0857935905809 16.0894498204031 16.0926736694575 16.0954648083538 16.0978229534516 16.0997478665363 16.1012393545474 16.1022972693528 16.1029215075755 16.1031120104670];
E2rms5 = [11.7267490804422 11.7269827617059 11.7279004335762 11.7295015514763 11.7317851861764 11.7347500249004 11.7383943728819 11.7427161553602 11.7477129200173 11.7533818398396 11.7597197164117 11.7667229836244 11.7743877117956 11.7827096121824 11.7916840419032 11.8013060092159 11.8115701791858 11.8224708796985 11.8340021078180 11.8461575364807 11.8589305214959 11.8723141088592 11.8863010423464 11.9008837713751 11.9160544591296 11.9318049909161 11.9481269827516 11.9650117901416 11.9824505170651 12.0004340251103 12.0189529427869 12.0379976749576 12.0575584123914 12.0776251414303 12.0981876537265 12.1192355560610 12.1407582801958 12.1627450927828 12.1851851052766 12.2080672838470 12.2313804592933 12.2551133369112 12.2792545063362 12.3037924513086 12.3287155593840 12.3540121315476 12.3796703917356 12.4056784962454 12.4320245430159 12.4586965807906 12.4856826181169 12.5129706322040 12.5405485776098 12.5684043947611 12.5965260182806 12.6249013851410 12.6535184426053 12.6823651559832 12.7114295161577 12.7406995469150 12.7701633120495 12.7998089222386 12.8296245417087 12.8595983946505 12.8897187714181 12.9199740344897 12.9503526241885 12.9808430641706 13.0114339666844 13.0421140375771 13.0728720810742 13.1036970043318 13.1345778217336 13.1655036589786 13.1964637569247 13.2274474752077 13.2584442956372 13.2894438253656 13.3204357998495 13.3514100855838 13.3823566826252 13.4132657269271 13.4441274924412 13.4749323930488 13.5056709842863 13.5363339648771 13.5669121780881 13.5973966129108 13.6277784050528 13.6580488377826 13.6881993425914 13.7182214997072 13.7481070384610 13.7778478374987 13.8074359248545 13.8368634778866 13.8661228230871 13.8952064357694 13.9241069396234 13.9528171061728 13.9813298541148 14.0096382485568 14.0377355001649 14.0656149641999 14.0932701394902 14.1206946673030 14.1478823301496 14.1748270505070 14.2015228894831 14.2279640454150 14.2541448524014 14.2800597787849 14.3057034255945 14.3310705249215 14.3561559382759 14.3809546548832 14.4054617899696 14.4296725829960 14.4535823958805 14.4771867111914 14.5004811303217 14.5234613716439 14.5461232686630 14.5684627681402 14.5904759282307 14.6121589166022 14.6335080085522 14.6545195851356 14.6751901312964 14.6955162339838 14.7154945803032 14.7351219556653 14.7543952419462 14.7733114156647 14.7918675461740 14.8100607938761 14.8278884084541 14.8453477271214 14.8624361728907 14.8791512528942 14.8954905566905 14.9114517546233 14.9270325962123 14.9422309085518 14.9570445947605 14.9714716324547 14.9855100722574 14.9991580363324 15.0124137169719 15.0252753751953 15.0377413394025 15.0498100040606 15.0614798284217 15.0727493352868 15.0836171097996 15.0940817982945 15.1041421071646 15.1137968017931 15.1230447054977 15.1318846985485 15.1403157171975 15.1483367527728 15.1559468507978 15.1631451101687 15.1699306823621 15.1763027706924 15.1822606296067 15.1878035640377 15.1929309287799 15.1976421279231 15.2019366143269 15.2058138891405 15.2092735013615 15.2123150474481 15.2149381709684 15.2171425622938 15.2189279583436 15.2202941423691 15.2212409437850 15.2217682380443 15.2218759465534 15.2215640366438 15.2208325215816 15.2196814606151 15.2181109590788 15.2161211685380 15.2137122869763 15.2108845590303 15.2076382762649 15.2039737774986 15.1998914491728 15.1953917257595 15.1904750902205 15.1851420745074 15.1793932601052 15.1732292786184 15.1666508124125 15.1596585952808 15.1522534131671 15.1444361049279 15.1362075631419 15.1275687349524 15.1185206229646 15.1090642861711 15.0992008409393 15.0889314620121 15.0782573835810 15.0671799003725 15.0557003687977 15.0438202081214 15.0315409016862 15.0188639981653 15.0057911128618 14.9923239290392 14.9784641992845 14.9642137469258 14.9495744674618 14.9345483300433 14.9191373789804 14.9033437352853 14.8871695982381 14.8706172470037 14.8536890422565 14.8363874278365 14.8187149324577 14.8006741714023 14.7822678482666 14.7634987567260 14.7443697823167 14.7248839042356 14.7050441971615 14.6848538330923 14.6643160831972 14.6434343196810 14.6222120176461 14.6006527569879 14.5787602242768 14.5565382146413 14.5339906336631 14.5111214992698 14.4879349436086 14.4644352149330 14.4406266794612 14.4165138232310 14.3921012539356 14.3673937027406 14.3423960260710 14.3171132073858 14.2915503589030 14.2657127233100 14.2396056754165 14.2132347237864 14.1866055123027 14.1597238216928 14.1325955710035 14.1052268190089 14.0776237655520 14.0497927528291 14.0217402665901 13.9934729372732 13.9649975410367 13.9363210007282 13.9074503867399 13.8783929177765 13.8491559615194 13.8197470351722 13.7901738059093 13.7604440911907 13.7305658589542 13.7005472276829 13.6703964663321 13.6401219941195 13.6097323801597 13.5792363429507 13.5486427497123 13.5179606155416 13.4871991024201 13.4563675180325 13.4254753144087 13.3945320863905 13.3635475699038 13.3325316400254 13.3014943088784 13.2704457232944 13.2393961622866 13.2083560343119 13.1773358743055 13.1463463405071 13.1153982110603 13.0845023803886 13.0536698553459 13.0229117511231 12.9922392869446 12.9616637815093 12.9311966482045 12.9008493900917 12.8706335946393 12.8405609282271 12.8106431304167 12.7808920079823 12.7513194287019 12.7219373149333 12.6927576369378 12.6637924059955 12.6350536672900 12.6065534925670 12.5783039725978 12.5503172094101 12.5226053083371 12.4951803698506 12.4680544812221 12.4412397079893 12.4147480852570 12.3885916088252 12.3627822261764 12.3373318272990 12.3122522353922 12.2875551974418 12.2632523746837 12.2393553329748 12.2158755330673 12.1928243208279 12.1702129173786 12.1480524092137 12.1263537382695 12.1051276919948 12.0843848934215 12.0641357912385 12.0443906499223 12.0251595399010 12.0064523277994 11.9882786667536 11.9706479868411 11.9535694856290 11.9370521188450 11.9211045912216 11.9057353474929 11.8909525635977 11.8767641380711 11.8631776836700 11.8502005192287 11.8378396617775 11.8261018189287 11.8149933815437 11.8045204167117 11.7946886610360 11.7855035142551 11.7769700332095 11.7690929261532 11.7618765474546 11.7553248926539 11.7494415939299 11.7442299159583 11.7396927521825 11.7358326214991 11.7326516653745 11.7301516453872 11.7283339412107 11.7271995490355 11.7267490804422];
E2rms6 = [14.9972806805146 14.9969727671015 14.9962714089320 14.9951766960253 14.9936887555316 14.9918077518361 14.9895338867175 14.9868673995527 14.9838085675735 14.9803577061645 14.9765151692155 14.9722813495190 14.9676566792177 14.9626416302884 14.9572367150985 14.9514424869774 14.9452595408633 14.9386885139816 14.9317300865723 14.9243849826708 14.9166539709237 14.9085378654627 14.9000375268198 14.8911538628817 14.8818878299005 14.8722404335391 14.8622127299744 14.8518058270255 14.8410208853472 14.8298591196434 14.8183217999493 14.8064102529356 14.7941258632560 14.7814700749503 14.7684443928700 14.7550503841596 14.7412896797539 14.7271639759422 14.7126750359506 14.6978246915573 14.6826148447616 14.6670474694635 14.6511246132005 14.6348483988925 14.6182210266355 14.6012447755129 14.5839220054414 14.5662551590406 14.5482467635193 14.5298994326054 14.5112158684719 14.4921988637009 14.4728513032546 14.4531761664704 14.4331765290549 14.4128555651109 14.3922165491472 14.3712628581240 14.3499979734701 14.3284254831292 14.3065490835899 14.2843725819112 14.2618998977529 14.2391350653773 14.2160822356526 14.1927456780343 14.1691297825237 14.1452390616034 14.1210781521544 14.0966518173272 14.0719649483848 14.0470225665159 14.0218298245820 13.9963920088438 13.9707145406158 13.9448029778735 13.9186630167975 13.8923004932490 13.8657213841869 13.8389318089947 13.8119380307279 13.7847464572986 13.7573636425327 13.7297962871646 13.7020512397158 13.6741354972630 13.6460562061022 13.6178206623001 13.5894363121042 13.5609107522479 13.5322517301014 13.5034671436900 13.4745650415757 13.4455536225766 13.4164412353324 13.3872363777097 13.3579476960385 13.3285839841834 13.2991541824180 13.2696673761374 13.2401327943685 13.2105598080890 13.1809579283565 13.1513368042153 13.1217062204205 13.0920760949294 13.0624564761910 13.0328575402012 13.0032895873513 12.9737630390408 12.9442884340576 12.9148764247291 12.8855377728471 12.8562833453387 12.8271241097186 12.7980711292809 12.7691355580743 12.7403286356157 12.7116616813761 12.6831460890241 12.6547933204312 12.6266148994341 12.5986224053775 12.5708274664043 12.5432417525373 12.5158769685246 12.4887448464615 12.4618571382081 12.4352256075973 12.4088620224197 12.3827781462374 12.3569857299936 12.3314965034453 12.3063221664228 12.2814743799240 12.2569647570620 12.2328048538653 12.2090061599452 12.1855800890386 12.1625379694661 12.1398910344629 12.1176504124480 12.0958271172279 12.0744320381310 12.0534759301170 12.0329694038584 12.0129229158155 11.9933467583118 11.9742510496500 11.9556457242488 11.9375405228522 11.9199449828104 11.9028684284481 11.8863199615483 11.8703084519553 11.8548425283335 11.8399305690736 11.8255806933943 11.8118007526190 11.7985983216917 11.7859806908995 11.7739548578584 11.7625275197511 11.7517050658542 11.7414935703510 11.7318987854579 11.7229261348717 11.7145807075649 11.7068672519170 11.6997901702185 11.6933535135508 11.6875609770517 11.6824158955742 11.6779212397596 11.6740796125205 11.6708932459461 11.6683639986427 11.6664933535063 11.6652824159394 11.6647319125098 11.6648421900558 11.6656132152495 11.6670445746067 11.6691354749422 11.6718847442843 11.6752908332341 11.6793518167688 11.6840653964875 11.6894289032903 11.6954393004916 11.7020931873562 11.7093868030455 11.7173160309772 11.7258764035790 11.7350631074273 11.7448709887625 11.7552945593759 11.7663280028336 11.7779651810476 11.7901996411679 11.8030246227887 11.8164330654391 11.8304176163640 11.8449706385568 11.8600842190546 11.8757501774410 11.8919600745872 11.9087052215679 11.9259766887738 11.9437653151698 11.9620617177099 11.9808563008710 12.0001392663018 12.0199006225593 12.0401301949165 12.0608176352430 12.0819524319052 12.1035239197079 12.1255212898365 12.1479335997906 12.1707497832917 12.1939586601650 12.2175489461508 12.2415092626516 12.2658281464166 12.2904940590988 12.3154953967275 12.3408204990555 12.3664576587739 12.3923951305789 12.4186211400919 12.4451238926136 12.4718915817056 12.4989123975945 12.5261745353707 12.5536662030173 12.5813756292144 12.6092910709351 12.6374008208344 12.6656932144172 12.6941566369677 12.7227795302707 12.7515503990807 12.7804578173697 12.8094904343330 12.8386369801560 12.8678862715387 12.8972272169909 12.9266488218677 12.9561401931862 12.9856905441818 13.0152891986464 13.0449255950114 13.0745892902043 13.1042699632766 13.1339574187976 13.1636415900132 13.1933125417968 13.2229604733669 13.2525757208003 13.2821487593145 13.3116702053656 13.3411308185241 13.3705215031607 13.3998333099390 13.4290574371083 13.4581852316313 13.4872081901165 13.5161179595841 13.5449063380633 13.5735652750282 13.6020868716806 13.6304633810687 13.6586872080703 13.6867509092422 13.7146471925150 13.7423689167823 13.7699090913494 13.7972608752677 13.8244175765646 13.8513726513640 13.8781197028889 13.9046524804013 13.9309648780153 13.9570509334458 13.9829048266778 14.0085208785458 14.0338935492532 14.0590174368235 14.0838872754880 14.1084979340254 14.1328444140293 14.1569218481550 14.1807254982996 14.2042507537514 14.2274931293121 14.2504482633708 14.2731119159614 14.2954799667947 14.3175484132686 14.3393133684552 14.3607710590891 14.3819178235280 14.4027501097154 14.4232644731414 14.4434575747843 14.4633261790827 14.4828671518776 14.5020774583914 14.5209541611880 14.5394944181665 14.5576954805525 14.5755546909140 14.5930694811841 14.6102373707205 14.6270559643616 14.6435229505256 14.6596360993258 14.6753932607098 14.6907923626326 14.7058314092472 14.7205084791435 14.7348217235955 14.7487693648649 14.7623496945155 14.7755610717804 14.7884019219628 14.8008707348538 14.8129660632179 14.8246865212930 14.8360307833466 14.8469975822547 14.8575857081350 14.8677940070240 14.8776213795758 14.8870667798289 14.8961292139938 14.9048077393045 14.9131014628961 14.9210095407384 14.9285311766059 14.9356656211024 14.9424121707258 14.9487701669730 14.9547389955033 14.9603180853374 14.9655069081083 14.9703049773615 14.9747118478878 14.9787271151299 14.9823504146043 14.9855814213964 14.9884198496910 14.9908654523534 14.9929180205546 14.9945773834525 14.9958434079127 14.9967159982821 14.9971950962085 14.9972806805146];

h = figure();
plot(thetad,E2rms6,'LineWidth',3.5);
ylabel('Campo elétrico (kV/cm)')
xlabel('Ângulo (°)')
grid
hold on
plot(thetad,Erms,':','LineWidth',3.5);
hold on
plot(thetad,E_crit,'--','LineWidth',3.5);
xlim([0 360])
% ylim([13.5 15.5])
legend('Original', 'Raio menor', 'Campo elétrico crítico')
% title('Campo elétrico superficial condutor um')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'solo_caso1_dist','-dpdf','-r0')