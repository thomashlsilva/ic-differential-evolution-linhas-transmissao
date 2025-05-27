% Raio maior
clear all; close all; clc
% Cálculo do campo elétrico superficial do caso 3 do sistema trifásico

e_0 = 8.854*(10^(-12));
r = 0.01815; % raio do condutor
n = 9; % número de condutores do sistema
nc = n*2; % número de condutores total (com imagens)
ci = (2*n-1); % número de cargas imagens por condutor

cond = 9;
%% matriz P

% posição dos condutores reais do sistema
xr = [-10.4785 -10.25 -10.0215 -0.2285 0 0.2285 10.0215 10.25 10.4785];
yr = [16.532 16.728 16.532 16.532 16.728 16.532 16.532 16.728 16.532];
yi = -yr;

% cálculo matriz P
P = zeros(n,n);
for i = 1:n
     for j = 1:n
        if(i==j)
           P(i,j) = (1/(2*pi*e_0))*log((4*yr(i))/(2*r));
        else
           P(i,j) = (1/(2*pi*e_0))*log(sqrt((xr(j)-xr(i))^2+(yi(j)-yr(i))^2)/(sqrt((xr(j)-xr(i))^2+(yr(j)-yr(i))^2)));
        end
     end
end

%% Cálculo tensão por fase

%tensão de operação máxima do sistema
V = 525*10^3;

%tensão condutores 1 e 2 fase a
V_ra = V/sqrt(3);
V_ia = 0;

%tensão condutores 3 e 4 fase b
V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

%tensão condutores 5 e 6 fase c
V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ];

%% Cálculo densidade de carga

rho = P\Vf;

%busco cada uma das posicoes de rho real e imaginario (ver indice)

%cabo 1 fase a
rho_r1 = real(rho(1)); 
rho_i1 = imag(rho(1));

%cabo 2 fase a
rho_r2 = real(rho(2));
rho_i2 = imag(rho(2));

%cabo 3 fase a
rho_r3 = real(rho(3));
rho_i3 = imag(rho(3));

%cabo 4 fase b
rho_r4 = real(rho(4));
rho_i4 = imag(rho(4));

%cabo 5 fase b
rho_r5 = real(rho(5));
rho_i5 = imag(rho(5));

%cabo 6 fase b
rho_r6 = real(rho(6));
rho_i6 = imag(rho(6));

%cabo 7 fase c
rho_r7 = real(rho(7));
rho_i7 = imag(rho(7));

%cabo 8 fase c
rho_r8 = real(rho(8));
rho_i8 = imag(rho(8));

%cabo 9 fase c
rho_r9 = real(rho(9));
rho_i9 = imag(rho(9));


%% Distância entre condutores

x = [xr xr]; % posição x dos condutores (com imagens)
y = [yr -yr]; % posição y dos condutores (com imagens)

%sendo: i -> o número do primeiro condutor e j -> o número do segundo condutor
%exemplo: 1,2 é a distância entre condutor 1 e o condutor 2, os próximos
%cálculos seguem esta mesma lógica
%sendo 10 a imagem do condutor 1, 12 a imagem do condutor 3 e 11 a imagem do
%condutor 2

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

% plot(x,y,'o')
% grid

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
            if x(j) > x(i)
                posx(i,j) = x(i) + delta(i,j);
            elseif x(i) > x(j)
                posx(i,j) = x(i) - delta(i,j);
            end
        elseif y(j) > y(i)
            if x(j) > x(i)
                phix(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posx(i,j) = x(i) + delta(i,j)*cos(phix(i,j));
            elseif x(i) > x(j)
                phix(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posx(i,j) = x(i) - delta(i,j)*cos(phix(i,j));
            end
        elseif y(i) > y(j)
            if x(j) > x(i)
                phix(i,j) = acos((y(i)-y(j))/(distance(i,j)));
                posx(i,j) = x(i) + delta(i,j)*sin(phix(i,j));
            elseif x(i) > x(j)
                phix(i,j) = acos((y(i)-y(j))/(distance(i,j)));
                posx(i,j) = x(i) - delta(i,j)*sin(phix(i,j));
            end
        end
    end
end

%% Cálculo das posições em y das cargas imagens

phiy = zeros(nc,nc);
posy = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j || y(j)==y(i)
            posy(i,j) = y(i);
        elseif x(i)==x(j)
            posy(i,j) = y(i) - delta(i,j);
        elseif y(j) > y(i)
            if  y(i) > 0 || y(i) < 0 && y(j) < 0
                phiy(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posy(i,j) = y(i) + delta(i,j)*sin(phiy(i,j));
            elseif y(i) < 0 && y(j) > 0
                phiy(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posy(i,j) = y(i) - delta(i,j)*sin(phiy(i,j));
            end
        elseif y(i) > y(j)
            phiy(i,j) = acos((y(i)-y(j))/(distance(i,j)));
            posy(i,j) = y(i) - delta(i,j)*cos(phiy(i,j));
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
ycj = y(cond); % posição y do centro condutor 1 fase a (para trocar o condutor avaliado deve-se mudar este valor)

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
E_xr17 = (-rho_r7/(2*pi*e_0)).*((xf - posx(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_xr18 = (-rho_r8/(2*pi*e_0)).*((xf - posx(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_xr19 = (-rho_r9/(2*pi*e_0)).*((xf - posx(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_xr110 = (rho_r1/(2*pi*e_0)).*((xf - posx(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_xr1_11 = (rho_r2/(2*pi*e_0)).*((xf - posx(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_xr1_12 = (rho_r3/(2*pi*e_0)).*((xf - posx(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_xr1_13 = (rho_r4/(2*pi*e_0)).*((xf - posx(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_xr1_14 = (rho_r5/(2*pi*e_0)).*((xf - posx(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_xr1_15 = (rho_r6/(2*pi*e_0)).*((xf - posx(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_xr1_16 = (rho_r7/(2*pi*e_0)).*((xf - posx(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_xr1_17 = (rho_r8/(2*pi*e_0)).*((xf - posx(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_xr1_18 = (rho_r9/(2*pi*e_0)).*((xf - posx(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));

E_xr21 = (-rho_r1/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xr23 = (-rho_r3/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xr24 = (-rho_r4/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xr25 = (-rho_r5/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xr26 = (-rho_r6/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_xr27 = (-rho_r7/(2*pi*e_0)).*((xf - posx(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_xr28 = (-rho_r8/(2*pi*e_0)).*((xf - posx(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_xr29 = (-rho_r9/(2*pi*e_0)).*((xf - posx(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_xr210 = (rho_r1/(2*pi*e_0)).*((xf - posx(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_xr2_11 = (rho_r2/(2*pi*e_0)).*((xf - posx(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_xr2_12 = (rho_r3/(2*pi*e_0)).*((xf - posx(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_xr2_13 = (rho_r4/(2*pi*e_0)).*((xf - posx(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_xr2_14 = (rho_r5/(2*pi*e_0)).*((xf - posx(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_xr2_15 = (rho_r6/(2*pi*e_0)).*((xf - posx(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_xr2_16 = (rho_r7/(2*pi*e_0)).*((xf - posx(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_xr2_17 = (rho_r8/(2*pi*e_0)).*((xf - posx(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_xr2_18 = (rho_r9/(2*pi*e_0)).*((xf - posx(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));

E_xr31 = (-rho_r1/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xr32 = (-rho_r2/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xr34 = (-rho_r4/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xr35 = (-rho_r5/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xr36 = (-rho_r6/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_xr37 = (-rho_r7/(2*pi*e_0)).*((xf - posx(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_xr38 = (-rho_r8/(2*pi*e_0)).*((xf - posx(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_xr39 = (-rho_r9/(2*pi*e_0)).*((xf - posx(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_xr310 = (rho_r1/(2*pi*e_0)).*((xf - posx(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_xr311 = (rho_r2/(2*pi*e_0)).*((xf - posx(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_xr312 = (rho_r3/(2*pi*e_0)).*((xf - posx(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_xr313 = (rho_r4/(2*pi*e_0)).*((xf - posx(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_xr314 = (rho_r5/(2*pi*e_0)).*((xf - posx(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_xr315 = (rho_r6/(2*pi*e_0)).*((xf - posx(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_xr316 = (rho_r7/(2*pi*e_0)).*((xf - posx(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_xr317 = (rho_r8/(2*pi*e_0)).*((xf - posx(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_xr318 = (rho_r9/(2*pi*e_0)).*((xf - posx(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));

E_xr41 = (-rho_r1/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xr42 = (-rho_r2/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xr43 = (-rho_r3/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xr45 = (-rho_r5/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xr46 = (-rho_r6/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_xr47 = (-rho_r7/(2*pi*e_0)).*((xf - posx(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_xr48 = (-rho_r8/(2*pi*e_0)).*((xf - posx(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_xr49 = (-rho_r9/(2*pi*e_0)).*((xf - posx(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_xr410 = (rho_r1/(2*pi*e_0)).*((xf - posx(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_xr411 = (rho_r2/(2*pi*e_0)).*((xf - posx(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_xr412 = (rho_r3/(2*pi*e_0)).*((xf - posx(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_xr413 = (rho_r4/(2*pi*e_0)).*((xf - posx(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_xr414 = (rho_r5/(2*pi*e_0)).*((xf - posx(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_xr415 = (rho_r6/(2*pi*e_0)).*((xf - posx(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_xr416 = (rho_r7/(2*pi*e_0)).*((xf - posx(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_xr417 = (rho_r8/(2*pi*e_0)).*((xf - posx(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_xr418 = (rho_r9/(2*pi*e_0)).*((xf - posx(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));

E_xr51 = (-rho_r1/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xr52 = (-rho_r2/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xr53 = (-rho_r3/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xr54 = (-rho_r4/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xr56 = (-rho_r6/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_xr57 = (-rho_r7/(2*pi*e_0)).*((xf - posx(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_xr58 = (-rho_r8/(2*pi*e_0)).*((xf - posx(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_xr59 = (-rho_r9/(2*pi*e_0)).*((xf - posx(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_xr510 = (rho_r1/(2*pi*e_0)).*((xf - posx(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_xr511 = (rho_r2/(2*pi*e_0)).*((xf - posx(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_xr512 = (rho_r3/(2*pi*e_0)).*((xf - posx(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_xr513 = (rho_r4/(2*pi*e_0)).*((xf - posx(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_xr514 = (rho_r5/(2*pi*e_0)).*((xf - posx(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_xr515 = (rho_r6/(2*pi*e_0)).*((xf - posx(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_xr516 = (rho_r7/(2*pi*e_0)).*((xf - posx(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_xr517 = (rho_r8/(2*pi*e_0)).*((xf - posx(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_xr518 = (rho_r9/(2*pi*e_0)).*((xf - posx(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));

E_xr61 = (-rho_r1/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xr62 = (-rho_r2/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xr63 = (-rho_r3/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xr64 = (-rho_r4/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xr65 = (-rho_r5/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_xr67 = (-rho_r7/(2*pi*e_0)).*((xf - posx(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_xr68 = (-rho_r8/(2*pi*e_0)).*((xf - posx(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_xr69 = (-rho_r9/(2*pi*e_0)).*((xf - posx(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_xr610 = (rho_r1/(2*pi*e_0)).*((xf - posx(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_xr611 = (rho_r2/(2*pi*e_0)).*((xf - posx(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_xr612 = (rho_r3/(2*pi*e_0)).*((xf - posx(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_xr613 = (rho_r4/(2*pi*e_0)).*((xf - posx(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_xr614 = (rho_r5/(2*pi*e_0)).*((xf - posx(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_xr615 = (rho_r6/(2*pi*e_0)).*((xf - posx(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_xr616 = (rho_r7/(2*pi*e_0)).*((xf - posx(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_xr617 = (rho_r8/(2*pi*e_0)).*((xf - posx(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_xr618 = (rho_r9/(2*pi*e_0)).*((xf - posx(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));

E_xr71 = (-rho_r1/(2*pi*e_0)).*((xf - posx(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_xr72 = (-rho_r2/(2*pi*e_0)).*((xf - posx(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_xr73 = (-rho_r3/(2*pi*e_0)).*((xf - posx(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_xr74 = (-rho_r4/(2*pi*e_0)).*((xf - posx(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_xr75 = (-rho_r5/(2*pi*e_0)).*((xf - posx(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_xr76 = (-rho_r6/(2*pi*e_0)).*((xf - posx(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_xr78 = (-rho_r8/(2*pi*e_0)).*((xf - posx(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_xr79 = (-rho_r9/(2*pi*e_0)).*((xf - posx(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_xr710 = (rho_r1/(2*pi*e_0)).*((xf - posx(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_xr711 = (rho_r2/(2*pi*e_0)).*((xf - posx(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_xr712 = (rho_r3/(2*pi*e_0)).*((xf - posx(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_xr713 = (rho_r4/(2*pi*e_0)).*((xf - posx(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_xr714 = (rho_r5/(2*pi*e_0)).*((xf - posx(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_xr715 = (rho_r6/(2*pi*e_0)).*((xf - posx(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_xr716 = (rho_r7/(2*pi*e_0)).*((xf - posx(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_xr717 = (rho_r8/(2*pi*e_0)).*((xf - posx(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_xr718 = (rho_r9/(2*pi*e_0)).*((xf - posx(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));

E_xr81 = (-rho_r1/(2*pi*e_0)).*((xf - posx(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_xr82 = (-rho_r2/(2*pi*e_0)).*((xf - posx(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_xr83 = (-rho_r3/(2*pi*e_0)).*((xf - posx(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_xr84 = (-rho_r4/(2*pi*e_0)).*((xf - posx(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_xr85 = (-rho_r5/(2*pi*e_0)).*((xf - posx(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_xr86 = (-rho_r6/(2*pi*e_0)).*((xf - posx(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_xr87 = (-rho_r7/(2*pi*e_0)).*((xf - posx(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_xr89 = (-rho_r9/(2*pi*e_0)).*((xf - posx(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_xr810 = (rho_r1/(2*pi*e_0)).*((xf - posx(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_xr811 = (rho_r2/(2*pi*e_0)).*((xf - posx(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_xr812 = (rho_r3/(2*pi*e_0)).*((xf - posx(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_xr813 = (rho_r4/(2*pi*e_0)).*((xf - posx(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_xr814 = (rho_r5/(2*pi*e_0)).*((xf - posx(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_xr815 = (rho_r6/(2*pi*e_0)).*((xf - posx(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_xr816 = (rho_r7/(2*pi*e_0)).*((xf - posx(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_xr817 = (rho_r8/(2*pi*e_0)).*((xf - posx(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_xr818 = (rho_r9/(2*pi*e_0)).*((xf - posx(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));

E_xr91 = (-rho_r1/(2*pi*e_0)).*((xf - posx(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_xr92 = (-rho_r2/(2*pi*e_0)).*((xf - posx(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_xr93 = (-rho_r3/(2*pi*e_0)).*((xf - posx(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_xr94 = (-rho_r4/(2*pi*e_0)).*((xf - posx(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_xr95 = (-rho_r5/(2*pi*e_0)).*((xf - posx(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_xr96 = (-rho_r6/(2*pi*e_0)).*((xf - posx(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_xr97 = (-rho_r7/(2*pi*e_0)).*((xf - posx(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_xr98 = (-rho_r8/(2*pi*e_0)).*((xf - posx(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_xr910 = (rho_r1/(2*pi*e_0)).*((xf - posx(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_xr911 = (rho_r2/(2*pi*e_0)).*((xf - posx(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_xr912 = (rho_r3/(2*pi*e_0)).*((xf - posx(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_xr913 = (rho_r4/(2*pi*e_0)).*((xf - posx(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_xr914 = (rho_r5/(2*pi*e_0)).*((xf - posx(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_xr915 = (rho_r6/(2*pi*e_0)).*((xf - posx(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_xr916 = (rho_r7/(2*pi*e_0)).*((xf - posx(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_xr917 = (rho_r8/(2*pi*e_0)).*((xf - posx(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_xr918 = (rho_r9/(2*pi*e_0)).*((xf - posx(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));

E_xr101 = (-rho_r1/(2*pi*e_0)).*((xf - posx(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_xr102 = (-rho_r2/(2*pi*e_0)).*((xf - posx(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_xr103 = (-rho_r3/(2*pi*e_0)).*((xf - posx(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_xr104 = (-rho_r4/(2*pi*e_0)).*((xf - posx(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_xr105 = (-rho_r5/(2*pi*e_0)).*((xf - posx(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_xr106 = (-rho_r6/(2*pi*e_0)).*((xf - posx(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_xr107 = (-rho_r7/(2*pi*e_0)).*((xf - posx(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_xr108 = (-rho_r8/(2*pi*e_0)).*((xf - posx(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_xr109 = (-rho_r9/(2*pi*e_0)).*((xf - posx(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_xr1011 = (rho_r2/(2*pi*e_0)).*((xf - posx(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_xr1012 = (rho_r3/(2*pi*e_0)).*((xf - posx(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_xr1013 = (rho_r4/(2*pi*e_0)).*((xf - posx(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_xr1014 = (rho_r5/(2*pi*e_0)).*((xf - posx(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_xr1015 = (rho_r6/(2*pi*e_0)).*((xf - posx(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_xr1016 = (rho_r7/(2*pi*e_0)).*((xf - posx(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_xr1017 = (rho_r8/(2*pi*e_0)).*((xf - posx(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_xr1018 = (rho_r9/(2*pi*e_0)).*((xf - posx(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));

E_xr11_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_xr11_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_xr113 = (-rho_r3/(2*pi*e_0)).*((xf - posx(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_xr114 = (-rho_r4/(2*pi*e_0)).*((xf - posx(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_xr115 = (-rho_r5/(2*pi*e_0)).*((xf - posx(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_xr116 = (-rho_r6/(2*pi*e_0)).*((xf - posx(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_xr117 = (-rho_r7/(2*pi*e_0)).*((xf - posx(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_xr118 = (-rho_r8/(2*pi*e_0)).*((xf - posx(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_xr119 = (-rho_r9/(2*pi*e_0)).*((xf - posx(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_xr1110 = (rho_r1/(2*pi*e_0)).*((xf - posx(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_xr1112 = (rho_r3/(2*pi*e_0)).*((xf - posx(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_xr1113 = (rho_r4/(2*pi*e_0)).*((xf - posx(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_xr1114 = (rho_r5/(2*pi*e_0)).*((xf - posx(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_xr1115 = (rho_r6/(2*pi*e_0)).*((xf - posx(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_xr1116 = (rho_r7/(2*pi*e_0)).*((xf - posx(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_xr1117 = (rho_r8/(2*pi*e_0)).*((xf - posx(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_xr1118 = (rho_r9/(2*pi*e_0)).*((xf - posx(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));

E_xr12_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_xr12_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_xr123 = (-rho_r3/(2*pi*e_0)).*((xf - posx(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_xr124 = (-rho_r4/(2*pi*e_0)).*((xf - posx(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_xr125 = (-rho_r5/(2*pi*e_0)).*((xf - posx(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_xr126 = (-rho_r6/(2*pi*e_0)).*((xf - posx(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_xr127 = (-rho_r7/(2*pi*e_0)).*((xf - posx(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_xr128 = (-rho_r8/(2*pi*e_0)).*((xf - posx(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_xr129 = (-rho_r9/(2*pi*e_0)).*((xf - posx(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_xr1210 = (rho_r1/(2*pi*e_0)).*((xf - posx(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_xr1211 = (rho_r2/(2*pi*e_0)).*((xf - posx(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_xr1213 = (rho_r4/(2*pi*e_0)).*((xf - posx(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_xr1214 = (rho_r5/(2*pi*e_0)).*((xf - posx(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_xr1215 = (rho_r6/(2*pi*e_0)).*((xf - posx(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_xr1216 = (rho_r7/(2*pi*e_0)).*((xf - posx(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_xr1217 = (rho_r8/(2*pi*e_0)).*((xf - posx(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_xr1218 = (rho_r9/(2*pi*e_0)).*((xf - posx(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));

E_xr131 = (-rho_r1/(2*pi*e_0)).*((xf - posx(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_xr132 = (-rho_r2/(2*pi*e_0)).*((xf - posx(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_xr133 = (-rho_r3/(2*pi*e_0)).*((xf - posx(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_xr134 = (-rho_r4/(2*pi*e_0)).*((xf - posx(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_xr135 = (-rho_r5/(2*pi*e_0)).*((xf - posx(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_xr136 = (-rho_r6/(2*pi*e_0)).*((xf - posx(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_xr137 = (-rho_r7/(2*pi*e_0)).*((xf - posx(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_xr138 = (-rho_r8/(2*pi*e_0)).*((xf - posx(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_xr139 = (-rho_r9/(2*pi*e_0)).*((xf - posx(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_xr1310 = (rho_r1/(2*pi*e_0)).*((xf - posx(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_xr1311 = (rho_r2/(2*pi*e_0)).*((xf - posx(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_xr1312 = (rho_r3/(2*pi*e_0)).*((xf - posx(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_xr1314 = (rho_r5/(2*pi*e_0)).*((xf - posx(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_xr1315 = (rho_r6/(2*pi*e_0)).*((xf - posx(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_xr1316 = (rho_r7/(2*pi*e_0)).*((xf - posx(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_xr1317 = (rho_r8/(2*pi*e_0)).*((xf - posx(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_xr1318 = (rho_r9/(2*pi*e_0)).*((xf - posx(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));

E_xr141 = (-rho_r1/(2*pi*e_0)).*((xf - posx(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_xr142 = (-rho_r2/(2*pi*e_0)).*((xf - posx(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_xr143 = (-rho_r3/(2*pi*e_0)).*((xf - posx(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_xr144 = (-rho_r4/(2*pi*e_0)).*((xf - posx(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_xr145 = (-rho_r5/(2*pi*e_0)).*((xf - posx(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_xr146 = (-rho_r6/(2*pi*e_0)).*((xf - posx(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_xr147 = (-rho_r7/(2*pi*e_0)).*((xf - posx(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_xr148 = (-rho_r8/(2*pi*e_0)).*((xf - posx(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_xr149 = (-rho_r9/(2*pi*e_0)).*((xf - posx(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_xr1410 = (rho_r1/(2*pi*e_0)).*((xf - posx(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_xr1411 = (rho_r2/(2*pi*e_0)).*((xf - posx(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_xr1412 = (rho_r3/(2*pi*e_0)).*((xf - posx(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_xr1413 = (rho_r4/(2*pi*e_0)).*((xf - posx(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_xr1415 = (rho_r6/(2*pi*e_0)).*((xf - posx(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_xr1416 = (rho_r7/(2*pi*e_0)).*((xf - posx(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_xr1417 = (rho_r8/(2*pi*e_0)).*((xf - posx(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_xr1418 = (rho_r9/(2*pi*e_0)).*((xf - posx(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));

E_xr151 = (-rho_r1/(2*pi*e_0)).*((xf - posx(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_xr152 = (-rho_r2/(2*pi*e_0)).*((xf - posx(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_xr153 = (-rho_r3/(2*pi*e_0)).*((xf - posx(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_xr154 = (-rho_r4/(2*pi*e_0)).*((xf - posx(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_xr155 = (-rho_r5/(2*pi*e_0)).*((xf - posx(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_xr156 = (-rho_r6/(2*pi*e_0)).*((xf - posx(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_xr157 = (-rho_r7/(2*pi*e_0)).*((xf - posx(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_xr158 = (-rho_r8/(2*pi*e_0)).*((xf - posx(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_xr159 = (-rho_r9/(2*pi*e_0)).*((xf - posx(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_xr1510 = (rho_r1/(2*pi*e_0)).*((xf - posx(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_xr1511 = (rho_r2/(2*pi*e_0)).*((xf - posx(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_xr1512 = (rho_r3/(2*pi*e_0)).*((xf - posx(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_xr1513 = (rho_r4/(2*pi*e_0)).*((xf - posx(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_xr1514 = (rho_r5/(2*pi*e_0)).*((xf - posx(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_xr1516 = (rho_r7/(2*pi*e_0)).*((xf - posx(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_xr1517 = (rho_r8/(2*pi*e_0)).*((xf - posx(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_xr1518 = (rho_r9/(2*pi*e_0)).*((xf - posx(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));

E_xr161 = (-rho_r1/(2*pi*e_0)).*((xf - posx(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_xr162 = (-rho_r2/(2*pi*e_0)).*((xf - posx(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_xr163 = (-rho_r3/(2*pi*e_0)).*((xf - posx(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_xr164 = (-rho_r4/(2*pi*e_0)).*((xf - posx(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_xr165 = (-rho_r5/(2*pi*e_0)).*((xf - posx(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_xr166 = (-rho_r6/(2*pi*e_0)).*((xf - posx(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_xr167 = (-rho_r7/(2*pi*e_0)).*((xf - posx(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_xr168 = (-rho_r8/(2*pi*e_0)).*((xf - posx(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_xr169 = (-rho_r9/(2*pi*e_0)).*((xf - posx(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_xr1610 = (rho_r1/(2*pi*e_0)).*((xf - posx(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_xr1611 = (rho_r2/(2*pi*e_0)).*((xf - posx(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_xr1612 = (rho_r3/(2*pi*e_0)).*((xf - posx(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_xr1613 = (rho_r4/(2*pi*e_0)).*((xf - posx(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_xr1614 = (rho_r5/(2*pi*e_0)).*((xf - posx(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_xr1615 = (rho_r6/(2*pi*e_0)).*((xf - posx(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_xr1617 = (rho_r8/(2*pi*e_0)).*((xf - posx(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_xr1618 = (rho_r9/(2*pi*e_0)).*((xf - posx(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));

E_xr171 = (-rho_r1/(2*pi*e_0)).*((xf - posx(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_xr172 = (-rho_r2/(2*pi*e_0)).*((xf - posx(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_xr173 = (-rho_r3/(2*pi*e_0)).*((xf - posx(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_xr174 = (-rho_r4/(2*pi*e_0)).*((xf - posx(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_xr175 = (-rho_r5/(2*pi*e_0)).*((xf - posx(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_xr176 = (-rho_r6/(2*pi*e_0)).*((xf - posx(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_xr177 = (-rho_r7/(2*pi*e_0)).*((xf - posx(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_xr178 = (-rho_r8/(2*pi*e_0)).*((xf - posx(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_xr179 = (-rho_r9/(2*pi*e_0)).*((xf - posx(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_xr1710 = (rho_r1/(2*pi*e_0)).*((xf - posx(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_xr1711 = (rho_r2/(2*pi*e_0)).*((xf - posx(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_xr1712 = (rho_r3/(2*pi*e_0)).*((xf - posx(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_xr1713 = (rho_r4/(2*pi*e_0)).*((xf - posx(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_xr1714 = (rho_r5/(2*pi*e_0)).*((xf - posx(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_xr1715 = (rho_r6/(2*pi*e_0)).*((xf - posx(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_xr1716 = (rho_r7/(2*pi*e_0)).*((xf - posx(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_xr1718 = (rho_r9/(2*pi*e_0)).*((xf - posx(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));

E_xr181 = (-rho_r1/(2*pi*e_0)).*((xf - posx(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_xr182 = (-rho_r2/(2*pi*e_0)).*((xf - posx(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_xr183 = (-rho_r3/(2*pi*e_0)).*((xf - posx(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_xr184 = (-rho_r4/(2*pi*e_0)).*((xf - posx(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_xr185 = (-rho_r5/(2*pi*e_0)).*((xf - posx(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_xr186 = (-rho_r6/(2*pi*e_0)).*((xf - posx(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_xr187 = (-rho_r7/(2*pi*e_0)).*((xf - posx(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_xr188 = (-rho_r8/(2*pi*e_0)).*((xf - posx(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_xr189 = (-rho_r9/(2*pi*e_0)).*((xf - posx(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_xr1810 = (rho_r1/(2*pi*e_0)).*((xf - posx(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_xr1811 = (rho_r2/(2*pi*e_0)).*((xf - posx(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_xr1812 = (rho_r3/(2*pi*e_0)).*((xf - posx(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_xr1813 = (rho_r4/(2*pi*e_0)).*((xf - posx(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_xr1814 = (rho_r5/(2*pi*e_0)).*((xf - posx(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_xr1815 = (rho_r6/(2*pi*e_0)).*((xf - posx(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_xr1816 = (rho_r7/(2*pi*e_0)).*((xf - posx(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_xr1817 = (rho_r8/(2*pi*e_0)).*((xf - posx(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));

%% E_xi componente x imaginario campo elétrico condutor 2 fase b assim segue:

E_xi12 = (-rho_i2/(2*pi*e_0)).*((xf - posx(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_xi13 = (-rho_i3/(2*pi*e_0)).*((xf - posx(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_xi14 = (-rho_i4/(2*pi*e_0)).*((xf - posx(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_xi15 = (-rho_i5/(2*pi*e_0)).*((xf - posx(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_xi16 = (-rho_i6/(2*pi*e_0)).*((xf - posx(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_xi17 = (-rho_i7/(2*pi*e_0)).*((xf - posx(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_xi18 = (-rho_i8/(2*pi*e_0)).*((xf - posx(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_xi19 = (-rho_i9/(2*pi*e_0)).*((xf - posx(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_xi110 = (rho_i1/(2*pi*e_0)).*((xf - posx(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_xi1_11 = (rho_i2/(2*pi*e_0)).*((xf - posx(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_xi1_12 = (rho_i3/(2*pi*e_0)).*((xf - posx(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_xi1_13 = (rho_i4/(2*pi*e_0)).*((xf - posx(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_xi1_14 = (rho_i5/(2*pi*e_0)).*((xf - posx(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_xi1_15 = (rho_i6/(2*pi*e_0)).*((xf - posx(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_xi1_16 = (rho_i7/(2*pi*e_0)).*((xf - posx(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_xi1_17 = (rho_i8/(2*pi*e_0)).*((xf - posx(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_xi1_18 = (rho_i9/(2*pi*e_0)).*((xf - posx(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));

E_xi21 = (-rho_i1/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xi23 = (-rho_i3/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xi24 = (-rho_i4/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xi25 = (-rho_i5/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xi26 = (-rho_i6/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_xi27 = (-rho_i7/(2*pi*e_0)).*((xf - posx(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_xi28 = (-rho_i8/(2*pi*e_0)).*((xf - posx(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_xi29 = (-rho_i9/(2*pi*e_0)).*((xf - posx(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_xi210 = (rho_i1/(2*pi*e_0)).*((xf - posx(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_xi2_11 = (rho_i2/(2*pi*e_0)).*((xf - posx(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_xi2_12 = (rho_i3/(2*pi*e_0)).*((xf - posx(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_xi2_13 = (rho_i4/(2*pi*e_0)).*((xf - posx(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_xi2_14 = (rho_i5/(2*pi*e_0)).*((xf - posx(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_xi2_15 = (rho_i6/(2*pi*e_0)).*((xf - posx(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_xi2_16 = (rho_i7/(2*pi*e_0)).*((xf - posx(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_xi2_17 = (rho_i8/(2*pi*e_0)).*((xf - posx(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_xi2_18 = (rho_i9/(2*pi*e_0)).*((xf - posx(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));

E_xi31 = (-rho_i1/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xi32 = (-rho_i2/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xi34 = (-rho_i4/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xi35 = (-rho_i5/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xi36 = (-rho_i6/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_xi37 = (-rho_i7/(2*pi*e_0)).*((xf - posx(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_xi38 = (-rho_i8/(2*pi*e_0)).*((xf - posx(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_xi39 = (-rho_i9/(2*pi*e_0)).*((xf - posx(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_xi310 = (rho_i1/(2*pi*e_0)).*((xf - posx(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_xi311 = (rho_i2/(2*pi*e_0)).*((xf - posx(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_xi312 = (rho_i3/(2*pi*e_0)).*((xf - posx(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_xi313 = (rho_i4/(2*pi*e_0)).*((xf - posx(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_xi314 = (rho_i5/(2*pi*e_0)).*((xf - posx(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_xi315 = (rho_i6/(2*pi*e_0)).*((xf - posx(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_xi316 = (rho_i7/(2*pi*e_0)).*((xf - posx(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_xi317 = (rho_i8/(2*pi*e_0)).*((xf - posx(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_xi318 = (rho_i9/(2*pi*e_0)).*((xf - posx(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));

E_xi41 = (-rho_i1/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xi42 = (-rho_i2/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xi43 = (-rho_i3/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xi45 = (-rho_i5/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xi46 = (-rho_i6/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_xi47 = (-rho_i7/(2*pi*e_0)).*((xf - posx(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_xi48 = (-rho_i8/(2*pi*e_0)).*((xf - posx(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_xi49 = (-rho_i9/(2*pi*e_0)).*((xf - posx(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_xi410 = (rho_i1/(2*pi*e_0)).*((xf - posx(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_xi411 = (rho_i2/(2*pi*e_0)).*((xf - posx(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_xi412 = (rho_i3/(2*pi*e_0)).*((xf - posx(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_xi413 = (rho_i4/(2*pi*e_0)).*((xf - posx(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_xi414 = (rho_i5/(2*pi*e_0)).*((xf - posx(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_xi415 = (rho_i6/(2*pi*e_0)).*((xf - posx(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_xi416 = (rho_i7/(2*pi*e_0)).*((xf - posx(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_xi417 = (rho_i8/(2*pi*e_0)).*((xf - posx(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_xi418 = (rho_i9/(2*pi*e_0)).*((xf - posx(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));

E_xi51 = (-rho_i1/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xi52 = (-rho_i2/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xi53 = (-rho_i3/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xi54 = (-rho_i4/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xi56 = (-rho_i6/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_xi57 = (-rho_i7/(2*pi*e_0)).*((xf - posx(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_xi58 = (-rho_i8/(2*pi*e_0)).*((xf - posx(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_xi59 = (-rho_i9/(2*pi*e_0)).*((xf - posx(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_xi510 = (rho_i1/(2*pi*e_0)).*((xf - posx(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_xi511 = (rho_i2/(2*pi*e_0)).*((xf - posx(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_xi512 = (rho_i3/(2*pi*e_0)).*((xf - posx(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_xi513 = (rho_i4/(2*pi*e_0)).*((xf - posx(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_xi514 = (rho_i5/(2*pi*e_0)).*((xf - posx(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_xi515 = (rho_i6/(2*pi*e_0)).*((xf - posx(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_xi516 = (rho_i7/(2*pi*e_0)).*((xf - posx(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_xi517 = (rho_i8/(2*pi*e_0)).*((xf - posx(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_xi518 = (rho_i9/(2*pi*e_0)).*((xf - posx(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));

E_xi61 = (-rho_i1/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xi62 = (-rho_i2/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xi63 = (-rho_i3/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xi64 = (-rho_i4/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xi65 = (-rho_i5/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_xi67 = (-rho_i7/(2*pi*e_0)).*((xf - posx(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_xi68 = (-rho_i8/(2*pi*e_0)).*((xf - posx(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_xi69 = (-rho_i9/(2*pi*e_0)).*((xf - posx(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_xi610 = (rho_i1/(2*pi*e_0)).*((xf - posx(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_xi611 = (rho_i2/(2*pi*e_0)).*((xf - posx(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_xi612 = (rho_i3/(2*pi*e_0)).*((xf - posx(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_xi613 = (rho_i4/(2*pi*e_0)).*((xf - posx(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_xi614 = (rho_i5/(2*pi*e_0)).*((xf - posx(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_xi615 = (rho_i6/(2*pi*e_0)).*((xf - posx(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_xi616 = (rho_i7/(2*pi*e_0)).*((xf - posx(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_xi617 = (rho_i8/(2*pi*e_0)).*((xf - posx(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_xi618 = (rho_i9/(2*pi*e_0)).*((xf - posx(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));

E_xi71 = (-rho_i1/(2*pi*e_0)).*((xf - posx(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_xi72 = (-rho_i2/(2*pi*e_0)).*((xf - posx(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_xi73 = (-rho_i3/(2*pi*e_0)).*((xf - posx(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_xi74 = (-rho_i4/(2*pi*e_0)).*((xf - posx(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_xi75 = (-rho_i5/(2*pi*e_0)).*((xf - posx(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_xi76 = (-rho_i6/(2*pi*e_0)).*((xf - posx(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_xi78 = (-rho_i8/(2*pi*e_0)).*((xf - posx(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_xi79 = (-rho_i9/(2*pi*e_0)).*((xf - posx(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_xi710 = (rho_i1/(2*pi*e_0)).*((xf - posx(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_xi711 = (rho_i2/(2*pi*e_0)).*((xf - posx(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_xi712 = (rho_i3/(2*pi*e_0)).*((xf - posx(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_xi713 = (rho_i4/(2*pi*e_0)).*((xf - posx(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_xi714 = (rho_i5/(2*pi*e_0)).*((xf - posx(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_xi715 = (rho_i6/(2*pi*e_0)).*((xf - posx(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_xi716 = (rho_i7/(2*pi*e_0)).*((xf - posx(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_xi717 = (rho_i8/(2*pi*e_0)).*((xf - posx(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_xi718 = (rho_i9/(2*pi*e_0)).*((xf - posx(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));

E_xi81 = (-rho_i1/(2*pi*e_0)).*((xf - posx(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_xi82 = (-rho_i2/(2*pi*e_0)).*((xf - posx(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_xi83 = (-rho_i3/(2*pi*e_0)).*((xf - posx(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_xi84 = (-rho_i4/(2*pi*e_0)).*((xf - posx(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_xi85 = (-rho_i5/(2*pi*e_0)).*((xf - posx(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_xi86 = (-rho_i6/(2*pi*e_0)).*((xf - posx(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_xi87 = (-rho_i7/(2*pi*e_0)).*((xf - posx(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_xi89 = (-rho_i9/(2*pi*e_0)).*((xf - posx(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_xi810 = (rho_i1/(2*pi*e_0)).*((xf - posx(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_xi811 = (rho_i2/(2*pi*e_0)).*((xf - posx(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_xi812 = (rho_i3/(2*pi*e_0)).*((xf - posx(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_xi813 = (rho_i4/(2*pi*e_0)).*((xf - posx(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_xi814 = (rho_i5/(2*pi*e_0)).*((xf - posx(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_xi815 = (rho_i6/(2*pi*e_0)).*((xf - posx(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_xi816 = (rho_i7/(2*pi*e_0)).*((xf - posx(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_xi817 = (rho_i8/(2*pi*e_0)).*((xf - posx(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_xi818 = (rho_i9/(2*pi*e_0)).*((xf - posx(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));

E_xi91 = (-rho_i1/(2*pi*e_0)).*((xf - posx(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_xi92 = (-rho_i2/(2*pi*e_0)).*((xf - posx(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_xi93 = (-rho_i3/(2*pi*e_0)).*((xf - posx(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_xi94 = (-rho_i4/(2*pi*e_0)).*((xf - posx(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_xi95 = (-rho_i5/(2*pi*e_0)).*((xf - posx(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_xi96 = (-rho_i6/(2*pi*e_0)).*((xf - posx(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_xi97 = (-rho_i7/(2*pi*e_0)).*((xf - posx(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_xi98 = (-rho_i8/(2*pi*e_0)).*((xf - posx(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_xi910 = (rho_i1/(2*pi*e_0)).*((xf - posx(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_xi911 = (rho_i2/(2*pi*e_0)).*((xf - posx(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_xi912 = (rho_i3/(2*pi*e_0)).*((xf - posx(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_xi913 = (rho_i4/(2*pi*e_0)).*((xf - posx(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_xi914 = (rho_i5/(2*pi*e_0)).*((xf - posx(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_xi915 = (rho_i6/(2*pi*e_0)).*((xf - posx(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_xi916 = (rho_i7/(2*pi*e_0)).*((xf - posx(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_xi917 = (rho_i8/(2*pi*e_0)).*((xf - posx(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_xi918 = (rho_i9/(2*pi*e_0)).*((xf - posx(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));

E_xi101 = (-rho_i1/(2*pi*e_0)).*((xf - posx(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_xi102 = (-rho_i2/(2*pi*e_0)).*((xf - posx(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_xi103 = (-rho_i3/(2*pi*e_0)).*((xf - posx(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_xi104 = (-rho_i4/(2*pi*e_0)).*((xf - posx(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_xi105 = (-rho_i5/(2*pi*e_0)).*((xf - posx(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_xi106 = (-rho_i6/(2*pi*e_0)).*((xf - posx(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_xi107 = (-rho_i7/(2*pi*e_0)).*((xf - posx(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_xi108 = (-rho_i8/(2*pi*e_0)).*((xf - posx(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_xi109 = (-rho_i9/(2*pi*e_0)).*((xf - posx(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_xi1011 = (rho_i2/(2*pi*e_0)).*((xf - posx(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_xi1012 = (rho_i3/(2*pi*e_0)).*((xf - posx(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_xi1013 = (rho_i4/(2*pi*e_0)).*((xf - posx(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_xi1014 = (rho_i5/(2*pi*e_0)).*((xf - posx(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_xi1015 = (rho_i6/(2*pi*e_0)).*((xf - posx(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_xi1016 = (rho_i7/(2*pi*e_0)).*((xf - posx(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_xi1017 = (rho_i8/(2*pi*e_0)).*((xf - posx(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_xi1018 = (rho_i9/(2*pi*e_0)).*((xf - posx(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));

E_xi11_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_xi11_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_xi113 = (-rho_i3/(2*pi*e_0)).*((xf - posx(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_xi114 = (-rho_i4/(2*pi*e_0)).*((xf - posx(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_xi115 = (-rho_i5/(2*pi*e_0)).*((xf - posx(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_xi116 = (-rho_i6/(2*pi*e_0)).*((xf - posx(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_xi117 = (-rho_i7/(2*pi*e_0)).*((xf - posx(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_xi118 = (-rho_i8/(2*pi*e_0)).*((xf - posx(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_xi119 = (-rho_i9/(2*pi*e_0)).*((xf - posx(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_xi1110 = (rho_i1/(2*pi*e_0)).*((xf - posx(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_xi1112 = (rho_i3/(2*pi*e_0)).*((xf - posx(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_xi1113 = (rho_i4/(2*pi*e_0)).*((xf - posx(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_xi1114 = (rho_i5/(2*pi*e_0)).*((xf - posx(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_xi1115 = (rho_i6/(2*pi*e_0)).*((xf - posx(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_xi1116 = (rho_i7/(2*pi*e_0)).*((xf - posx(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_xi1117 = (rho_i8/(2*pi*e_0)).*((xf - posx(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_xi1118 = (rho_i9/(2*pi*e_0)).*((xf - posx(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));

E_xi12_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_xi12_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_xi123 = (-rho_i3/(2*pi*e_0)).*((xf - posx(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_xi124 = (-rho_i4/(2*pi*e_0)).*((xf - posx(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_xi125 = (-rho_i5/(2*pi*e_0)).*((xf - posx(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_xi126 = (-rho_i6/(2*pi*e_0)).*((xf - posx(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_xi127 = (-rho_i7/(2*pi*e_0)).*((xf - posx(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_xi128 = (-rho_i8/(2*pi*e_0)).*((xf - posx(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_xi129 = (-rho_i9/(2*pi*e_0)).*((xf - posx(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_xi1210 = (rho_i1/(2*pi*e_0)).*((xf - posx(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_xi1211 = (rho_i2/(2*pi*e_0)).*((xf - posx(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_xi1213 = (rho_i4/(2*pi*e_0)).*((xf - posx(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_xi1214 = (rho_i5/(2*pi*e_0)).*((xf - posx(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_xi1215 = (rho_i6/(2*pi*e_0)).*((xf - posx(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_xi1216 = (rho_i7/(2*pi*e_0)).*((xf - posx(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_xi1217 = (rho_i8/(2*pi*e_0)).*((xf - posx(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_xi1218 = (rho_i9/(2*pi*e_0)).*((xf - posx(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));

E_xi131 = (-rho_i1/(2*pi*e_0)).*((xf - posx(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_xi132 = (-rho_i2/(2*pi*e_0)).*((xf - posx(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_xi133 = (-rho_i3/(2*pi*e_0)).*((xf - posx(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_xi134 = (-rho_i4/(2*pi*e_0)).*((xf - posx(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_xi135 = (-rho_i5/(2*pi*e_0)).*((xf - posx(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_xi136 = (-rho_i6/(2*pi*e_0)).*((xf - posx(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_xi137 = (-rho_i7/(2*pi*e_0)).*((xf - posx(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_xi138 = (-rho_i8/(2*pi*e_0)).*((xf - posx(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_xi139 = (-rho_i9/(2*pi*e_0)).*((xf - posx(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_xi1310 = (rho_i1/(2*pi*e_0)).*((xf - posx(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_xi1311 = (rho_i2/(2*pi*e_0)).*((xf - posx(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_xi1312 = (rho_i3/(2*pi*e_0)).*((xf - posx(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_xi1314 = (rho_i5/(2*pi*e_0)).*((xf - posx(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_xi1315 = (rho_i6/(2*pi*e_0)).*((xf - posx(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_xi1316 = (rho_i7/(2*pi*e_0)).*((xf - posx(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_xi1317 = (rho_i8/(2*pi*e_0)).*((xf - posx(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_xi1318 = (rho_i9/(2*pi*e_0)).*((xf - posx(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));

E_xi141 = (-rho_i1/(2*pi*e_0)).*((xf - posx(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_xi142 = (-rho_i2/(2*pi*e_0)).*((xf - posx(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_xi143 = (-rho_i3/(2*pi*e_0)).*((xf - posx(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_xi144 = (-rho_i4/(2*pi*e_0)).*((xf - posx(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_xi145 = (-rho_i5/(2*pi*e_0)).*((xf - posx(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_xi146 = (-rho_i6/(2*pi*e_0)).*((xf - posx(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_xi147 = (-rho_i7/(2*pi*e_0)).*((xf - posx(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_xi148 = (-rho_i8/(2*pi*e_0)).*((xf - posx(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_xi149 = (-rho_i9/(2*pi*e_0)).*((xf - posx(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_xi1410 = (rho_i1/(2*pi*e_0)).*((xf - posx(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_xi1411 = (rho_i2/(2*pi*e_0)).*((xf - posx(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_xi1412 = (rho_i3/(2*pi*e_0)).*((xf - posx(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_xi1413 = (rho_i4/(2*pi*e_0)).*((xf - posx(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_xi1415 = (rho_i6/(2*pi*e_0)).*((xf - posx(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_xi1416 = (rho_i7/(2*pi*e_0)).*((xf - posx(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_xi1417 = (rho_i8/(2*pi*e_0)).*((xf - posx(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_xi1418 = (rho_i9/(2*pi*e_0)).*((xf - posx(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));

E_xi151 = (-rho_i1/(2*pi*e_0)).*((xf - posx(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_xi152 = (-rho_i2/(2*pi*e_0)).*((xf - posx(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_xi153 = (-rho_i3/(2*pi*e_0)).*((xf - posx(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_xi154 = (-rho_i4/(2*pi*e_0)).*((xf - posx(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_xi155 = (-rho_i5/(2*pi*e_0)).*((xf - posx(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_xi156 = (-rho_i6/(2*pi*e_0)).*((xf - posx(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_xi157 = (-rho_i7/(2*pi*e_0)).*((xf - posx(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_xi158 = (-rho_i8/(2*pi*e_0)).*((xf - posx(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_xi159 = (-rho_i9/(2*pi*e_0)).*((xf - posx(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_xi1510 = (rho_i1/(2*pi*e_0)).*((xf - posx(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_xi1511 = (rho_i2/(2*pi*e_0)).*((xf - posx(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_xi1512 = (rho_i3/(2*pi*e_0)).*((xf - posx(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_xi1513 = (rho_i4/(2*pi*e_0)).*((xf - posx(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_xi1514 = (rho_i5/(2*pi*e_0)).*((xf - posx(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_xi1516 = (rho_i7/(2*pi*e_0)).*((xf - posx(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_xi1517 = (rho_i8/(2*pi*e_0)).*((xf - posx(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_xi1518 = (rho_i9/(2*pi*e_0)).*((xf - posx(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));

E_xi161 = (-rho_i1/(2*pi*e_0)).*((xf - posx(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_xi162 = (-rho_i2/(2*pi*e_0)).*((xf - posx(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_xi163 = (-rho_i3/(2*pi*e_0)).*((xf - posx(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_xi164 = (-rho_i4/(2*pi*e_0)).*((xf - posx(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_xi165 = (-rho_i5/(2*pi*e_0)).*((xf - posx(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_xi166 = (-rho_i6/(2*pi*e_0)).*((xf - posx(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_xi167 = (-rho_i7/(2*pi*e_0)).*((xf - posx(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_xi168 = (-rho_i8/(2*pi*e_0)).*((xf - posx(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_xi169 = (-rho_i9/(2*pi*e_0)).*((xf - posx(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_xi1610 = (rho_i1/(2*pi*e_0)).*((xf - posx(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_xi1611 = (rho_i2/(2*pi*e_0)).*((xf - posx(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_xi1612 = (rho_i3/(2*pi*e_0)).*((xf - posx(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_xi1613 = (rho_i4/(2*pi*e_0)).*((xf - posx(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_xi1614 = (rho_i5/(2*pi*e_0)).*((xf - posx(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_xi1615 = (rho_i6/(2*pi*e_0)).*((xf - posx(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_xi1617 = (rho_i8/(2*pi*e_0)).*((xf - posx(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_xi1618 = (rho_i9/(2*pi*e_0)).*((xf - posx(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));

E_xi171 = (-rho_i1/(2*pi*e_0)).*((xf - posx(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_xi172 = (-rho_i2/(2*pi*e_0)).*((xf - posx(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_xi173 = (-rho_i3/(2*pi*e_0)).*((xf - posx(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_xi174 = (-rho_i4/(2*pi*e_0)).*((xf - posx(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_xi175 = (-rho_i5/(2*pi*e_0)).*((xf - posx(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_xi176 = (-rho_i6/(2*pi*e_0)).*((xf - posx(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_xi177 = (-rho_i7/(2*pi*e_0)).*((xf - posx(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_xi178 = (-rho_i8/(2*pi*e_0)).*((xf - posx(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_xi179 = (-rho_i9/(2*pi*e_0)).*((xf - posx(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_xi1710 = (rho_i1/(2*pi*e_0)).*((xf - posx(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_xi1711 = (rho_i2/(2*pi*e_0)).*((xf - posx(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_xi1712 = (rho_i3/(2*pi*e_0)).*((xf - posx(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_xi1713 = (rho_i4/(2*pi*e_0)).*((xf - posx(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_xi1714 = (rho_i5/(2*pi*e_0)).*((xf - posx(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_xi1715 = (rho_i6/(2*pi*e_0)).*((xf - posx(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_xi1716 = (rho_i7/(2*pi*e_0)).*((xf - posx(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_xi1718 = (rho_i9/(2*pi*e_0)).*((xf - posx(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));

E_xi181 = (-rho_i1/(2*pi*e_0)).*((xf - posx(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_xi182 = (-rho_i2/(2*pi*e_0)).*((xf - posx(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_xi183 = (-rho_i3/(2*pi*e_0)).*((xf - posx(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_xi184 = (-rho_i4/(2*pi*e_0)).*((xf - posx(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_xi185 = (-rho_i5/(2*pi*e_0)).*((xf - posx(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_xi186 = (-rho_i6/(2*pi*e_0)).*((xf - posx(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_xi187 = (-rho_i7/(2*pi*e_0)).*((xf - posx(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_xi188 = (-rho_i8/(2*pi*e_0)).*((xf - posx(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_xi189 = (-rho_i9/(2*pi*e_0)).*((xf - posx(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_xi1810 = (rho_i1/(2*pi*e_0)).*((xf - posx(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_xi1811 = (rho_i2/(2*pi*e_0)).*((xf - posx(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_xi1812 = (rho_i3/(2*pi*e_0)).*((xf - posx(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_xi1813 = (rho_i4/(2*pi*e_0)).*((xf - posx(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_xi1814 = (rho_i5/(2*pi*e_0)).*((xf - posx(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_xi1815 = (rho_i6/(2*pi*e_0)).*((xf - posx(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_xi1816 = (rho_i7/(2*pi*e_0)).*((xf - posx(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_xi1817 = (rho_i8/(2*pi*e_0)).*((xf - posx(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));

%% E_yr componente y real campo elétrico condutor 2 fase b

E_yr12 = (-rho_r2/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yr13 = (-rho_r3/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yr14 = (-rho_r4/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yr15 = (-rho_r5/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yr16 = (-rho_r6/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_yr17 = (-rho_r7/(2*pi*e_0)).*((yf - posy(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_yr18 = (-rho_r8/(2*pi*e_0)).*((yf - posy(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_yr19 = (-rho_r9/(2*pi*e_0)).*((yf - posy(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_yr110 = (rho_r1/(2*pi*e_0)).*((yf - posy(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_yr1_11 = (rho_r2/(2*pi*e_0)).*((yf - posy(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_yr1_12 = (rho_r3/(2*pi*e_0)).*((yf - posy(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_yr1_13 = (rho_r4/(2*pi*e_0)).*((yf - posy(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_yr1_14 = (rho_r5/(2*pi*e_0)).*((yf - posy(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_yr1_15 = (rho_r6/(2*pi*e_0)).*((yf - posy(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_yr1_16 = (rho_r7/(2*pi*e_0)).*((yf - posy(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_yr1_17 = (rho_r8/(2*pi*e_0)).*((yf - posy(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_yr1_18 = (rho_r9/(2*pi*e_0)).*((yf - posy(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));

E_yr21 = (-rho_r1/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yr23 = (-rho_r3/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yr24 = (-rho_r4/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yr25 = (-rho_r5/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yr26 = (-rho_r6/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_yr27 = (-rho_r7/(2*pi*e_0)).*((yf - posy(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_yr28 = (-rho_r8/(2*pi*e_0)).*((yf - posy(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_yr29 = (-rho_r9/(2*pi*e_0)).*((yf - posy(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_yr210 = (rho_r1/(2*pi*e_0)).*((yf - posy(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_yr2_11 = (rho_r2/(2*pi*e_0)).*((yf - posy(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_yr2_12 = (rho_r3/(2*pi*e_0)).*((yf - posy(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_yr2_13 = (rho_r4/(2*pi*e_0)).*((yf - posy(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_yr2_14 = (rho_r5/(2*pi*e_0)).*((yf - posy(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_yr2_15 = (rho_r6/(2*pi*e_0)).*((yf - posy(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_yr2_16 = (rho_r7/(2*pi*e_0)).*((yf - posy(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_yr2_17 = (rho_r8/(2*pi*e_0)).*((yf - posy(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_yr2_18 = (rho_r9/(2*pi*e_0)).*((yf - posy(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));

E_yr31 = (-rho_r1/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yr32 = (-rho_r2/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yr34 = (-rho_r4/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yr35 = (-rho_r5/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yr36 = (-rho_r6/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_yr37 = (-rho_r7/(2*pi*e_0)).*((yf - posy(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_yr38 = (-rho_r8/(2*pi*e_0)).*((yf - posy(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_yr39 = (-rho_r9/(2*pi*e_0)).*((yf - posy(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_yr310 = (rho_r1/(2*pi*e_0)).*((yf - posy(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_yr311 = (rho_r2/(2*pi*e_0)).*((yf - posy(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_yr312 = (rho_r3/(2*pi*e_0)).*((yf - posy(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_yr313 = (rho_r4/(2*pi*e_0)).*((yf - posy(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_yr314 = (rho_r5/(2*pi*e_0)).*((yf - posy(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_yr315 = (rho_r6/(2*pi*e_0)).*((yf - posy(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_yr316 = (rho_r7/(2*pi*e_0)).*((yf - posy(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_yr317 = (rho_r8/(2*pi*e_0)).*((yf - posy(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_yr318 = (rho_r9/(2*pi*e_0)).*((yf - posy(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));

E_yr41 = (-rho_r1/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yr42 = (-rho_r2/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yr43 = (-rho_r3/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yr45 = (-rho_r5/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yr46 = (-rho_r6/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_yr47 = (-rho_r7/(2*pi*e_0)).*((yf - posy(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_yr48 = (-rho_r8/(2*pi*e_0)).*((yf - posy(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_yr49 = (-rho_r9/(2*pi*e_0)).*((yf - posy(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_yr410 = (rho_r1/(2*pi*e_0)).*((yf - posy(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_yr411 = (rho_r2/(2*pi*e_0)).*((yf - posy(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_yr412 = (rho_r3/(2*pi*e_0)).*((yf - posy(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_yr413 = (rho_r4/(2*pi*e_0)).*((yf - posy(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_yr414 = (rho_r5/(2*pi*e_0)).*((yf - posy(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_yr415 = (rho_r6/(2*pi*e_0)).*((yf - posy(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_yr416 = (rho_r7/(2*pi*e_0)).*((yf - posy(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_yr417 = (rho_r8/(2*pi*e_0)).*((yf - posy(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_yr418 = (rho_r9/(2*pi*e_0)).*((yf - posy(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));

E_yr51 = (-rho_r1/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yr52 = (-rho_r2/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yr53 = (-rho_r3/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yr54 = (-rho_r4/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yr56 = (-rho_r6/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_yr57 = (-rho_r7/(2*pi*e_0)).*((yf - posy(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_yr58 = (-rho_r8/(2*pi*e_0)).*((yf - posy(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_yr59 = (-rho_r9/(2*pi*e_0)).*((yf - posy(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_yr510 = (rho_r1/(2*pi*e_0)).*((yf - posy(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_yr511 = (rho_r2/(2*pi*e_0)).*((yf - posy(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_yr512 = (rho_r3/(2*pi*e_0)).*((yf - posy(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_yr513 = (rho_r4/(2*pi*e_0)).*((yf - posy(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_yr514 = (rho_r5/(2*pi*e_0)).*((yf - posy(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_yr515 = (rho_r6/(2*pi*e_0)).*((yf - posy(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_yr516 = (rho_r7/(2*pi*e_0)).*((yf - posy(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_yr517 = (rho_r8/(2*pi*e_0)).*((yf - posy(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_yr518 = (rho_r9/(2*pi*e_0)).*((yf - posy(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));

E_yr61 = (-rho_r1/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yr62 = (-rho_r2/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yr63 = (-rho_r3/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yr64 = (-rho_r4/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yr65 = (-rho_r5/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_yr67 = (-rho_r7/(2*pi*e_0)).*((yf - posy(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_yr68 = (-rho_r8/(2*pi*e_0)).*((yf - posy(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_yr69 = (-rho_r9/(2*pi*e_0)).*((yf - posy(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_yr610 = (rho_r1/(2*pi*e_0)).*((yf - posy(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_yr611 = (rho_r2/(2*pi*e_0)).*((yf - posy(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_yr612 = (rho_r3/(2*pi*e_0)).*((yf - posy(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_yr613 = (rho_r4/(2*pi*e_0)).*((yf - posy(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_yr614 = (rho_r5/(2*pi*e_0)).*((yf - posy(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_yr615 = (rho_r6/(2*pi*e_0)).*((yf - posy(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_yr616 = (rho_r7/(2*pi*e_0)).*((yf - posy(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_yr617 = (rho_r8/(2*pi*e_0)).*((yf - posy(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_yr618 = (rho_r9/(2*pi*e_0)).*((yf - posy(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));

E_yr71 = (-rho_r1/(2*pi*e_0)).*((yf - posy(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_yr72 = (-rho_r2/(2*pi*e_0)).*((yf - posy(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_yr73 = (-rho_r3/(2*pi*e_0)).*((yf - posy(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_yr74 = (-rho_r4/(2*pi*e_0)).*((yf - posy(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_yr75 = (-rho_r5/(2*pi*e_0)).*((yf - posy(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_yr76 = (-rho_r6/(2*pi*e_0)).*((yf - posy(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_yr78 = (-rho_r8/(2*pi*e_0)).*((yf - posy(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_yr79 = (-rho_r9/(2*pi*e_0)).*((yf - posy(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_yr710 = (rho_r1/(2*pi*e_0)).*((yf - posy(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_yr711 = (rho_r2/(2*pi*e_0)).*((yf - posy(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_yr712 = (rho_r3/(2*pi*e_0)).*((yf - posy(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_yr713 = (rho_r4/(2*pi*e_0)).*((yf - posy(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_yr714 = (rho_r5/(2*pi*e_0)).*((yf - posy(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_yr715 = (rho_r6/(2*pi*e_0)).*((yf - posy(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_yr716 = (rho_r7/(2*pi*e_0)).*((yf - posy(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_yr717 = (rho_r8/(2*pi*e_0)).*((yf - posy(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_yr718 = (rho_r9/(2*pi*e_0)).*((yf - posy(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));

E_yr81 = (-rho_r1/(2*pi*e_0)).*((yf - posy(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_yr82 = (-rho_r2/(2*pi*e_0)).*((yf - posy(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_yr83 = (-rho_r3/(2*pi*e_0)).*((yf - posy(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_yr84 = (-rho_r4/(2*pi*e_0)).*((yf - posy(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_yr85 = (-rho_r5/(2*pi*e_0)).*((yf - posy(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_yr86 = (-rho_r6/(2*pi*e_0)).*((yf - posy(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_yr87 = (-rho_r7/(2*pi*e_0)).*((yf - posy(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_yr89 = (-rho_r9/(2*pi*e_0)).*((yf - posy(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_yr810 = (rho_r1/(2*pi*e_0)).*((yf - posy(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_yr811 = (rho_r2/(2*pi*e_0)).*((yf - posy(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_yr812 = (rho_r3/(2*pi*e_0)).*((yf - posy(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_yr813 = (rho_r4/(2*pi*e_0)).*((yf - posy(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_yr814 = (rho_r5/(2*pi*e_0)).*((yf - posy(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_yr815 = (rho_r6/(2*pi*e_0)).*((yf - posy(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_yr816 = (rho_r7/(2*pi*e_0)).*((yf - posy(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_yr817 = (rho_r8/(2*pi*e_0)).*((yf - posy(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_yr818 = (rho_r9/(2*pi*e_0)).*((yf - posy(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));

E_yr91 = (-rho_r1/(2*pi*e_0)).*((yf - posy(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_yr92 = (-rho_r2/(2*pi*e_0)).*((yf - posy(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_yr93 = (-rho_r3/(2*pi*e_0)).*((yf - posy(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_yr94 = (-rho_r4/(2*pi*e_0)).*((yf - posy(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_yr95 = (-rho_r5/(2*pi*e_0)).*((yf - posy(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_yr96 = (-rho_r6/(2*pi*e_0)).*((yf - posy(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_yr97 = (-rho_r7/(2*pi*e_0)).*((yf - posy(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_yr98 = (-rho_r8/(2*pi*e_0)).*((yf - posy(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_yr910 = (rho_r1/(2*pi*e_0)).*((yf - posy(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_yr911 = (rho_r2/(2*pi*e_0)).*((yf - posy(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_yr912 = (rho_r3/(2*pi*e_0)).*((yf - posy(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_yr913 = (rho_r4/(2*pi*e_0)).*((yf - posy(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_yr914 = (rho_r5/(2*pi*e_0)).*((yf - posy(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_yr915 = (rho_r6/(2*pi*e_0)).*((yf - posy(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_yr916 = (rho_r7/(2*pi*e_0)).*((yf - posy(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_yr917 = (rho_r8/(2*pi*e_0)).*((yf - posy(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_yr918 = (rho_r9/(2*pi*e_0)).*((yf - posy(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));

E_yr101 = (-rho_r1/(2*pi*e_0)).*((yf - posy(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_yr102 = (-rho_r2/(2*pi*e_0)).*((yf - posy(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_yr103 = (-rho_r3/(2*pi*e_0)).*((yf - posy(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_yr104 = (-rho_r4/(2*pi*e_0)).*((yf - posy(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_yr105 = (-rho_r5/(2*pi*e_0)).*((yf - posy(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_yr106 = (-rho_r6/(2*pi*e_0)).*((yf - posy(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_yr107 = (-rho_r7/(2*pi*e_0)).*((yf - posy(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_yr108 = (-rho_r8/(2*pi*e_0)).*((yf - posy(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_yr109 = (-rho_r9/(2*pi*e_0)).*((yf - posy(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_yr1011 = (rho_r2/(2*pi*e_0)).*((yf - posy(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_yr1012 = (rho_r3/(2*pi*e_0)).*((yf - posy(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_yr1013 = (rho_r4/(2*pi*e_0)).*((yf - posy(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_yr1014 = (rho_r5/(2*pi*e_0)).*((yf - posy(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_yr1015 = (rho_r6/(2*pi*e_0)).*((yf - posy(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_yr1016 = (rho_r7/(2*pi*e_0)).*((yf - posy(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_yr1017 = (rho_r8/(2*pi*e_0)).*((yf - posy(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_yr1018 = (rho_r9/(2*pi*e_0)).*((yf - posy(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));

E_yr11_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_yr11_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_yr113 = (-rho_r3/(2*pi*e_0)).*((yf - posy(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_yr114 = (-rho_r4/(2*pi*e_0)).*((yf - posy(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_yr115 = (-rho_r5/(2*pi*e_0)).*((yf - posy(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_yr116 = (-rho_r6/(2*pi*e_0)).*((yf - posy(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_yr117 = (-rho_r7/(2*pi*e_0)).*((yf - posy(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_yr118 = (-rho_r8/(2*pi*e_0)).*((yf - posy(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_yr119 = (-rho_r9/(2*pi*e_0)).*((yf - posy(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_yr1110 = (rho_r1/(2*pi*e_0)).*((yf - posy(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_yr1112 = (rho_r3/(2*pi*e_0)).*((yf - posy(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_yr1113 = (rho_r4/(2*pi*e_0)).*((yf - posy(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_yr1114 = (rho_r5/(2*pi*e_0)).*((yf - posy(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_yr1115 = (rho_r6/(2*pi*e_0)).*((yf - posy(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_yr1116 = (rho_r7/(2*pi*e_0)).*((yf - posy(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_yr1117 = (rho_r8/(2*pi*e_0)).*((yf - posy(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_yr1118 = (rho_r9/(2*pi*e_0)).*((yf - posy(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));

E_yr12_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_yr12_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_yr123 = (-rho_r3/(2*pi*e_0)).*((yf - posy(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_yr124 = (-rho_r4/(2*pi*e_0)).*((yf - posy(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_yr125 = (-rho_r5/(2*pi*e_0)).*((yf - posy(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_yr126 = (-rho_r6/(2*pi*e_0)).*((yf - posy(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_yr127 = (-rho_r7/(2*pi*e_0)).*((yf - posy(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_yr128 = (-rho_r8/(2*pi*e_0)).*((yf - posy(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_yr129 = (-rho_r9/(2*pi*e_0)).*((yf - posy(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_yr1210 = (rho_r1/(2*pi*e_0)).*((yf - posy(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_yr1211 = (rho_r2/(2*pi*e_0)).*((yf - posy(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_yr1213 = (rho_r4/(2*pi*e_0)).*((yf - posy(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_yr1214 = (rho_r5/(2*pi*e_0)).*((yf - posy(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_yr1215 = (rho_r6/(2*pi*e_0)).*((yf - posy(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_yr1216 = (rho_r7/(2*pi*e_0)).*((yf - posy(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_yr1217 = (rho_r8/(2*pi*e_0)).*((yf - posy(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_yr1218 = (rho_r9/(2*pi*e_0)).*((yf - posy(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));

E_yr131 = (-rho_r1/(2*pi*e_0)).*((yf - posy(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_yr132 = (-rho_r2/(2*pi*e_0)).*((yf - posy(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_yr133 = (-rho_r3/(2*pi*e_0)).*((yf - posy(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_yr134 = (-rho_r4/(2*pi*e_0)).*((yf - posy(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_yr135 = (-rho_r5/(2*pi*e_0)).*((yf - posy(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_yr136 = (-rho_r6/(2*pi*e_0)).*((yf - posy(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_yr137 = (-rho_r7/(2*pi*e_0)).*((yf - posy(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_yr138 = (-rho_r8/(2*pi*e_0)).*((yf - posy(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_yr139 = (-rho_r9/(2*pi*e_0)).*((yf - posy(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_yr1310 = (rho_r1/(2*pi*e_0)).*((yf - posy(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_yr1311 = (rho_r2/(2*pi*e_0)).*((yf - posy(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_yr1312 = (rho_r3/(2*pi*e_0)).*((yf - posy(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_yr1314 = (rho_r5/(2*pi*e_0)).*((yf - posy(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_yr1315 = (rho_r6/(2*pi*e_0)).*((yf - posy(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_yr1316 = (rho_r7/(2*pi*e_0)).*((yf - posy(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_yr1317 = (rho_r8/(2*pi*e_0)).*((yf - posy(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_yr1318 = (rho_r9/(2*pi*e_0)).*((yf - posy(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));

E_yr141 = (-rho_r1/(2*pi*e_0)).*((yf - posy(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_yr142 = (-rho_r2/(2*pi*e_0)).*((yf - posy(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_yr143 = (-rho_r3/(2*pi*e_0)).*((yf - posy(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_yr144 = (-rho_r4/(2*pi*e_0)).*((yf - posy(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_yr145 = (-rho_r5/(2*pi*e_0)).*((yf - posy(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_yr146 = (-rho_r6/(2*pi*e_0)).*((yf - posy(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_yr147 = (-rho_r7/(2*pi*e_0)).*((yf - posy(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_yr148 = (-rho_r8/(2*pi*e_0)).*((yf - posy(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_yr149 = (-rho_r9/(2*pi*e_0)).*((yf - posy(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_yr1410 = (rho_r1/(2*pi*e_0)).*((yf - posy(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_yr1411 = (rho_r2/(2*pi*e_0)).*((yf - posy(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_yr1412 = (rho_r3/(2*pi*e_0)).*((yf - posy(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_yr1413 = (rho_r4/(2*pi*e_0)).*((yf - posy(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_yr1415 = (rho_r6/(2*pi*e_0)).*((yf - posy(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_yr1416 = (rho_r7/(2*pi*e_0)).*((yf - posy(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_yr1417 = (rho_r8/(2*pi*e_0)).*((yf - posy(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_yr1418 = (rho_r9/(2*pi*e_0)).*((yf - posy(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));

E_yr151 = (-rho_r1/(2*pi*e_0)).*((yf - posy(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_yr152 = (-rho_r2/(2*pi*e_0)).*((yf - posy(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_yr153 = (-rho_r3/(2*pi*e_0)).*((yf - posy(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_yr154 = (-rho_r4/(2*pi*e_0)).*((yf - posy(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_yr155 = (-rho_r5/(2*pi*e_0)).*((yf - posy(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_yr156 = (-rho_r6/(2*pi*e_0)).*((yf - posy(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_yr157 = (-rho_r7/(2*pi*e_0)).*((yf - posy(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_yr158 = (-rho_r8/(2*pi*e_0)).*((yf - posy(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_yr159 = (-rho_r9/(2*pi*e_0)).*((yf - posy(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_yr1510 = (rho_r1/(2*pi*e_0)).*((yf - posy(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_yr1511 = (rho_r2/(2*pi*e_0)).*((yf - posy(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_yr1512 = (rho_r3/(2*pi*e_0)).*((yf - posy(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_yr1513 = (rho_r4/(2*pi*e_0)).*((yf - posy(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_yr1514 = (rho_r5/(2*pi*e_0)).*((yf - posy(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_yr1516 = (rho_r7/(2*pi*e_0)).*((yf - posy(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_yr1517 = (rho_r8/(2*pi*e_0)).*((yf - posy(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_yr1518 = (rho_r9/(2*pi*e_0)).*((yf - posy(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));

E_yr161 = (-rho_r1/(2*pi*e_0)).*((yf - posy(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_yr162 = (-rho_r2/(2*pi*e_0)).*((yf - posy(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_yr163 = (-rho_r3/(2*pi*e_0)).*((yf - posy(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_yr164 = (-rho_r4/(2*pi*e_0)).*((yf - posy(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_yr165 = (-rho_r5/(2*pi*e_0)).*((yf - posy(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_yr166 = (-rho_r6/(2*pi*e_0)).*((yf - posy(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_yr167 = (-rho_r7/(2*pi*e_0)).*((yf - posy(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_yr168 = (-rho_r8/(2*pi*e_0)).*((yf - posy(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_yr169 = (-rho_r9/(2*pi*e_0)).*((yf - posy(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_yr1610 = (rho_r1/(2*pi*e_0)).*((yf - posy(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_yr1611 = (rho_r2/(2*pi*e_0)).*((yf - posy(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_yr1612 = (rho_r3/(2*pi*e_0)).*((yf - posy(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_yr1613 = (rho_r4/(2*pi*e_0)).*((yf - posy(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_yr1614 = (rho_r5/(2*pi*e_0)).*((yf - posy(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_yr1615 = (rho_r6/(2*pi*e_0)).*((yf - posy(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_yr1617 = (rho_r8/(2*pi*e_0)).*((yf - posy(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_yr1618 = (rho_r9/(2*pi*e_0)).*((yf - posy(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));

E_yr171 = (-rho_r1/(2*pi*e_0)).*((yf - posy(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_yr172 = (-rho_r2/(2*pi*e_0)).*((yf - posy(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_yr173 = (-rho_r3/(2*pi*e_0)).*((yf - posy(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_yr174 = (-rho_r4/(2*pi*e_0)).*((yf - posy(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_yr175 = (-rho_r5/(2*pi*e_0)).*((yf - posy(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_yr176 = (-rho_r6/(2*pi*e_0)).*((yf - posy(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_yr177 = (-rho_r7/(2*pi*e_0)).*((yf - posy(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_yr178 = (-rho_r8/(2*pi*e_0)).*((yf - posy(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_yr179 = (-rho_r9/(2*pi*e_0)).*((yf - posy(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_yr1710 = (rho_r1/(2*pi*e_0)).*((yf - posy(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_yr1711 = (rho_r2/(2*pi*e_0)).*((yf - posy(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_yr1712 = (rho_r3/(2*pi*e_0)).*((yf - posy(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_yr1713 = (rho_r4/(2*pi*e_0)).*((yf - posy(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_yr1714 = (rho_r5/(2*pi*e_0)).*((yf - posy(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_yr1715 = (rho_r6/(2*pi*e_0)).*((yf - posy(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_yr1716 = (rho_r7/(2*pi*e_0)).*((yf - posy(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_yr1718 = (rho_r9/(2*pi*e_0)).*((yf - posy(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));

E_yr181 = (-rho_r1/(2*pi*e_0)).*((yf - posy(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_yr182 = (-rho_r2/(2*pi*e_0)).*((yf - posy(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_yr183 = (-rho_r3/(2*pi*e_0)).*((yf - posy(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_yr184 = (-rho_r4/(2*pi*e_0)).*((yf - posy(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_yr185 = (-rho_r5/(2*pi*e_0)).*((yf - posy(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_yr186 = (-rho_r6/(2*pi*e_0)).*((yf - posy(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_yr187 = (-rho_r7/(2*pi*e_0)).*((yf - posy(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_yr188 = (-rho_r8/(2*pi*e_0)).*((yf - posy(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_yr189 = (-rho_r9/(2*pi*e_0)).*((yf - posy(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_yr1810 = (rho_r1/(2*pi*e_0)).*((yf - posy(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_yr1811 = (rho_r2/(2*pi*e_0)).*((yf - posy(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_yr1812 = (rho_r3/(2*pi*e_0)).*((yf - posy(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_yr1813 = (rho_r4/(2*pi*e_0)).*((yf - posy(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_yr1814 = (rho_r5/(2*pi*e_0)).*((yf - posy(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_yr1815 = (rho_r6/(2*pi*e_0)).*((yf - posy(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_yr1816 = (rho_r7/(2*pi*e_0)).*((yf - posy(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_yr1817 = (rho_r8/(2*pi*e_0)).*((yf - posy(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));

%% E_yi21 componente y imaginario campo elétrico condutor 2 fase b

E_yi12 = (-rho_i2/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yi13 = (-rho_i3/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yi14 = (-rho_i4/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yi15 = (-rho_i5/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yi16 = (-rho_i6/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_yi17 = (-rho_i7/(2*pi*e_0)).*((yf - posy(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_yi18 = (-rho_i8/(2*pi*e_0)).*((yf - posy(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_yi19 = (-rho_i9/(2*pi*e_0)).*((yf - posy(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_yi110 = (rho_i1/(2*pi*e_0)).*((yf - posy(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_yi1_11 = (rho_i2/(2*pi*e_0)).*((yf - posy(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_yi1_12 = (rho_i3/(2*pi*e_0)).*((yf - posy(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_yi1_13 = (rho_i4/(2*pi*e_0)).*((yf - posy(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_yi1_14 = (rho_i5/(2*pi*e_0)).*((yf - posy(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_yi1_15 = (rho_i6/(2*pi*e_0)).*((yf - posy(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_yi1_16 = (rho_i7/(2*pi*e_0)).*((yf - posy(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_yi1_17 = (rho_i8/(2*pi*e_0)).*((yf - posy(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_yi1_18 = (rho_i9/(2*pi*e_0)).*((yf - posy(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));

E_yi21 = (-rho_i1/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yi23 = (-rho_i3/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yi24 = (-rho_i4/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yi25 = (-rho_i5/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yi26 = (-rho_i6/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_yi27 = (-rho_i7/(2*pi*e_0)).*((yf - posy(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_yi28 = (-rho_i8/(2*pi*e_0)).*((yf - posy(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_yi29 = (-rho_i9/(2*pi*e_0)).*((yf - posy(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_yi210 = (rho_i1/(2*pi*e_0)).*((yf - posy(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_yi2_11 = (rho_i2/(2*pi*e_0)).*((yf - posy(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_yi2_12 = (rho_i3/(2*pi*e_0)).*((yf - posy(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_yi2_13 = (rho_i4/(2*pi*e_0)).*((yf - posy(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_yi2_14 = (rho_i5/(2*pi*e_0)).*((yf - posy(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_yi2_15 = (rho_i6/(2*pi*e_0)).*((yf - posy(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_yi2_16 = (rho_i7/(2*pi*e_0)).*((yf - posy(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_yi2_17 = (rho_i8/(2*pi*e_0)).*((yf - posy(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_yi2_18 = (rho_i9/(2*pi*e_0)).*((yf - posy(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));

E_yi31 = (-rho_i1/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yi32 = (-rho_i2/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yi34 = (-rho_i4/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yi35 = (-rho_i5/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yi36 = (-rho_i6/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_yi37 = (-rho_i7/(2*pi*e_0)).*((yf - posy(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_yi38 = (-rho_i8/(2*pi*e_0)).*((yf - posy(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_yi39 = (-rho_i9/(2*pi*e_0)).*((yf - posy(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_yi310 = (rho_i1/(2*pi*e_0)).*((yf - posy(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_yi311 = (rho_i2/(2*pi*e_0)).*((yf - posy(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_yi312 = (rho_i3/(2*pi*e_0)).*((yf - posy(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_yi313 = (rho_i4/(2*pi*e_0)).*((yf - posy(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_yi314 = (rho_i5/(2*pi*e_0)).*((yf - posy(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_yi315 = (rho_i6/(2*pi*e_0)).*((yf - posy(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_yi316 = (rho_i7/(2*pi*e_0)).*((yf - posy(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_yi317 = (rho_i8/(2*pi*e_0)).*((yf - posy(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_yi318 = (rho_i9/(2*pi*e_0)).*((yf - posy(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));

E_yi41 = (-rho_i1/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yi42 = (-rho_i2/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yi43 = (-rho_i3/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yi45 = (-rho_i5/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yi46 = (-rho_i6/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_yi47 = (-rho_i7/(2*pi*e_0)).*((yf - posy(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_yi48 = (-rho_i8/(2*pi*e_0)).*((yf - posy(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_yi49 = (-rho_i9/(2*pi*e_0)).*((yf - posy(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_yi410 = (rho_i1/(2*pi*e_0)).*((yf - posy(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_yi411 = (rho_i2/(2*pi*e_0)).*((yf - posy(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_yi412 = (rho_i3/(2*pi*e_0)).*((yf - posy(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_yi413 = (rho_i4/(2*pi*e_0)).*((yf - posy(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_yi414 = (rho_i5/(2*pi*e_0)).*((yf - posy(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_yi415 = (rho_i6/(2*pi*e_0)).*((yf - posy(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_yi416 = (rho_i7/(2*pi*e_0)).*((yf - posy(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_yi417 = (rho_i8/(2*pi*e_0)).*((yf - posy(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_yi418 = (rho_i9/(2*pi*e_0)).*((yf - posy(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));

E_yi51 = (-rho_i1/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yi52 = (-rho_i2/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yi53 = (-rho_i3/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yi54 = (-rho_i4/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yi56 = (-rho_i6/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_yi57 = (-rho_i7/(2*pi*e_0)).*((yf - posy(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_yi58 = (-rho_i8/(2*pi*e_0)).*((yf - posy(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_yi59 = (-rho_i9/(2*pi*e_0)).*((yf - posy(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_yi510 = (rho_i1/(2*pi*e_0)).*((yf - posy(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_yi511 = (rho_i2/(2*pi*e_0)).*((yf - posy(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_yi512 = (rho_i3/(2*pi*e_0)).*((yf - posy(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_yi513 = (rho_i4/(2*pi*e_0)).*((yf - posy(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_yi514 = (rho_i5/(2*pi*e_0)).*((yf - posy(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_yi515 = (rho_i6/(2*pi*e_0)).*((yf - posy(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_yi516 = (rho_i7/(2*pi*e_0)).*((yf - posy(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_yi517 = (rho_i8/(2*pi*e_0)).*((yf - posy(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_yi518 = (rho_i9/(2*pi*e_0)).*((yf - posy(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));

E_yi61 = (-rho_i1/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yi62 = (-rho_i2/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yi63 = (-rho_i3/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yi64 = (-rho_i4/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yi65 = (-rho_i5/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_yi67 = (-rho_i7/(2*pi*e_0)).*((yf - posy(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_yi68 = (-rho_i8/(2*pi*e_0)).*((yf - posy(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_yi69 = (-rho_i9/(2*pi*e_0)).*((yf - posy(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_yi610 = (rho_i1/(2*pi*e_0)).*((yf - posy(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_yi611 = (rho_i2/(2*pi*e_0)).*((yf - posy(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_yi612 = (rho_i3/(2*pi*e_0)).*((yf - posy(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_yi613 = (rho_i4/(2*pi*e_0)).*((yf - posy(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_yi614 = (rho_i5/(2*pi*e_0)).*((yf - posy(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_yi615 = (rho_i6/(2*pi*e_0)).*((yf - posy(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_yi616 = (rho_i7/(2*pi*e_0)).*((yf - posy(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_yi617 = (rho_i8/(2*pi*e_0)).*((yf - posy(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_yi618 = (rho_i9/(2*pi*e_0)).*((yf - posy(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));

E_yi71 = (-rho_i1/(2*pi*e_0)).*((yf - posy(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_yi72 = (-rho_i2/(2*pi*e_0)).*((yf - posy(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_yi73 = (-rho_i3/(2*pi*e_0)).*((yf - posy(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_yi74 = (-rho_i4/(2*pi*e_0)).*((yf - posy(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_yi75 = (-rho_i5/(2*pi*e_0)).*((yf - posy(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_yi76 = (-rho_i6/(2*pi*e_0)).*((yf - posy(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_yi78 = (-rho_i8/(2*pi*e_0)).*((yf - posy(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_yi79 = (-rho_i9/(2*pi*e_0)).*((yf - posy(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_yi710 = (rho_i1/(2*pi*e_0)).*((yf - posy(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_yi711 = (rho_i2/(2*pi*e_0)).*((yf - posy(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_yi712 = (rho_i3/(2*pi*e_0)).*((yf - posy(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_yi713 = (rho_i4/(2*pi*e_0)).*((yf - posy(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_yi714 = (rho_i5/(2*pi*e_0)).*((yf - posy(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_yi715 = (rho_i6/(2*pi*e_0)).*((yf - posy(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_yi716 = (rho_i7/(2*pi*e_0)).*((yf - posy(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_yi717 = (rho_i8/(2*pi*e_0)).*((yf - posy(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_yi718 = (rho_i9/(2*pi*e_0)).*((yf - posy(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));

E_yi81 = (-rho_i1/(2*pi*e_0)).*((yf - posy(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_yi82 = (-rho_i2/(2*pi*e_0)).*((yf - posy(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_yi83 = (-rho_i3/(2*pi*e_0)).*((yf - posy(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_yi84 = (-rho_i4/(2*pi*e_0)).*((yf - posy(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_yi85 = (-rho_i5/(2*pi*e_0)).*((yf - posy(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_yi86 = (-rho_i6/(2*pi*e_0)).*((yf - posy(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_yi87 = (-rho_i7/(2*pi*e_0)).*((yf - posy(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_yi89 = (-rho_i9/(2*pi*e_0)).*((yf - posy(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_yi810 = (rho_i1/(2*pi*e_0)).*((yf - posy(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_yi811 = (rho_i2/(2*pi*e_0)).*((yf - posy(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_yi812 = (rho_i3/(2*pi*e_0)).*((yf - posy(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_yi813 = (rho_i4/(2*pi*e_0)).*((yf - posy(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_yi814 = (rho_i5/(2*pi*e_0)).*((yf - posy(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_yi815 = (rho_i6/(2*pi*e_0)).*((yf - posy(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_yi816 = (rho_i7/(2*pi*e_0)).*((yf - posy(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_yi817 = (rho_i8/(2*pi*e_0)).*((yf - posy(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_yi818 = (rho_i9/(2*pi*e_0)).*((yf - posy(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));

E_yi91 = (-rho_i1/(2*pi*e_0)).*((yf - posy(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_yi92 = (-rho_i2/(2*pi*e_0)).*((yf - posy(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_yi93 = (-rho_i3/(2*pi*e_0)).*((yf - posy(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_yi94 = (-rho_i4/(2*pi*e_0)).*((yf - posy(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_yi95 = (-rho_i5/(2*pi*e_0)).*((yf - posy(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_yi96 = (-rho_i6/(2*pi*e_0)).*((yf - posy(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_yi97 = (-rho_i7/(2*pi*e_0)).*((yf - posy(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_yi98 = (-rho_i8/(2*pi*e_0)).*((yf - posy(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_yi910 = (rho_i1/(2*pi*e_0)).*((yf - posy(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_yi911 = (rho_i2/(2*pi*e_0)).*((yf - posy(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_yi912 = (rho_i3/(2*pi*e_0)).*((yf - posy(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_yi913 = (rho_i4/(2*pi*e_0)).*((yf - posy(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_yi914 = (rho_i5/(2*pi*e_0)).*((yf - posy(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_yi915 = (rho_i6/(2*pi*e_0)).*((yf - posy(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_yi916 = (rho_i7/(2*pi*e_0)).*((yf - posy(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_yi917 = (rho_i8/(2*pi*e_0)).*((yf - posy(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_yi918 = (rho_i9/(2*pi*e_0)).*((yf - posy(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));

E_yi101 = (-rho_i1/(2*pi*e_0)).*((yf - posy(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_yi102 = (-rho_i2/(2*pi*e_0)).*((yf - posy(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_yi103 = (-rho_i3/(2*pi*e_0)).*((yf - posy(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_yi104 = (-rho_i4/(2*pi*e_0)).*((yf - posy(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_yi105 = (-rho_i5/(2*pi*e_0)).*((yf - posy(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_yi106 = (-rho_i6/(2*pi*e_0)).*((yf - posy(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_yi107 = (-rho_i7/(2*pi*e_0)).*((yf - posy(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_yi108 = (-rho_i8/(2*pi*e_0)).*((yf - posy(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_yi109 = (-rho_i9/(2*pi*e_0)).*((yf - posy(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_yi1011 = (rho_i2/(2*pi*e_0)).*((yf - posy(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_yi1012 = (rho_i3/(2*pi*e_0)).*((yf - posy(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_yi1013 = (rho_i4/(2*pi*e_0)).*((yf - posy(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_yi1014 = (rho_i5/(2*pi*e_0)).*((yf - posy(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_yi1015 = (rho_i6/(2*pi*e_0)).*((yf - posy(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_yi1016 = (rho_i7/(2*pi*e_0)).*((yf - posy(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_yi1017 = (rho_i8/(2*pi*e_0)).*((yf - posy(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_yi1018 = (rho_i9/(2*pi*e_0)).*((yf - posy(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));

E_yi11_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_yi11_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_yi113 = (-rho_i3/(2*pi*e_0)).*((yf - posy(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_yi114 = (-rho_i4/(2*pi*e_0)).*((yf - posy(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_yi115 = (-rho_i5/(2*pi*e_0)).*((yf - posy(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_yi116 = (-rho_i6/(2*pi*e_0)).*((yf - posy(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_yi117 = (-rho_i7/(2*pi*e_0)).*((yf - posy(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_yi118 = (-rho_i8/(2*pi*e_0)).*((yf - posy(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_yi119 = (-rho_i9/(2*pi*e_0)).*((yf - posy(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_yi1110 = (rho_i1/(2*pi*e_0)).*((yf - posy(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_yi1112 = (rho_i3/(2*pi*e_0)).*((yf - posy(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_yi1113 = (rho_i4/(2*pi*e_0)).*((yf - posy(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_yi1114 = (rho_i5/(2*pi*e_0)).*((yf - posy(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_yi1115 = (rho_i6/(2*pi*e_0)).*((yf - posy(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_yi1116 = (rho_i7/(2*pi*e_0)).*((yf - posy(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_yi1117 = (rho_i8/(2*pi*e_0)).*((yf - posy(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_yi1118 = (rho_i9/(2*pi*e_0)).*((yf - posy(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));

E_yi12_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_yi12_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_yi123 = (-rho_i3/(2*pi*e_0)).*((yf - posy(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_yi124 = (-rho_i4/(2*pi*e_0)).*((yf - posy(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_yi125 = (-rho_i5/(2*pi*e_0)).*((yf - posy(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_yi126 = (-rho_i6/(2*pi*e_0)).*((yf - posy(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_yi127 = (-rho_i7/(2*pi*e_0)).*((yf - posy(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_yi128 = (-rho_i8/(2*pi*e_0)).*((yf - posy(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_yi129 = (-rho_i9/(2*pi*e_0)).*((yf - posy(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_yi1210 = (rho_i1/(2*pi*e_0)).*((yf - posy(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_yi1211 = (rho_i2/(2*pi*e_0)).*((yf - posy(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_yi1213 = (rho_i4/(2*pi*e_0)).*((yf - posy(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_yi1214 = (rho_i5/(2*pi*e_0)).*((yf - posy(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_yi1215 = (rho_i6/(2*pi*e_0)).*((yf - posy(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_yi1216 = (rho_i7/(2*pi*e_0)).*((yf - posy(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_yi1217 = (rho_i8/(2*pi*e_0)).*((yf - posy(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_yi1218 = (rho_i9/(2*pi*e_0)).*((yf - posy(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));

E_yi131 = (-rho_i1/(2*pi*e_0)).*((yf - posy(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_yi132 = (-rho_i2/(2*pi*e_0)).*((yf - posy(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_yi133 = (-rho_i3/(2*pi*e_0)).*((yf - posy(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_yi134 = (-rho_i4/(2*pi*e_0)).*((yf - posy(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_yi135 = (-rho_i5/(2*pi*e_0)).*((yf - posy(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_yi136 = (-rho_i6/(2*pi*e_0)).*((yf - posy(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_yi137 = (-rho_i7/(2*pi*e_0)).*((yf - posy(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_yi138 = (-rho_i8/(2*pi*e_0)).*((yf - posy(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_yi139 = (-rho_i9/(2*pi*e_0)).*((yf - posy(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_yi1310 = (rho_i1/(2*pi*e_0)).*((yf - posy(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_yi1311 = (rho_i2/(2*pi*e_0)).*((yf - posy(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_yi1312 = (rho_i3/(2*pi*e_0)).*((yf - posy(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_yi1314 = (rho_i5/(2*pi*e_0)).*((yf - posy(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_yi1315 = (rho_i6/(2*pi*e_0)).*((yf - posy(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_yi1316 = (rho_i7/(2*pi*e_0)).*((yf - posy(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_yi1317 = (rho_i8/(2*pi*e_0)).*((yf - posy(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_yi1318 = (rho_i9/(2*pi*e_0)).*((yf - posy(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));

E_yi141 = (-rho_i1/(2*pi*e_0)).*((yf - posy(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_yi142 = (-rho_i2/(2*pi*e_0)).*((yf - posy(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_yi143 = (-rho_i3/(2*pi*e_0)).*((yf - posy(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_yi144 = (-rho_i4/(2*pi*e_0)).*((yf - posy(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_yi145 = (-rho_i5/(2*pi*e_0)).*((yf - posy(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_yi146 = (-rho_i6/(2*pi*e_0)).*((yf - posy(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_yi147 = (-rho_i7/(2*pi*e_0)).*((yf - posy(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_yi148 = (-rho_i8/(2*pi*e_0)).*((yf - posy(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_yi149 = (-rho_i9/(2*pi*e_0)).*((yf - posy(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_yi1410 = (rho_i1/(2*pi*e_0)).*((yf - posy(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_yi1411 = (rho_i2/(2*pi*e_0)).*((yf - posy(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_yi1412 = (rho_i3/(2*pi*e_0)).*((yf - posy(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_yi1413 = (rho_i4/(2*pi*e_0)).*((yf - posy(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_yi1415 = (rho_i6/(2*pi*e_0)).*((yf - posy(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_yi1416 = (rho_i7/(2*pi*e_0)).*((yf - posy(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_yi1417 = (rho_i8/(2*pi*e_0)).*((yf - posy(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_yi1418 = (rho_i9/(2*pi*e_0)).*((yf - posy(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));

E_yi151 = (-rho_i1/(2*pi*e_0)).*((yf - posy(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_yi152 = (-rho_i2/(2*pi*e_0)).*((yf - posy(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_yi153 = (-rho_i3/(2*pi*e_0)).*((yf - posy(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_yi154 = (-rho_i4/(2*pi*e_0)).*((yf - posy(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_yi155 = (-rho_i5/(2*pi*e_0)).*((yf - posy(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_yi156 = (-rho_i6/(2*pi*e_0)).*((yf - posy(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_yi157 = (-rho_i7/(2*pi*e_0)).*((yf - posy(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_yi158 = (-rho_i8/(2*pi*e_0)).*((yf - posy(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_yi159 = (-rho_i9/(2*pi*e_0)).*((yf - posy(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_yi1510 = (rho_i1/(2*pi*e_0)).*((yf - posy(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_yi1511 = (rho_i2/(2*pi*e_0)).*((yf - posy(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_yi1512 = (rho_i3/(2*pi*e_0)).*((yf - posy(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_yi1513 = (rho_i4/(2*pi*e_0)).*((yf - posy(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_yi1514 = (rho_i5/(2*pi*e_0)).*((yf - posy(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_yi1516 = (rho_i7/(2*pi*e_0)).*((yf - posy(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_yi1517 = (rho_i8/(2*pi*e_0)).*((yf - posy(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_yi1518 = (rho_i9/(2*pi*e_0)).*((yf - posy(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));

E_yi161 = (-rho_i1/(2*pi*e_0)).*((yf - posy(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_yi162 = (-rho_i2/(2*pi*e_0)).*((yf - posy(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_yi163 = (-rho_i3/(2*pi*e_0)).*((yf - posy(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_yi164 = (-rho_i4/(2*pi*e_0)).*((yf - posy(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_yi165 = (-rho_i5/(2*pi*e_0)).*((yf - posy(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_yi166 = (-rho_i6/(2*pi*e_0)).*((yf - posy(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_yi167 = (-rho_i7/(2*pi*e_0)).*((yf - posy(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_yi168 = (-rho_i8/(2*pi*e_0)).*((yf - posy(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_yi169 = (-rho_i9/(2*pi*e_0)).*((yf - posy(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_yi1610 = (rho_i1/(2*pi*e_0)).*((yf - posy(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_yi1611 = (rho_i2/(2*pi*e_0)).*((yf - posy(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_yi1612 = (rho_i3/(2*pi*e_0)).*((yf - posy(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_yi1613 = (rho_i4/(2*pi*e_0)).*((yf - posy(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_yi1614 = (rho_i5/(2*pi*e_0)).*((yf - posy(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_yi1615 = (rho_i6/(2*pi*e_0)).*((yf - posy(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_yi1617 = (rho_i8/(2*pi*e_0)).*((yf - posy(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_yi1618 = (rho_i9/(2*pi*e_0)).*((yf - posy(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));

E_yi171 = (-rho_i1/(2*pi*e_0)).*((yf - posy(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_yi172 = (-rho_i2/(2*pi*e_0)).*((yf - posy(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_yi173 = (-rho_i3/(2*pi*e_0)).*((yf - posy(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_yi174 = (-rho_i4/(2*pi*e_0)).*((yf - posy(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_yi175 = (-rho_i5/(2*pi*e_0)).*((yf - posy(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_yi176 = (-rho_i6/(2*pi*e_0)).*((yf - posy(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_yi177 = (-rho_i7/(2*pi*e_0)).*((yf - posy(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_yi178 = (-rho_i8/(2*pi*e_0)).*((yf - posy(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_yi179 = (-rho_i9/(2*pi*e_0)).*((yf - posy(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_yi1710 = (rho_i1/(2*pi*e_0)).*((yf - posy(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_yi1711 = (rho_i2/(2*pi*e_0)).*((yf - posy(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_yi1712 = (rho_i3/(2*pi*e_0)).*((yf - posy(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_yi1713 = (rho_i4/(2*pi*e_0)).*((yf - posy(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_yi1714 = (rho_i5/(2*pi*e_0)).*((yf - posy(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_yi1715 = (rho_i6/(2*pi*e_0)).*((yf - posy(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_yi1716 = (rho_i7/(2*pi*e_0)).*((yf - posy(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_yi1718 = (rho_i9/(2*pi*e_0)).*((yf - posy(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));

E_yi181 = (-rho_i1/(2*pi*e_0)).*((yf - posy(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_yi182 = (-rho_i2/(2*pi*e_0)).*((yf - posy(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_yi183 = (-rho_i3/(2*pi*e_0)).*((yf - posy(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_yi184 = (-rho_i4/(2*pi*e_0)).*((yf - posy(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_yi185 = (-rho_i5/(2*pi*e_0)).*((yf - posy(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_yi186 = (-rho_i6/(2*pi*e_0)).*((yf - posy(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_yi187 = (-rho_i7/(2*pi*e_0)).*((yf - posy(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_yi188 = (-rho_i8/(2*pi*e_0)).*((yf - posy(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_yi189 = (-rho_i9/(2*pi*e_0)).*((yf - posy(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_yi1810 = (rho_i1/(2*pi*e_0)).*((yf - posy(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_yi1811 = (rho_i2/(2*pi*e_0)).*((yf - posy(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_yi1812 = (rho_i3/(2*pi*e_0)).*((yf - posy(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_yi1813 = (rho_i4/(2*pi*e_0)).*((yf - posy(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_yi1814 = (rho_i5/(2*pi*e_0)).*((yf - posy(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_yi1815 = (rho_i6/(2*pi*e_0)).*((yf - posy(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_yi1816 = (rho_i7/(2*pi*e_0)).*((yf - posy(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_yi1817 = (rho_i8/(2*pi*e_0)).*((yf - posy(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));

%% aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr12 ;	E_xr13 	;	E_xr14 	;	E_xr15 	;	E_xr16 	;	E_xr17 	;	E_xr18 	;	E_xr19 	;	E_xr110 	;	E_xr1_11 	;	E_xr1_12 	;	E_xr1_13 	;	E_xr1_14 	;	E_xr1_15 	;	E_xr1_16 	;	E_xr1_17 	;	E_xr1_18 	;	E_xr21 	;	E_xr23 	;	E_xr24 	;	E_xr25 	;	E_xr26 	;	E_xr27 	;	E_xr28 	;	E_xr29 	;	E_xr210 	;	E_xr2_11 	;	E_xr2_12 	;	E_xr2_13 	;	E_xr2_14 	;	E_xr2_15 	;	E_xr2_16 	;	E_xr2_17 	;	E_xr2_18 	;	E_xr31 	;	E_xr32 	;	E_xr34 	;	E_xr35 	;	E_xr36 	;	E_xr37 	;	E_xr38 	;	E_xr39 	;	E_xr310 	;	E_xr311 	;	E_xr312 	;	E_xr313 	;	E_xr314 	;	E_xr315 	;	E_xr316 	;	E_xr317 	;	E_xr318 	;	E_xr41 	;	E_xr42 	;	E_xr43 	;	E_xr45 	;	E_xr46 	;	E_xr47 	;	E_xr48 	;	E_xr49 	;	E_xr410 	;	E_xr411 	;	E_xr412 	;	E_xr413 	;	E_xr414 	;	E_xr415 	;	E_xr416 	;	E_xr417 	;	E_xr418 	;	E_xr51 	;	E_xr52 	;	E_xr53 	;	E_xr54 	;	E_xr56 	;	E_xr57 	;	E_xr58 	;	E_xr59 	;	E_xr510 	;	E_xr511 	;	E_xr512 	;	E_xr513 	;	E_xr514 	;	E_xr515 	;	E_xr516 	;	E_xr517 	;	E_xr518 	;	E_xr61 	;	E_xr62 	;	E_xr63 	;	E_xr64 	;	E_xr65 	;	E_xr67 	;	E_xr68 	;	E_xr69 	;	E_xr610 	;	E_xr611 	;	E_xr612 	;	E_xr613 	;	E_xr614 	;	E_xr615 	;	E_xr616 	;	E_xr617 	;	E_xr618 	;	E_xr71 	;	E_xr72 	;	E_xr73 	;	E_xr74 	;	E_xr75 	;	E_xr76 	;	E_xr78 	;	E_xr79 	;	E_xr710 	;	E_xr711 	;	E_xr712 	;	E_xr713 	;	E_xr714 	;	E_xr715 	;	E_xr716 	;	E_xr717 	;	E_xr718 	;	E_xr81 	;	E_xr82 	;	E_xr83 	;	E_xr84 	;	E_xr85 	;	E_xr86 	;	E_xr87 	;	E_xr89 	;	E_xr810 	;	E_xr811 	;	E_xr812 	;	E_xr813 	;	E_xr814 	;	E_xr815 	;	E_xr816 	;	E_xr817 	;	E_xr818 	;	E_xr91 	;	E_xr92 	;	E_xr93 	;	E_xr94 	;	E_xr95 	;	E_xr96 	;	E_xr97 	;	E_xr98 	;	E_xr910 	;	E_xr911 	;	E_xr912 	;	E_xr913 	;	E_xr914 	;	E_xr915 	;	E_xr916 	;	E_xr917 	;	E_xr918 	;	E_xr101 	;	E_xr102 	;	E_xr103 	;	E_xr104 	;	E_xr105 	;	E_xr106 	;	E_xr107 	;	E_xr108 	;	E_xr109 	;	E_xr1011 	;	E_xr1012 	;	E_xr1013 	;	E_xr1014 	;	E_xr1015 	;	E_xr1016 	;	E_xr1017 	;	E_xr1018 	;	E_xr11_1 	;	E_xr11_2 	;	E_xr113 	;	E_xr114 	;	E_xr115 	;	E_xr116 	;	E_xr117 	;	E_xr118 	;	E_xr119 	;	E_xr1110 	;	E_xr1112 	;	E_xr1113 	;	E_xr1114 	;	E_xr1115 	;	E_xr1116 	;	E_xr1117 	;	E_xr1118 	;	E_xr12_1 	;	E_xr12_2 	;	E_xr123 	;	E_xr124 	;	E_xr125 	;	E_xr126 	;	E_xr127 	;	E_xr128 	;	E_xr129 	;	E_xr1210 	;	E_xr1211 	;	E_xr1213 	;	E_xr1214 	;	E_xr1215 	;	E_xr1216 	;	E_xr1217 	;	E_xr1218 	;	E_xr131 	;	E_xr132 	;	E_xr133 	;	E_xr134 	;	E_xr135 	;	E_xr136 	;	E_xr137 	;	E_xr138 	;	E_xr139 	;	E_xr1310 	;	E_xr1311 	;	E_xr1312 	;	E_xr1314 	;	E_xr1315 	;	E_xr1316 	;	E_xr1317 	;	E_xr1318 	;	E_xr141 	;	E_xr142 	;	E_xr143 	;	E_xr144 	;	E_xr145 	;	E_xr146 	;	E_xr147 	;	E_xr148 	;	E_xr149 	;	E_xr1410 	;	E_xr1411 	;	E_xr1412 	;	E_xr1413 	;	E_xr1415 	;	E_xr1416 	;	E_xr1417 	;	E_xr1418 	;	E_xr151 	;	E_xr152 	;	E_xr153 	;	E_xr154 	;	E_xr155 	;	E_xr156 	;	E_xr157 	;	E_xr158 	;	E_xr159 	;	E_xr1510 	;	E_xr1511 	;	E_xr1512 	;	E_xr1513 	;	E_xr1514 	;	E_xr1516 	;	E_xr1517 	;	E_xr1518 	;	E_xr161 	;	E_xr162 	;	E_xr163 	;	E_xr164 	;	E_xr165 	;	E_xr166 	;	E_xr167 	;	E_xr168 	;	E_xr169 	;	E_xr1610 	;	E_xr1611 	;	E_xr1612 	;	E_xr1613 	;	E_xr1614 	;	E_xr1615 	;	E_xr1617 	;	E_xr1618 	;	E_xr171 	;	E_xr172 	;	E_xr173 	;	E_xr174 	;	E_xr175 	;	E_xr176 	;	E_xr177 	;	E_xr178 	;	E_xr179 	;	E_xr1710 	;	E_xr1711 	;	E_xr1712 	;	E_xr1713 	;	E_xr1714 	;	E_xr1715 	;	E_xr1716 	;	E_xr1718 	;	E_xr181 	;	E_xr182 	;	E_xr183 	;	E_xr184 	;	E_xr185 	;	E_xr186 	;	E_xr187 	;	E_xr188 	;	E_xr189 	;	E_xr1810 	;	E_xr1811 	;	E_xr1812 	;	E_xr1813 	;	E_xr1814 	;	E_xr1815 	;	E_xr1816 	;	E_xr1817 ];

E_xi = [ E_xi12 ;	E_xi13 	;	E_xi14 	;	E_xi15 	;	E_xi16 	;	E_xi17 	;	E_xi18 	;	E_xi19 	;	E_xi110 	;	E_xi1_11 	;	E_xi1_12 	;	E_xi1_13 	;	E_xi1_14 	;	E_xi1_15 	;	E_xi1_16 	;	E_xi1_17 	;	E_xi1_18 	;	E_xi21 	;	E_xi23 	;	E_xi24 	;	E_xi25 	;	E_xi26 	;	E_xi27 	;	E_xi28 	;	E_xi29 	;	E_xi210 	;	E_xi2_11 	;	E_xi2_12 	;	E_xi2_13 	;	E_xi2_14 	;	E_xi2_15 	;	E_xi2_16 	;	E_xi2_17 	;	E_xi2_18 	;	E_xi31 	;	E_xi32 	;	E_xi34 	;	E_xi35 	;	E_xi36 	;	E_xi37 	;	E_xi38 	;	E_xi39 	;	E_xi310 	;	E_xi311 	;	E_xi312 	;	E_xi313 	;	E_xi314 	;	E_xi315 	;	E_xi316 	;	E_xi317 	;	E_xi318 	;	E_xi41 	;	E_xi42 	;	E_xi43 	;	E_xi45 	;	E_xi46 	;	E_xi47 	;	E_xi48 	;	E_xi49 	;	E_xi410 	;	E_xi411 	;	E_xi412 	;	E_xi413 	;	E_xi414 	;	E_xi415 	;	E_xi416 	;	E_xi417 	;	E_xi418 	;	E_xi51 	;	E_xi52 	;	E_xi53 	;	E_xi54 	;	E_xi56 	;	E_xi57 	;	E_xi58 	;	E_xi59 	;	E_xi510 	;	E_xi511 	;	E_xi512 	;	E_xi513 	;	E_xi514 	;	E_xi515 	;	E_xi516 	;	E_xi517 	;	E_xi518 	;	E_xi61 	;	E_xi62 	;	E_xi63 	;	E_xi64 	;	E_xi65 	;	E_xi67 	;	E_xi68 	;	E_xi69 	;	E_xi610 	;	E_xi611 	;	E_xi612 	;	E_xi613 	;	E_xi614 	;	E_xi615 	;	E_xi616 	;	E_xi617 	;	E_xi618 	;	E_xi71 	;	E_xi72 	;	E_xi73 	;	E_xi74 	;	E_xi75 	;	E_xi76 	;	E_xi78 	;	E_xi79 	;	E_xi710 	;	E_xi711 	;	E_xi712 	;	E_xi713 	;	E_xi714 	;	E_xi715 	;	E_xi716 	;	E_xi717 	;	E_xi718 	;	E_xi81 	;	E_xi82 	;	E_xi83 	;	E_xi84 	;	E_xi85 	;	E_xi86 	;	E_xi87 	;	E_xi89 	;	E_xi810 	;	E_xi811 	;	E_xi812 	;	E_xi813 	;	E_xi814 	;	E_xi815 	;	E_xi816 	;	E_xi817 	;	E_xi818 	;	E_xi91 	;	E_xi92 	;	E_xi93 	;	E_xi94 	;	E_xi95 	;	E_xi96 	;	E_xi97 	;	E_xi98 	;	E_xi910 	;	E_xi911 	;	E_xi912 	;	E_xi913 	;	E_xi914 	;	E_xi915 	;	E_xi916 	;	E_xi917 	;	E_xi918 	;	E_xi101 	;	E_xi102 	;	E_xi103 	;	E_xi104 	;	E_xi105 	;	E_xi106 	;	E_xi107 	;	E_xi108 	;	E_xi109 	;	E_xi1011 	;	E_xi1012 	;	E_xi1013 	;	E_xi1014 	;	E_xi1015 	;	E_xi1016 	;	E_xi1017 	;	E_xi1018 	;	E_xi11_1 	;	E_xi11_2 	;	E_xi113 	;	E_xi114 	;	E_xi115 	;	E_xi116 	;	E_xi117 	;	E_xi118 	;	E_xi119 	;	E_xi1110 	;	E_xi1112 	;	E_xi1113 	;	E_xi1114 	;	E_xi1115 	;	E_xi1116 	;	E_xi1117 	;	E_xi1118 	;	E_xi12_1 	;	E_xi12_2 	;	E_xi123 	;	E_xi124 	;	E_xi125 	;	E_xi126 	;	E_xi127 	;	E_xi128 	;	E_xi129 	;	E_xi1210 	;	E_xi1211 	;	E_xi1213 	;	E_xi1214 	;	E_xi1215 	;	E_xi1216 	;	E_xi1217 	;	E_xi1218 	;	E_xi131 	;	E_xi132 	;	E_xi133 	;	E_xi134 	;	E_xi135 	;	E_xi136 	;	E_xi137 	;	E_xi138 	;	E_xi139 	;	E_xi1310 	;	E_xi1311 	;	E_xi1312 	;	E_xi1314 	;	E_xi1315 	;	E_xi1316 	;	E_xi1317 	;	E_xi1318 	;	E_xi141 	;	E_xi142 	;	E_xi143 	;	E_xi144 	;	E_xi145 	;	E_xi146 	;	E_xi147 	;	E_xi148 	;	E_xi149 	;	E_xi1410 	;	E_xi1411 	;	E_xi1412 	;	E_xi1413 	;	E_xi1415 	;	E_xi1416 	;	E_xi1417 	;	E_xi1418 	;	E_xi151 	;	E_xi152 	;	E_xi153 	;	E_xi154 	;	E_xi155 	;	E_xi156 	;	E_xi157 	;	E_xi158 	;	E_xi159 	;	E_xi1510 	;	E_xi1511 	;	E_xi1512 	;	E_xi1513 	;	E_xi1514 	;	E_xi1516 	;	E_xi1517 	;	E_xi1518 	;	E_xi161 	;	E_xi162 	;	E_xi163 	;	E_xi164 	;	E_xi165 	;	E_xi166 	;	E_xi167 	;	E_xi168 	;	E_xi169 	;	E_xi1610 	;	E_xi1611 	;	E_xi1612 	;	E_xi1613 	;	E_xi1614 	;	E_xi1615 	;	E_xi1617 	;	E_xi1618 	;	E_xi171 	;	E_xi172 	;	E_xi173 	;	E_xi174 	;	E_xi175 	;	E_xi176 	;	E_xi177 	;	E_xi178 	;	E_xi179 	;	E_xi1710 	;	E_xi1711 	;	E_xi1712 	;	E_xi1713 	;	E_xi1714 	;	E_xi1715 	;	E_xi1716 	;	E_xi1718 	;	E_xi181 	;	E_xi182 	;	E_xi183 	;	E_xi184 	;	E_xi185 	;	E_xi186 	;	E_xi187 	;	E_xi188 	;	E_xi189 	;	E_xi1810 	;	E_xi1811 	;	E_xi1812 	;	E_xi1813 	;	E_xi1814 	;	E_xi1815 	;	E_xi1816 	;	E_xi1817 ]; 

E_yr = [ E_yr12 ;	E_yr13 	;	E_yr14 	;	E_yr15 	;	E_yr16 	;	E_yr17 	;	E_yr18 	;	E_yr19 	;	E_yr110 	;	E_yr1_11 	;	E_yr1_12 	;	E_yr1_13 	;	E_yr1_14 	;	E_yr1_15 	;	E_yr1_16 	;	E_yr1_17 	;	E_yr1_18 	;	E_yr21 	;	E_yr23 	;	E_yr24 	;	E_yr25 	;	E_yr26 	;	E_yr27 	;	E_yr28 	;	E_yr29 	;	E_yr210 	;	E_yr2_11 	;	E_yr2_12 	;	E_yr2_13 	;	E_yr2_14 	;	E_yr2_15 	;	E_yr2_16 	;	E_yr2_17 	;	E_yr2_18 	;	E_yr31 	;	E_yr32 	;	E_yr34 	;	E_yr35 	;	E_yr36 	;	E_yr37 	;	E_yr38 	;	E_yr39 	;	E_yr310 	;	E_yr311 	;	E_yr312 	;	E_yr313 	;	E_yr314 	;	E_yr315 	;	E_yr316 	;	E_yr317 	;	E_yr318 	;	E_yr41 	;	E_yr42 	;	E_yr43 	;	E_yr45 	;	E_yr46 	;	E_yr47 	;	E_yr48 	;	E_yr49 	;	E_yr410 	;	E_yr411 	;	E_yr412 	;	E_yr413 	;	E_yr414 	;	E_yr415 	;	E_yr416 	;	E_yr417 	;	E_yr418 	;	E_yr51 	;	E_yr52 	;	E_yr53 	;	E_yr54 	;	E_yr56 	;	E_yr57 	;	E_yr58 	;	E_yr59 	;	E_yr510 	;	E_yr511 	;	E_yr512 	;	E_yr513 	;	E_yr514 	;	E_yr515 	;	E_yr516 	;	E_yr517 	;	E_yr518 	;	E_yr61 	;	E_yr62 	;	E_yr63 	;	E_yr64 	;	E_yr65 	;	E_yr67 	;	E_yr68 	;	E_yr69 	;	E_yr610 	;	E_yr611 	;	E_yr612 	;	E_yr613 	;	E_yr614 	;	E_yr615 	;	E_yr616 	;	E_yr617 	;	E_yr618 	;	E_yr71 	;	E_yr72 	;	E_yr73 	;	E_yr74 	;	E_yr75 	;	E_yr76 	;	E_yr78 	;	E_yr79 	;	E_yr710 	;	E_yr711 	;	E_yr712 	;	E_yr713 	;	E_yr714 	;	E_yr715 	;	E_yr716 	;	E_yr717 	;	E_yr718 	;	E_yr81 	;	E_yr82 	;	E_yr83 	;	E_yr84 	;	E_yr85 	;	E_yr86 	;	E_yr87 	;	E_yr89 	;	E_yr810 	;	E_yr811 	;	E_yr812 	;	E_yr813 	;	E_yr814 	;	E_yr815 	;	E_yr816 	;	E_yr817 	;	E_yr818 	;	E_yr91 	;	E_yr92 	;	E_yr93 	;	E_yr94 	;	E_yr95 	;	E_yr96 	;	E_yr97 	;	E_yr98 	;	E_yr910 	;	E_yr911 	;	E_yr912 	;	E_yr913 	;	E_yr914 	;	E_yr915 	;	E_yr916 	;	E_yr917 	;	E_yr918 	;	E_yr101 	;	E_yr102 	;	E_yr103 	;	E_yr104 	;	E_yr105 	;	E_yr106 	;	E_yr107 	;	E_yr108 	;	E_yr109 	;	E_yr1011 	;	E_yr1012 	;	E_yr1013 	;	E_yr1014 	;	E_yr1015 	;	E_yr1016 	;	E_yr1017 	;	E_yr1018 	;	E_yr11_1 	;	E_yr11_2 	;	E_yr113 	;	E_yr114 	;	E_yr115 	;	E_yr116 	;	E_yr117 	;	E_yr118 	;	E_yr119 	;	E_yr1110 	;	E_yr1112 	;	E_yr1113 	;	E_yr1114 	;	E_yr1115 	;	E_yr1116 	;	E_yr1117 	;	E_yr1118 	;	E_yr12_1 	;	E_yr12_2 	;	E_yr123 	;	E_yr124 	;	E_yr125 	;	E_yr126 	;	E_yr127 	;	E_yr128 	;	E_yr129 	;	E_yr1210 	;	E_yr1211 	;	E_yr1213 	;	E_yr1214 	;	E_yr1215 	;	E_yr1216 	;	E_yr1217 	;	E_yr1218 	;	E_yr131 	;	E_yr132 	;	E_yr133 	;	E_yr134 	;	E_yr135 	;	E_yr136 	;	E_yr137 	;	E_yr138 	;	E_yr139 	;	E_yr1310 	;	E_yr1311 	;	E_yr1312 	;	E_yr1314 	;	E_yr1315 	;	E_yr1316 	;	E_yr1317 	;	E_yr1318 	;	E_yr141 	;	E_yr142 	;	E_yr143 	;	E_yr144 	;	E_yr145 	;	E_yr146 	;	E_yr147 	;	E_yr148 	;	E_yr149 	;	E_yr1410 	;	E_yr1411 	;	E_yr1412 	;	E_yr1413 	;	E_yr1415 	;	E_yr1416 	;	E_yr1417 	;	E_yr1418 	;	E_yr151 	;	E_yr152 	;	E_yr153 	;	E_yr154 	;	E_yr155 	;	E_yr156 	;	E_yr157 	;	E_yr158 	;	E_yr159 	;	E_yr1510 	;	E_yr1511 	;	E_yr1512 	;	E_yr1513 	;	E_yr1514 	;	E_yr1516 	;	E_yr1517 	;	E_yr1518 	;	E_yr161 	;	E_yr162 	;	E_yr163 	;	E_yr164 	;	E_yr165 	;	E_yr166 	;	E_yr167 	;	E_yr168 	;	E_yr169 	;	E_yr1610 	;	E_yr1611 	;	E_yr1612 	;	E_yr1613 	;	E_yr1614 	;	E_yr1615 	;	E_yr1617 	;	E_yr1618 	;	E_yr171 	;	E_yr172 	;	E_yr173 	;	E_yr174 	;	E_yr175 	;	E_yr176 	;	E_yr177 	;	E_yr178 	;	E_yr179 	;	E_yr1710 	;	E_yr1711 	;	E_yr1712 	;	E_yr1713 	;	E_yr1714 	;	E_yr1715 	;	E_yr1716 	;	E_yr1718 	;	E_yr181 	;	E_yr182 	;	E_yr183 	;	E_yr184 	;	E_yr185 	;	E_yr186 	;	E_yr187 	;	E_yr188 	;	E_yr189 	;	E_yr1810 	;	E_yr1811 	;	E_yr1812 	;	E_yr1813 	;	E_yr1814 	;	E_yr1815 	;	E_yr1816 	;	E_yr1817 ];

E_yi = [ E_yi12 ;	E_yi13 	;	E_yi14 	;	E_yi15 	;	E_yi16 	;	E_yi17 	;	E_yi18 	;	E_yi19 	;	E_yi110 	;	E_yi1_11 	;	E_yi1_12 	;	E_yi1_13 	;	E_yi1_14 	;	E_yi1_15 	;	E_yi1_16 	;	E_yi1_17 	;	E_yi1_18 	;	E_yi21 	;	E_yi23 	;	E_yi24 	;	E_yi25 	;	E_yi26 	;	E_yi27 	;	E_yi28 	;	E_yi29 	;	E_yi210 	;	E_yi2_11 	;	E_yi2_12 	;	E_yi2_13 	;	E_yi2_14 	;	E_yi2_15 	;	E_yi2_16 	;	E_yi2_17 	;	E_yi2_18 	;	E_yi31 	;	E_yi32 	;	E_yi34 	;	E_yi35 	;	E_yi36 	;	E_yi37 	;	E_yi38 	;	E_yi39 	;	E_yi310 	;	E_yi311 	;	E_yi312 	;	E_yi313 	;	E_yi314 	;	E_yi315 	;	E_yi316 	;	E_yi317 	;	E_yi318 	;	E_yi41 	;	E_yi42 	;	E_yi43 	;	E_yi45 	;	E_yi46 	;	E_yi47 	;	E_yi48 	;	E_yi49 	;	E_yi410 	;	E_yi411 	;	E_yi412 	;	E_yi413 	;	E_yi414 	;	E_yi415 	;	E_yi416 	;	E_yi417 	;	E_yi418 	;	E_yi51 	;	E_yi52 	;	E_yi53 	;	E_yi54 	;	E_yi56 	;	E_yi57 	;	E_yi58 	;	E_yi59 	;	E_yi510 	;	E_yi511 	;	E_yi512 	;	E_yi513 	;	E_yi514 	;	E_yi515 	;	E_yi516 	;	E_yi517 	;	E_yi518 	;	E_yi61 	;	E_yi62 	;	E_yi63 	;	E_yi64 	;	E_yi65 	;	E_yi67 	;	E_yi68 	;	E_yi69 	;	E_yi610 	;	E_yi611 	;	E_yi612 	;	E_yi613 	;	E_yi614 	;	E_yi615 	;	E_yi616 	;	E_yi617 	;	E_yi618 	;	E_yi71 	;	E_yi72 	;	E_yi73 	;	E_yi74 	;	E_yi75 	;	E_yi76 	;	E_yi78 	;	E_yi79 	;	E_yi710 	;	E_yi711 	;	E_yi712 	;	E_yi713 	;	E_yi714 	;	E_yi715 	;	E_yi716 	;	E_yi717 	;	E_yi718 	;	E_yi81 	;	E_yi82 	;	E_yi83 	;	E_yi84 	;	E_yi85 	;	E_yi86 	;	E_yi87 	;	E_yi89 	;	E_yi810 	;	E_yi811 	;	E_yi812 	;	E_yi813 	;	E_yi814 	;	E_yi815 	;	E_yi816 	;	E_yi817 	;	E_yi818 	;	E_yi91 	;	E_yi92 	;	E_yi93 	;	E_yi94 	;	E_yi95 	;	E_yi96 	;	E_yi97 	;	E_yi98 	;	E_yi910 	;	E_yi911 	;	E_yi912 	;	E_yi913 	;	E_yi914 	;	E_yi915 	;	E_yi916 	;	E_yi917 	;	E_yi918 	;	E_yi101 	;	E_yi102 	;	E_yi103 	;	E_yi104 	;	E_yi105 	;	E_yi106 	;	E_yi107 	;	E_yi108 	;	E_yi109 	;	E_yi1011 	;	E_yi1012 	;	E_yi1013 	;	E_yi1014 	;	E_yi1015 	;	E_yi1016 	;	E_yi1017 	;	E_yi1018 	;	E_yi11_1 	;	E_yi11_2 	;	E_yi113 	;	E_yi114 	;	E_yi115 	;	E_yi116 	;	E_yi117 	;	E_yi118 	;	E_yi119 	;	E_yi1110 	;	E_yi1112 	;	E_yi1113 	;	E_yi1114 	;	E_yi1115 	;	E_yi1116 	;	E_yi1117 	;	E_yi1118 	;	E_yi12_1 	;	E_yi12_2 	;	E_yi123 	;	E_yi124 	;	E_yi125 	;	E_yi126 	;	E_yi127 	;	E_yi128 	;	E_yi129 	;	E_yi1210 	;	E_yi1211 	;	E_yi1213 	;	E_yi1214 	;	E_yi1215 	;	E_yi1216 	;	E_yi1217 	;	E_yi1218 	;	E_yi131 	;	E_yi132 	;	E_yi133 	;	E_yi134 	;	E_yi135 	;	E_yi136 	;	E_yi137 	;	E_yi138 	;	E_yi139 	;	E_yi1310 	;	E_yi1311 	;	E_yi1312 	;	E_yi1314 	;	E_yi1315 	;	E_yi1316 	;	E_yi1317 	;	E_yi1318 	;	E_yi141 	;	E_yi142 	;	E_yi143 	;	E_yi144 	;	E_yi145 	;	E_yi146 	;	E_yi147 	;	E_yi148 	;	E_yi149 	;	E_yi1410 	;	E_yi1411 	;	E_yi1412 	;	E_yi1413 	;	E_yi1415 	;	E_yi1416 	;	E_yi1417 	;	E_yi1418 	;	E_yi151 	;	E_yi152 	;	E_yi153 	;	E_yi154 	;	E_yi155 	;	E_yi156 	;	E_yi157 	;	E_yi158 	;	E_yi159 	;	E_yi1510 	;	E_yi1511 	;	E_yi1512 	;	E_yi1513 	;	E_yi1514 	;	E_yi1516 	;	E_yi1517 	;	E_yi1518 	;	E_yi161 	;	E_yi162 	;	E_yi163 	;	E_yi164 	;	E_yi165 	;	E_yi166 	;	E_yi167 	;	E_yi168 	;	E_yi169 	;	E_yi1610 	;	E_yi1611 	;	E_yi1612 	;	E_yi1613 	;	E_yi1614 	;	E_yi1615 	;	E_yi1617 	;	E_yi1618 	;	E_yi171 	;	E_yi172 	;	E_yi173 	;	E_yi174 	;	E_yi175 	;	E_yi176 	;	E_yi177 	;	E_yi178 	;	E_yi179 	;	E_yi1710 	;	E_yi1711 	;	E_yi1712 	;	E_yi1713 	;	E_yi1714 	;	E_yi1715 	;	E_yi1716 	;	E_yi1718 	;	E_yi181 	;	E_yi182 	;	E_yi183 	;	E_yi184 	;	E_yi185 	;	E_yi186 	;	E_yi187 	;	E_yi188 	;	E_yi189 	;	E_yi1810 	;	E_yi1811 	;	E_yi1812 	;	E_yi1813 	;	E_yi1814 	;	E_yi1815 	;	E_yi1816 	;	E_yi1817 ];

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

E3rms1 = [13.6557657488035 13.6388759914355 13.6226246253573 13.6070188653673 13.5920656825778 13.5777717990280 13.5641436824283 13.5511875410312 13.5389093186505 13.5273146898132 13.5164090550811 13.5061975365174 13.4966849733357 13.4878759177020 13.4797746307443 13.4723850787210 13.4657109294062 13.4597555486574 13.4545219971974 13.4500130276006 13.4462310814990 13.4431782870158 13.4408564564066 13.4392670839487 13.4384113440640 13.4382900896585 13.4389038507486 13.4402528332789 13.4423369182383 13.4451556609909 13.4487082908863 13.4529937111149 13.4580104988128 13.4637569054513 13.4702308574496 13.4774299570812 13.4853514836207 13.4939923947566 13.5033493282572 13.5134186039111 13.5241962256922 13.5356778842188 13.5478589594226 13.5607345235026 13.5742993440978 13.5885478877214 13.6034743234069 13.6190725266278 13.6353360833975 13.6522582946360 13.6698321807231 13.6880504862817 13.7069056851580 13.7263899856073 13.7464953356683 13.7672134287149 13.7885357092113 13.8104533786063 13.8329574014206 13.8560385114530 13.8796872181738 13.9038938132227 13.9286483770557 13.9539407857046 13.9797607176481 14.0060976607959 14.0329409195466 14.0602796219479 14.0881027269278 14.1163990315878 14.1451571785663 14.1743656634164 14.2040128420782 14.2340869383184 14.2645760512427 14.2954681627880 14.3267511452291 14.3584127686832 14.3904407085861 14.4228225531652 14.4555458108679 14.4885979177443 14.5219662448160 14.5556381053477 14.5896007620776 14.6238414343954 14.6583473054046 14.6931055289493 14.7281032364988 14.7633275439877 14.7987655585233 14.8344043849979 14.8702311325813 14.9062329211167 14.9423968873489 14.9787101910981 15.0151600212282 15.0517336015307 15.0884181964513 15.1252011166870 15.1620697246213 15.1990114396298 15.2360137432284 15.2730641840791 15.3101503828332 15.3472600368216 15.3843809245974 15.4215009103217 15.4586079479884 15.4956900854780 15.5327354684873 15.5697323442839 15.6066690652756 15.6435340924938 15.6803159988405 15.7170034722401 15.7535853186061 15.7900504646874 15.8263879607040 15.8625869829161 15.8986368359816 15.9345269551895 15.9702469085677 16.0057863988168 16.0411352651487 16.0762834849470 16.1112211753300 16.1459385945671 16.1804261433725 16.2146743660640 16.2486739516441 16.2824157346995 16.3158906962453 16.3490899644251 16.3820048151153 16.4146266724275 16.4469471091121 16.4789578468485 16.5106507564711 16.5420178580819 16.5730513210865 16.6037434641450 16.6340867550495 16.6640738105234 16.6936973959335 16.7229504249667 16.7518259592128 16.7803172076901 16.8084175263241 16.8361204173618 16.8634195287450 16.8903086534142 16.9167817285869 16.9428328349915 16.9684561960526 16.9936461770469 17.0183972842238 17.0427041639093 17.0665616015490 17.0899645207708 17.1129079823936 17.1353871834249 17.1573974560555 17.1789342666205 17.1999932145702 17.2205700314067 17.2406605796381 17.2602608517198 17.2793669689893 17.2979751806030 17.3160818624845 17.3336835162645 17.3507767682328 17.3673583683010 17.3834251889728 17.3989742243188 17.4140025889777 17.4285075171612 17.4424863616836 17.4559365930031 17.4688557982792 17.4812416804631 17.4930920574021 17.5044048609594 17.5151781361694 17.5254100404162 17.5350988426322 17.5442429225300 17.5528407698531 17.5608909836674 17.5683922716749 17.5753434495541 17.5817434403374 17.5875912738197 17.5928860859906 17.5976271185056 17.6018137181934 17.6054453365887 17.6085215294921 17.6110419565851 17.6130063810560 17.6144146692757 17.6152667904931 17.6155628165729 17.6153029217777 17.6144873825497 17.6131165773711 17.6111909866227 17.6087111924921 17.6056778789084 17.6020918315194 17.5979539376897 17.5932651865303 17.5880266689771 17.5822395778816 17.5759052081419 17.5690249568544 17.5616003235080 17.5536329102084 17.5451244219056 17.5360766666787 17.5264915560352 17.5163711052307 17.5057174336215 17.4945327650384 17.4828194281912 17.4705798570744 17.4578165914209 17.4445322771612 17.4307296669003 17.4164116204163 17.4015811051861 17.3862411968930 17.3703950800071 17.3540460483135 17.3371975054955 17.3198529657166 17.3020160542133 17.2836905078803 17.2648801758968 17.2455890203105 17.2258211166719 17.2055806546321 17.1848719385597 17.1636993881775 17.1420675391259 17.1199810436099 17.0974446709752 17.0744633083081 17.0510419609933 17.0271857533211 17.0028999289957 16.9781898517022 16.9530610056260 16.9275189959334 16.9015695492660 16.8752185141971 16.8484718616580 16.8213356853354 16.7938162020485 16.7659197520944 16.7376527995494 16.7090219325443 16.6800338634903 16.6506954292858 16.6210135914589 16.5909954362797 16.5606481748026 16.5299791429175 16.4989958012638 16.4677057351734 16.4361166545064 16.4042363934525 16.3720729102639 16.3396342869396 16.3069287288111 16.2739645641184 16.2407502434578 16.2072943392140 16.1736055448741 16.1396926742856 16.1055646608628 16.0712305566646 16.0366995314311 16.0019808715277 15.9670839787923 15.9320183693154 15.8967936721183 15.8614196277447 15.8259060867734 15.7902630082023 15.7545004577952 15.7186286062627 15.6826577274070 15.6465981961289 15.6104604863504 15.5742551688281 15.5379929088584 15.5016844638993 15.4653406810622 15.4289724945089 15.3925909227414 15.3562070657727 15.3198321022114 15.2834772861995 15.2471539442832 15.2108734721353 15.1746473311988 15.1384870451855 15.1024041964896 15.0664104224833 15.0305174116944 14.9947368998820 14.9590806659917 14.9235605280189 14.8881883387328 14.8529759813346 14.8179353649610 14.7830784201210 14.7484170940045 14.7139633457015 14.6797291413100 14.6457264489455 14.6119672336557 14.5784634522446 14.5452270479876 14.5122699452702 14.4796040441291 14.4472412147176 14.4151932916675 14.3834720684043 14.3520892913504 14.3210566540810 14.2903857914031 14.2600882733574 14.2301755991846 14.2006591911984 14.1715503886524 14.1428604415075 14.1146005041976 14.0867816293331 14.0594147613744 14.0325107302889 14.0060802451683 13.9801338878379 13.9546821064524 13.9297352090846 13.9053033573196 13.8813965598409 13.8580246660516 13.8351973596877 13.8129241524805 13.7912143778366 13.7700771845636 13.7495215306344 13.7295561770143 13.7101896815374 13.6914303928571 13.6732864444658 13.6557657488035];
E3rms2 = [13.6501618118740 13.6846854427050 13.7192048031420 13.7537084932569 13.7881851048765 13.8226232271449 13.8570114520718 13.8913383800518 13.9255926253544 13.9597628215654 13.9938376269920 14.0278057300087 14.0616558543526 14.0953767643425 14.1289572700528 14.1623862323873 14.1956525681005 14.2287452547209 14.2616533353951 14.2943659236419 14.3268722080101 14.3591614566515 14.3912230217772 14.4230463440252 14.4546209567253 14.4859364900374 14.5169826750212 14.5477493475456 14.5782264521345 14.6084040456635 14.6382723009700 14.6678215103377 14.6970420888591 14.7259245777108 14.7544596472774 14.7826381001957 14.8104508742650 14.8378890452541 14.8649438295873 14.8916065869388 14.9178688226869 14.9437221902975 14.9691584935631 14.9941696887700 15.0187478867383 15.0428853547790 15.0665745185265 15.0898079637154 15.1125784378128 15.1348788516029 15.1567022806500 15.1780419666916 15.1988913189375 15.2192439152942 15.2390935035041 15.2584340021986 15.2772595019060 15.2955642659517 15.3133427313231 15.3305895094299 15.3472993868442 15.3634673259449 15.3790884655237 15.3941581213284 15.4086717865524 15.4226251322861 15.4360140078977 15.4488344413897 15.4610826397069 15.4727549890007 15.4838480548697 15.4943585825265 15.5042834969956 15.5136199032095 15.5223650861333 15.5305165108284 15.5380718225027 15.5450288465418 15.5513855885047 15.5571402341192 15.5622911492444 15.5668368798107 15.5707761517784 15.5741078710362 15.5768311233089 15.5789451740722 15.5804494684112 15.5813436309217 15.5816274655424 15.5813009554350 15.5803642628199 15.5788177288163 15.5766618732673 15.5738973945769 15.5705251694913 15.5665462529490 15.5619618778408 15.5567734548168 15.5509825720598 15.5445909950652 15.5376006663883 15.5300137053998 15.5218324080202 15.5130592464484 15.5036968688686 15.4937480991401 15.4832159364815 15.4721035551314 15.4604143039907 15.4481517062269 15.4353194588884 15.4219214324779 15.4079616704723 15.3934443888792 15.3783739756902 15.3627549903524 15.3465921631817 15.3298903947634 15.3126547552601 15.2948904837537 15.2766029874851 15.2577978410619 15.2384807856391 15.2186577280111 15.1983347396925 15.1775180559036 15.1562140745227 15.1344293549673 15.1121706170131 15.0894447395332 15.0662587592071 15.0426198691026 15.0185354172315 14.9940129050021 14.9690599855984 14.9436844622841 14.9178942866183 14.8916975565766 14.8651025146013 14.8381175455496 14.8107511745509 14.7830120647697 14.7549090150790 14.7264509576311 14.6976469553195 14.6685061991645 14.6390380055783 14.6092518135252 14.5791571815935 14.5487637849507 14.5180814122031 14.4871199621374 14.4558894403688 14.4243999558824 14.3926617174617 14.3606850300216 14.3284802908300 14.2960579856406 14.2634286846967 14.2306030386702 14.1975917744739 14.1644056909890 14.1310556546992 14.0975525952301 14.0639075008027 14.0301314135925 13.9962354250144 13.9622306709311 13.9281283267742 13.8939396026002 13.8596757380826 13.8253479974385 13.7909676642932 13.7565460365001 13.7220944209095 13.6876241280906 13.6531464670256 13.6186727397689 13.5842142360856 13.5497822280674 13.5153879647379 13.4810426666611 13.4467575205467 13.4125436738588 13.3784122294523 13.3443742402292 13.3104407038206 13.2766225573169 13.2429306720322 13.2093758483340 13.1759688105271 13.1427202018053 13.1096405792799 13.0767404090971 13.0440300616345 13.0115198068044 12.9792198094605 12.9471401249165 12.9152906945725 12.8836813416894 12.8523217672705 12.8212215461015 12.7903901229157 12.7598368087207 12.7295707772841 12.6996010617559 12.6699365514879 12.6405859889960 12.6115579671060 12.5828609262745 12.5545031520929 12.5264927729681 12.4988377579852 12.4715459149723 12.4446248887380 12.4180821595045 12.3919250415252 12.3661606818999 12.3407960595819 12.3158379845532 12.2912930972102 12.2671678679238 12.2434685967797 12.2202014135024 12.1973722775544 12.1749869784135 12.1530511360015 12.1315702012972 12.1105494570956 12.0899940189192 12.0699088360785 12.0502986928790 12.0311682099378 12.0125218456676 11.9943638978329 11.9766985052461 11.9595296495594 11.9428611571508 11.9266967010848 11.9110398031738 11.8958938360750 11.8812620254771 11.8671474523098 11.8535530550059 11.8404816318088 11.8279358430506 11.8159182135030 11.8044311346858 11.7934768671882 11.7830575429516 11.7731751675678 11.7638316224828 11.7550286672078 11.7467679414549 11.7390509671990 11.7318791506947 11.7252537844028 11.7191760488316 11.7136470142828 11.7086676425048 11.7042387882406 11.7003612006571 11.6970355246633 11.6942623020990 11.6920419728112 11.6903748755838 11.6892612489430 11.6887012318087 11.6886948640464 11.6892420868165 11.6903427428371 11.6919965764676 11.6942032336612 11.6969622617697 11.7002731092131 11.7041351249802 11.7085475580330 11.7135095565271 11.7190201669388 11.7250783330262 11.7316828946782 11.7388325866576 11.7465260371917 11.7547617664775 11.7635381850756 11.7728535921944 11.7827061738955 11.7930940012031 11.8040150281390 11.8154670896947 11.8274478997176 11.8399550487864 11.8529860019871 11.8665380967037 11.8806085403527 11.8951944081090 11.9102926406240 11.9259000417388 11.9420132762254 11.9586288675304 11.9757431955597 11.9933524945007 12.0114528506910 12.0300402005595 12.0491103286098 12.0686588655154 12.0886812862716 12.1091729084686 12.1301288906419 12.1515442307573 12.1734137648091 12.1957321655430 12.2184939413206 12.2416934351185 12.2653248236920 12.2893821168660 12.3138591570293 12.3387496187468 12.3640470085840 12.3897446650811 12.4158357589270 12.4423132932984 12.4691701043951 12.4963988621598 12.5239920712032 12.5519420718967 12.5802410416766 12.6088809965389 12.6378537927293 12.6671511286100 12.6967645467533 12.7266854361889 12.7569050348657 12.7874144322929 12.8182045723546 12.8492662563297 12.8805901460505 12.9121667672884 12.9439865132458 12.9760396482629 13.0083163116576 13.0408065217155 13.0735001798469 13.1063870748673 13.1394568874135 13.1726991944990 13.2061034741788 13.2396591103393 13.2733553975797 13.3071815462151 13.3411266873439 13.3751798780284 13.4093301065295 13.4435662976280 13.4778773179883 13.5122519815974 13.5466790552307 13.5811472639687 13.6156452967366 13.6501618118740];
E3rms3 = [18.0176384861983 18.0044362506906 17.9906571469549 17.9763036656042 17.9613784295952 17.9458841951581 17.9298238527490 17.9132004280061 17.8960170827291 17.8782771158482 17.8599839644212 17.8411412046186 17.8217525527287 17.8018218661398 17.7813531443568 17.7603505299806 17.7388183097124 17.7167609153320 17.6941829246808 17.6710890626222 17.6474842020003 17.6233733645872 17.5987617219904 17.5736545965642 17.5480574622964 17.5219759456380 17.4954158263692 17.4683830383567 17.4408836703466 17.4129239666714 17.3845103279516 17.3556493117400 17.3263476331115 17.2966121652393 17.2664499398691 17.2358681477855 17.2048741391978 17.1734754240664 17.1416796723620 17.1094947142829 17.0769285403507 17.0439893014961 17.0106853090008 16.9770250344177 16.9430171093616 16.9086703252445 16.8739936328825 16.8389961420632 16.8036871209493 16.7680759954383 16.7321723483776 16.6959859186936 16.6595266003975 16.6228044414849 16.5858296427125 16.5486125562423 16.5111636841948 16.4734936770346 16.4356133318634 16.3975335905375 16.3592655376959 16.3208203986164 16.2822095369458 16.2434444522839 16.2045367776200 16.1654982766331 16.1263408408210 16.0870764864983 16.0477173516346 16.0082756925386 15.9687638803990 15.9291943976273 15.8895798341156 15.8499328832595 15.8102663378831 15.7705930859725 15.7309261062616 15.6912784636670 15.6516633045490 15.6120938518383 15.5725833999931 15.5331453097967 15.4937930030383 15.4545399569986 15.4153996988084 15.3763857996858 15.3375118689819 15.2987915481458 15.2602385044907 15.2218664248887 15.1836890093028 15.1457199642051 15.1079729958756 15.0704618036029 15.0332000727312 14.9962014676817 14.9594796247960 14.9230481451436 14.8869205872218 14.8511104595940 14.8156312134349 14.7804962350329 14.7457188382262 14.7113122567970 14.6772896368205 14.6436640289769 14.6104483808517 14.5776555292131 14.5452981922848 14.5133889620036 14.4819402963183 14.4509645114914 14.4204737744063 14.3904800949689 14.3609953184919 14.3320311181779 14.3035989876496 14.2757102335730 14.2483759683173 14.2216071027754 14.1954143392336 14.1698081643636 14.1447988423565 14.1203964081491 14.0966106608272 14.0734511571397 14.0509272051893 14.0290478582772 14.0078219089208 13.9872578830354 13.9673640343390 13.9481483388979 13.9296184899170 13.9117818927108 13.8946456598935 13.8782166067927 13.8625012470845 13.8475057886479 13.8332361296728 13.8196978549935 13.8068962326685 13.7948362108066 13.7835224146496 13.7729591439035 13.7631503703172 13.7540997355422 13.7458105492337 13.7382857874132 13.7315280911052 13.7255397652282 13.7203227777569 13.7158787591334 13.7122090019539 13.7093144609151 13.7071957530136 13.7058531580091 13.7052866191396 13.7054957441018 13.7064798062554 13.7082377461090 13.7107681730227 13.7140693671596 13.7181392816789 13.7229755451487 13.7285754641981 13.7349360263728 13.7420539032230 13.7499254535973 13.7585467271354 13.7679134679615 13.7780211185775 13.7888648239323 13.8004394356744 13.8127395165801 13.8257593451458 13.8394929203342 13.8539339664821 13.8690759383456 13.8849120262883 13.9014351615979 13.9186380219203 13.9365130368233 13.9550523934601 13.9742480423307 13.9940917031448 14.0145748707713 14.0356888212629 14.0574246179571 14.0797731176374 14.1027249767603 14.1262706577301 14.1504004352140 14.1751044025010 14.2003724778915 14.2261944111041 14.2525597897069 14.2794580455643 14.3068784612794 14.3348101766370 14.3632421950504 14.3921633899855 14.4215625113726 14.4514281919902 14.4817489538222 14.5125132143924 14.5437092930334 14.5753254171457 14.6073497283838 14.6397702888021 14.6725750869372 14.7057520438395 14.7392890190307 14.7731738163936 14.8073941900011 14.8419378498565 14.8767924675671 14.9119456819220 14.9473851044022 14.9830983246025 15.0190729155436 15.0552964389182 15.0917564502337 15.1284405038533 15.1653361579467 15.2024309793419 15.2397125482765 15.2771684630362 15.3147863445074 15.3525538406175 15.3904586306635 15.4284884295436 15.4666309918867 15.5048741160481 15.5432056480467 15.5816134853378 15.6200855805194 15.6586099449145 15.6971746520496 15.7357678410166 15.7743777197509 15.8129925681663 15.8516007412256 15.8901906718667 15.9287508738485 15.9672699445059 16.0057365673422 16.0441395145968 16.0824676496588 16.1207099294023 16.1588554064030 16.1968932311065 16.2348126538280 16.2726030267230 16.3102538056388 16.3477545518564 16.3850949337785 16.4222647285090 16.4592538233509 16.4960522172118 16.5326500219385 16.5690374635686 16.6052048834962 16.6411427395666 16.6768416070871 16.7122921797847 16.7474852706686 16.7824118128411 16.8170628602142 16.8514295882227 16.8855032943823 16.9192753988732 16.9527374450121 16.9858810996862 17.0186981537284 17.0511805222514 17.0833202448914 17.1151094860677 17.1465405351283 17.1776058065109 17.2082978398101 17.2386092998268 17.2685329766031 17.2980617853694 17.3271887664959 17.3559070854039 17.3842100324314 17.4120910226921 17.4395435958899 17.4665614161161 17.4931382716338 17.5192680746050 17.5449448608636 17.5701627895871 17.5949161430306 17.6191993262010 17.6430068665324 17.6663334135513 17.6891737385234 17.7115227341175 17.7333754140370 17.7547269126625 17.7755724846874 17.7959075047494 17.8157274670798 17.8350279851170 17.8538047911741 17.8720537360653 17.8897707887772 17.9069520361077 17.9235936823453 17.9396920489487 17.9552435742315 17.9702448130678 17.9846924366003 17.9985832319842 18.0119141021078 18.0246820653868 18.0368842555104 18.0485179212614 18.0595804263196 18.0700692491120 18.0799819826607 18.0893163344642 18.0980701263970 18.1062412946494 18.1138278896610 18.1208280761011 18.1272401328685 18.1330624531219 18.1382935443118 18.1429320282863 18.1469766413723 18.1504262345196 18.1532797734634 18.1555363388988 18.1571951267183 18.1582554482276 18.1587167304578 18.1585785164315 18.1578404655182 18.1565023537882 18.1545640743949 18.1520256380070 18.1488871732476 18.1451489271688 18.1408112657603 18.1358746744785 18.1303397588145 18.1242072448688 18.1174779799838 18.1101529333671 18.1022331967743 18.0937199851943 18.0846146375741 18.0749186175529 18.0646335142382 18.0537610429881 18.0423030462262 18.0302614942715 18.0176384861983];
E3rms4 = [14.9336025219877 14.9154370317619 14.8979888734257 14.8812658215758 14.8652753798514 14.8500247751069 14.8355209517207 14.8217705660458 14.8087799810162 14.7965552609055 14.7851021662620 14.7744261490075 14.7645323477297 14.7554255831520 14.7471103538082 14.7395908319179 14.7328708594716 14.7269539445310 14.7218432577600 14.7175416291692 14.7140515451149 14.7113751455287 14.7095142213860 14.7084702124428 14.7082442052190 14.7088369312317 14.7102487655211 14.7124797254100 14.7155294695634 14.7193972973045 14.7240821482165 14.7295826020223 14.7358968787384 14.7430228391278 14.7509579854092 14.7596994622751 14.7692440581870 14.7795882069423 14.7907279895309 14.8026591362914 14.8153770293020 14.8288767051074 14.8431528576586 14.8581998415815 14.8740116756690 14.8905820466679 14.9079043132981 14.9259715105640 14.9447763542722 14.9643112458308 14.9845682772637 15.0055392364692 15.0272156126988 15.0495886022561 15.0726491144085 15.0963877774924 15.1207949452350 15.1458607032387 15.1715748756680 15.1979270320793 15.2249064944385 15.2525023442775 15.2807034299896 15.3094983742690 15.3388755816639 15.3688232462542 15.3993293594177 15.4303817177042 15.4619679307857 15.4940754294860 15.5266914738814 15.5598031614164 15.5933974351324 15.6274610918580 15.6619807904768 15.6969430601833 15.7323343087491 15.7681408307928 15.8043488160202 15.8409443574604 15.8779134596573 15.9152420468138 15.9529159709205 15.9909210197818 16.0292429249972 16.0678673698905 16.1067799973100 16.1459664173965 16.1854122151934 16.2251029582139 16.2650242038635 16.3051615067586 16.3455004259204 16.3860265318610 16.4267254134894 16.4675826849617 16.5085839923050 16.5497150199538 16.5909614971145 16.6323092039911 16.6737439778344 16.7152517188522 16.7568183959477 16.7984300523059 16.8400728108017 16.8817328792423 16.9233965554504 16.9650502321780 17.0066804018424 17.0482736610800 17.0898167151626 17.1312963822251 17.1726995972905 17.2140134162032 17.2552250193010 17.2963217149912 17.3372909431154 17.3781202781892 17.4187974324136 17.4593102586084 17.4996467529153 17.5397950573749 17.5797434623498 17.6194804087750 17.6589944902912 17.6982744551870 17.7373092082392 17.7760878123883 17.8145994902725 17.8528336256314 17.8907797646041 17.9284276168440 17.9657670565625 18.0027881234155 18.0394810232934 18.0758361289828 18.1118439807329 18.1474952866883 18.1827809232600 18.2176919353528 18.2522195365280 18.2863551090605 18.3200902039127 18.3534165406255 18.3863260071157 18.4188106594212 18.4508627213471 18.4824745840585 18.5136388055960 18.5443481103425 18.5745953884276 18.6043736950597 18.6336762498265 18.6624964359410 18.6908277994431 18.7186640483492 18.7459990517766 18.7728268390310 18.7991415986492 18.8249376774254 18.8502095794088 18.8749519648653 18.8991596492444 18.9228276020941 18.9459509459953 18.9685249554517 18.9905450557961 19.0120068220679 19.0329059779007 19.0532383943969 19.0730000890090 19.0921872244118 19.1107961073900 19.1288231877220 19.1462650570763 19.1631184479130 19.1793802324006 19.1950474213403 19.2101171631084 19.2245867426137 19.2384535802706 19.2517152309924 19.2643693832037 19.2764138578749 19.2878466075785 19.2986657155694 19.3088693948893 19.3184559874983 19.3274239634312 19.3357719199827 19.3434985809207 19.3506027957298 19.3570835388794 19.3629399091316 19.3681711288711 19.3727765434735 19.3767556206973 19.3801079501213 19.3828332425968 19.3849313297557 19.3864021635249 19.3872458157040 19.3874624775471 19.3870524594012 19.3860161903671 19.3843542179928 19.3820672080113 19.3791559441008 19.3756213276787 19.3714643777381 19.3666862307102 19.3612881403599 19.3552714777060 19.3486377309948 19.3413885056871 19.3335255244735 19.3250506273333 19.3159657716113 19.3062730321399 19.2959746013646 19.2850727895211 19.2735700248265 19.2614688536981 19.2487719410020 19.2354820703130 19.2216021442279 19.2071351846572 19.1920843331805 19.1764528513936 19.1602441212914 19.1434616456542 19.1261090484730 19.1081900753533 19.0897085939880 19.0706685945889 19.0510741903552 19.0309296179572 19.0102392380275 18.9890075356338 18.9672391208027 18.9449387290022 18.9221112216622 18.8987615866746 18.8748949388956 18.8505165206786 18.8256317023288 18.8002459826424 18.7743649893684 18.7479944797054 18.7211403407388 18.6938085899442 18.6660053755748 18.6377369771178 18.6090098056969 18.5798304044342 18.5502054488257 18.5201417470791 18.4896462404163 18.4587260033466 18.4273882439185 18.3956403039384 18.3634896591408 18.3309439193323 18.2980108284825 18.2646982647957 18.2310142407134 18.1969669028849 18.1625645320644 18.1278155430179 18.0927284842827 18.0573120379601 18.0215750193942 17.9855263768115 17.9491751908926 17.9125306742899 17.8756021710385 17.8383991559630 17.8009312339424 17.7632081391637 17.7252397342447 17.6870360093033 17.6486070809678 17.6099631912548 17.5711147063990 17.5320721155883 17.4928460295960 17.4534471793415 17.4138864143411 17.3741747010715 17.3343231212426 17.2943428699423 17.2542452537374 17.2140416885950 17.1737436977728 17.1333629095622 17.0929110549323 17.0523999650704 17.0118415688045 16.9712478899297 16.9306310444050 16.8900032374492 16.8493767605118 16.8087639881393 16.7681773747231 16.7276294511148 16.6871328211549 16.6467001580526 16.6063442006748 16.5660777496883 16.5259136636081 16.4858648547160 16.4459442848618 16.4061649611499 16.3665399315020 16.3270822801161 16.2878051227844 16.2487216021356 16.2098448827163 16.1711881459999 16.1327645852538 16.0945874003260 16.0566697922923 16.0190249580206 15.9816660846136 15.9446063437747 15.9078588860427 15.8714368349463 15.8353532810707 15.7996212760209 15.7642538262928 15.7292638870803 15.6946643559788 15.6604680666248 15.6266877822562 15.5933361892026 15.5604258903168 15.5279693983319 15.4959791291841 15.4644673952584 15.4334463986087 15.4029282241302 15.3729248326888 15.3434480542335 15.3145095808786 15.2861209599625 15.2582935871059 15.2310386992549 15.2043673677295 15.1782904912787 15.1528187891495 15.1279627941755 15.1037328458984 15.0801390837190 15.0571914400955 15.0348996337892 15.0132731631688 14.9923212995807 14.9720530807901 14.9524773045050 14.9336025219877];
E3rms5 = [14.7413069329094 14.7787709122928 14.8162495640611 14.8537305767865 14.8912016240262 14.9286503702653 14.9660644768408 15.0034316078397 15.0407394359626 15.0779756483473 15.1151279523466 15.1521840812515 15.1891317999579 15.2259589105664 15.2626532579137 15.2992027350285 15.3355952885076 15.3718189238066 15.4078617104448 15.4437117871103 15.4793573666774 15.5147867411200 15.5499882863166 15.5849504667593 15.6196618401507 15.6541110618791 15.6882868894049 15.7221781865045 15.7557739274210 15.7890632008838 15.8220352140170 15.8546792961280 15.8869849023728 15.9189416173141 15.9505391583392 15.9817673789794 16.0126162721073 16.0430759730039 16.0731367623166 16.1027890689197 16.1320234726160 16.1608307067850 16.1892016608562 16.2171273827328 16.2445990810607 16.2716081274251 16.2981460584130 16.3242045776157 16.3497755574855 16.3748510411352 16.3994232440207 16.4234845555435 16.4470275405606 16.4700449408102 16.4925296762572 16.5144748463458 16.5358737312010 16.5567197927264 16.5770066756585 16.5967282085212 16.6158784045457 16.6344514625120 16.6524417675281 16.6698438917621 16.6866525951119 16.7028628258335 16.7184697211044 16.7334686075564 16.7478550017559 16.7616246106498 16.7747733319798 16.7872972546195 16.7991926589515 16.8104560171349 16.8210839934017 16.8310734442958 16.8404214188942 16.8491251590123 16.8571820993719 16.8645898677679 16.8713462852051 16.8774493660097 16.8828973179615 16.8876885423654 16.8918216341301 16.8952953818586 16.8981087678772 16.9002609683126 16.9017513530941 16.9025794860084 16.9027451247083 16.9022482207246 16.9010889194637 16.8992675602177 16.8967846761147 16.8936409941430 16.8898374350790 16.8853751134641 16.8802553375476 16.8744796092332 16.8680496239919 16.8609672707831 16.8532346319508 16.8448539831157 16.8358277930390 16.8261587234704 16.8158496289850 16.8049035568008 16.7933237465645 16.7811136301046 16.7682768311886 16.7548171652353 16.7407386389673 16.7260454501036 16.7107419869327 16.6948328279141 16.6783227412010 16.6612166841567 16.6435198027688 16.6252374310906 16.6063750905709 16.5869384893586 16.5669335215533 16.5463662663802 16.5252429873321 16.5035701312111 16.4813543271387 16.4586023854762 16.4353212966701 16.4115182300274 16.3872005324305 16.3623757269296 16.3370515112923 16.3112357564388 16.2849365048031 16.2581619685922 16.2309205279613 16.2032207290716 16.1750712820792 16.1464810589883 16.1174590914242 16.0880145682887 16.0581568333101 16.0278953824856 15.9972398614010 15.9662000624546 15.9347859219481 15.9030075170732 15.8708750627729 15.8383989084936 15.8055895348163 15.7724575499586 15.7390136861702 15.7052687960046 15.6712338484730 15.6369199250746 15.6023382157160 15.5675000145135 15.5324167154740 15.4970998080688 15.4615608726939 15.4258115760158 15.3898636662210 15.3537289681463 15.3174193783223 15.2809468599054 15.2443234375262 15.2075611920380 15.1706722551855 15.1336688041864 15.0965630562391 15.0593672629538 15.0220937047208 14.9847546850139 14.9473625246388 14.9099295559328 14.8724681169186 14.8349905454244 14.7975091731713 14.7600363198382 14.7225842871115 14.6851653527260 14.6477917645044 14.6104757344030 14.5732294325751 14.5360649814533 14.4989944498624 14.4620298471757 14.4251831175129 14.3884661339972 14.3518906930726 14.3154685088983 14.2792112078132 14.2431303229025 14.2072372886480 14.1715434356924 14.1360599857046 14.1007980463812 14.0657686065544 14.0309825314604 13.9964505581205 13.9621832909020 13.9281911972053 13.8944846033348 13.8610736905273 13.8279684911519 13.7951788850991 13.7627145963438 13.7305851896958 13.6988000677514 13.6673684680352 13.6362994603398 13.6056019442632 13.5752846469664 13.5453561211219 13.5158247430685 13.4866987111820 13.4579860444473 13.4296945812476 13.4018319783428 13.3744057100705 13.3474230677402 13.3208911592295 13.2948169087818 13.2692070569871 13.2440681609757 13.2194065947648 13.1952285498160 13.1715400357465 13.1483468812214 13.1256547349970 13.1034690671349 13.0817951703379 13.0606381614644 13.0400029831364 13.0198944054875 13.0003170280310 12.9812752816303 12.9627734305484 12.9448155746155 12.9274056514402 12.9105474387134 12.8942445565526 12.8785004698943 12.8633184909458 12.8487017816102 12.8346533559899 12.8211760828452 12.8082726880795 12.7959457571751 12.7841977376561 12.7730309414494 12.7624475472563 12.7524496028467 12.7430390272762 12.7342176130537 12.7259870282208 12.7183488183335 12.7113044083482 12.7048551044081 12.6990020955182 12.6937464550926 12.6890891423850 12.6850310037800 12.6815727739644 12.6787150769419 12.6764584269144 12.6748032289972 12.6737497798289 12.6732982679619 12.6734487741608 12.6742012715077 12.6755556253657 12.6775115931836 12.6800688241543 12.6832268586906 12.6869851277972 12.6913429522433 12.6962995416346 12.7018539933040 12.7080052910828 12.7147523039546 12.7220937845451 12.7300283675179 12.7385545678491 12.7476707789820 12.7573752708969 12.7676661880732 12.7785415473724 12.7899992358485 12.8020370084674 12.8146524858112 12.8278431516716 12.8416063506608 12.8559392857583 12.8708390158437 12.8863024532226 12.9023263611457 12.9189073513502 12.9360418816081 12.9537262533143 12.9719566091075 12.9907289305518 13.0100390358776 13.0298825777806 13.0502550413287 13.0711517419332 13.0925678234453 13.1144982563365 13.1369378360203 13.1598811812911 13.1833227329006 13.2072567522796 13.2316773204068 13.2565783368476 13.2819535189344 13.3077964011576 13.3341003346875 13.3608584871192 13.3880638423742 13.4157092008268 13.4437871795895 13.4722902130264 13.5012105534469 13.5305402720337 13.5602712599402 13.5903952296134 13.6209037163330 13.6517880799447 13.6830395067931 13.7146490118901 13.7466074412575 13.7789054744911 13.8115336275182 13.8444822555508 13.8777415562432 13.9113015730185 13.9451521986156 13.9792831787753 14.0136841161371 14.0483444742945 14.0832535820024 14.1184006375651 14.1537747133680 14.1893647605409 14.2251596137843 14.2611479963083 14.2973185249027 14.3336597151227 14.3701599865828 14.4068076683455 14.4435910044153 14.4804981593009 14.5175172236688 14.5546362200491 14.5918431086153 14.6291257930045 14.6664721261867 14.7038699163664 14.7413069329094];
E3rms6 = [19.2315956601318 19.2174271375321 19.2026571460017 19.1872883327523 19.1713234810122 19.1547655110603 19.1376174812795 19.1198825892196 19.1015641726813 19.0826657108018 19.0631908251624 19.0431432808916 19.0225269877881 19.0013460014337 18.9796045243194 18.9573069069684 18.9344576490555 18.9110614005220 18.8871229626926 18.8626472893622 18.8376394878981 18.8121048203144 18.7860487043182 18.7594767143619 18.7323945826600 18.7048082001640 18.6767236175535 18.6481470461437 18.6190848588008 18.5895435907961 18.5595299406301 18.5290507708133 18.4981131085894 18.4667241466326 18.4348912436532 18.4026219249808 18.3699238830809 18.3368049779860 18.3032732376726 18.2693368583929 18.2350042048658 18.2002838104747 18.1651843772927 18.1297147761058 18.0938840462796 18.0577013955783 18.0211761998447 17.9843180026339 17.9471365146779 17.9096416132960 17.8718433416589 17.8337519079531 17.7953776844199 17.7567312062719 17.7178231704869 17.6786644344554 17.6392660145263 17.5996390843801 17.5597949732953 17.5197451642345 17.4795012918211 17.4390751401487 17.3984786404351 17.3577238685348 17.3168230422841 17.2757885186993 17.2346327909952 17.1933684854582 17.1520083581440 17.1105652914168 17.0690522903276 17.0274824787785 16.9858690955987 16.9442254903635 16.9025651191073 16.8609015398314 16.8192484078493 16.7776194709706 16.7360285644962 16.6944896060661 16.6530165903256 16.6116235834178 16.5703247173483 16.5291341841399 16.4880662298484 16.4471351484456 16.4063552754974 16.3657409817496 16.3253066665024 16.2850667509036 16.2450356710659 16.2052278710617 16.1656577957850 16.1263398837080 16.0872885594726 16.0485182264473 16.0100432590969 15.9718779953081 15.9340367286004 15.8965337002662 15.8593830914125 15.8225990149523 15.7861955075238 15.7501865213667 15.7145859161407 15.6794074507115 15.6446647749144 15.6103714213013 15.5765407968718 15.5431861747941 15.5103206861649 15.4779573117737 15.4461088738689 15.4147880280292 15.3840072550070 15.3537788526915 15.3241149281030 15.2950273895008 15.2665279385275 15.2386280625188 15.2113390268677 15.1846718675313 15.1586373836677 15.1332461303954 15.1085084117249 15.0844342736148 15.0610334972248 15.0383155923258 15.0162897908909 14.9949650408821 14.9743500002541 14.9544530311273 14.9352821942134 14.9168452434378 14.8991496208060 14.8822024514950 14.8660105391979 14.8505803616956 14.8359180667113 14.8220294679907 14.8089200416663 14.7965949228786 14.7850589026698 14.7743164251528 14.7643715849468 14.7552281249075 14.7468894341217 14.7393585461936 14.7326381378067 14.7267305275785 14.7216376751967 14.7173611808272 14.7139022848184 14.7112618676837 14.7094404503637 14.7084381947607 14.7082549045571 14.7088900263062 14.7103426507848 14.7126115146205 14.7156950021797 14.7195911477059 14.7242976377288 14.7298118136994 14.7361306748918 14.7432508815221 14.7511687581169 14.7598802970942 14.7693811625752 14.7796666944013 14.7907319123657 14.8025715206349 14.8151799123734 14.8285511745470 14.8426790929076 14.8575571571467 14.8731785662151 14.8895362337933 14.9066227939140 14.9244306067231 14.9429517643711 14.9621780970302 14.9821011790272 15.0027123350845 15.0240026466610 15.0459629583862 15.0685838845780 15.0918558158383 15.1157689257166 15.1403131774360 15.1654783306719 15.1912539483789 15.2176294036537 15.2445938866337 15.2721364114168 15.3002458230018 15.3289108042352 15.3581198827693 15.3878614380081 15.4181237080569 15.4488947966416 15.4801626800239 15.5119152138732 15.5441401401180 15.5768250937533 15.6099576096034 15.6435251290425 15.6775150066555 15.7119145168391 15.7467108603474 15.7818911707650 15.8174425209080 15.8533519291466 15.8896063656598 15.9261927585971 15.9630980001512 16.0003089525500 16.0378124539479 16.0755953242316 16.1136443707124 16.1519463937317 16.1904881921560 16.2292565687669 16.2682383355477 16.3074203188505 16.3467893644746 16.3863323426018 16.4260361526466 16.4658877279736 16.5058740405094 16.5459821052281 16.5861989845357 16.6265117925096 16.6669076990606 16.7073739339358 16.7478977906220 16.7884666301350 16.8290678846886 16.8696890612269 16.9103177448735 16.9509416022273 16.9915483845706 17.0321259309443 17.0726621711093 17.1131451284217 17.1535629225343 17.1939037720595 17.2341559970652 17.2743080214970 17.3143483754535 17.3542656974202 17.3940487363147 17.4336863535001 17.4731675246688 17.5124813416080 17.5516170138990 17.5905638705109 17.6293113612900 17.6678490583601 17.7061666574381 17.7442539790644 17.7821009697341 17.8196977029552 17.8570343802125 17.8941013318740 17.9308890179959 17.9673880290675 18.0035890866578 18.0394830440524 18.0750608867296 18.1103137328572 18.1452328336716 18.1798095738152 18.2140354716102 18.2479021792862 18.2814014831129 18.3145253035474 18.3472656952592 18.3796148471731 18.4115650824082 18.4431088582013 18.4742387658117 18.5049475303339 18.5352280105138 18.5650731985224 18.5944762196814 18.6234303321812 18.6519289267517 18.6799655263157 18.7075337856270 18.7346274908513 18.7612405591942 18.7873670384211 18.8130011064464 18.8381370708604 18.8627693684527 18.8868925647322 18.9105013534284 18.9335905560042 18.9561551211446 18.9781901242551 18.9996907669474 19.0206523765386 19.0410704055489 19.0609404311834 19.0802581548575 19.0990194016902 19.1172201200395 19.1348563810114 19.1519243780117 19.1684204262947 19.1843409625281 19.1996825443751 19.2144418500870 19.2286156781277 19.2422009467856 19.2551946938580 19.2675940762921 19.2793963699042 19.2905989690790 19.3011993865359 19.3111952530732 19.3205843173736 19.3293644458107 19.3375336223190 19.3450899482475 19.3520316422639 19.3583570403002 19.3640645955104 19.3691528782504 19.3736205761232 19.3774664940202 19.3806895542155 19.3832887964845 19.3852633782543 19.3866125747967 19.3873357794329 19.3874325038085 19.3869023781563 19.3857451516303 19.3839606926600 19.3815489893261 19.3785101497918 19.3748444027583 19.3705520979417 19.3656337066054 19.3600898221087 19.3539211604943 19.3471285611078 19.3397129872488 19.3316755268491 19.3230173931916 19.3137399256473 19.3038445904522 19.2933329815044 19.2822068211965 19.2704679612696 19.2581183836966 19.2451602015896 19.2315956601318];
E3rms7 = [13.9099520935357 13.8930881647886 13.8769076535932 13.8614177561667 13.8466254123577 13.8325373002626 13.8191598309756 13.8064991434670 13.7945610996135 13.7833512793611 13.7728749760620 13.7631371919580 13.7541426338483 13.7458957089130 13.7384005207470 13.7316608655541 13.7256802285612 13.7204617806170 13.7160083750078 13.7123225444793 13.7094064984821 13.7072621206460 13.7058909664657 13.7052942612363 13.7054728982260 13.7064274370675 13.7081581024395 13.7106647829428 13.7139470302764 13.7180040586326 13.7228347443746 13.7284376259609 13.7348109041233 13.7419524423330 13.7498597674924 13.7585300709232 13.7679602096004 13.7781467076531 13.7890857581235 13.8007732250024 13.8132046454950 13.8263752325807 13.8402798777902 13.8549131542676 13.8702693200576 13.8863423216555 13.9031257977747 13.9206130833895 13.9387972139635 13.9576709299470 13.9772266814655 13.9974566332422 14.0183526697189 14.0399064003899 14.0621091653244 14.0849520408726 14.1084258455772 14.1325211462291 14.1572282641185 14.1825372814096 14.2084380477051 14.2349201867196 14.2619731031056 14.2895859893948 14.3177478330554 14.3464474236656 14.3756733601638 14.4054140582055 14.4356577575938 14.4663925297755 14.4976062854103 14.5292867819583 14.5614216313645 14.5939983077138 14.6270041549620 14.6604263946496 14.6942521336330 14.7284683718151 14.7630620098530 14.7980198568622 14.8333286380812 14.8689750024904 14.9049455304229 14.9412267410767 14.9778050999903 15.0146670264679 15.0517989008899 15.0891870719909 15.1268178639961 15.1646775837169 15.2027525275136 15.2410289881651 15.2794932616212 15.3181316536585 15.3569304863697 15.3958761046012 15.4349548821844 15.4741532280847 15.5134575923932 15.5528544721910 15.5923304172499 15.6318720356061 15.6714659989758 15.7110990480291 15.7507579975011 15.7904297411500 15.8301012565678 15.8697596098340 15.9093919600077 15.9489855634481 15.9885277780110 16.0280060670691 16.0674080033454 16.1067212726624 16.1459336774503 16.1850331401621 16.2240077065000 16.2628455485224 16.3015349675354 16.3400643969203 16.3784224047448 16.4165976962493 16.4545791162057 16.4923556511001 16.5299164312187 16.5672507325525 16.6043479786024 16.6411977420364 16.6777897462197 16.7141138666090 16.7501601320655 16.7859187259836 16.8213799873616 16.8565344117252 16.8913726519517 16.9258855189900 16.9600639824781 16.9938991712430 17.0273823737345 17.0605050383409 17.0932587736251 17.1256353484713 17.1576266921551 17.1892248943314 17.2204222049311 17.2512110340179 17.2815839515484 17.3115336870698 17.3410531293658 17.3701353260362 17.3987734830288 17.4269609640984 17.4546912902362 17.4819581390490 17.5087553440864 17.5350768941340 17.5609169324648 17.5862697560715 17.6111298148320 17.6354917106890 17.6593501967722 17.6827001765035 17.7055367026947 17.7278549766055 17.7496503470086 17.7709183092127 17.7916545040999 17.8118547171436 17.8315148774168 17.8506310565945 17.8691994679642 17.8872164654225 17.9046785424780 17.9215823312597 17.9379246015277 17.9537022596840 17.9689123478024 17.9835520426587 17.9976186547786 18.0111096274937 18.0240225360063 18.0363550864832 18.0481051151578 18.0592705874414 18.0698495970661 18.0798403652443 18.0892412398433 18.0980506945884 18.1062673282782 18.1138898640370 18.1209171485826 18.1273481515167 18.1331819646457 18.1384178013296 18.1430549958463 18.1470930027926 18.1505313965143 18.1533698705575 18.1556082371416 18.1572464266822 18.1582844873201 18.1587225844958 18.1585610005373 18.1578001342863 18.1564405007640 18.1544827308310 18.1519275709217 18.1487758827726 18.1450286431935 18.1406869438639 18.1357519911636 18.1302251060235 18.1241077238003 18.1174013941977 18.1101077811944 18.1022286630095 18.0937659320824 18.0847215950929 18.0750977730072 18.0648967011218 18.0541207291623 18.0427723213923 18.0308540567392 18.0183686289483 18.0053188467561 17.9917076340902 17.9775380302660 17.9628131902299 17.9475363848033 17.9317110009418 17.9153405420130 17.8984286280991 17.8809789962764 17.8629955009683 17.8444821142397 17.8254429261466 17.8058821450810 17.7858040981240 17.7652132313896 17.7441141104070 17.7225114204586 17.7004099669671 17.6778146758426 17.6547305938487 17.6311628889863 17.6071168508016 17.5825978907779 17.5576115426538 17.5321634627651 17.5062594303416 17.4799053478557 17.4531072412673 17.4258712603296 17.3982036788452 17.3701108948828 17.3415994310087 17.3126759344775 17.2833471773957 17.2536200568526 17.2235015950343 17.1929989393039 17.1621193622392 17.1308702616432 17.0992591605073 17.0672937069553 17.0349816741231 17.0023309600102 16.9693495872609 16.9360457029564 16.9024275782702 16.8685036081533 16.8342823109152 16.7997723277689 16.7649824223123 16.7299214799607 16.6945985072818 16.6590226313245 16.6232030988204 16.5871492753718 16.5508706445232 16.5143768067797 16.4776774785757 16.4407824911172 16.4037017891850 16.3664454298493 16.3290235810874 16.2914465203398 16.2537246329651 16.2158684106115 16.1778884495129 16.1397954486586 16.1016002079287 16.0633136260670 16.0249466986200 15.9865105157432 15.9480162599236 15.9094752035979 15.8708987066655 15.8322982139207 15.7936852523594 15.7550714283931 15.7164684249601 15.6778879985204 15.6393419759688 15.6008422514028 15.5624007828296 15.5240295887262 15.4857407445247 15.4475463789564 15.4094586703142 15.3714898426009 15.3336521615625 15.2959579306233 15.2584194867041 15.2210491959523 15.1838594493353 15.1468626581728 15.1100712495166 15.0734976614683 15.0371543383715 15.0010537259179 14.9652082661466 14.9296303923465 14.8943325238687 14.8593270608516 14.8246263788389 14.7902428233260 14.7561887042108 14.7224762901727 14.6891178029507 14.6561254115786 14.6235112265102 14.5912872937004 14.5594655886139 14.5280580101573 14.4970763745800 14.4665324092841 14.4364377466314 14.4068039176535 14.3776423457618 14.3489643404002 14.3207810906681 14.2931036589311 14.2659429743934 14.2393098266659 14.2132148593223 14.1876685634505 14.1626812712132 14.1382631494004 14.1144241930222 14.0911742188938 14.0685228592730 14.0464795555167 14.0250535517869 14.0042538887956 13.9840893976159 13.9645686935448 13.9457001700414 13.9274919927368 13.9099520935357];
E3rms8 = [13.5325796763007 13.5669942164497 13.6014408799985 13.6359084057103 13.6703855119852 13.7048609022424 13.7393232702944 13.7737613056983 13.8081636990844 13.8425191474441 13.8768163593863 13.9110440603443 13.9451909977377 13.9792459460667 14.0131977119682 14.0470351391797 14.0807471134600 14.1143225674173 14.1477504852695 14.1810199075203 14.2141199355495 14.2470397361247 14.2797685458075 14.3122956752775 14.3446105135636 14.3767025321573 14.4085612890651 14.4401764327118 14.4715377057861 14.5026349489486 14.5334581044605 14.5639972196960 14.5942424505417 14.6241840647173 14.6538124449543 14.6831180921029 14.7120916281143 14.7407237989280 14.7690054772475 14.7969276652313 14.8244814970536 14.8516582414023 14.8784493038404 14.9048462291038 14.9308407032802 14.9564245559132 14.9815897619891 15.0063284438762 15.0306328731271 15.0544954722405 15.0779088163087 15.1008656346035 15.1233588120735 15.1453813907747 15.1669265712225 15.1879877136643 15.2085583393128 15.2286321314787 15.2482029366736 15.2672647656152 15.2858117942172 15.3038383644893 15.3213389854052 15.3383083337120 15.3547412546921 15.3706327628912 15.3859780427809 15.4007724494018 15.4150115089621 15.4286909193971 15.4418065509094 15.4543544464389 15.4663308221674 15.4777320679260 15.4885547476332 15.4987955996762 15.5084515372814 15.5175196488690 15.5259971983746 15.5338816255725 15.5411705463702 15.5478617530800 15.5539532147134 15.5594430772162 15.5643296637143 15.5686114747679 15.5722871885715 15.5753556611972 15.5778159267695 15.5796671976887 15.5809088648117 15.5815404976352 15.5815618444662 15.5809728326032 15.5797735684617 15.5779643377683 15.5755456056637 15.5725180168525 15.5688823957177 15.5646397464415 15.5597912530908 15.5543382797111 15.5482823703965 15.5416252493515 15.5343688209309 15.5265151696586 15.5180665602392 15.5090254375431 15.4993944265712 15.4891763323772 15.4783741399986 15.4669910143419 15.4550303000112 15.4424955211707 15.4293903812982 15.4157187629551 15.4014847274919 15.3866925147378 15.3713465425924 15.3554514066477 15.3390118797015 15.3220329112349 15.3045196268566 15.2864773276542 15.2679114895276 15.2488277624214 15.2292319695198 15.2091301063616 15.1885283398891 15.1674330074097 15.1458506155210 15.1237878389018 15.1012515190708 15.0782486630341 15.0547864418528 15.0308721891286 15.0065133993927 14.9817177263931 14.9564929813032 14.9308471308213 14.9047882951697 14.8783247459934 14.8514649041595 14.8242173374468 14.7965907581176 14.7685940204056 14.7402361178737 14.7115261806593 14.6824734726225 14.6530873883683 14.6233774501654 14.5933533047378 14.5630247199537 14.5324015813996 14.5014938888322 14.4703117525259 14.4388653895004 14.4071651196501 14.3752213617363 14.3430446293049 14.3106455264695 14.2780347436003 14.2452230529115 14.2122213039446 14.1790404189578 14.1456913882123 14.1121852651763 14.0785331616434 14.0447462427567 14.0108357219587 13.9768128558701 13.9426889390949 13.9084752989556 13.8741832901774 13.8398242895124 13.8054096903113 13.7709508970594 13.7364593198716 13.7019463689592 13.6674234490694 13.6329019539064 13.5983932605510 13.5639087238714 13.5294596709312 13.4950573954214 13.4607131521058 13.4264381512892 13.3922435533288 13.3581404631754 13.3241399249770 13.2902529167331 13.2564903450136 13.2228630397510 13.1893817491207 13.1560571344963 13.1228997655116 13.0899201152232 13.0571285553843 13.0245353518252 12.9921506599817 12.9599845205297 12.9280468551789 12.8963474625910 12.8648960144615 12.8337020517602 12.8027749811102 12.7721240713650 12.7417584503314 12.7116871016763 12.6819188620126 12.6524624181696 12.6233263046438 12.5945189012345 12.5660484308852 12.5379229577008 12.5101503851645 12.4827384545438 12.4556947434976 12.4290266648803 12.4027414657179 12.3768462264006 12.3513478600538 12.3262531120947 12.3015685599785 12.2773006131242 12.2534555130250 12.2300393335153 12.2070579812278 12.1845171962025 12.1624225526543 12.1407794598945 12.1195931634065 12.0988687460356 12.0786111293531 12.0588250750984 12.0395151867670 12.0206859113004 12.0023415408762 11.9844862147786 11.9671239213761 11.9502585001392 11.9338936437556 11.9180329002722 11.9026796752969 11.8878372342488 11.8735087045876 11.8596970781189 11.8464052132622 11.8336358373381 11.8213915488208 11.8096748196107 11.7984879972177 11.7878333069555 11.7777128540651 11.7681286257690 11.7590824932826 11.7505762137409 11.7426114320422 11.7351896825994 11.7283123910030 11.7219808755800 11.7161963488383 11.7109599188007 11.7062725902132 11.7021352656420 11.6985487464251 11.6955137335027 11.6930308280930 11.6911005322707 11.6897232493428 11.6888992841397 11.6886288431256 11.6889120343798 11.6897488674287 11.6911392529427 11.6930830022634 11.6955798268338 11.6986293374433 11.7022310433782 11.7063843514042 11.7110885646352 11.7163428812946 11.7221463933210 11.7284980848870 11.7353968308029 11.7428413948083 11.7508304277824 11.7593624658541 11.7684359284338 11.7780491161804 11.7882002088761 11.7988872632882 11.8101082109263 11.8218608558098 11.8341428721750 11.8469518021683 11.8602850535269 11.8741398972498 11.8885134652909 11.9034027482503 11.9188045930998 11.9347157009419 11.9511326248082 11.9680517675258 11.9854693796185 12.0033815573163 12.0217842406127 12.0406732114403 12.0600440919151 12.0798923427110 12.1002132615406 12.1210019817552 12.1422534710800 12.1639625304778 12.1861237931718 12.2087317237893 12.2317806177029 12.2552646004839 12.2791776275630 12.3035134840331 12.3282657846489 12.3534279739845 12.3789933267824 12.4049549484825 12.4313057759506 12.4580385783705 12.4851459583423 12.5126203531674 12.5404540363246 12.5686391191185 12.5971675525516 12.6260311293461 12.6552214861766 12.6847301060806 12.7145483210377 12.7446673147528 12.7750781255762 12.8057716496470 12.8367386441431 12.8679697307437 12.8994553992237 12.9311860111976 12.9631518040303 12.9953428948708 13.0277492848216 13.0603608632479 13.0931674121976 13.1261586109460 13.1593240406352 13.1926531890354 13.2261354553732 13.2597601552756 13.2935165257700 13.3273937303744 13.3613808642331 13.3954669593307 13.4296409897377 13.4638918769098 13.4982084950162 13.5325796763007];
E3rms9 = [17.4624627201649 17.4492777075197 17.4355628607234 17.4213206508466 17.4065536661368 17.3912646129602 17.3754563167696 17.3591317230831 17.3422938984901 17.3249460316532 17.3070914343429 17.2887335424658 17.2698759171184 17.2505222456230 17.2306763426035 17.2103421510287 17.1895237432866 17.1682253222373 17.1464512222759 17.1242059103782 17.1014939871492 17.0783201878661 17.0546893834857 17.0306065816615 17.0060769277404 16.9811057057079 16.9556983391736 16.9298603922545 16.9035975704965 16.8769157217167 16.8498208368435 16.8223190507052 16.7944166427679 16.7661200378643 16.7374358068255 16.7083706671091 16.6789314833506 16.6491252678604 16.6189591810578 16.5884405318655 16.5575767779897 16.5263755261909 16.4948445324174 16.4629917019188 16.4308250892359 16.3983528981362 16.3655834814277 16.3325253407334 16.2991871261087 16.2655776356150 16.2317058147548 16.1975807558216 16.1632116971315 16.1286080221499 16.0937792584992 16.0587350768370 16.0234852896443 15.9880398498509 15.9524088493683 15.9166025174525 15.8806312189761 15.8445054525334 15.8082358484202 15.7718331664700 15.7353082937439 15.6986722420865 15.6619361455123 15.6251112574622 15.5882089479015 15.5512407002657 15.5142181082635 15.4771528724838 15.4400567969173 15.4029417852454 15.3658198370307 15.3287030437177 15.2916035844857 15.2545337219506 15.2175057976952 15.1805322276621 15.1436254973804 15.1067981570306 15.0700628163927 15.0334321396002 14.9969188397656 14.9605356734829 14.9242954351414 14.8882109511534 14.8522950739867 14.8165606761211 14.7810206438436 14.7456878709272 14.7105752521895 14.6756956769517 14.6410620223425 14.6066871465748 14.5725838820422 14.5387650283789 14.5052433454084 14.4720315460304 14.4391422890114 14.4065881717306 14.3743817228570 14.3425353949830 14.3110615572090 14.2799724876900 14.2492803661668 14.2189972664716 14.1891351490271 14.1597058533270 14.1307210904543 14.1021924356008 14.0741313205877 14.0465490264783 14.0194566761697 13.9928652270811 13.9667854638841 13.9412279913246 13.9162032270768 13.8917213947551 13.8677925169706 13.8444264085085 13.8216326696412 13.7994206795329 13.7777995898155 13.7567783182729 13.7363655426979 13.7165696948991 13.6973989548780 13.6788612451656 13.6609642253788 13.6437152869139 13.6271215478810 13.6111898482158 13.5959267450058 13.5813385080383 13.5674311155636 13.5542102502767 13.5416812955463 13.5298493318681 13.5187191335583 13.5082951656927 13.4985815812959 13.4895822187780 13.4813005996140 13.4737399263008 13.4669030805527 13.4607926217545 13.4554107856862 13.4507594834974 13.4468403009509 13.4436544979116 13.4412030081075 13.4394864391476 13.4385050727897 13.4382588654687 13.4387474490744 13.4399701319909 13.4419259003580 13.4446134196096 13.4480310362286 13.4521767797449 13.4570483649753 13.4626431944778 13.4689583612442 13.4759906515930 13.4837365482926 13.4921922338869 13.5013535942181 13.5112162221481 13.5217754214757 13.5330262110293 13.5449633289397 13.5575812370860 13.5708741257037 13.5848359181446 13.5994602757952 13.6147406031312 13.6306700529126 13.6472415315038 13.6644477043119 13.6822810013481 13.7007336228913 13.7197975452440 13.7394645265892 13.7597261129313 13.7805736441101 13.8019982598885 13.8239909060988 13.8465423408524 13.8696431407947 13.8932837073985 13.9174542732995 13.9421449086583 13.9673455275395 13.9930458943118 14.0192356300610 14.0459042189972 14.0730410148594 14.1006352473216 14.1286760283741 14.1571523586895 14.1860531339580 14.2153671511943 14.2450831150187 14.2751896438693 14.3056752762033 14.3365284766202 14.3677376419445 14.3992911072405 14.4311771517731 14.4633840048935 14.4958998518521 14.5287128395485 14.5618110821906 14.5951826668820 14.6288156591125 14.6626981081764 14.6968180524992 14.7311635248530 14.7657225575008 14.8004831872346 14.8354334603095 14.8705614372824 14.9058551977501 14.9413028449820 14.9768925104379 15.0126123581952 15.0484505892614 15.0843954457707 15.1204352150803 15.1565582337612 15.1927528914517 15.2290076346484 15.2653109703302 15.3016514695110 15.3380177706664 15.3743985830519 15.4107826899008 15.4471589515366 15.4835163083364 15.5198437836281 15.5561304864393 15.5923656141625 15.6285384551254 15.6646383909970 15.7006548991634 15.7365775549479 15.7723960337494 15.8081001130580 15.8436796744196 15.8791247052296 15.9144253004933 15.9495716644646 15.9845541121730 16.0193630708923 16.0539890814993 16.0884227997474 16.1226549974452 16.1566765635630 16.1904785052509 16.2240519487713 16.2573881403572 16.2904784469834 16.3233143570838 16.3558874811701 16.3881895523984 16.4202124270381 16.4519480849357 16.4833886298234 16.5145262896496 16.5453534167969 16.5758624882613 16.6060461057707 16.6358969958606 16.6654080098595 16.6945721238879 16.7233824387432 16.7518321797985 16.7799146968075 16.8076234636914 16.8349520783086 16.8618942621388 16.8884438599705 16.9145948395468 16.9403412911670 16.9656774272808 16.9905975820387 17.0150962108248 17.0391678897777 17.0628073152544 17.0860093033398 17.1087687892582 17.1310808268447 17.1529405879592 17.1743433619040 17.1952845548294 17.2157596891255 17.2357644028288 17.2552944490019 17.2743456951218 17.2929141224666 17.3109958254973 17.3285870112602 17.3456839987563 17.3622832183652 17.3783812112314 17.3939746286974 17.4090602317067 17.4236348902498 17.4376955828132 17.4512393958357 17.4642635231895 17.4767652656655 17.4887420304979 17.5001913308691 17.5111107854901 17.5214981181378 17.5313511572677 17.5406678356152 17.5494461898453 17.5576843602057 17.5653805902147 17.5725332263738 17.5791407179212 17.5852016165857 17.5907145763903 17.5956783534788 17.6000918059779 17.6039538938652 17.6072636789139 17.6100203246190 17.6122230961901 17.6138713605664 17.6149645864493 17.6155023444014 17.6154843069303 17.6149102486708 17.6137800465283 17.6120936799148 17.6098512309883 17.6070528849263 17.6036989302535 17.5997897591810 17.5953258679873 17.5903078574368 17.5847364332251 17.5786124064673 17.5719366941993 17.5647103199408 17.5569344142555 17.5486102153751 17.5397390698315 17.5303224331338 17.5203618704599 17.5098590573948 17.4988157806819 17.4872339390101 17.4751155438237 17.4624627201649];

h = figure();
plot(thetad,E3rms9,'-','LineWidth',3.5);
ylabel('Campo elétrico (kV/cm)')
xlabel('Ângulo (°)')
grid
hold on
plot(thetad,Erms,':','LineWidth',3.5);
hold on
plot(thetad,E_crit,'--','LineWidth',3.5);
xlim([0 360])
% ylim([13.5 15.5])
legend('Original','Raio maior', 'Campo elétrico crítico')
% title('Campo elétrico superficial condutor um')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'solo_caso1_dist','-dpdf','-r0')