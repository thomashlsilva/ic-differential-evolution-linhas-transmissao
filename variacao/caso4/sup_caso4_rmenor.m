% Raio menor
clear all; close all; clc
% Cálculo do campo elétrico superficial do caso 4 do sistema trifásico

e_0 = 8.854*(10^(-12));
r = 0.00788; % raio do condutor
n = 12; % número de condutores do sistema
nc = n*2; % número de condutores total (com imagens)
ci = (2*n-1); % número de cargas imagens por condutor

cond = 12;
%% matriz P

% posição dos condutores reais do sistema
xr = [-7.975 -7.025 -7.025 -7.975 -0.475 0.475 0.475 -0.475 7.025 7.975 7.975 7.025];
yr = [18.45 18.45 17.5 17.5 25.95 25.95 25 25 18.45 18.45 17.5 17.5];
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

%tensão condutores fase a
V_ra = V/sqrt(3);
V_ia = 0;

%tensão condutores fase b
V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

%tensão condutores fase c
V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ];

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

%cabo 4 fase a
rho_r4 = real(rho(4));
rho_i4 = imag(rho(4));

%cabo 5 fase b
rho_r5 = real(rho(5));
rho_i5 = imag(rho(5));

%cabo 6 fase b
rho_r6 = real(rho(6));
rho_i6 = imag(rho(6));

%cabo 7 fase b
rho_r7 = real(rho(7));
rho_i7 = imag(rho(7));

%cabo 8 fase b
rho_r8 = real(rho(8));
rho_i8 = imag(rho(8));

%cabo 9 fase c
rho_r9 = real(rho(9));
rho_i9 = imag(rho(9));

%cabo 10 fase c
rho_r10 = real(rho(10));
rho_i10 = imag(rho(10));

%cabo 11 fase c
rho_r11 = real(rho(11));
rho_i11 = imag(rho(11));

%cabo 12 fase c
rho_r12 = real(rho(12));
rho_i12 = imag(rho(12));


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
        elseif x(i)~=x(j) && y(j) > y(i)
            if x(j) > x(i)
                phix(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posx(i,j) = x(i) + delta(i,j)*cos(phix(i,j));
            elseif x(i) > x(j)
                phix(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posx(i,j) = x(i) - delta(i,j)*cos(phix(i,j));
            end
        elseif x(i)~=x(j) &&  y(i) > y(j)
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
            if y(j) > y(i) || y(i) < 0 && y(j) < 0
                posy(i,j) = y(i) + delta(i,j);
            elseif y(i) > y(j)
                posy(i,j) = y(i) - delta(i,j);
            end
        elseif x(i)~=x(j) && y(j) > y(i)
            if  y(i) > 0 || y(i) < 0 && y(j) < 0
                phiy(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posy(i,j) = y(i) + delta(i,j)*sin(phiy(i,j));
            elseif y(i) < 0 && y(j) > 0
                phiy(i,j) = asin((y(j)-y(i))/(distance(i,j)));
                posy(i,j) = y(i) - delta(i,j)*sin(phiy(i,j));
            end
        elseif x(i)~=x(j) && y(i) > y(j)
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
E_xr110 = (-rho_r10/(2*pi*e_0)).*((xf - posx(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_xr1_11 = (-rho_r11/(2*pi*e_0)).*((xf - posx(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_xr1_12 = (-rho_r12/(2*pi*e_0)).*((xf - posx(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_xr1_13 = (rho_r1/(2*pi*e_0)).*((xf - posx(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_xr1_14 = (rho_r2/(2*pi*e_0)).*((xf - posx(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_xr1_15 = (rho_r3/(2*pi*e_0)).*((xf - posx(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_xr1_16 = (rho_r4/(2*pi*e_0)).*((xf - posx(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_xr1_17 = (rho_r5/(2*pi*e_0)).*((xf - posx(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_xr1_18 = (rho_r6/(2*pi*e_0)).*((xf - posx(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));
E_xr1_19 = (rho_r7/(2*pi*e_0)).*((xf - posx(1,19))./((xf - posx(1,19)).^2 + (yf - posy(1,19)).^2));
E_xr1_20 = (rho_r8/(2*pi*e_0)).*((xf - posx(1,20))./((xf - posx(1,20)).^2 + (yf - posy(1,20)).^2));
E_xr1_21 = (rho_r9/(2*pi*e_0)).*((xf - posx(1,21))./((xf - posx(1,21)).^2 + (yf - posy(1,21)).^2));
E_xr1_22 = (rho_r10/(2*pi*e_0)).*((xf - posx(1,22))./((xf - posx(1,22)).^2 + (yf - posy(1,22)).^2));
E_xr1_23 = (rho_r11/(2*pi*e_0)).*((xf - posx(1,23))./((xf - posx(1,23)).^2 + (yf - posy(1,23)).^2));
E_xr1_24 = (rho_r12/(2*pi*e_0)).*((xf - posx(1,24))./((xf - posx(1,24)).^2 + (yf - posy(1,24)).^2));

E_xr21 = (-rho_r1/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xr23 = (-rho_r3/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xr24 = (-rho_r4/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xr25 = (-rho_r5/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xr26 = (-rho_r6/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_xr27 = (-rho_r7/(2*pi*e_0)).*((xf - posx(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_xr28 = (-rho_r8/(2*pi*e_0)).*((xf - posx(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_xr29 = (-rho_r9/(2*pi*e_0)).*((xf - posx(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_xr210 = (-rho_r10/(2*pi*e_0)).*((xf - posx(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_xr2_11 = (-rho_r11/(2*pi*e_0)).*((xf - posx(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_xr2_12 = (-rho_r12/(2*pi*e_0)).*((xf - posx(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_xr2_13 = (rho_r1/(2*pi*e_0)).*((xf - posx(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_xr2_14 = (rho_r2/(2*pi*e_0)).*((xf - posx(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_xr2_15 = (rho_r3/(2*pi*e_0)).*((xf - posx(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_xr2_16 = (rho_r4/(2*pi*e_0)).*((xf - posx(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_xr2_17 = (rho_r5/(2*pi*e_0)).*((xf - posx(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_xr2_18 = (rho_r6/(2*pi*e_0)).*((xf - posx(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));
E_xr2_19 = (rho_r7/(2*pi*e_0)).*((xf - posx(2,19))./((xf - posx(2,19)).^2 + (yf - posy(2,19)).^2));
E_xr2_20 = (rho_r8/(2*pi*e_0)).*((xf - posx(2,20))./((xf - posx(2,20)).^2 + (yf - posy(2,20)).^2));
E_xr2_21 = (rho_r9/(2*pi*e_0)).*((xf - posx(2,21))./((xf - posx(2,21)).^2 + (yf - posy(2,21)).^2));
E_xr2_22 = (rho_r10/(2*pi*e_0)).*((xf - posx(2,22))./((xf - posx(2,22)).^2 + (yf - posy(2,22)).^2));
E_xr2_23 = (rho_r11/(2*pi*e_0)).*((xf - posx(2,23))./((xf - posx(2,23)).^2 + (yf - posy(2,23)).^2));
E_xr2_24 = (rho_r12/(2*pi*e_0)).*((xf - posx(2,24))./((xf - posx(2,24)).^2 + (yf - posy(2,24)).^2));

E_xr31 = (-rho_r1/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xr32 = (-rho_r2/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xr34 = (-rho_r4/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xr35 = (-rho_r5/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xr36 = (-rho_r6/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_xr37 = (-rho_r7/(2*pi*e_0)).*((xf - posx(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_xr38 = (-rho_r8/(2*pi*e_0)).*((xf - posx(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_xr39 = (-rho_r9/(2*pi*e_0)).*((xf - posx(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_xr310 = (-rho_r10/(2*pi*e_0)).*((xf - posx(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_xr311 = (-rho_r11/(2*pi*e_0)).*((xf - posx(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_xr312 = (-rho_r12/(2*pi*e_0)).*((xf - posx(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_xr313 = (rho_r1/(2*pi*e_0)).*((xf - posx(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_xr314 = (rho_r2/(2*pi*e_0)).*((xf - posx(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_xr315 = (rho_r3/(2*pi*e_0)).*((xf - posx(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_xr316 = (rho_r4/(2*pi*e_0)).*((xf - posx(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_xr317 = (rho_r5/(2*pi*e_0)).*((xf - posx(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_xr318 = (rho_r6/(2*pi*e_0)).*((xf - posx(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));
E_xr319 = (rho_r7/(2*pi*e_0)).*((xf - posx(3,19))./((xf - posx(3,19)).^2 + (yf - posy(3,19)).^2));
E_xr320 = (rho_r8/(2*pi*e_0)).*((xf - posx(3,20))./((xf - posx(3,20)).^2 + (yf - posy(3,20)).^2));
E_xr321 = (rho_r9/(2*pi*e_0)).*((xf - posx(3,21))./((xf - posx(3,21)).^2 + (yf - posy(3,21)).^2));
E_xr322 = (rho_r10/(2*pi*e_0)).*((xf - posx(3,22))./((xf - posx(3,22)).^2 + (yf - posy(3,22)).^2));
E_xr323 = (rho_r11/(2*pi*e_0)).*((xf - posx(3,23))./((xf - posx(3,23)).^2 + (yf - posy(3,23)).^2));
E_xr324 = (rho_r12/(2*pi*e_0)).*((xf - posx(3,24))./((xf - posx(3,24)).^2 + (yf - posy(3,24)).^2));

E_xr41 = (-rho_r1/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xr42 = (-rho_r2/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xr43 = (-rho_r3/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xr45 = (-rho_r5/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xr46 = (-rho_r6/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_xr47 = (-rho_r7/(2*pi*e_0)).*((xf - posx(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_xr48 = (-rho_r8/(2*pi*e_0)).*((xf - posx(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_xr49 = (-rho_r9/(2*pi*e_0)).*((xf - posx(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_xr410 = (-rho_r10/(2*pi*e_0)).*((xf - posx(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_xr411 = (-rho_r11/(2*pi*e_0)).*((xf - posx(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_xr412 = (-rho_r12/(2*pi*e_0)).*((xf - posx(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_xr413 = (rho_r1/(2*pi*e_0)).*((xf - posx(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_xr414 = (rho_r2/(2*pi*e_0)).*((xf - posx(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_xr415 = (rho_r3/(2*pi*e_0)).*((xf - posx(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_xr416 = (rho_r4/(2*pi*e_0)).*((xf - posx(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_xr417 = (rho_r5/(2*pi*e_0)).*((xf - posx(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_xr418 = (rho_r6/(2*pi*e_0)).*((xf - posx(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));
E_xr419 = (rho_r7/(2*pi*e_0)).*((xf - posx(4,19))./((xf - posx(4,19)).^2 + (yf - posy(4,19)).^2));
E_xr420 = (rho_r8/(2*pi*e_0)).*((xf - posx(4,20))./((xf - posx(4,20)).^2 + (yf - posy(4,20)).^2));
E_xr421 = (rho_r9/(2*pi*e_0)).*((xf - posx(4,21))./((xf - posx(4,21)).^2 + (yf - posy(4,21)).^2));
E_xr422 = (rho_r10/(2*pi*e_0)).*((xf - posx(4,22))./((xf - posx(4,22)).^2 + (yf - posy(4,22)).^2));
E_xr423 = (rho_r11/(2*pi*e_0)).*((xf - posx(4,23))./((xf - posx(4,23)).^2 + (yf - posy(4,23)).^2));
E_xr424 = (rho_r12/(2*pi*e_0)).*((xf - posx(4,24))./((xf - posx(4,24)).^2 + (yf - posy(4,24)).^2));

E_xr51 = (-rho_r1/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xr52 = (-rho_r2/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xr53 = (-rho_r3/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xr54 = (-rho_r4/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xr56 = (-rho_r6/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_xr57 = (-rho_r7/(2*pi*e_0)).*((xf - posx(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_xr58 = (-rho_r8/(2*pi*e_0)).*((xf - posx(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_xr59 = (-rho_r9/(2*pi*e_0)).*((xf - posx(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_xr510 = (-rho_r10/(2*pi*e_0)).*((xf - posx(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_xr511 = (-rho_r11/(2*pi*e_0)).*((xf - posx(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_xr512 = (-rho_r12/(2*pi*e_0)).*((xf - posx(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_xr513 = (rho_r1/(2*pi*e_0)).*((xf - posx(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_xr514 = (rho_r2/(2*pi*e_0)).*((xf - posx(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_xr515 = (rho_r3/(2*pi*e_0)).*((xf - posx(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_xr516 = (rho_r4/(2*pi*e_0)).*((xf - posx(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_xr517 = (rho_r5/(2*pi*e_0)).*((xf - posx(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_xr518 = (rho_r6/(2*pi*e_0)).*((xf - posx(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));
E_xr519 = (rho_r7/(2*pi*e_0)).*((xf - posx(5,19))./((xf - posx(5,19)).^2 + (yf - posy(5,19)).^2));
E_xr520 = (rho_r8/(2*pi*e_0)).*((xf - posx(5,20))./((xf - posx(5,20)).^2 + (yf - posy(5,20)).^2));
E_xr521 = (rho_r9/(2*pi*e_0)).*((xf - posx(5,21))./((xf - posx(5,21)).^2 + (yf - posy(5,21)).^2));
E_xr522 = (rho_r10/(2*pi*e_0)).*((xf - posx(5,22))./((xf - posx(5,22)).^2 + (yf - posy(5,22)).^2));
E_xr523 = (rho_r11/(2*pi*e_0)).*((xf - posx(5,23))./((xf - posx(5,23)).^2 + (yf - posy(5,23)).^2));
E_xr524 = (rho_r12/(2*pi*e_0)).*((xf - posx(5,24))./((xf - posx(5,24)).^2 + (yf - posy(5,24)).^2));

E_xr61 = (-rho_r1/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xr62 = (-rho_r2/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xr63 = (-rho_r3/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xr64 = (-rho_r4/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xr65 = (-rho_r5/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_xr67 = (-rho_r7/(2*pi*e_0)).*((xf - posx(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_xr68 = (-rho_r8/(2*pi*e_0)).*((xf - posx(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_xr69 = (-rho_r9/(2*pi*e_0)).*((xf - posx(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_xr610 = (-rho_r10/(2*pi*e_0)).*((xf - posx(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_xr611 = (-rho_r11/(2*pi*e_0)).*((xf - posx(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_xr612 = (-rho_r12/(2*pi*e_0)).*((xf - posx(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_xr613 = (rho_r1/(2*pi*e_0)).*((xf - posx(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_xr614 = (rho_r2/(2*pi*e_0)).*((xf - posx(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_xr615 = (rho_r3/(2*pi*e_0)).*((xf - posx(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_xr616 = (rho_r4/(2*pi*e_0)).*((xf - posx(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_xr617 = (rho_r5/(2*pi*e_0)).*((xf - posx(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_xr618 = (rho_r6/(2*pi*e_0)).*((xf - posx(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));
E_xr619 = (rho_r7/(2*pi*e_0)).*((xf - posx(6,19))./((xf - posx(6,19)).^2 + (yf - posy(6,19)).^2));
E_xr620 = (rho_r8/(2*pi*e_0)).*((xf - posx(6,20))./((xf - posx(6,20)).^2 + (yf - posy(6,20)).^2));
E_xr621 = (rho_r9/(2*pi*e_0)).*((xf - posx(6,21))./((xf - posx(6,21)).^2 + (yf - posy(6,21)).^2));
E_xr622 = (rho_r10/(2*pi*e_0)).*((xf - posx(6,22))./((xf - posx(6,22)).^2 + (yf - posy(6,22)).^2));
E_xr623 = (rho_r11/(2*pi*e_0)).*((xf - posx(6,23))./((xf - posx(6,23)).^2 + (yf - posy(6,23)).^2));
E_xr624 = (rho_r12/(2*pi*e_0)).*((xf - posx(6,24))./((xf - posx(6,24)).^2 + (yf - posy(6,24)).^2));

E_xr71 = (-rho_r1/(2*pi*e_0)).*((xf - posx(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_xr72 = (-rho_r2/(2*pi*e_0)).*((xf - posx(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_xr73 = (-rho_r3/(2*pi*e_0)).*((xf - posx(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_xr74 = (-rho_r4/(2*pi*e_0)).*((xf - posx(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_xr75 = (-rho_r5/(2*pi*e_0)).*((xf - posx(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_xr76 = (-rho_r6/(2*pi*e_0)).*((xf - posx(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_xr78 = (-rho_r8/(2*pi*e_0)).*((xf - posx(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_xr79 = (-rho_r9/(2*pi*e_0)).*((xf - posx(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_xr710 = (-rho_r10/(2*pi*e_0)).*((xf - posx(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_xr711 = (-rho_r11/(2*pi*e_0)).*((xf - posx(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_xr712 = (-rho_r12/(2*pi*e_0)).*((xf - posx(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_xr713 = (rho_r1/(2*pi*e_0)).*((xf - posx(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_xr714 = (rho_r2/(2*pi*e_0)).*((xf - posx(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_xr715 = (rho_r3/(2*pi*e_0)).*((xf - posx(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_xr716 = (rho_r4/(2*pi*e_0)).*((xf - posx(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_xr717 = (rho_r5/(2*pi*e_0)).*((xf - posx(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_xr718 = (rho_r6/(2*pi*e_0)).*((xf - posx(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));
E_xr719 = (rho_r7/(2*pi*e_0)).*((xf - posx(7,19))./((xf - posx(7,19)).^2 + (yf - posy(7,19)).^2));
E_xr720 = (rho_r8/(2*pi*e_0)).*((xf - posx(7,20))./((xf - posx(7,20)).^2 + (yf - posy(7,20)).^2));
E_xr721 = (rho_r9/(2*pi*e_0)).*((xf - posx(7,21))./((xf - posx(7,21)).^2 + (yf - posy(7,21)).^2));
E_xr722 = (rho_r10/(2*pi*e_0)).*((xf - posx(7,22))./((xf - posx(7,22)).^2 + (yf - posy(7,22)).^2));
E_xr723 = (rho_r11/(2*pi*e_0)).*((xf - posx(7,23))./((xf - posx(7,23)).^2 + (yf - posy(7,23)).^2));
E_xr724 = (rho_r12/(2*pi*e_0)).*((xf - posx(7,24))./((xf - posx(7,24)).^2 + (yf - posy(7,24)).^2));

E_xr81 = (-rho_r1/(2*pi*e_0)).*((xf - posx(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_xr82 = (-rho_r2/(2*pi*e_0)).*((xf - posx(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_xr83 = (-rho_r3/(2*pi*e_0)).*((xf - posx(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_xr84 = (-rho_r4/(2*pi*e_0)).*((xf - posx(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_xr85 = (-rho_r5/(2*pi*e_0)).*((xf - posx(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_xr86 = (-rho_r6/(2*pi*e_0)).*((xf - posx(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_xr87 = (-rho_r7/(2*pi*e_0)).*((xf - posx(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_xr89 = (-rho_r9/(2*pi*e_0)).*((xf - posx(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_xr810 = (-rho_r10/(2*pi*e_0)).*((xf - posx(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_xr811 = (-rho_r11/(2*pi*e_0)).*((xf - posx(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_xr812 = (-rho_r12/(2*pi*e_0)).*((xf - posx(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_xr813 = (rho_r1/(2*pi*e_0)).*((xf - posx(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_xr814 = (rho_r2/(2*pi*e_0)).*((xf - posx(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_xr815 = (rho_r3/(2*pi*e_0)).*((xf - posx(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_xr816 = (rho_r4/(2*pi*e_0)).*((xf - posx(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_xr817 = (rho_r5/(2*pi*e_0)).*((xf - posx(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_xr818 = (rho_r6/(2*pi*e_0)).*((xf - posx(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));
E_xr819 = (rho_r7/(2*pi*e_0)).*((xf - posx(8,19))./((xf - posx(8,19)).^2 + (yf - posy(8,19)).^2));
E_xr820 = (rho_r8/(2*pi*e_0)).*((xf - posx(8,20))./((xf - posx(8,20)).^2 + (yf - posy(8,20)).^2));
E_xr821 = (rho_r9/(2*pi*e_0)).*((xf - posx(8,21))./((xf - posx(8,21)).^2 + (yf - posy(8,21)).^2));
E_xr822 = (rho_r10/(2*pi*e_0)).*((xf - posx(8,22))./((xf - posx(8,22)).^2 + (yf - posy(8,22)).^2));
E_xr823 = (rho_r11/(2*pi*e_0)).*((xf - posx(8,23))./((xf - posx(8,23)).^2 + (yf - posy(8,23)).^2));
E_xr824 = (rho_r12/(2*pi*e_0)).*((xf - posx(8,24))./((xf - posx(8,24)).^2 + (yf - posy(8,24)).^2));

E_xr91 = (-rho_r1/(2*pi*e_0)).*((xf - posx(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_xr92 = (-rho_r2/(2*pi*e_0)).*((xf - posx(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_xr93 = (-rho_r3/(2*pi*e_0)).*((xf - posx(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_xr94 = (-rho_r4/(2*pi*e_0)).*((xf - posx(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_xr95 = (-rho_r5/(2*pi*e_0)).*((xf - posx(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_xr96 = (-rho_r6/(2*pi*e_0)).*((xf - posx(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_xr97 = (-rho_r7/(2*pi*e_0)).*((xf - posx(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_xr98 = (-rho_r8/(2*pi*e_0)).*((xf - posx(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_xr910 = (-rho_r10/(2*pi*e_0)).*((xf - posx(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_xr911 = (-rho_r11/(2*pi*e_0)).*((xf - posx(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_xr912 = (-rho_r12/(2*pi*e_0)).*((xf - posx(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_xr913 = (rho_r1/(2*pi*e_0)).*((xf - posx(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_xr914 = (rho_r2/(2*pi*e_0)).*((xf - posx(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_xr915 = (rho_r3/(2*pi*e_0)).*((xf - posx(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_xr916 = (rho_r4/(2*pi*e_0)).*((xf - posx(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_xr917 = (rho_r5/(2*pi*e_0)).*((xf - posx(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_xr918 = (rho_r6/(2*pi*e_0)).*((xf - posx(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));
E_xr919 = (rho_r7/(2*pi*e_0)).*((xf - posx(9,19))./((xf - posx(9,19)).^2 + (yf - posy(9,19)).^2));
E_xr920 = (rho_r8/(2*pi*e_0)).*((xf - posx(9,20))./((xf - posx(9,20)).^2 + (yf - posy(9,20)).^2));
E_xr921 = (rho_r9/(2*pi*e_0)).*((xf - posx(9,21))./((xf - posx(9,21)).^2 + (yf - posy(9,21)).^2));
E_xr922 = (rho_r10/(2*pi*e_0)).*((xf - posx(9,22))./((xf - posx(9,22)).^2 + (yf - posy(9,22)).^2));
E_xr923 = (rho_r11/(2*pi*e_0)).*((xf - posx(9,23))./((xf - posx(9,23)).^2 + (yf - posy(9,23)).^2));
E_xr924 = (rho_r12/(2*pi*e_0)).*((xf - posx(9,24))./((xf - posx(9,24)).^2 + (yf - posy(9,24)).^2));

E_xr101 = (-rho_r1/(2*pi*e_0)).*((xf - posx(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_xr102 = (-rho_r2/(2*pi*e_0)).*((xf - posx(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_xr103 = (-rho_r3/(2*pi*e_0)).*((xf - posx(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_xr104 = (-rho_r4/(2*pi*e_0)).*((xf - posx(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_xr105 = (-rho_r5/(2*pi*e_0)).*((xf - posx(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_xr106 = (-rho_r6/(2*pi*e_0)).*((xf - posx(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_xr107 = (-rho_r7/(2*pi*e_0)).*((xf - posx(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_xr108 = (-rho_r8/(2*pi*e_0)).*((xf - posx(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_xr109 = (-rho_r9/(2*pi*e_0)).*((xf - posx(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_xr1011 = (-rho_r11/(2*pi*e_0)).*((xf - posx(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_xr1012 = (-rho_r12/(2*pi*e_0)).*((xf - posx(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_xr1013 = (rho_r1/(2*pi*e_0)).*((xf - posx(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_xr1014 = (rho_r2/(2*pi*e_0)).*((xf - posx(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_xr1015 = (rho_r3/(2*pi*e_0)).*((xf - posx(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_xr1016 = (rho_r4/(2*pi*e_0)).*((xf - posx(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_xr1017 = (rho_r5/(2*pi*e_0)).*((xf - posx(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_xr1018 = (rho_r6/(2*pi*e_0)).*((xf - posx(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));
E_xr10_19 = (rho_r7/(2*pi*e_0)).*((xf - posx(10,19))./((xf - posx(10,19)).^2 + (yf - posy(10,19)).^2));
E_xr10_20 = (rho_r8/(2*pi*e_0)).*((xf - posx(10,20))./((xf - posx(10,20)).^2 + (yf - posy(10,20)).^2));
E_xr10_21 = (rho_r9/(2*pi*e_0)).*((xf - posx(10,21))./((xf - posx(10,21)).^2 + (yf - posy(10,21)).^2));
E_xr10_22 = (rho_r10/(2*pi*e_0)).*((xf - posx(10,22))./((xf - posx(10,22)).^2 + (yf - posy(10,22)).^2));
E_xr10_23 = (rho_r11/(2*pi*e_0)).*((xf - posx(10,23))./((xf - posx(10,23)).^2 + (yf - posy(10,23)).^2));
E_xr10_24 = (rho_r12/(2*pi*e_0)).*((xf - posx(10,24))./((xf - posx(10,24)).^2 + (yf - posy(10,24)).^2));

E_xr11_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_xr11_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_xr113 = (-rho_r3/(2*pi*e_0)).*((xf - posx(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_xr114 = (-rho_r4/(2*pi*e_0)).*((xf - posx(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_xr115 = (-rho_r5/(2*pi*e_0)).*((xf - posx(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_xr116 = (-rho_r6/(2*pi*e_0)).*((xf - posx(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_xr117 = (-rho_r7/(2*pi*e_0)).*((xf - posx(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_xr118 = (-rho_r8/(2*pi*e_0)).*((xf - posx(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_xr119 = (-rho_r9/(2*pi*e_0)).*((xf - posx(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_xr1110 = (-rho_r10/(2*pi*e_0)).*((xf - posx(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_xr1112 = (-rho_r12/(2*pi*e_0)).*((xf - posx(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_xr1113 = (rho_r1/(2*pi*e_0)).*((xf - posx(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_xr1114 = (rho_r2/(2*pi*e_0)).*((xf - posx(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_xr1115 = (rho_r3/(2*pi*e_0)).*((xf - posx(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_xr1116 = (rho_r4/(2*pi*e_0)).*((xf - posx(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_xr1117 = (rho_r5/(2*pi*e_0)).*((xf - posx(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_xr1118 = (rho_r6/(2*pi*e_0)).*((xf - posx(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));
E_xr1019 = (rho_r7/(2*pi*e_0)).*((xf - posx(11,19))./((xf - posx(11,19)).^2 + (yf - posy(11,19)).^2));
E_xr1020 = (rho_r8/(2*pi*e_0)).*((xf - posx(11,20))./((xf - posx(11,20)).^2 + (yf - posy(11,20)).^2));
E_xr1021 = (rho_r9/(2*pi*e_0)).*((xf - posx(11,21))./((xf - posx(11,21)).^2 + (yf - posy(11,21)).^2));
E_xr1022 = (rho_r10/(2*pi*e_0)).*((xf - posx(11,22))./((xf - posx(11,22)).^2 + (yf - posy(11,22)).^2));
E_xr1023 = (rho_r11/(2*pi*e_0)).*((xf - posx(11,23))./((xf - posx(11,23)).^2 + (yf - posy(11,23)).^2));
E_xr1024 = (rho_r12/(2*pi*e_0)).*((xf - posx(11,24))./((xf - posx(11,24)).^2 + (yf - posy(11,24)).^2));

E_xr12_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_xr12_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_xr123 = (-rho_r3/(2*pi*e_0)).*((xf - posx(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_xr124 = (-rho_r4/(2*pi*e_0)).*((xf - posx(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_xr125 = (-rho_r5/(2*pi*e_0)).*((xf - posx(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_xr126 = (-rho_r6/(2*pi*e_0)).*((xf - posx(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_xr127 = (-rho_r7/(2*pi*e_0)).*((xf - posx(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_xr128 = (-rho_r8/(2*pi*e_0)).*((xf - posx(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_xr129 = (-rho_r9/(2*pi*e_0)).*((xf - posx(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_xr1210 = (-rho_r10/(2*pi*e_0)).*((xf - posx(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_xr1211 = (-rho_r11/(2*pi*e_0)).*((xf - posx(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_xr1213 = (rho_r1/(2*pi*e_0)).*((xf - posx(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_xr1214 = (rho_r2/(2*pi*e_0)).*((xf - posx(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_xr1215 = (rho_r3/(2*pi*e_0)).*((xf - posx(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_xr1216 = (rho_r4/(2*pi*e_0)).*((xf - posx(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_xr1217 = (rho_r5/(2*pi*e_0)).*((xf - posx(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_xr1218 = (rho_r6/(2*pi*e_0)).*((xf - posx(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));
E_xr1219 = (rho_r7/(2*pi*e_0)).*((xf - posx(12,19))./((xf - posx(12,19)).^2 + (yf - posy(12,19)).^2));
E_xr1220 = (rho_r8/(2*pi*e_0)).*((xf - posx(12,20))./((xf - posx(12,20)).^2 + (yf - posy(12,20)).^2));
E_xr1221 = (rho_r9/(2*pi*e_0)).*((xf - posx(12,21))./((xf - posx(12,21)).^2 + (yf - posy(12,21)).^2));
E_xr1222 = (rho_r10/(2*pi*e_0)).*((xf - posx(12,22))./((xf - posx(12,22)).^2 + (yf - posy(12,22)).^2));
E_xr1223 = (rho_r11/(2*pi*e_0)).*((xf - posx(12,23))./((xf - posx(12,23)).^2 + (yf - posy(12,23)).^2));
E_xr1224 = (rho_r12/(2*pi*e_0)).*((xf - posx(12,24))./((xf - posx(12,24)).^2 + (yf - posy(12,24)).^2));

E_xr131 = (-rho_r1/(2*pi*e_0)).*((xf - posx(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_xr132 = (-rho_r2/(2*pi*e_0)).*((xf - posx(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_xr133 = (-rho_r3/(2*pi*e_0)).*((xf - posx(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_xr134 = (-rho_r4/(2*pi*e_0)).*((xf - posx(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_xr135 = (-rho_r5/(2*pi*e_0)).*((xf - posx(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_xr136 = (-rho_r6/(2*pi*e_0)).*((xf - posx(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_xr137 = (-rho_r7/(2*pi*e_0)).*((xf - posx(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_xr138 = (-rho_r8/(2*pi*e_0)).*((xf - posx(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_xr139 = (-rho_r9/(2*pi*e_0)).*((xf - posx(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_xr1310 = (-rho_r10/(2*pi*e_0)).*((xf - posx(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_xr1311 = (-rho_r11/(2*pi*e_0)).*((xf - posx(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_xr1312 = (-rho_r12/(2*pi*e_0)).*((xf - posx(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_xr1314 = (rho_r2/(2*pi*e_0)).*((xf - posx(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_xr1315 = (rho_r3/(2*pi*e_0)).*((xf - posx(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_xr1316 = (rho_r4/(2*pi*e_0)).*((xf - posx(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_xr1317 = (rho_r5/(2*pi*e_0)).*((xf - posx(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_xr1318 = (rho_r6/(2*pi*e_0)).*((xf - posx(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));
E_xr1319 = (rho_r7/(2*pi*e_0)).*((xf - posx(13,19))./((xf - posx(13,19)).^2 + (yf - posy(13,19)).^2));
E_xr1320 = (rho_r8/(2*pi*e_0)).*((xf - posx(13,20))./((xf - posx(13,20)).^2 + (yf - posy(13,20)).^2));
E_xr1321 = (rho_r9/(2*pi*e_0)).*((xf - posx(13,21))./((xf - posx(13,21)).^2 + (yf - posy(13,21)).^2));
E_xr1322 = (rho_r10/(2*pi*e_0)).*((xf - posx(13,22))./((xf - posx(13,22)).^2 + (yf - posy(13,22)).^2));
E_xr1323 = (rho_r11/(2*pi*e_0)).*((xf - posx(13,23))./((xf - posx(13,23)).^2 + (yf - posy(13,23)).^2));
E_xr1324 = (rho_r12/(2*pi*e_0)).*((xf - posx(13,24))./((xf - posx(13,24)).^2 + (yf - posy(13,24)).^2));

E_xr141 = (-rho_r1/(2*pi*e_0)).*((xf - posx(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_xr142 = (-rho_r2/(2*pi*e_0)).*((xf - posx(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_xr143 = (-rho_r3/(2*pi*e_0)).*((xf - posx(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_xr144 = (-rho_r4/(2*pi*e_0)).*((xf - posx(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_xr145 = (-rho_r5/(2*pi*e_0)).*((xf - posx(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_xr146 = (-rho_r6/(2*pi*e_0)).*((xf - posx(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_xr147 = (-rho_r7/(2*pi*e_0)).*((xf - posx(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_xr148 = (-rho_r8/(2*pi*e_0)).*((xf - posx(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_xr149 = (-rho_r9/(2*pi*e_0)).*((xf - posx(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_xr1410 = (-rho_r10/(2*pi*e_0)).*((xf - posx(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_xr1411 = (-rho_r11/(2*pi*e_0)).*((xf - posx(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_xr1412 = (-rho_r12/(2*pi*e_0)).*((xf - posx(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_xr1413 = (rho_r1/(2*pi*e_0)).*((xf - posx(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_xr1415 = (rho_r3/(2*pi*e_0)).*((xf - posx(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_xr1416 = (rho_r4/(2*pi*e_0)).*((xf - posx(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_xr1417 = (rho_r5/(2*pi*e_0)).*((xf - posx(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_xr1418 = (rho_r6/(2*pi*e_0)).*((xf - posx(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));
E_xr1419 = (rho_r7/(2*pi*e_0)).*((xf - posx(14,19))./((xf - posx(14,19)).^2 + (yf - posy(14,19)).^2));
E_xr1420 = (rho_r8/(2*pi*e_0)).*((xf - posx(14,20))./((xf - posx(14,20)).^2 + (yf - posy(14,20)).^2));
E_xr1421 = (rho_r9/(2*pi*e_0)).*((xf - posx(14,21))./((xf - posx(14,21)).^2 + (yf - posy(14,21)).^2));
E_xr1422 = (rho_r10/(2*pi*e_0)).*((xf - posx(14,22))./((xf - posx(14,22)).^2 + (yf - posy(14,22)).^2));
E_xr1423 = (rho_r11/(2*pi*e_0)).*((xf - posx(14,23))./((xf - posx(14,23)).^2 + (yf - posy(14,23)).^2));
E_xr1424 = (rho_r12/(2*pi*e_0)).*((xf - posx(14,24))./((xf - posx(14,24)).^2 + (yf - posy(14,24)).^2));

E_xr151 = (-rho_r1/(2*pi*e_0)).*((xf - posx(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_xr152 = (-rho_r2/(2*pi*e_0)).*((xf - posx(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_xr153 = (-rho_r3/(2*pi*e_0)).*((xf - posx(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_xr154 = (-rho_r4/(2*pi*e_0)).*((xf - posx(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_xr155 = (-rho_r5/(2*pi*e_0)).*((xf - posx(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_xr156 = (-rho_r6/(2*pi*e_0)).*((xf - posx(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_xr157 = (-rho_r7/(2*pi*e_0)).*((xf - posx(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_xr158 = (-rho_r8/(2*pi*e_0)).*((xf - posx(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_xr159 = (-rho_r9/(2*pi*e_0)).*((xf - posx(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_xr1510 = (-rho_r10/(2*pi*e_0)).*((xf - posx(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_xr1511 = (-rho_r11/(2*pi*e_0)).*((xf - posx(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_xr1512 = (-rho_r12/(2*pi*e_0)).*((xf - posx(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_xr1513 = (rho_r1/(2*pi*e_0)).*((xf - posx(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_xr1514 = (rho_r2/(2*pi*e_0)).*((xf - posx(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_xr1516 = (rho_r4/(2*pi*e_0)).*((xf - posx(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_xr1517 = (rho_r5/(2*pi*e_0)).*((xf - posx(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_xr1518 = (rho_r6/(2*pi*e_0)).*((xf - posx(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));
E_xr1519 = (rho_r7/(2*pi*e_0)).*((xf - posx(15,19))./((xf - posx(15,19)).^2 + (yf - posy(15,19)).^2));
E_xr1520 = (rho_r8/(2*pi*e_0)).*((xf - posx(15,20))./((xf - posx(15,20)).^2 + (yf - posy(15,20)).^2));
E_xr1521 = (rho_r9/(2*pi*e_0)).*((xf - posx(15,21))./((xf - posx(15,21)).^2 + (yf - posy(15,21)).^2));
E_xr1522 = (rho_r10/(2*pi*e_0)).*((xf - posx(15,22))./((xf - posx(15,22)).^2 + (yf - posy(15,22)).^2));
E_xr1523 = (rho_r11/(2*pi*e_0)).*((xf - posx(15,23))./((xf - posx(15,23)).^2 + (yf - posy(15,23)).^2));
E_xr1524 = (rho_r12/(2*pi*e_0)).*((xf - posx(15,24))./((xf - posx(15,24)).^2 + (yf - posy(15,24)).^2));

E_xr161 = (-rho_r1/(2*pi*e_0)).*((xf - posx(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_xr162 = (-rho_r2/(2*pi*e_0)).*((xf - posx(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_xr163 = (-rho_r3/(2*pi*e_0)).*((xf - posx(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_xr164 = (-rho_r4/(2*pi*e_0)).*((xf - posx(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_xr165 = (-rho_r5/(2*pi*e_0)).*((xf - posx(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_xr166 = (-rho_r6/(2*pi*e_0)).*((xf - posx(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_xr167 = (-rho_r7/(2*pi*e_0)).*((xf - posx(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_xr168 = (-rho_r8/(2*pi*e_0)).*((xf - posx(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_xr169 = (-rho_r9/(2*pi*e_0)).*((xf - posx(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_xr1610 = (-rho_r10/(2*pi*e_0)).*((xf - posx(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_xr1611 = (-rho_r11/(2*pi*e_0)).*((xf - posx(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_xr1612 = (-rho_r12/(2*pi*e_0)).*((xf - posx(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_xr1613 = (rho_r1/(2*pi*e_0)).*((xf - posx(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_xr1614 = (rho_r2/(2*pi*e_0)).*((xf - posx(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_xr1615 = (rho_r3/(2*pi*e_0)).*((xf - posx(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_xr1617 = (rho_r5/(2*pi*e_0)).*((xf - posx(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_xr1618 = (rho_r6/(2*pi*e_0)).*((xf - posx(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));
E_xr1619 = (rho_r7/(2*pi*e_0)).*((xf - posx(16,19))./((xf - posx(16,19)).^2 + (yf - posy(16,19)).^2));
E_xr1620 = (rho_r8/(2*pi*e_0)).*((xf - posx(16,20))./((xf - posx(16,20)).^2 + (yf - posy(16,20)).^2));
E_xr1621 = (rho_r9/(2*pi*e_0)).*((xf - posx(16,21))./((xf - posx(16,21)).^2 + (yf - posy(16,21)).^2));
E_xr1622 = (rho_r10/(2*pi*e_0)).*((xf - posx(16,22))./((xf - posx(16,22)).^2 + (yf - posy(16,22)).^2));
E_xr1623 = (rho_r11/(2*pi*e_0)).*((xf - posx(16,23))./((xf - posx(16,23)).^2 + (yf - posy(16,23)).^2));
E_xr1624 = (rho_r12/(2*pi*e_0)).*((xf - posx(16,24))./((xf - posx(16,24)).^2 + (yf - posy(16,24)).^2));

E_xr171 = (-rho_r1/(2*pi*e_0)).*((xf - posx(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_xr172 = (-rho_r2/(2*pi*e_0)).*((xf - posx(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_xr173 = (-rho_r3/(2*pi*e_0)).*((xf - posx(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_xr174 = (-rho_r4/(2*pi*e_0)).*((xf - posx(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_xr175 = (-rho_r5/(2*pi*e_0)).*((xf - posx(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_xr176 = (-rho_r6/(2*pi*e_0)).*((xf - posx(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_xr177 = (-rho_r7/(2*pi*e_0)).*((xf - posx(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_xr178 = (-rho_r8/(2*pi*e_0)).*((xf - posx(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_xr179 = (-rho_r9/(2*pi*e_0)).*((xf - posx(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_xr1710 = (-rho_r10/(2*pi*e_0)).*((xf - posx(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_xr1711 = (-rho_r11/(2*pi*e_0)).*((xf - posx(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_xr1712 = (-rho_r12/(2*pi*e_0)).*((xf - posx(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_xr1713 = (rho_r1/(2*pi*e_0)).*((xf - posx(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_xr1714 = (rho_r2/(2*pi*e_0)).*((xf - posx(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_xr1715 = (rho_r3/(2*pi*e_0)).*((xf - posx(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_xr1716 = (rho_r4/(2*pi*e_0)).*((xf - posx(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_xr1718 = (rho_r6/(2*pi*e_0)).*((xf - posx(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));
E_xr1719 = (rho_r7/(2*pi*e_0)).*((xf - posx(17,19))./((xf - posx(17,19)).^2 + (yf - posy(17,19)).^2));
E_xr1720 = (rho_r8/(2*pi*e_0)).*((xf - posx(17,20))./((xf - posx(17,20)).^2 + (yf - posy(17,20)).^2));
E_xr1721 = (rho_r9/(2*pi*e_0)).*((xf - posx(17,21))./((xf - posx(17,21)).^2 + (yf - posy(17,21)).^2));
E_xr1722 = (rho_r10/(2*pi*e_0)).*((xf - posx(17,22))./((xf - posx(17,22)).^2 + (yf - posy(17,22)).^2));
E_xr1723 = (rho_r11/(2*pi*e_0)).*((xf - posx(17,23))./((xf - posx(17,23)).^2 + (yf - posy(17,23)).^2));
E_xr1724 = (rho_r12/(2*pi*e_0)).*((xf - posx(17,24))./((xf - posx(17,24)).^2 + (yf - posy(17,24)).^2));

E_xr181 = (-rho_r1/(2*pi*e_0)).*((xf - posx(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_xr182 = (-rho_r2/(2*pi*e_0)).*((xf - posx(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_xr183 = (-rho_r3/(2*pi*e_0)).*((xf - posx(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_xr184 = (-rho_r4/(2*pi*e_0)).*((xf - posx(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_xr185 = (-rho_r5/(2*pi*e_0)).*((xf - posx(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_xr186 = (-rho_r6/(2*pi*e_0)).*((xf - posx(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_xr187 = (-rho_r7/(2*pi*e_0)).*((xf - posx(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_xr188 = (-rho_r8/(2*pi*e_0)).*((xf - posx(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_xr189 = (-rho_r9/(2*pi*e_0)).*((xf - posx(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_xr1810 = (-rho_r10/(2*pi*e_0)).*((xf - posx(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_xr1811 = (-rho_r11/(2*pi*e_0)).*((xf - posx(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_xr1812 = (-rho_r12/(2*pi*e_0)).*((xf - posx(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_xr1813 = (rho_r1/(2*pi*e_0)).*((xf - posx(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_xr1814 = (rho_r2/(2*pi*e_0)).*((xf - posx(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_xr1815 = (rho_r3/(2*pi*e_0)).*((xf - posx(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_xr1816 = (rho_r4/(2*pi*e_0)).*((xf - posx(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_xr1817 = (rho_r5/(2*pi*e_0)).*((xf - posx(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));
E_xr1819 = (rho_r7/(2*pi*e_0)).*((xf - posx(18,19))./((xf - posx(18,19)).^2 + (yf - posy(18,19)).^2));
E_xr1820 = (rho_r8/(2*pi*e_0)).*((xf - posx(18,20))./((xf - posx(18,20)).^2 + (yf - posy(18,20)).^2));
E_xr1821 = (rho_r9/(2*pi*e_0)).*((xf - posx(18,21))./((xf - posx(18,21)).^2 + (yf - posy(18,21)).^2));
E_xr1822 = (rho_r10/(2*pi*e_0)).*((xf - posx(18,22))./((xf - posx(18,22)).^2 + (yf - posy(18,22)).^2));
E_xr1823 = (rho_r11/(2*pi*e_0)).*((xf - posx(18,23))./((xf - posx(18,23)).^2 + (yf - posy(18,23)).^2));
E_xr1824 = (rho_r12/(2*pi*e_0)).*((xf - posx(18,24))./((xf - posx(18,24)).^2 + (yf - posy(18,24)).^2));

E_xr191 = (-rho_r1/(2*pi*e_0)).*((xf - posx(19,1))./((xf - posx(19,1)).^2 + (yf - posy(19,1)).^2));
E_xr192 = (-rho_r2/(2*pi*e_0)).*((xf - posx(19,2))./((xf - posx(19,2)).^2 + (yf - posy(19,2)).^2));
E_xr193 = (-rho_r3/(2*pi*e_0)).*((xf - posx(19,3))./((xf - posx(19,3)).^2 + (yf - posy(19,3)).^2));
E_xr194 = (-rho_r4/(2*pi*e_0)).*((xf - posx(19,4))./((xf - posx(19,4)).^2 + (yf - posy(19,4)).^2));
E_xr195 = (-rho_r5/(2*pi*e_0)).*((xf - posx(19,5))./((xf - posx(19,5)).^2 + (yf - posy(19,5)).^2));
E_xr196 = (-rho_r6/(2*pi*e_0)).*((xf - posx(19,6))./((xf - posx(19,6)).^2 + (yf - posy(19,6)).^2));
E_xr197 = (-rho_r7/(2*pi*e_0)).*((xf - posx(19,7))./((xf - posx(19,7)).^2 + (yf - posy(19,7)).^2));
E_xr198 = (-rho_r8/(2*pi*e_0)).*((xf - posx(19,8))./((xf - posx(19,8)).^2 + (yf - posy(19,8)).^2));
E_xr199 = (-rho_r9/(2*pi*e_0)).*((xf - posx(19,9))./((xf - posx(19,9)).^2 + (yf - posy(19,9)).^2));
E_xr1910 = (-rho_r10/(2*pi*e_0)).*((xf - posx(19,10))./((xf - posx(19,10)).^2 + (yf - posy(19,10)).^2));
E_xr1911 = (-rho_r11/(2*pi*e_0)).*((xf - posx(19,11))./((xf - posx(19,11)).^2 + (yf - posy(19,11)).^2));
E_xr1912 = (-rho_r12/(2*pi*e_0)).*((xf - posx(19,12))./((xf - posx(19,12)).^2 + (yf - posy(19,12)).^2));
E_xr1913 = (rho_r1/(2*pi*e_0)).*((xf - posx(19,13))./((xf - posx(19,13)).^2 + (yf - posy(19,13)).^2));
E_xr1914 = (rho_r2/(2*pi*e_0)).*((xf - posx(19,14))./((xf - posx(19,14)).^2 + (yf - posy(19,14)).^2));
E_xr1915 = (rho_r3/(2*pi*e_0)).*((xf - posx(19,15))./((xf - posx(19,15)).^2 + (yf - posy(19,15)).^2));
E_xr1916 = (rho_r4/(2*pi*e_0)).*((xf - posx(19,16))./((xf - posx(19,16)).^2 + (yf - posy(19,16)).^2));
E_xr1917 = (rho_r5/(2*pi*e_0)).*((xf - posx(19,17))./((xf - posx(19,17)).^2 + (yf - posy(19,17)).^2));
E_xr1918 = (rho_r6/(2*pi*e_0)).*((xf - posx(19,18))./((xf - posx(19,18)).^2 + (yf - posy(19,18)).^2));
E_xr1920 = (rho_r8/(2*pi*e_0)).*((xf - posx(19,20))./((xf - posx(19,20)).^2 + (yf - posy(19,20)).^2));
E_xr1921 = (rho_r9/(2*pi*e_0)).*((xf - posx(19,21))./((xf - posx(19,21)).^2 + (yf - posy(19,21)).^2));
E_xr1922 = (rho_r10/(2*pi*e_0)).*((xf - posx(19,22))./((xf - posx(19,22)).^2 + (yf - posy(19,22)).^2));
E_xr1923 = (rho_r11/(2*pi*e_0)).*((xf - posx(19,23))./((xf - posx(19,23)).^2 + (yf - posy(19,23)).^2));
E_xr1924 = (rho_r12/(2*pi*e_0)).*((xf - posx(19,24))./((xf - posx(19,24)).^2 + (yf - posy(19,24)).^2));

E_xr201 = (-rho_r1/(2*pi*e_0)).*((xf - posx(20,1))./((xf - posx(20,1)).^2 + (yf - posy(20,1)).^2));
E_xr202 = (-rho_r2/(2*pi*e_0)).*((xf - posx(20,2))./((xf - posx(20,2)).^2 + (yf - posy(20,2)).^2));
E_xr203 = (-rho_r3/(2*pi*e_0)).*((xf - posx(20,3))./((xf - posx(20,3)).^2 + (yf - posy(20,3)).^2));
E_xr204 = (-rho_r4/(2*pi*e_0)).*((xf - posx(20,4))./((xf - posx(20,4)).^2 + (yf - posy(20,4)).^2));
E_xr205 = (-rho_r5/(2*pi*e_0)).*((xf - posx(20,5))./((xf - posx(20,5)).^2 + (yf - posy(20,5)).^2));
E_xr206 = (-rho_r6/(2*pi*e_0)).*((xf - posx(20,6))./((xf - posx(20,6)).^2 + (yf - posy(20,6)).^2));
E_xr207 = (-rho_r7/(2*pi*e_0)).*((xf - posx(20,7))./((xf - posx(20,7)).^2 + (yf - posy(20,7)).^2));
E_xr208 = (-rho_r8/(2*pi*e_0)).*((xf - posx(20,8))./((xf - posx(20,8)).^2 + (yf - posy(20,8)).^2));
E_xr209 = (-rho_r9/(2*pi*e_0)).*((xf - posx(20,9))./((xf - posx(20,9)).^2 + (yf - posy(20,9)).^2));
E_xr2010 = (-rho_r10/(2*pi*e_0)).*((xf - posx(20,10))./((xf - posx(20,10)).^2 + (yf - posy(20,10)).^2));
E_xr2011 = (-rho_r11/(2*pi*e_0)).*((xf - posx(20,11))./((xf - posx(20,11)).^2 + (yf - posy(20,11)).^2));
E_xr2012 = (-rho_r12/(2*pi*e_0)).*((xf - posx(20,12))./((xf - posx(20,12)).^2 + (yf - posy(20,12)).^2));
E_xr2013 = (rho_r1/(2*pi*e_0)).*((xf - posx(20,13))./((xf - posx(20,13)).^2 + (yf - posy(20,13)).^2));
E_xr2014 = (rho_r2/(2*pi*e_0)).*((xf - posx(20,14))./((xf - posx(20,14)).^2 + (yf - posy(20,14)).^2));
E_xr2015 = (rho_r3/(2*pi*e_0)).*((xf - posx(20,15))./((xf - posx(20,15)).^2 + (yf - posy(20,15)).^2));
E_xr2016 = (rho_r4/(2*pi*e_0)).*((xf - posx(20,16))./((xf - posx(20,16)).^2 + (yf - posy(20,16)).^2));
E_xr2017 = (rho_r5/(2*pi*e_0)).*((xf - posx(20,17))./((xf - posx(20,17)).^2 + (yf - posy(20,17)).^2));
E_xr2018 = (rho_r6/(2*pi*e_0)).*((xf - posx(20,18))./((xf - posx(20,18)).^2 + (yf - posy(20,18)).^2));
E_xr2019 = (rho_r7/(2*pi*e_0)).*((xf - posx(20,19))./((xf - posx(20,19)).^2 + (yf - posy(20,19)).^2));
E_xr2021 = (rho_r9/(2*pi*e_0)).*((xf - posx(20,21))./((xf - posx(20,21)).^2 + (yf - posy(20,21)).^2));
E_xr2022 = (rho_r10/(2*pi*e_0)).*((xf - posx(20,22))./((xf - posx(20,22)).^2 + (yf - posy(20,22)).^2));
E_xr2023 = (rho_r11/(2*pi*e_0)).*((xf - posx(20,23))./((xf - posx(20,23)).^2 + (yf - posy(20,23)).^2));
E_xr2024 = (rho_r12/(2*pi*e_0)).*((xf - posx(20,24))./((xf - posx(20,24)).^2 + (yf - posy(20,24)).^2));

E_xr21_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(21,1))./((xf - posx(21,1)).^2 + (yf - posy(21,1)).^2));
E_xr21_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(21,2))./((xf - posx(21,2)).^2 + (yf - posy(21,2)).^2));
E_xr21_3 = (-rho_r3/(2*pi*e_0)).*((xf - posx(21,3))./((xf - posx(21,3)).^2 + (yf - posy(21,3)).^2));
E_xr214 = (-rho_r4/(2*pi*e_0)).*((xf - posx(21,4))./((xf - posx(21,4)).^2 + (yf - posy(21,4)).^2));
E_xr215 = (-rho_r5/(2*pi*e_0)).*((xf - posx(21,5))./((xf - posx(21,5)).^2 + (yf - posy(21,5)).^2));
E_xr216 = (-rho_r6/(2*pi*e_0)).*((xf - posx(21,6))./((xf - posx(21,6)).^2 + (yf - posy(21,6)).^2));
E_xr217 = (-rho_r7/(2*pi*e_0)).*((xf - posx(21,7))./((xf - posx(21,7)).^2 + (yf - posy(21,7)).^2));
E_xr218 = (-rho_r8/(2*pi*e_0)).*((xf - posx(21,8))./((xf - posx(21,8)).^2 + (yf - posy(21,8)).^2));
E_xr219 = (-rho_r9/(2*pi*e_0)).*((xf - posx(21,9))./((xf - posx(21,9)).^2 + (yf - posy(21,9)).^2));
E_xr2110 = (-rho_r10/(2*pi*e_0)).*((xf - posx(21,10))./((xf - posx(21,10)).^2 + (yf - posy(21,10)).^2));
E_xr2111 = (-rho_r11/(2*pi*e_0)).*((xf - posx(21,11))./((xf - posx(21,11)).^2 + (yf - posy(21,11)).^2));
E_xr2112 = (-rho_r12/(2*pi*e_0)).*((xf - posx(21,12))./((xf - posx(21,12)).^2 + (yf - posy(21,12)).^2));
E_xr2113 = (rho_r1/(2*pi*e_0)).*((xf - posx(21,13))./((xf - posx(21,13)).^2 + (yf - posy(21,13)).^2));
E_xr2114 = (rho_r2/(2*pi*e_0)).*((xf - posx(21,14))./((xf - posx(21,14)).^2 + (yf - posy(21,14)).^2));
E_xr2115 = (rho_r3/(2*pi*e_0)).*((xf - posx(21,15))./((xf - posx(21,15)).^2 + (yf - posy(21,15)).^2));
E_xr2116 = (rho_r4/(2*pi*e_0)).*((xf - posx(21,16))./((xf - posx(21,16)).^2 + (yf - posy(21,16)).^2));
E_xr2117 = (rho_r5/(2*pi*e_0)).*((xf - posx(21,17))./((xf - posx(21,17)).^2 + (yf - posy(21,17)).^2));
E_xr2118 = (rho_r6/(2*pi*e_0)).*((xf - posx(21,18))./((xf - posx(21,18)).^2 + (yf - posy(21,18)).^2));
E_xr2119 = (rho_r7/(2*pi*e_0)).*((xf - posx(21,19))./((xf - posx(21,19)).^2 + (yf - posy(21,19)).^2));
E_xr2120 = (rho_r8/(2*pi*e_0)).*((xf - posx(21,20))./((xf - posx(21,20)).^2 + (yf - posy(21,20)).^2));
E_xr2122 = (rho_r10/(2*pi*e_0)).*((xf - posx(21,22))./((xf - posx(21,22)).^2 + (yf - posy(21,22)).^2));
E_xr2123 = (rho_r11/(2*pi*e_0)).*((xf - posx(21,23))./((xf - posx(21,23)).^2 + (yf - posy(21,23)).^2));
E_xr2124 = (rho_r12/(2*pi*e_0)).*((xf - posx(21,24))./((xf - posx(21,24)).^2 + (yf - posy(21,24)).^2));

E_xr22_1 = (-rho_r1/(2*pi*e_0)).*((xf - posx(22,1))./((xf - posx(22,1)).^2 + (yf - posy(22,1)).^2));
E_xr22_2 = (-rho_r2/(2*pi*e_0)).*((xf - posx(22,2))./((xf - posx(22,2)).^2 + (yf - posy(22,2)).^2));
E_xr22_3 = (-rho_r3/(2*pi*e_0)).*((xf - posx(22,3))./((xf - posx(22,3)).^2 + (yf - posy(22,3)).^2));
E_xr22_4 = (-rho_r4/(2*pi*e_0)).*((xf - posx(22,4))./((xf - posx(22,4)).^2 + (yf - posy(22,4)).^2));
E_xr225 = (-rho_r5/(2*pi*e_0)).*((xf - posx(22,5))./((xf - posx(22,5)).^2 + (yf - posy(22,5)).^2));
E_xr226 = (-rho_r6/(2*pi*e_0)).*((xf - posx(22,6))./((xf - posx(22,6)).^2 + (yf - posy(22,6)).^2));
E_xr227 = (-rho_r7/(2*pi*e_0)).*((xf - posx(22,7))./((xf - posx(22,7)).^2 + (yf - posy(22,7)).^2));
E_xr228 = (-rho_r8/(2*pi*e_0)).*((xf - posx(22,8))./((xf - posx(22,8)).^2 + (yf - posy(22,8)).^2));
E_xr229 = (-rho_r9/(2*pi*e_0)).*((xf - posx(22,9))./((xf - posx(22,9)).^2 + (yf - posy(22,9)).^2));
E_xr2210 = (-rho_r10/(2*pi*e_0)).*((xf - posx(22,10))./((xf - posx(22,10)).^2 + (yf - posy(22,10)).^2));
E_xr2211 = (-rho_r11/(2*pi*e_0)).*((xf - posx(22,11))./((xf - posx(22,11)).^2 + (yf - posy(22,11)).^2));
E_xr2212 = (-rho_r12/(2*pi*e_0)).*((xf - posx(22,12))./((xf - posx(22,12)).^2 + (yf - posy(22,12)).^2));
E_xr2213 = (rho_r1/(2*pi*e_0)).*((xf - posx(22,13))./((xf - posx(22,13)).^2 + (yf - posy(22,13)).^2));
E_xr2214 = (rho_r2/(2*pi*e_0)).*((xf - posx(22,14))./((xf - posx(22,14)).^2 + (yf - posy(22,14)).^2));
E_xr2215 = (rho_r3/(2*pi*e_0)).*((xf - posx(22,15))./((xf - posx(22,15)).^2 + (yf - posy(22,15)).^2));
E_xr2216 = (rho_r4/(2*pi*e_0)).*((xf - posx(22,16))./((xf - posx(22,16)).^2 + (yf - posy(22,16)).^2));
E_xr2217 = (rho_r5/(2*pi*e_0)).*((xf - posx(22,17))./((xf - posx(22,17)).^2 + (yf - posy(22,17)).^2));
E_xr2218 = (rho_r6/(2*pi*e_0)).*((xf - posx(22,18))./((xf - posx(22,18)).^2 + (yf - posy(22,18)).^2));
E_xr2219 = (rho_r7/(2*pi*e_0)).*((xf - posx(22,19))./((xf - posx(22,19)).^2 + (yf - posy(22,19)).^2));
E_xr2220 = (rho_r8/(2*pi*e_0)).*((xf - posx(22,20))./((xf - posx(22,20)).^2 + (yf - posy(22,20)).^2));
E_xr2221 = (rho_r9/(2*pi*e_0)).*((xf - posx(22,21))./((xf - posx(22,21)).^2 + (yf - posy(22,21)).^2));
E_xr2223 = (rho_r11/(2*pi*e_0)).*((xf - posx(22,23))./((xf - posx(22,23)).^2 + (yf - posy(22,23)).^2));
E_xr2224 = (rho_r12/(2*pi*e_0)).*((xf - posx(22,24))./((xf - posx(22,24)).^2 + (yf - posy(22,24)).^2));

E_xr231 = (-rho_r1/(2*pi*e_0)).*((xf - posx(23,1))./((xf - posx(23,1)).^2 + (yf - posy(23,1)).^2));
E_xr232 = (-rho_r2/(2*pi*e_0)).*((xf - posx(23,2))./((xf - posx(23,2)).^2 + (yf - posy(23,2)).^2));
E_xr233 = (-rho_r3/(2*pi*e_0)).*((xf - posx(23,3))./((xf - posx(23,3)).^2 + (yf - posy(23,3)).^2));
E_xr234 = (-rho_r4/(2*pi*e_0)).*((xf - posx(23,4))./((xf - posx(23,4)).^2 + (yf - posy(23,4)).^2));
E_xr235 = (-rho_r5/(2*pi*e_0)).*((xf - posx(23,5))./((xf - posx(23,5)).^2 + (yf - posy(23,5)).^2));
E_xr236 = (-rho_r6/(2*pi*e_0)).*((xf - posx(23,6))./((xf - posx(23,6)).^2 + (yf - posy(23,6)).^2));
E_xr237 = (-rho_r7/(2*pi*e_0)).*((xf - posx(23,7))./((xf - posx(23,7)).^2 + (yf - posy(23,7)).^2));
E_xr238 = (-rho_r8/(2*pi*e_0)).*((xf - posx(23,8))./((xf - posx(23,8)).^2 + (yf - posy(23,8)).^2));
E_xr239 = (-rho_r9/(2*pi*e_0)).*((xf - posx(23,9))./((xf - posx(23,9)).^2 + (yf - posy(23,9)).^2));
E_xr2310 = (-rho_r10/(2*pi*e_0)).*((xf - posx(23,10))./((xf - posx(23,10)).^2 + (yf - posy(23,10)).^2));
E_xr2311 = (-rho_r11/(2*pi*e_0)).*((xf - posx(23,11))./((xf - posx(23,11)).^2 + (yf - posy(23,11)).^2));
E_xr2312 = (-rho_r12/(2*pi*e_0)).*((xf - posx(23,12))./((xf - posx(23,12)).^2 + (yf - posy(23,12)).^2));
E_xr2313 = (rho_r1/(2*pi*e_0)).*((xf - posx(23,13))./((xf - posx(23,13)).^2 + (yf - posy(23,13)).^2));
E_xr2314 = (rho_r2/(2*pi*e_0)).*((xf - posx(23,14))./((xf - posx(23,14)).^2 + (yf - posy(23,14)).^2));
E_xr2315 = (rho_r3/(2*pi*e_0)).*((xf - posx(23,15))./((xf - posx(23,15)).^2 + (yf - posy(23,15)).^2));
E_xr2316 = (rho_r4/(2*pi*e_0)).*((xf - posx(23,16))./((xf - posx(23,16)).^2 + (yf - posy(23,16)).^2));
E_xr2317 = (rho_r5/(2*pi*e_0)).*((xf - posx(23,17))./((xf - posx(23,17)).^2 + (yf - posy(23,17)).^2));
E_xr2318 = (rho_r6/(2*pi*e_0)).*((xf - posx(23,18))./((xf - posx(23,18)).^2 + (yf - posy(23,18)).^2));
E_xr2319 = (rho_r7/(2*pi*e_0)).*((xf - posx(23,19))./((xf - posx(23,19)).^2 + (yf - posy(23,19)).^2));
E_xr2320 = (rho_r8/(2*pi*e_0)).*((xf - posx(23,20))./((xf - posx(23,20)).^2 + (yf - posy(23,20)).^2));
E_xr2321 = (rho_r9/(2*pi*e_0)).*((xf - posx(23,21))./((xf - posx(23,21)).^2 + (yf - posy(23,21)).^2));
E_xr2322 = (rho_r10/(2*pi*e_0)).*((xf - posx(23,22))./((xf - posx(23,22)).^2 + (yf - posy(23,22)).^2));
E_xr2324 = (rho_r12/(2*pi*e_0)).*((xf - posx(23,24))./((xf - posx(23,24)).^2 + (yf - posy(23,24)).^2));

E_xr241 = (-rho_r1/(2*pi*e_0)).*((xf - posx(24,1))./((xf - posx(24,1)).^2 + (yf - posy(24,1)).^2));
E_xr242 = (-rho_r2/(2*pi*e_0)).*((xf - posx(24,2))./((xf - posx(24,2)).^2 + (yf - posy(24,2)).^2));
E_xr243 = (-rho_r3/(2*pi*e_0)).*((xf - posx(24,3))./((xf - posx(24,3)).^2 + (yf - posy(24,3)).^2));
E_xr244 = (-rho_r4/(2*pi*e_0)).*((xf - posx(24,4))./((xf - posx(24,4)).^2 + (yf - posy(24,4)).^2));
E_xr245 = (-rho_r5/(2*pi*e_0)).*((xf - posx(24,5))./((xf - posx(24,5)).^2 + (yf - posy(24,5)).^2));
E_xr246 = (-rho_r6/(2*pi*e_0)).*((xf - posx(24,6))./((xf - posx(24,6)).^2 + (yf - posy(24,6)).^2));
E_xr247 = (-rho_r7/(2*pi*e_0)).*((xf - posx(24,7))./((xf - posx(24,7)).^2 + (yf - posy(24,7)).^2));
E_xr248 = (-rho_r8/(2*pi*e_0)).*((xf - posx(24,8))./((xf - posx(24,8)).^2 + (yf - posy(24,8)).^2));
E_xr249 = (-rho_r9/(2*pi*e_0)).*((xf - posx(24,9))./((xf - posx(24,9)).^2 + (yf - posy(24,9)).^2));
E_xr2410 = (-rho_r10/(2*pi*e_0)).*((xf - posx(24,10))./((xf - posx(24,10)).^2 + (yf - posy(24,10)).^2));
E_xr2411 = (-rho_r11/(2*pi*e_0)).*((xf - posx(24,11))./((xf - posx(24,11)).^2 + (yf - posy(24,11)).^2));
E_xr2412 = (-rho_r12/(2*pi*e_0)).*((xf - posx(24,12))./((xf - posx(24,12)).^2 + (yf - posy(24,12)).^2));
E_xr2413 = (rho_r1/(2*pi*e_0)).*((xf - posx(24,13))./((xf - posx(24,13)).^2 + (yf - posy(24,13)).^2));
E_xr2414 = (rho_r2/(2*pi*e_0)).*((xf - posx(24,14))./((xf - posx(24,14)).^2 + (yf - posy(24,14)).^2));
E_xr2415 = (rho_r3/(2*pi*e_0)).*((xf - posx(24,15))./((xf - posx(24,15)).^2 + (yf - posy(24,15)).^2));
E_xr2416 = (rho_r4/(2*pi*e_0)).*((xf - posx(24,16))./((xf - posx(24,16)).^2 + (yf - posy(24,16)).^2));
E_xr2417 = (rho_r5/(2*pi*e_0)).*((xf - posx(24,17))./((xf - posx(24,17)).^2 + (yf - posy(24,17)).^2));
E_xr2418 = (rho_r6/(2*pi*e_0)).*((xf - posx(24,18))./((xf - posx(24,18)).^2 + (yf - posy(24,18)).^2));
E_xr2419 = (rho_r7/(2*pi*e_0)).*((xf - posx(24,19))./((xf - posx(24,19)).^2 + (yf - posy(24,19)).^2));
E_xr2420 = (rho_r8/(2*pi*e_0)).*((xf - posx(24,20))./((xf - posx(24,20)).^2 + (yf - posy(24,20)).^2));
E_xr2421 = (rho_r9/(2*pi*e_0)).*((xf - posx(24,21))./((xf - posx(24,21)).^2 + (yf - posy(24,21)).^2));
E_xr2422 = (rho_r10/(2*pi*e_0)).*((xf - posx(24,22))./((xf - posx(24,22)).^2 + (yf - posy(24,22)).^2));
E_xr2423 = (rho_r11/(2*pi*e_0)).*((xf - posx(24,23))./((xf - posx(24,23)).^2 + (yf - posy(24,23)).^2));

%% E_xi componente x imaginario campo elétrico condutor 2 fase b assim segue:

E_xi12 = (-rho_i2/(2*pi*e_0)).*((xf - posx(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_xi13 = (-rho_i3/(2*pi*e_0)).*((xf - posx(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_xi14 = (-rho_i4/(2*pi*e_0)).*((xf - posx(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_xi15 = (-rho_i5/(2*pi*e_0)).*((xf - posx(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_xi16 = (-rho_i6/(2*pi*e_0)).*((xf - posx(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_xi17 = (-rho_i7/(2*pi*e_0)).*((xf - posx(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_xi18 = (-rho_i8/(2*pi*e_0)).*((xf - posx(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_xi19 = (-rho_i9/(2*pi*e_0)).*((xf - posx(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_xi110 = (-rho_i10/(2*pi*e_0)).*((xf - posx(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_xi1_11 = (-rho_i11/(2*pi*e_0)).*((xf - posx(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_xi1_12 = (-rho_i12/(2*pi*e_0)).*((xf - posx(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_xi1_13 = (rho_i1/(2*pi*e_0)).*((xf - posx(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_xi1_14 = (rho_i2/(2*pi*e_0)).*((xf - posx(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_xi1_15 = (rho_i3/(2*pi*e_0)).*((xf - posx(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_xi1_16 = (rho_i4/(2*pi*e_0)).*((xf - posx(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_xi1_17 = (rho_i5/(2*pi*e_0)).*((xf - posx(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_xi1_18 = (rho_i6/(2*pi*e_0)).*((xf - posx(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));
E_xi1_19 = (rho_i7/(2*pi*e_0)).*((xf - posx(1,19))./((xf - posx(1,19)).^2 + (yf - posy(1,19)).^2));
E_xi1_20 = (rho_i8/(2*pi*e_0)).*((xf - posx(1,20))./((xf - posx(1,20)).^2 + (yf - posy(1,20)).^2));
E_xi1_21 = (rho_i9/(2*pi*e_0)).*((xf - posx(1,21))./((xf - posx(1,21)).^2 + (yf - posy(1,21)).^2));
E_xi1_22 = (rho_i10/(2*pi*e_0)).*((xf - posx(1,22))./((xf - posx(1,22)).^2 + (yf - posy(1,22)).^2));
E_xi1_23 = (rho_i11/(2*pi*e_0)).*((xf - posx(1,23))./((xf - posx(1,23)).^2 + (yf - posy(1,23)).^2));
E_xi1_24 = (rho_i12/(2*pi*e_0)).*((xf - posx(1,24))./((xf - posx(1,24)).^2 + (yf - posy(1,24)).^2));

E_xi21 = (-rho_i1/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xi23 = (-rho_i3/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xi24 = (-rho_i4/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xi25 = (-rho_i5/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xi26 = (-rho_i6/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_xi27 = (-rho_i7/(2*pi*e_0)).*((xf - posx(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_xi28 = (-rho_i8/(2*pi*e_0)).*((xf - posx(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_xi29 = (-rho_i9/(2*pi*e_0)).*((xf - posx(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_xi210 = (-rho_i10/(2*pi*e_0)).*((xf - posx(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_xi2_11 = (-rho_i11/(2*pi*e_0)).*((xf - posx(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_xi2_12 = (-rho_i12/(2*pi*e_0)).*((xf - posx(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_xi2_13 = (rho_i1/(2*pi*e_0)).*((xf - posx(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_xi2_14 = (rho_i2/(2*pi*e_0)).*((xf - posx(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_xi2_15 = (rho_i3/(2*pi*e_0)).*((xf - posx(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_xi2_16 = (rho_i4/(2*pi*e_0)).*((xf - posx(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_xi2_17 = (rho_i5/(2*pi*e_0)).*((xf - posx(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_xi2_18 = (rho_i6/(2*pi*e_0)).*((xf - posx(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));
E_xi2_19 = (rho_i7/(2*pi*e_0)).*((xf - posx(2,19))./((xf - posx(2,19)).^2 + (yf - posy(2,19)).^2));
E_xi2_20 = (rho_i8/(2*pi*e_0)).*((xf - posx(2,20))./((xf - posx(2,20)).^2 + (yf - posy(2,20)).^2));
E_xi2_21 = (rho_i9/(2*pi*e_0)).*((xf - posx(2,21))./((xf - posx(2,21)).^2 + (yf - posy(2,21)).^2));
E_xi2_22 = (rho_i10/(2*pi*e_0)).*((xf - posx(2,22))./((xf - posx(2,22)).^2 + (yf - posy(2,22)).^2));
E_xi2_23 = (rho_i11/(2*pi*e_0)).*((xf - posx(2,23))./((xf - posx(2,23)).^2 + (yf - posy(2,23)).^2));
E_xi2_24 = (rho_i12/(2*pi*e_0)).*((xf - posx(2,24))./((xf - posx(2,24)).^2 + (yf - posy(2,24)).^2));

E_xi31 = (-rho_i1/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xi32 = (-rho_i2/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xi34 = (-rho_i4/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xi35 = (-rho_i5/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xi36 = (-rho_i6/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_xi37 = (-rho_i7/(2*pi*e_0)).*((xf - posx(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_xi38 = (-rho_i8/(2*pi*e_0)).*((xf - posx(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_xi39 = (-rho_i9/(2*pi*e_0)).*((xf - posx(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_xi310 = (-rho_i10/(2*pi*e_0)).*((xf - posx(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_xi311 = (-rho_i11/(2*pi*e_0)).*((xf - posx(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_xi312 = (-rho_i12/(2*pi*e_0)).*((xf - posx(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_xi313 = (rho_i1/(2*pi*e_0)).*((xf - posx(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_xi314 = (rho_i2/(2*pi*e_0)).*((xf - posx(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_xi315 = (rho_i3/(2*pi*e_0)).*((xf - posx(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_xi316 = (rho_i4/(2*pi*e_0)).*((xf - posx(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_xi317 = (rho_i5/(2*pi*e_0)).*((xf - posx(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_xi318 = (rho_i6/(2*pi*e_0)).*((xf - posx(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));
E_xi319 = (rho_i7/(2*pi*e_0)).*((xf - posx(3,19))./((xf - posx(3,19)).^2 + (yf - posy(3,19)).^2));
E_xi320 = (rho_i8/(2*pi*e_0)).*((xf - posx(3,20))./((xf - posx(3,20)).^2 + (yf - posy(3,20)).^2));
E_xi321 = (rho_i9/(2*pi*e_0)).*((xf - posx(3,21))./((xf - posx(3,21)).^2 + (yf - posy(3,21)).^2));
E_xi322 = (rho_i10/(2*pi*e_0)).*((xf - posx(3,22))./((xf - posx(3,22)).^2 + (yf - posy(3,22)).^2));
E_xi323 = (rho_i11/(2*pi*e_0)).*((xf - posx(3,23))./((xf - posx(3,23)).^2 + (yf - posy(3,23)).^2));
E_xi324 = (rho_i12/(2*pi*e_0)).*((xf - posx(3,24))./((xf - posx(3,24)).^2 + (yf - posy(3,24)).^2));

E_xi41 = (-rho_i1/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xi42 = (-rho_i2/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xi43 = (-rho_i3/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xi45 = (-rho_i5/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xi46 = (-rho_i6/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_xi47 = (-rho_i7/(2*pi*e_0)).*((xf - posx(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_xi48 = (-rho_i8/(2*pi*e_0)).*((xf - posx(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_xi49 = (-rho_i9/(2*pi*e_0)).*((xf - posx(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_xi410 = (-rho_i10/(2*pi*e_0)).*((xf - posx(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_xi411 = (-rho_i11/(2*pi*e_0)).*((xf - posx(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_xi412 = (-rho_i12/(2*pi*e_0)).*((xf - posx(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_xi413 = (rho_i1/(2*pi*e_0)).*((xf - posx(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_xi414 = (rho_i2/(2*pi*e_0)).*((xf - posx(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_xi415 = (rho_i3/(2*pi*e_0)).*((xf - posx(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_xi416 = (rho_i4/(2*pi*e_0)).*((xf - posx(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_xi417 = (rho_i5/(2*pi*e_0)).*((xf - posx(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_xi418 = (rho_i6/(2*pi*e_0)).*((xf - posx(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));
E_xi419 = (rho_i7/(2*pi*e_0)).*((xf - posx(4,19))./((xf - posx(4,19)).^2 + (yf - posy(4,19)).^2));
E_xi420 = (rho_i8/(2*pi*e_0)).*((xf - posx(4,20))./((xf - posx(4,20)).^2 + (yf - posy(4,20)).^2));
E_xi421 = (rho_i9/(2*pi*e_0)).*((xf - posx(4,21))./((xf - posx(4,21)).^2 + (yf - posy(4,21)).^2));
E_xi422 = (rho_i10/(2*pi*e_0)).*((xf - posx(4,22))./((xf - posx(4,22)).^2 + (yf - posy(4,22)).^2));
E_xi423 = (rho_i11/(2*pi*e_0)).*((xf - posx(4,23))./((xf - posx(4,23)).^2 + (yf - posy(4,23)).^2));
E_xi424 = (rho_i12/(2*pi*e_0)).*((xf - posx(4,24))./((xf - posx(4,24)).^2 + (yf - posy(4,24)).^2));

E_xi51 = (-rho_i1/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xi52 = (-rho_i2/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xi53 = (-rho_i3/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xi54 = (-rho_i4/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xi56 = (-rho_i6/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_xi57 = (-rho_i7/(2*pi*e_0)).*((xf - posx(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_xi58 = (-rho_i8/(2*pi*e_0)).*((xf - posx(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_xi59 = (-rho_i9/(2*pi*e_0)).*((xf - posx(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_xi510 = (-rho_i10/(2*pi*e_0)).*((xf - posx(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_xi511 = (-rho_i11/(2*pi*e_0)).*((xf - posx(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_xi512 = (-rho_i12/(2*pi*e_0)).*((xf - posx(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_xi513 = (rho_i1/(2*pi*e_0)).*((xf - posx(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_xi514 = (rho_i2/(2*pi*e_0)).*((xf - posx(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_xi515 = (rho_i3/(2*pi*e_0)).*((xf - posx(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_xi516 = (rho_i4/(2*pi*e_0)).*((xf - posx(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_xi517 = (rho_i5/(2*pi*e_0)).*((xf - posx(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_xi518 = (rho_i6/(2*pi*e_0)).*((xf - posx(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));
E_xi519 = (rho_i7/(2*pi*e_0)).*((xf - posx(5,19))./((xf - posx(5,19)).^2 + (yf - posy(5,19)).^2));
E_xi520 = (rho_i8/(2*pi*e_0)).*((xf - posx(5,20))./((xf - posx(5,20)).^2 + (yf - posy(5,20)).^2));
E_xi521 = (rho_i9/(2*pi*e_0)).*((xf - posx(5,21))./((xf - posx(5,21)).^2 + (yf - posy(5,21)).^2));
E_xi522 = (rho_i10/(2*pi*e_0)).*((xf - posx(5,22))./((xf - posx(5,22)).^2 + (yf - posy(5,22)).^2));
E_xi523 = (rho_i11/(2*pi*e_0)).*((xf - posx(5,23))./((xf - posx(5,23)).^2 + (yf - posy(5,23)).^2));
E_xi524 = (rho_i12/(2*pi*e_0)).*((xf - posx(5,24))./((xf - posx(5,24)).^2 + (yf - posy(5,24)).^2));

E_xi61 = (-rho_i1/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xi62 = (-rho_i2/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xi63 = (-rho_i3/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xi64 = (-rho_i4/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xi65 = (-rho_i5/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_xi67 = (-rho_i7/(2*pi*e_0)).*((xf - posx(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_xi68 = (-rho_i8/(2*pi*e_0)).*((xf - posx(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_xi69 = (-rho_i9/(2*pi*e_0)).*((xf - posx(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_xi610 = (-rho_i10/(2*pi*e_0)).*((xf - posx(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_xi611 = (-rho_i11/(2*pi*e_0)).*((xf - posx(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_xi612 = (-rho_i12/(2*pi*e_0)).*((xf - posx(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_xi613 = (rho_i1/(2*pi*e_0)).*((xf - posx(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_xi614 = (rho_i2/(2*pi*e_0)).*((xf - posx(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_xi615 = (rho_i3/(2*pi*e_0)).*((xf - posx(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_xi616 = (rho_i4/(2*pi*e_0)).*((xf - posx(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_xi617 = (rho_i5/(2*pi*e_0)).*((xf - posx(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_xi618 = (rho_i6/(2*pi*e_0)).*((xf - posx(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));
E_xi619 = (rho_i7/(2*pi*e_0)).*((xf - posx(6,19))./((xf - posx(6,19)).^2 + (yf - posy(6,19)).^2));
E_xi620 = (rho_i8/(2*pi*e_0)).*((xf - posx(6,20))./((xf - posx(6,20)).^2 + (yf - posy(6,20)).^2));
E_xi621 = (rho_i9/(2*pi*e_0)).*((xf - posx(6,21))./((xf - posx(6,21)).^2 + (yf - posy(6,21)).^2));
E_xi622 = (rho_i10/(2*pi*e_0)).*((xf - posx(6,22))./((xf - posx(6,22)).^2 + (yf - posy(6,22)).^2));
E_xi623 = (rho_i11/(2*pi*e_0)).*((xf - posx(6,23))./((xf - posx(6,23)).^2 + (yf - posy(6,23)).^2));
E_xi624 = (rho_i12/(2*pi*e_0)).*((xf - posx(6,24))./((xf - posx(6,24)).^2 + (yf - posy(6,24)).^2));

E_xi71 = (-rho_i1/(2*pi*e_0)).*((xf - posx(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_xi72 = (-rho_i2/(2*pi*e_0)).*((xf - posx(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_xi73 = (-rho_i3/(2*pi*e_0)).*((xf - posx(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_xi74 = (-rho_i4/(2*pi*e_0)).*((xf - posx(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_xi75 = (-rho_i5/(2*pi*e_0)).*((xf - posx(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_xi76 = (-rho_i6/(2*pi*e_0)).*((xf - posx(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_xi78 = (-rho_i8/(2*pi*e_0)).*((xf - posx(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_xi79 = (-rho_i9/(2*pi*e_0)).*((xf - posx(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_xi710 = (-rho_i10/(2*pi*e_0)).*((xf - posx(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_xi711 = (-rho_i11/(2*pi*e_0)).*((xf - posx(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_xi712 = (-rho_i12/(2*pi*e_0)).*((xf - posx(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_xi713 = (rho_i1/(2*pi*e_0)).*((xf - posx(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_xi714 = (rho_i2/(2*pi*e_0)).*((xf - posx(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_xi715 = (rho_i3/(2*pi*e_0)).*((xf - posx(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_xi716 = (rho_i4/(2*pi*e_0)).*((xf - posx(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_xi717 = (rho_i5/(2*pi*e_0)).*((xf - posx(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_xi718 = (rho_i6/(2*pi*e_0)).*((xf - posx(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));
E_xi719 = (rho_i7/(2*pi*e_0)).*((xf - posx(7,19))./((xf - posx(7,19)).^2 + (yf - posy(7,19)).^2));
E_xi720 = (rho_i8/(2*pi*e_0)).*((xf - posx(7,20))./((xf - posx(7,20)).^2 + (yf - posy(7,20)).^2));
E_xi721 = (rho_i9/(2*pi*e_0)).*((xf - posx(7,21))./((xf - posx(7,21)).^2 + (yf - posy(7,21)).^2));
E_xi722 = (rho_i10/(2*pi*e_0)).*((xf - posx(7,22))./((xf - posx(7,22)).^2 + (yf - posy(7,22)).^2));
E_xi723 = (rho_i11/(2*pi*e_0)).*((xf - posx(7,23))./((xf - posx(7,23)).^2 + (yf - posy(7,23)).^2));
E_xi724 = (rho_i12/(2*pi*e_0)).*((xf - posx(7,24))./((xf - posx(7,24)).^2 + (yf - posy(7,24)).^2));

E_xi81 = (-rho_i1/(2*pi*e_0)).*((xf - posx(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_xi82 = (-rho_i2/(2*pi*e_0)).*((xf - posx(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_xi83 = (-rho_i3/(2*pi*e_0)).*((xf - posx(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_xi84 = (-rho_i4/(2*pi*e_0)).*((xf - posx(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_xi85 = (-rho_i5/(2*pi*e_0)).*((xf - posx(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_xi86 = (-rho_i6/(2*pi*e_0)).*((xf - posx(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_xi87 = (-rho_i7/(2*pi*e_0)).*((xf - posx(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_xi89 = (-rho_i9/(2*pi*e_0)).*((xf - posx(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_xi810 = (-rho_i10/(2*pi*e_0)).*((xf - posx(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_xi811 = (-rho_i11/(2*pi*e_0)).*((xf - posx(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_xi812 = (-rho_i12/(2*pi*e_0)).*((xf - posx(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_xi813 = (rho_i1/(2*pi*e_0)).*((xf - posx(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_xi814 = (rho_i2/(2*pi*e_0)).*((xf - posx(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_xi815 = (rho_i3/(2*pi*e_0)).*((xf - posx(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_xi816 = (rho_i4/(2*pi*e_0)).*((xf - posx(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_xi817 = (rho_i5/(2*pi*e_0)).*((xf - posx(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_xi818 = (rho_i6/(2*pi*e_0)).*((xf - posx(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));
E_xi819 = (rho_i7/(2*pi*e_0)).*((xf - posx(8,19))./((xf - posx(8,19)).^2 + (yf - posy(8,19)).^2));
E_xi820 = (rho_i8/(2*pi*e_0)).*((xf - posx(8,20))./((xf - posx(8,20)).^2 + (yf - posy(8,20)).^2));
E_xi821 = (rho_i9/(2*pi*e_0)).*((xf - posx(8,21))./((xf - posx(8,21)).^2 + (yf - posy(8,21)).^2));
E_xi822 = (rho_i10/(2*pi*e_0)).*((xf - posx(8,22))./((xf - posx(8,22)).^2 + (yf - posy(8,22)).^2));
E_xi823 = (rho_i11/(2*pi*e_0)).*((xf - posx(8,23))./((xf - posx(8,23)).^2 + (yf - posy(8,23)).^2));
E_xi824 = (rho_i12/(2*pi*e_0)).*((xf - posx(8,24))./((xf - posx(8,24)).^2 + (yf - posy(8,24)).^2));

E_xi91 = (-rho_i1/(2*pi*e_0)).*((xf - posx(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_xi92 = (-rho_i2/(2*pi*e_0)).*((xf - posx(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_xi93 = (-rho_i3/(2*pi*e_0)).*((xf - posx(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_xi94 = (-rho_i4/(2*pi*e_0)).*((xf - posx(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_xi95 = (-rho_i5/(2*pi*e_0)).*((xf - posx(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_xi96 = (-rho_i6/(2*pi*e_0)).*((xf - posx(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_xi97 = (-rho_i7/(2*pi*e_0)).*((xf - posx(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_xi98 = (-rho_i8/(2*pi*e_0)).*((xf - posx(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_xi910 = (-rho_i10/(2*pi*e_0)).*((xf - posx(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_xi911 = (-rho_i11/(2*pi*e_0)).*((xf - posx(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_xi912 = (-rho_i12/(2*pi*e_0)).*((xf - posx(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_xi913 = (rho_i1/(2*pi*e_0)).*((xf - posx(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_xi914 = (rho_i2/(2*pi*e_0)).*((xf - posx(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_xi915 = (rho_i3/(2*pi*e_0)).*((xf - posx(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_xi916 = (rho_i4/(2*pi*e_0)).*((xf - posx(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_xi917 = (rho_i5/(2*pi*e_0)).*((xf - posx(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_xi918 = (rho_i6/(2*pi*e_0)).*((xf - posx(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));
E_xi919 = (rho_i7/(2*pi*e_0)).*((xf - posx(9,19))./((xf - posx(9,19)).^2 + (yf - posy(9,19)).^2));
E_xi920 = (rho_i8/(2*pi*e_0)).*((xf - posx(9,20))./((xf - posx(9,20)).^2 + (yf - posy(9,20)).^2));
E_xi921 = (rho_i9/(2*pi*e_0)).*((xf - posx(9,21))./((xf - posx(9,21)).^2 + (yf - posy(9,21)).^2));
E_xi922 = (rho_i10/(2*pi*e_0)).*((xf - posx(9,22))./((xf - posx(9,22)).^2 + (yf - posy(9,22)).^2));
E_xi923 = (rho_i11/(2*pi*e_0)).*((xf - posx(9,23))./((xf - posx(9,23)).^2 + (yf - posy(9,23)).^2));
E_xi924 = (rho_i12/(2*pi*e_0)).*((xf - posx(9,24))./((xf - posx(9,24)).^2 + (yf - posy(9,24)).^2));

E_xi101 = (-rho_i1/(2*pi*e_0)).*((xf - posx(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_xi102 = (-rho_i2/(2*pi*e_0)).*((xf - posx(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_xi103 = (-rho_i3/(2*pi*e_0)).*((xf - posx(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_xi104 = (-rho_i4/(2*pi*e_0)).*((xf - posx(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_xi105 = (-rho_i5/(2*pi*e_0)).*((xf - posx(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_xi106 = (-rho_i6/(2*pi*e_0)).*((xf - posx(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_xi107 = (-rho_i7/(2*pi*e_0)).*((xf - posx(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_xi108 = (-rho_i8/(2*pi*e_0)).*((xf - posx(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_xi109 = (-rho_i9/(2*pi*e_0)).*((xf - posx(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_xi1011 = (-rho_i11/(2*pi*e_0)).*((xf - posx(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_xi1012 = (-rho_i12/(2*pi*e_0)).*((xf - posx(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_xi1013 = (rho_i1/(2*pi*e_0)).*((xf - posx(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_xi1014 = (rho_i2/(2*pi*e_0)).*((xf - posx(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_xi1015 = (rho_i3/(2*pi*e_0)).*((xf - posx(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_xi1016 = (rho_i4/(2*pi*e_0)).*((xf - posx(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_xi1017 = (rho_i5/(2*pi*e_0)).*((xf - posx(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_xi1018 = (rho_i6/(2*pi*e_0)).*((xf - posx(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));
E_xi10_19 = (rho_i7/(2*pi*e_0)).*((xf - posx(10,19))./((xf - posx(10,19)).^2 + (yf - posy(10,19)).^2));
E_xi10_20 = (rho_i8/(2*pi*e_0)).*((xf - posx(10,20))./((xf - posx(10,20)).^2 + (yf - posy(10,20)).^2));
E_xi10_21 = (rho_i9/(2*pi*e_0)).*((xf - posx(10,21))./((xf - posx(10,21)).^2 + (yf - posy(10,21)).^2));
E_xi10_22 = (rho_i10/(2*pi*e_0)).*((xf - posx(10,22))./((xf - posx(10,22)).^2 + (yf - posy(10,22)).^2));
E_xi10_23 = (rho_i11/(2*pi*e_0)).*((xf - posx(10,23))./((xf - posx(10,23)).^2 + (yf - posy(10,23)).^2));
E_xi10_24 = (rho_i12/(2*pi*e_0)).*((xf - posx(10,24))./((xf - posx(10,24)).^2 + (yf - posy(10,24)).^2));

E_xi11_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_xi11_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_xi113 = (-rho_i3/(2*pi*e_0)).*((xf - posx(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_xi114 = (-rho_i4/(2*pi*e_0)).*((xf - posx(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_xi115 = (-rho_i5/(2*pi*e_0)).*((xf - posx(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_xi116 = (-rho_i6/(2*pi*e_0)).*((xf - posx(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_xi117 = (-rho_i7/(2*pi*e_0)).*((xf - posx(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_xi118 = (-rho_i8/(2*pi*e_0)).*((xf - posx(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_xi119 = (-rho_i9/(2*pi*e_0)).*((xf - posx(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_xi1110 = (-rho_i10/(2*pi*e_0)).*((xf - posx(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_xi1112 = (-rho_i12/(2*pi*e_0)).*((xf - posx(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_xi1113 = (rho_i1/(2*pi*e_0)).*((xf - posx(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_xi1114 = (rho_i2/(2*pi*e_0)).*((xf - posx(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_xi1115 = (rho_i3/(2*pi*e_0)).*((xf - posx(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_xi1116 = (rho_i4/(2*pi*e_0)).*((xf - posx(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_xi1117 = (rho_i5/(2*pi*e_0)).*((xf - posx(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_xi1118 = (rho_i6/(2*pi*e_0)).*((xf - posx(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));
E_xi1019 = (rho_i7/(2*pi*e_0)).*((xf - posx(11,19))./((xf - posx(11,19)).^2 + (yf - posy(11,19)).^2));
E_xi1020 = (rho_i8/(2*pi*e_0)).*((xf - posx(11,20))./((xf - posx(11,20)).^2 + (yf - posy(11,20)).^2));
E_xi1021 = (rho_i9/(2*pi*e_0)).*((xf - posx(11,21))./((xf - posx(11,21)).^2 + (yf - posy(11,21)).^2));
E_xi1022 = (rho_i10/(2*pi*e_0)).*((xf - posx(11,22))./((xf - posx(11,22)).^2 + (yf - posy(11,22)).^2));
E_xi1023 = (rho_i11/(2*pi*e_0)).*((xf - posx(11,23))./((xf - posx(11,23)).^2 + (yf - posy(11,23)).^2));
E_xi1024 = (rho_i12/(2*pi*e_0)).*((xf - posx(11,24))./((xf - posx(11,24)).^2 + (yf - posy(11,24)).^2));

E_xi12_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_xi12_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_xi123 = (-rho_i3/(2*pi*e_0)).*((xf - posx(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_xi124 = (-rho_i4/(2*pi*e_0)).*((xf - posx(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_xi125 = (-rho_i5/(2*pi*e_0)).*((xf - posx(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_xi126 = (-rho_i6/(2*pi*e_0)).*((xf - posx(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_xi127 = (-rho_i7/(2*pi*e_0)).*((xf - posx(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_xi128 = (-rho_i8/(2*pi*e_0)).*((xf - posx(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_xi129 = (-rho_i9/(2*pi*e_0)).*((xf - posx(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_xi1210 = (-rho_i10/(2*pi*e_0)).*((xf - posx(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_xi1211 = (-rho_i11/(2*pi*e_0)).*((xf - posx(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_xi1213 = (rho_i1/(2*pi*e_0)).*((xf - posx(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_xi1214 = (rho_i2/(2*pi*e_0)).*((xf - posx(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_xi1215 = (rho_i3/(2*pi*e_0)).*((xf - posx(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_xi1216 = (rho_i4/(2*pi*e_0)).*((xf - posx(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_xi1217 = (rho_i5/(2*pi*e_0)).*((xf - posx(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_xi1218 = (rho_i6/(2*pi*e_0)).*((xf - posx(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));
E_xi1219 = (rho_i7/(2*pi*e_0)).*((xf - posx(12,19))./((xf - posx(12,19)).^2 + (yf - posy(12,19)).^2));
E_xi1220 = (rho_i8/(2*pi*e_0)).*((xf - posx(12,20))./((xf - posx(12,20)).^2 + (yf - posy(12,20)).^2));
E_xi1221 = (rho_i9/(2*pi*e_0)).*((xf - posx(12,21))./((xf - posx(12,21)).^2 + (yf - posy(12,21)).^2));
E_xi1222 = (rho_i10/(2*pi*e_0)).*((xf - posx(12,22))./((xf - posx(12,22)).^2 + (yf - posy(12,22)).^2));
E_xi1223 = (rho_i11/(2*pi*e_0)).*((xf - posx(12,23))./((xf - posx(12,23)).^2 + (yf - posy(12,23)).^2));
E_xi1224 = (rho_i12/(2*pi*e_0)).*((xf - posx(12,24))./((xf - posx(12,24)).^2 + (yf - posy(12,24)).^2));

E_xi131 = (-rho_i1/(2*pi*e_0)).*((xf - posx(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_xi132 = (-rho_i2/(2*pi*e_0)).*((xf - posx(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_xi133 = (-rho_i3/(2*pi*e_0)).*((xf - posx(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_xi134 = (-rho_i4/(2*pi*e_0)).*((xf - posx(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_xi135 = (-rho_i5/(2*pi*e_0)).*((xf - posx(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_xi136 = (-rho_i6/(2*pi*e_0)).*((xf - posx(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_xi137 = (-rho_i7/(2*pi*e_0)).*((xf - posx(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_xi138 = (-rho_i8/(2*pi*e_0)).*((xf - posx(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_xi139 = (-rho_i9/(2*pi*e_0)).*((xf - posx(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_xi1310 = (-rho_i10/(2*pi*e_0)).*((xf - posx(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_xi1311 = (-rho_i11/(2*pi*e_0)).*((xf - posx(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_xi1312 = (-rho_i12/(2*pi*e_0)).*((xf - posx(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_xi1314 = (rho_i2/(2*pi*e_0)).*((xf - posx(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_xi1315 = (rho_i3/(2*pi*e_0)).*((xf - posx(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_xi1316 = (rho_i4/(2*pi*e_0)).*((xf - posx(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_xi1317 = (rho_i5/(2*pi*e_0)).*((xf - posx(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_xi1318 = (rho_i6/(2*pi*e_0)).*((xf - posx(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));
E_xi1319 = (rho_i7/(2*pi*e_0)).*((xf - posx(13,19))./((xf - posx(13,19)).^2 + (yf - posy(13,19)).^2));
E_xi1320 = (rho_i8/(2*pi*e_0)).*((xf - posx(13,20))./((xf - posx(13,20)).^2 + (yf - posy(13,20)).^2));
E_xi1321 = (rho_i9/(2*pi*e_0)).*((xf - posx(13,21))./((xf - posx(13,21)).^2 + (yf - posy(13,21)).^2));
E_xi1322 = (rho_i10/(2*pi*e_0)).*((xf - posx(13,22))./((xf - posx(13,22)).^2 + (yf - posy(13,22)).^2));
E_xi1323 = (rho_i11/(2*pi*e_0)).*((xf - posx(13,23))./((xf - posx(13,23)).^2 + (yf - posy(13,23)).^2));
E_xi1324 = (rho_i12/(2*pi*e_0)).*((xf - posx(13,24))./((xf - posx(13,24)).^2 + (yf - posy(13,24)).^2));

E_xi141 = (-rho_i1/(2*pi*e_0)).*((xf - posx(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_xi142 = (-rho_i2/(2*pi*e_0)).*((xf - posx(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_xi143 = (-rho_i3/(2*pi*e_0)).*((xf - posx(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_xi144 = (-rho_i4/(2*pi*e_0)).*((xf - posx(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_xi145 = (-rho_i5/(2*pi*e_0)).*((xf - posx(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_xi146 = (-rho_i6/(2*pi*e_0)).*((xf - posx(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_xi147 = (-rho_i7/(2*pi*e_0)).*((xf - posx(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_xi148 = (-rho_i8/(2*pi*e_0)).*((xf - posx(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_xi149 = (-rho_i9/(2*pi*e_0)).*((xf - posx(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_xi1410 = (-rho_i10/(2*pi*e_0)).*((xf - posx(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_xi1411 = (-rho_i11/(2*pi*e_0)).*((xf - posx(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_xi1412 = (-rho_i12/(2*pi*e_0)).*((xf - posx(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_xi1413 = (rho_i1/(2*pi*e_0)).*((xf - posx(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_xi1415 = (rho_i3/(2*pi*e_0)).*((xf - posx(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_xi1416 = (rho_i4/(2*pi*e_0)).*((xf - posx(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_xi1417 = (rho_i5/(2*pi*e_0)).*((xf - posx(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_xi1418 = (rho_i6/(2*pi*e_0)).*((xf - posx(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));
E_xi1419 = (rho_i7/(2*pi*e_0)).*((xf - posx(14,19))./((xf - posx(14,19)).^2 + (yf - posy(14,19)).^2));
E_xi1420 = (rho_i8/(2*pi*e_0)).*((xf - posx(14,20))./((xf - posx(14,20)).^2 + (yf - posy(14,20)).^2));
E_xi1421 = (rho_i9/(2*pi*e_0)).*((xf - posx(14,21))./((xf - posx(14,21)).^2 + (yf - posy(14,21)).^2));
E_xi1422 = (rho_i10/(2*pi*e_0)).*((xf - posx(14,22))./((xf - posx(14,22)).^2 + (yf - posy(14,22)).^2));
E_xi1423 = (rho_i11/(2*pi*e_0)).*((xf - posx(14,23))./((xf - posx(14,23)).^2 + (yf - posy(14,23)).^2));
E_xi1424 = (rho_i12/(2*pi*e_0)).*((xf - posx(14,24))./((xf - posx(14,24)).^2 + (yf - posy(14,24)).^2));

E_xi151 = (-rho_i1/(2*pi*e_0)).*((xf - posx(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_xi152 = (-rho_i2/(2*pi*e_0)).*((xf - posx(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_xi153 = (-rho_i3/(2*pi*e_0)).*((xf - posx(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_xi154 = (-rho_i4/(2*pi*e_0)).*((xf - posx(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_xi155 = (-rho_i5/(2*pi*e_0)).*((xf - posx(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_xi156 = (-rho_i6/(2*pi*e_0)).*((xf - posx(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_xi157 = (-rho_i7/(2*pi*e_0)).*((xf - posx(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_xi158 = (-rho_i8/(2*pi*e_0)).*((xf - posx(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_xi159 = (-rho_i9/(2*pi*e_0)).*((xf - posx(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_xi1510 = (-rho_i10/(2*pi*e_0)).*((xf - posx(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_xi1511 = (-rho_i11/(2*pi*e_0)).*((xf - posx(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_xi1512 = (-rho_i12/(2*pi*e_0)).*((xf - posx(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_xi1513 = (rho_i1/(2*pi*e_0)).*((xf - posx(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_xi1514 = (rho_i2/(2*pi*e_0)).*((xf - posx(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_xi1516 = (rho_i4/(2*pi*e_0)).*((xf - posx(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_xi1517 = (rho_i5/(2*pi*e_0)).*((xf - posx(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_xi1518 = (rho_i6/(2*pi*e_0)).*((xf - posx(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));
E_xi1519 = (rho_i7/(2*pi*e_0)).*((xf - posx(15,19))./((xf - posx(15,19)).^2 + (yf - posy(15,19)).^2));
E_xi1520 = (rho_i8/(2*pi*e_0)).*((xf - posx(15,20))./((xf - posx(15,20)).^2 + (yf - posy(15,20)).^2));
E_xi1521 = (rho_i9/(2*pi*e_0)).*((xf - posx(15,21))./((xf - posx(15,21)).^2 + (yf - posy(15,21)).^2));
E_xi1522 = (rho_i10/(2*pi*e_0)).*((xf - posx(15,22))./((xf - posx(15,22)).^2 + (yf - posy(15,22)).^2));
E_xi1523 = (rho_i11/(2*pi*e_0)).*((xf - posx(15,23))./((xf - posx(15,23)).^2 + (yf - posy(15,23)).^2));
E_xi1524 = (rho_i12/(2*pi*e_0)).*((xf - posx(15,24))./((xf - posx(15,24)).^2 + (yf - posy(15,24)).^2));

E_xi161 = (-rho_i1/(2*pi*e_0)).*((xf - posx(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_xi162 = (-rho_i2/(2*pi*e_0)).*((xf - posx(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_xi163 = (-rho_i3/(2*pi*e_0)).*((xf - posx(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_xi164 = (-rho_i4/(2*pi*e_0)).*((xf - posx(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_xi165 = (-rho_i5/(2*pi*e_0)).*((xf - posx(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_xi166 = (-rho_i6/(2*pi*e_0)).*((xf - posx(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_xi167 = (-rho_i7/(2*pi*e_0)).*((xf - posx(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_xi168 = (-rho_i8/(2*pi*e_0)).*((xf - posx(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_xi169 = (-rho_i9/(2*pi*e_0)).*((xf - posx(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_xi1610 = (-rho_i10/(2*pi*e_0)).*((xf - posx(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_xi1611 = (-rho_i11/(2*pi*e_0)).*((xf - posx(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_xi1612 = (-rho_i12/(2*pi*e_0)).*((xf - posx(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_xi1613 = (rho_i1/(2*pi*e_0)).*((xf - posx(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_xi1614 = (rho_i2/(2*pi*e_0)).*((xf - posx(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_xi1615 = (rho_i3/(2*pi*e_0)).*((xf - posx(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_xi1617 = (rho_i5/(2*pi*e_0)).*((xf - posx(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_xi1618 = (rho_i6/(2*pi*e_0)).*((xf - posx(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));
E_xi1619 = (rho_i7/(2*pi*e_0)).*((xf - posx(16,19))./((xf - posx(16,19)).^2 + (yf - posy(16,19)).^2));
E_xi1620 = (rho_i8/(2*pi*e_0)).*((xf - posx(16,20))./((xf - posx(16,20)).^2 + (yf - posy(16,20)).^2));
E_xi1621 = (rho_i9/(2*pi*e_0)).*((xf - posx(16,21))./((xf - posx(16,21)).^2 + (yf - posy(16,21)).^2));
E_xi1622 = (rho_i10/(2*pi*e_0)).*((xf - posx(16,22))./((xf - posx(16,22)).^2 + (yf - posy(16,22)).^2));
E_xi1623 = (rho_i11/(2*pi*e_0)).*((xf - posx(16,23))./((xf - posx(16,23)).^2 + (yf - posy(16,23)).^2));
E_xi1624 = (rho_i12/(2*pi*e_0)).*((xf - posx(16,24))./((xf - posx(16,24)).^2 + (yf - posy(16,24)).^2));

E_xi171 = (-rho_i1/(2*pi*e_0)).*((xf - posx(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_xi172 = (-rho_i2/(2*pi*e_0)).*((xf - posx(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_xi173 = (-rho_i3/(2*pi*e_0)).*((xf - posx(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_xi174 = (-rho_i4/(2*pi*e_0)).*((xf - posx(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_xi175 = (-rho_i5/(2*pi*e_0)).*((xf - posx(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_xi176 = (-rho_i6/(2*pi*e_0)).*((xf - posx(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_xi177 = (-rho_i7/(2*pi*e_0)).*((xf - posx(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_xi178 = (-rho_i8/(2*pi*e_0)).*((xf - posx(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_xi179 = (-rho_i9/(2*pi*e_0)).*((xf - posx(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_xi1710 = (-rho_i10/(2*pi*e_0)).*((xf - posx(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_xi1711 = (-rho_i11/(2*pi*e_0)).*((xf - posx(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_xi1712 = (-rho_i12/(2*pi*e_0)).*((xf - posx(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_xi1713 = (rho_i1/(2*pi*e_0)).*((xf - posx(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_xi1714 = (rho_i2/(2*pi*e_0)).*((xf - posx(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_xi1715 = (rho_i3/(2*pi*e_0)).*((xf - posx(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_xi1716 = (rho_i4/(2*pi*e_0)).*((xf - posx(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_xi1718 = (rho_i6/(2*pi*e_0)).*((xf - posx(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));
E_xi1719 = (rho_i7/(2*pi*e_0)).*((xf - posx(17,19))./((xf - posx(17,19)).^2 + (yf - posy(17,19)).^2));
E_xi1720 = (rho_i8/(2*pi*e_0)).*((xf - posx(17,20))./((xf - posx(17,20)).^2 + (yf - posy(17,20)).^2));
E_xi1721 = (rho_i9/(2*pi*e_0)).*((xf - posx(17,21))./((xf - posx(17,21)).^2 + (yf - posy(17,21)).^2));
E_xi1722 = (rho_i10/(2*pi*e_0)).*((xf - posx(17,22))./((xf - posx(17,22)).^2 + (yf - posy(17,22)).^2));
E_xi1723 = (rho_i11/(2*pi*e_0)).*((xf - posx(17,23))./((xf - posx(17,23)).^2 + (yf - posy(17,23)).^2));
E_xi1724 = (rho_i12/(2*pi*e_0)).*((xf - posx(17,24))./((xf - posx(17,24)).^2 + (yf - posy(17,24)).^2));

E_xi181 = (-rho_i1/(2*pi*e_0)).*((xf - posx(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_xi182 = (-rho_i2/(2*pi*e_0)).*((xf - posx(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_xi183 = (-rho_i3/(2*pi*e_0)).*((xf - posx(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_xi184 = (-rho_i4/(2*pi*e_0)).*((xf - posx(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_xi185 = (-rho_i5/(2*pi*e_0)).*((xf - posx(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_xi186 = (-rho_i6/(2*pi*e_0)).*((xf - posx(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_xi187 = (-rho_i7/(2*pi*e_0)).*((xf - posx(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_xi188 = (-rho_i8/(2*pi*e_0)).*((xf - posx(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_xi189 = (-rho_i9/(2*pi*e_0)).*((xf - posx(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_xi1810 = (-rho_i10/(2*pi*e_0)).*((xf - posx(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_xi1811 = (-rho_i11/(2*pi*e_0)).*((xf - posx(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_xi1812 = (-rho_i12/(2*pi*e_0)).*((xf - posx(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_xi1813 = (rho_i1/(2*pi*e_0)).*((xf - posx(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_xi1814 = (rho_i2/(2*pi*e_0)).*((xf - posx(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_xi1815 = (rho_i3/(2*pi*e_0)).*((xf - posx(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_xi1816 = (rho_i4/(2*pi*e_0)).*((xf - posx(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_xi1817 = (rho_i5/(2*pi*e_0)).*((xf - posx(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));
E_xi1819 = (rho_i7/(2*pi*e_0)).*((xf - posx(18,19))./((xf - posx(18,19)).^2 + (yf - posy(18,19)).^2));
E_xi1820 = (rho_i8/(2*pi*e_0)).*((xf - posx(18,20))./((xf - posx(18,20)).^2 + (yf - posy(18,20)).^2));
E_xi1821 = (rho_i9/(2*pi*e_0)).*((xf - posx(18,21))./((xf - posx(18,21)).^2 + (yf - posy(18,21)).^2));
E_xi1822 = (rho_i10/(2*pi*e_0)).*((xf - posx(18,22))./((xf - posx(18,22)).^2 + (yf - posy(18,22)).^2));
E_xi1823 = (rho_i11/(2*pi*e_0)).*((xf - posx(18,23))./((xf - posx(18,23)).^2 + (yf - posy(18,23)).^2));
E_xi1824 = (rho_i12/(2*pi*e_0)).*((xf - posx(18,24))./((xf - posx(18,24)).^2 + (yf - posy(18,24)).^2));

E_xi191 = (-rho_i1/(2*pi*e_0)).*((xf - posx(19,1))./((xf - posx(19,1)).^2 + (yf - posy(19,1)).^2));
E_xi192 = (-rho_i2/(2*pi*e_0)).*((xf - posx(19,2))./((xf - posx(19,2)).^2 + (yf - posy(19,2)).^2));
E_xi193 = (-rho_i3/(2*pi*e_0)).*((xf - posx(19,3))./((xf - posx(19,3)).^2 + (yf - posy(19,3)).^2));
E_xi194 = (-rho_i4/(2*pi*e_0)).*((xf - posx(19,4))./((xf - posx(19,4)).^2 + (yf - posy(19,4)).^2));
E_xi195 = (-rho_i5/(2*pi*e_0)).*((xf - posx(19,5))./((xf - posx(19,5)).^2 + (yf - posy(19,5)).^2));
E_xi196 = (-rho_i6/(2*pi*e_0)).*((xf - posx(19,6))./((xf - posx(19,6)).^2 + (yf - posy(19,6)).^2));
E_xi197 = (-rho_i7/(2*pi*e_0)).*((xf - posx(19,7))./((xf - posx(19,7)).^2 + (yf - posy(19,7)).^2));
E_xi198 = (-rho_i8/(2*pi*e_0)).*((xf - posx(19,8))./((xf - posx(19,8)).^2 + (yf - posy(19,8)).^2));
E_xi199 = (-rho_i9/(2*pi*e_0)).*((xf - posx(19,9))./((xf - posx(19,9)).^2 + (yf - posy(19,9)).^2));
E_xi1910 = (-rho_i10/(2*pi*e_0)).*((xf - posx(19,10))./((xf - posx(19,10)).^2 + (yf - posy(19,10)).^2));
E_xi1911 = (-rho_i11/(2*pi*e_0)).*((xf - posx(19,11))./((xf - posx(19,11)).^2 + (yf - posy(19,11)).^2));
E_xi1912 = (-rho_i12/(2*pi*e_0)).*((xf - posx(19,12))./((xf - posx(19,12)).^2 + (yf - posy(19,12)).^2));
E_xi1913 = (rho_i1/(2*pi*e_0)).*((xf - posx(19,13))./((xf - posx(19,13)).^2 + (yf - posy(19,13)).^2));
E_xi1914 = (rho_i2/(2*pi*e_0)).*((xf - posx(19,14))./((xf - posx(19,14)).^2 + (yf - posy(19,14)).^2));
E_xi1915 = (rho_i3/(2*pi*e_0)).*((xf - posx(19,15))./((xf - posx(19,15)).^2 + (yf - posy(19,15)).^2));
E_xi1916 = (rho_i4/(2*pi*e_0)).*((xf - posx(19,16))./((xf - posx(19,16)).^2 + (yf - posy(19,16)).^2));
E_xi1917 = (rho_i5/(2*pi*e_0)).*((xf - posx(19,17))./((xf - posx(19,17)).^2 + (yf - posy(19,17)).^2));
E_xi1918 = (rho_i6/(2*pi*e_0)).*((xf - posx(19,18))./((xf - posx(19,18)).^2 + (yf - posy(19,18)).^2));
E_xi1920 = (rho_i8/(2*pi*e_0)).*((xf - posx(19,20))./((xf - posx(19,20)).^2 + (yf - posy(19,20)).^2));
E_xi1921 = (rho_i9/(2*pi*e_0)).*((xf - posx(19,21))./((xf - posx(19,21)).^2 + (yf - posy(19,21)).^2));
E_xi1922 = (rho_i10/(2*pi*e_0)).*((xf - posx(19,22))./((xf - posx(19,22)).^2 + (yf - posy(19,22)).^2));
E_xi1923 = (rho_i11/(2*pi*e_0)).*((xf - posx(19,23))./((xf - posx(19,23)).^2 + (yf - posy(19,23)).^2));
E_xi1924 = (rho_i12/(2*pi*e_0)).*((xf - posx(19,24))./((xf - posx(19,24)).^2 + (yf - posy(19,24)).^2));

E_xi201 = (-rho_i1/(2*pi*e_0)).*((xf - posx(20,1))./((xf - posx(20,1)).^2 + (yf - posy(20,1)).^2));
E_xi202 = (-rho_i2/(2*pi*e_0)).*((xf - posx(20,2))./((xf - posx(20,2)).^2 + (yf - posy(20,2)).^2));
E_xi203 = (-rho_i3/(2*pi*e_0)).*((xf - posx(20,3))./((xf - posx(20,3)).^2 + (yf - posy(20,3)).^2));
E_xi204 = (-rho_i4/(2*pi*e_0)).*((xf - posx(20,4))./((xf - posx(20,4)).^2 + (yf - posy(20,4)).^2));
E_xi205 = (-rho_i5/(2*pi*e_0)).*((xf - posx(20,5))./((xf - posx(20,5)).^2 + (yf - posy(20,5)).^2));
E_xi206 = (-rho_i6/(2*pi*e_0)).*((xf - posx(20,6))./((xf - posx(20,6)).^2 + (yf - posy(20,6)).^2));
E_xi207 = (-rho_i7/(2*pi*e_0)).*((xf - posx(20,7))./((xf - posx(20,7)).^2 + (yf - posy(20,7)).^2));
E_xi208 = (-rho_i8/(2*pi*e_0)).*((xf - posx(20,8))./((xf - posx(20,8)).^2 + (yf - posy(20,8)).^2));
E_xi209 = (-rho_i9/(2*pi*e_0)).*((xf - posx(20,9))./((xf - posx(20,9)).^2 + (yf - posy(20,9)).^2));
E_xi2010 = (-rho_i10/(2*pi*e_0)).*((xf - posx(20,10))./((xf - posx(20,10)).^2 + (yf - posy(20,10)).^2));
E_xi2011 = (-rho_i11/(2*pi*e_0)).*((xf - posx(20,11))./((xf - posx(20,11)).^2 + (yf - posy(20,11)).^2));
E_xi2012 = (-rho_i12/(2*pi*e_0)).*((xf - posx(20,12))./((xf - posx(20,12)).^2 + (yf - posy(20,12)).^2));
E_xi2013 = (rho_i1/(2*pi*e_0)).*((xf - posx(20,13))./((xf - posx(20,13)).^2 + (yf - posy(20,13)).^2));
E_xi2014 = (rho_i2/(2*pi*e_0)).*((xf - posx(20,14))./((xf - posx(20,14)).^2 + (yf - posy(20,14)).^2));
E_xi2015 = (rho_i3/(2*pi*e_0)).*((xf - posx(20,15))./((xf - posx(20,15)).^2 + (yf - posy(20,15)).^2));
E_xi2016 = (rho_i4/(2*pi*e_0)).*((xf - posx(20,16))./((xf - posx(20,16)).^2 + (yf - posy(20,16)).^2));
E_xi2017 = (rho_i5/(2*pi*e_0)).*((xf - posx(20,17))./((xf - posx(20,17)).^2 + (yf - posy(20,17)).^2));
E_xi2018 = (rho_i6/(2*pi*e_0)).*((xf - posx(20,18))./((xf - posx(20,18)).^2 + (yf - posy(20,18)).^2));
E_xi2019 = (rho_i7/(2*pi*e_0)).*((xf - posx(20,19))./((xf - posx(20,19)).^2 + (yf - posy(20,19)).^2));
E_xi2021 = (rho_i9/(2*pi*e_0)).*((xf - posx(20,21))./((xf - posx(20,21)).^2 + (yf - posy(20,21)).^2));
E_xi2022 = (rho_i10/(2*pi*e_0)).*((xf - posx(20,22))./((xf - posx(20,22)).^2 + (yf - posy(20,22)).^2));
E_xi2023 = (rho_i11/(2*pi*e_0)).*((xf - posx(20,23))./((xf - posx(20,23)).^2 + (yf - posy(20,23)).^2));
E_xi2024 = (rho_i12/(2*pi*e_0)).*((xf - posx(20,24))./((xf - posx(20,24)).^2 + (yf - posy(20,24)).^2));

E_xi21_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(21,1))./((xf - posx(21,1)).^2 + (yf - posy(21,1)).^2));
E_xi21_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(21,2))./((xf - posx(21,2)).^2 + (yf - posy(21,2)).^2));
E_xi21_3 = (-rho_i3/(2*pi*e_0)).*((xf - posx(21,3))./((xf - posx(21,3)).^2 + (yf - posy(21,3)).^2));
E_xi214 = (-rho_i4/(2*pi*e_0)).*((xf - posx(21,4))./((xf - posx(21,4)).^2 + (yf - posy(21,4)).^2));
E_xi215 = (-rho_i5/(2*pi*e_0)).*((xf - posx(21,5))./((xf - posx(21,5)).^2 + (yf - posy(21,5)).^2));
E_xi216 = (-rho_i6/(2*pi*e_0)).*((xf - posx(21,6))./((xf - posx(21,6)).^2 + (yf - posy(21,6)).^2));
E_xi217 = (-rho_i7/(2*pi*e_0)).*((xf - posx(21,7))./((xf - posx(21,7)).^2 + (yf - posy(21,7)).^2));
E_xi218 = (-rho_i8/(2*pi*e_0)).*((xf - posx(21,8))./((xf - posx(21,8)).^2 + (yf - posy(21,8)).^2));
E_xi219 = (-rho_i9/(2*pi*e_0)).*((xf - posx(21,9))./((xf - posx(21,9)).^2 + (yf - posy(21,9)).^2));
E_xi2110 = (-rho_i10/(2*pi*e_0)).*((xf - posx(21,10))./((xf - posx(21,10)).^2 + (yf - posy(21,10)).^2));
E_xi2111 = (-rho_i11/(2*pi*e_0)).*((xf - posx(21,11))./((xf - posx(21,11)).^2 + (yf - posy(21,11)).^2));
E_xi2112 = (-rho_i12/(2*pi*e_0)).*((xf - posx(21,12))./((xf - posx(21,12)).^2 + (yf - posy(21,12)).^2));
E_xi2113 = (rho_i1/(2*pi*e_0)).*((xf - posx(21,13))./((xf - posx(21,13)).^2 + (yf - posy(21,13)).^2));
E_xi2114 = (rho_i2/(2*pi*e_0)).*((xf - posx(21,14))./((xf - posx(21,14)).^2 + (yf - posy(21,14)).^2));
E_xi2115 = (rho_i3/(2*pi*e_0)).*((xf - posx(21,15))./((xf - posx(21,15)).^2 + (yf - posy(21,15)).^2));
E_xi2116 = (rho_i4/(2*pi*e_0)).*((xf - posx(21,16))./((xf - posx(21,16)).^2 + (yf - posy(21,16)).^2));
E_xi2117 = (rho_i5/(2*pi*e_0)).*((xf - posx(21,17))./((xf - posx(21,17)).^2 + (yf - posy(21,17)).^2));
E_xi2118 = (rho_i6/(2*pi*e_0)).*((xf - posx(21,18))./((xf - posx(21,18)).^2 + (yf - posy(21,18)).^2));
E_xi2119 = (rho_i7/(2*pi*e_0)).*((xf - posx(21,19))./((xf - posx(21,19)).^2 + (yf - posy(21,19)).^2));
E_xi2120 = (rho_i8/(2*pi*e_0)).*((xf - posx(21,20))./((xf - posx(21,20)).^2 + (yf - posy(21,20)).^2));
E_xi2122 = (rho_i10/(2*pi*e_0)).*((xf - posx(21,22))./((xf - posx(21,22)).^2 + (yf - posy(21,22)).^2));
E_xi2123 = (rho_i11/(2*pi*e_0)).*((xf - posx(21,23))./((xf - posx(21,23)).^2 + (yf - posy(21,23)).^2));
E_xi2124 = (rho_i12/(2*pi*e_0)).*((xf - posx(21,24))./((xf - posx(21,24)).^2 + (yf - posy(21,24)).^2));

E_xi22_1 = (-rho_i1/(2*pi*e_0)).*((xf - posx(22,1))./((xf - posx(22,1)).^2 + (yf - posy(22,1)).^2));
E_xi22_2 = (-rho_i2/(2*pi*e_0)).*((xf - posx(22,2))./((xf - posx(22,2)).^2 + (yf - posy(22,2)).^2));
E_xi22_3 = (-rho_i3/(2*pi*e_0)).*((xf - posx(22,3))./((xf - posx(22,3)).^2 + (yf - posy(22,3)).^2));
E_xi22_4 = (-rho_i4/(2*pi*e_0)).*((xf - posx(22,4))./((xf - posx(22,4)).^2 + (yf - posy(22,4)).^2));
E_xi225 = (-rho_i5/(2*pi*e_0)).*((xf - posx(22,5))./((xf - posx(22,5)).^2 + (yf - posy(22,5)).^2));
E_xi226 = (-rho_i6/(2*pi*e_0)).*((xf - posx(22,6))./((xf - posx(22,6)).^2 + (yf - posy(22,6)).^2));
E_xi227 = (-rho_i7/(2*pi*e_0)).*((xf - posx(22,7))./((xf - posx(22,7)).^2 + (yf - posy(22,7)).^2));
E_xi228 = (-rho_i8/(2*pi*e_0)).*((xf - posx(22,8))./((xf - posx(22,8)).^2 + (yf - posy(22,8)).^2));
E_xi229 = (-rho_i9/(2*pi*e_0)).*((xf - posx(22,9))./((xf - posx(22,9)).^2 + (yf - posy(22,9)).^2));
E_xi2210 = (-rho_i10/(2*pi*e_0)).*((xf - posx(22,10))./((xf - posx(22,10)).^2 + (yf - posy(22,10)).^2));
E_xi2211 = (-rho_i11/(2*pi*e_0)).*((xf - posx(22,11))./((xf - posx(22,11)).^2 + (yf - posy(22,11)).^2));
E_xi2212 = (-rho_i12/(2*pi*e_0)).*((xf - posx(22,12))./((xf - posx(22,12)).^2 + (yf - posy(22,12)).^2));
E_xi2213 = (rho_i1/(2*pi*e_0)).*((xf - posx(22,13))./((xf - posx(22,13)).^2 + (yf - posy(22,13)).^2));
E_xi2214 = (rho_i2/(2*pi*e_0)).*((xf - posx(22,14))./((xf - posx(22,14)).^2 + (yf - posy(22,14)).^2));
E_xi2215 = (rho_i3/(2*pi*e_0)).*((xf - posx(22,15))./((xf - posx(22,15)).^2 + (yf - posy(22,15)).^2));
E_xi2216 = (rho_i4/(2*pi*e_0)).*((xf - posx(22,16))./((xf - posx(22,16)).^2 + (yf - posy(22,16)).^2));
E_xi2217 = (rho_i5/(2*pi*e_0)).*((xf - posx(22,17))./((xf - posx(22,17)).^2 + (yf - posy(22,17)).^2));
E_xi2218 = (rho_i6/(2*pi*e_0)).*((xf - posx(22,18))./((xf - posx(22,18)).^2 + (yf - posy(22,18)).^2));
E_xi2219 = (rho_i7/(2*pi*e_0)).*((xf - posx(22,19))./((xf - posx(22,19)).^2 + (yf - posy(22,19)).^2));
E_xi2220 = (rho_i8/(2*pi*e_0)).*((xf - posx(22,20))./((xf - posx(22,20)).^2 + (yf - posy(22,20)).^2));
E_xi2221 = (rho_i9/(2*pi*e_0)).*((xf - posx(22,21))./((xf - posx(22,21)).^2 + (yf - posy(22,21)).^2));
E_xi2223 = (rho_i11/(2*pi*e_0)).*((xf - posx(22,23))./((xf - posx(22,23)).^2 + (yf - posy(22,23)).^2));
E_xi2224 = (rho_i12/(2*pi*e_0)).*((xf - posx(22,24))./((xf - posx(22,24)).^2 + (yf - posy(22,24)).^2));

E_xi231 = (-rho_i1/(2*pi*e_0)).*((xf - posx(23,1))./((xf - posx(23,1)).^2 + (yf - posy(23,1)).^2));
E_xi232 = (-rho_i2/(2*pi*e_0)).*((xf - posx(23,2))./((xf - posx(23,2)).^2 + (yf - posy(23,2)).^2));
E_xi233 = (-rho_i3/(2*pi*e_0)).*((xf - posx(23,3))./((xf - posx(23,3)).^2 + (yf - posy(23,3)).^2));
E_xi234 = (-rho_i4/(2*pi*e_0)).*((xf - posx(23,4))./((xf - posx(23,4)).^2 + (yf - posy(23,4)).^2));
E_xi235 = (-rho_i5/(2*pi*e_0)).*((xf - posx(23,5))./((xf - posx(23,5)).^2 + (yf - posy(23,5)).^2));
E_xi236 = (-rho_i6/(2*pi*e_0)).*((xf - posx(23,6))./((xf - posx(23,6)).^2 + (yf - posy(23,6)).^2));
E_xi237 = (-rho_i7/(2*pi*e_0)).*((xf - posx(23,7))./((xf - posx(23,7)).^2 + (yf - posy(23,7)).^2));
E_xi238 = (-rho_i8/(2*pi*e_0)).*((xf - posx(23,8))./((xf - posx(23,8)).^2 + (yf - posy(23,8)).^2));
E_xi239 = (-rho_i9/(2*pi*e_0)).*((xf - posx(23,9))./((xf - posx(23,9)).^2 + (yf - posy(23,9)).^2));
E_xi2310 = (-rho_i10/(2*pi*e_0)).*((xf - posx(23,10))./((xf - posx(23,10)).^2 + (yf - posy(23,10)).^2));
E_xi2311 = (-rho_i11/(2*pi*e_0)).*((xf - posx(23,11))./((xf - posx(23,11)).^2 + (yf - posy(23,11)).^2));
E_xi2312 = (-rho_i12/(2*pi*e_0)).*((xf - posx(23,12))./((xf - posx(23,12)).^2 + (yf - posy(23,12)).^2));
E_xi2313 = (rho_i1/(2*pi*e_0)).*((xf - posx(23,13))./((xf - posx(23,13)).^2 + (yf - posy(23,13)).^2));
E_xi2314 = (rho_i2/(2*pi*e_0)).*((xf - posx(23,14))./((xf - posx(23,14)).^2 + (yf - posy(23,14)).^2));
E_xi2315 = (rho_i3/(2*pi*e_0)).*((xf - posx(23,15))./((xf - posx(23,15)).^2 + (yf - posy(23,15)).^2));
E_xi2316 = (rho_i4/(2*pi*e_0)).*((xf - posx(23,16))./((xf - posx(23,16)).^2 + (yf - posy(23,16)).^2));
E_xi2317 = (rho_i5/(2*pi*e_0)).*((xf - posx(23,17))./((xf - posx(23,17)).^2 + (yf - posy(23,17)).^2));
E_xi2318 = (rho_i6/(2*pi*e_0)).*((xf - posx(23,18))./((xf - posx(23,18)).^2 + (yf - posy(23,18)).^2));
E_xi2319 = (rho_i7/(2*pi*e_0)).*((xf - posx(23,19))./((xf - posx(23,19)).^2 + (yf - posy(23,19)).^2));
E_xi2320 = (rho_i8/(2*pi*e_0)).*((xf - posx(23,20))./((xf - posx(23,20)).^2 + (yf - posy(23,20)).^2));
E_xi2321 = (rho_i9/(2*pi*e_0)).*((xf - posx(23,21))./((xf - posx(23,21)).^2 + (yf - posy(23,21)).^2));
E_xi2322 = (rho_i10/(2*pi*e_0)).*((xf - posx(23,22))./((xf - posx(23,22)).^2 + (yf - posy(23,22)).^2));
E_xi2324 = (rho_i12/(2*pi*e_0)).*((xf - posx(23,24))./((xf - posx(23,24)).^2 + (yf - posy(23,24)).^2));

E_xi241 = (-rho_i1/(2*pi*e_0)).*((xf - posx(24,1))./((xf - posx(24,1)).^2 + (yf - posy(24,1)).^2));
E_xi242 = (-rho_i2/(2*pi*e_0)).*((xf - posx(24,2))./((xf - posx(24,2)).^2 + (yf - posy(24,2)).^2));
E_xi243 = (-rho_i3/(2*pi*e_0)).*((xf - posx(24,3))./((xf - posx(24,3)).^2 + (yf - posy(24,3)).^2));
E_xi244 = (-rho_i4/(2*pi*e_0)).*((xf - posx(24,4))./((xf - posx(24,4)).^2 + (yf - posy(24,4)).^2));
E_xi245 = (-rho_i5/(2*pi*e_0)).*((xf - posx(24,5))./((xf - posx(24,5)).^2 + (yf - posy(24,5)).^2));
E_xi246 = (-rho_i6/(2*pi*e_0)).*((xf - posx(24,6))./((xf - posx(24,6)).^2 + (yf - posy(24,6)).^2));
E_xi247 = (-rho_i7/(2*pi*e_0)).*((xf - posx(24,7))./((xf - posx(24,7)).^2 + (yf - posy(24,7)).^2));
E_xi248 = (-rho_i8/(2*pi*e_0)).*((xf - posx(24,8))./((xf - posx(24,8)).^2 + (yf - posy(24,8)).^2));
E_xi249 = (-rho_i9/(2*pi*e_0)).*((xf - posx(24,9))./((xf - posx(24,9)).^2 + (yf - posy(24,9)).^2));
E_xi2410 = (-rho_i10/(2*pi*e_0)).*((xf - posx(24,10))./((xf - posx(24,10)).^2 + (yf - posy(24,10)).^2));
E_xi2411 = (-rho_i11/(2*pi*e_0)).*((xf - posx(24,11))./((xf - posx(24,11)).^2 + (yf - posy(24,11)).^2));
E_xi2412 = (-rho_i12/(2*pi*e_0)).*((xf - posx(24,12))./((xf - posx(24,12)).^2 + (yf - posy(24,12)).^2));
E_xi2413 = (rho_i1/(2*pi*e_0)).*((xf - posx(24,13))./((xf - posx(24,13)).^2 + (yf - posy(24,13)).^2));
E_xi2414 = (rho_i2/(2*pi*e_0)).*((xf - posx(24,14))./((xf - posx(24,14)).^2 + (yf - posy(24,14)).^2));
E_xi2415 = (rho_i3/(2*pi*e_0)).*((xf - posx(24,15))./((xf - posx(24,15)).^2 + (yf - posy(24,15)).^2));
E_xi2416 = (rho_i4/(2*pi*e_0)).*((xf - posx(24,16))./((xf - posx(24,16)).^2 + (yf - posy(24,16)).^2));
E_xi2417 = (rho_i5/(2*pi*e_0)).*((xf - posx(24,17))./((xf - posx(24,17)).^2 + (yf - posy(24,17)).^2));
E_xi2418 = (rho_i6/(2*pi*e_0)).*((xf - posx(24,18))./((xf - posx(24,18)).^2 + (yf - posy(24,18)).^2));
E_xi2419 = (rho_i7/(2*pi*e_0)).*((xf - posx(24,19))./((xf - posx(24,19)).^2 + (yf - posy(24,19)).^2));
E_xi2420 = (rho_i8/(2*pi*e_0)).*((xf - posx(24,20))./((xf - posx(24,20)).^2 + (yf - posy(24,20)).^2));
E_xi2421 = (rho_i9/(2*pi*e_0)).*((xf - posx(24,21))./((xf - posx(24,21)).^2 + (yf - posy(24,21)).^2));
E_xi2422 = (rho_i10/(2*pi*e_0)).*((xf - posx(24,22))./((xf - posx(24,22)).^2 + (yf - posy(24,22)).^2));
E_xi2423 = (rho_i11/(2*pi*e_0)).*((xf - posx(24,23))./((xf - posx(24,23)).^2 + (yf - posy(24,23)).^2));

%% E_yr componente y real campo elétrico condutor 2 fase b

E_yr12 = (-rho_r2/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yr13 = (-rho_r3/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yr14 = (-rho_r4/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yr15 = (-rho_r5/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yr16 = (-rho_r6/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_yr17 = (-rho_r7/(2*pi*e_0)).*((yf - posy(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_yr18 = (-rho_r8/(2*pi*e_0)).*((yf - posy(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_yr19 = (-rho_r9/(2*pi*e_0)).*((yf - posy(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_yr110 = (-rho_r10/(2*pi*e_0)).*((yf - posy(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_yr1_11 = (-rho_r11/(2*pi*e_0)).*((yf - posy(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_yr1_12 = (-rho_r12/(2*pi*e_0)).*((yf - posy(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_yr1_13 = (rho_r1/(2*pi*e_0)).*((yf - posy(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_yr1_14 = (rho_r2/(2*pi*e_0)).*((yf - posy(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_yr1_15 = (rho_r3/(2*pi*e_0)).*((yf - posy(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_yr1_16 = (rho_r4/(2*pi*e_0)).*((yf - posy(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_yr1_17 = (rho_r5/(2*pi*e_0)).*((yf - posy(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_yr1_18 = (rho_r6/(2*pi*e_0)).*((yf - posy(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));
E_yr1_19 = (rho_r7/(2*pi*e_0)).*((yf - posy(1,19))./((xf - posx(1,19)).^2 + (yf - posy(1,19)).^2));
E_yr1_20 = (rho_r8/(2*pi*e_0)).*((yf - posy(1,20))./((xf - posx(1,20)).^2 + (yf - posy(1,20)).^2));
E_yr1_21 = (rho_r9/(2*pi*e_0)).*((yf - posy(1,21))./((xf - posx(1,21)).^2 + (yf - posy(1,21)).^2));
E_yr1_22 = (rho_r10/(2*pi*e_0)).*((yf - posy(1,22))./((xf - posx(1,22)).^2 + (yf - posy(1,22)).^2));
E_yr1_23 = (rho_r11/(2*pi*e_0)).*((yf - posy(1,23))./((xf - posx(1,23)).^2 + (yf - posy(1,23)).^2));
E_yr1_24 = (rho_r12/(2*pi*e_0)).*((yf - posy(1,24))./((xf - posx(1,24)).^2 + (yf - posy(1,24)).^2));

E_yr21 = (-rho_r1/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yr23 = (-rho_r3/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yr24 = (-rho_r4/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yr25 = (-rho_r5/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yr26 = (-rho_r6/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_yr27 = (-rho_r7/(2*pi*e_0)).*((yf - posy(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_yr28 = (-rho_r8/(2*pi*e_0)).*((yf - posy(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_yr29 = (-rho_r9/(2*pi*e_0)).*((yf - posy(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_yr210 = (-rho_r10/(2*pi*e_0)).*((yf - posy(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_yr2_11 = (-rho_r11/(2*pi*e_0)).*((yf - posy(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_yr2_12 = (-rho_r12/(2*pi*e_0)).*((yf - posy(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_yr2_13 = (rho_r1/(2*pi*e_0)).*((yf - posy(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_yr2_14 = (rho_r2/(2*pi*e_0)).*((yf - posy(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_yr2_15 = (rho_r3/(2*pi*e_0)).*((yf - posy(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_yr2_16 = (rho_r4/(2*pi*e_0)).*((yf - posy(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_yr2_17 = (rho_r5/(2*pi*e_0)).*((yf - posy(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_yr2_18 = (rho_r6/(2*pi*e_0)).*((yf - posy(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));
E_yr2_19 = (rho_r7/(2*pi*e_0)).*((yf - posy(2,19))./((xf - posx(2,19)).^2 + (yf - posy(2,19)).^2));
E_yr2_20 = (rho_r8/(2*pi*e_0)).*((yf - posy(2,20))./((xf - posx(2,20)).^2 + (yf - posy(2,20)).^2));
E_yr2_21 = (rho_r9/(2*pi*e_0)).*((yf - posy(2,21))./((xf - posx(2,21)).^2 + (yf - posy(2,21)).^2));
E_yr2_22 = (rho_r10/(2*pi*e_0)).*((yf - posy(2,22))./((xf - posx(2,22)).^2 + (yf - posy(2,22)).^2));
E_yr2_23 = (rho_r11/(2*pi*e_0)).*((yf - posy(2,23))./((xf - posx(2,23)).^2 + (yf - posy(2,23)).^2));
E_yr2_24 = (rho_r12/(2*pi*e_0)).*((yf - posy(2,24))./((xf - posx(2,24)).^2 + (yf - posy(2,24)).^2));

E_yr31 = (-rho_r1/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yr32 = (-rho_r2/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yr34 = (-rho_r4/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yr35 = (-rho_r5/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yr36 = (-rho_r6/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_yr37 = (-rho_r7/(2*pi*e_0)).*((yf - posy(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_yr38 = (-rho_r8/(2*pi*e_0)).*((yf - posy(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_yr39 = (-rho_r9/(2*pi*e_0)).*((yf - posy(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_yr310 = (-rho_r10/(2*pi*e_0)).*((yf - posy(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_yr311 = (-rho_r11/(2*pi*e_0)).*((yf - posy(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_yr312 = (-rho_r12/(2*pi*e_0)).*((yf - posy(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_yr313 = (rho_r1/(2*pi*e_0)).*((yf - posy(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_yr314 = (rho_r2/(2*pi*e_0)).*((yf - posy(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_yr315 = (rho_r3/(2*pi*e_0)).*((yf - posy(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_yr316 = (rho_r4/(2*pi*e_0)).*((yf - posy(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_yr317 = (rho_r5/(2*pi*e_0)).*((yf - posy(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_yr318 = (rho_r6/(2*pi*e_0)).*((yf - posy(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));
E_yr319 = (rho_r7/(2*pi*e_0)).*((yf - posy(3,19))./((xf - posx(3,19)).^2 + (yf - posy(3,19)).^2));
E_yr320 = (rho_r8/(2*pi*e_0)).*((yf - posy(3,20))./((xf - posx(3,20)).^2 + (yf - posy(3,20)).^2));
E_yr321 = (rho_r9/(2*pi*e_0)).*((yf - posy(3,21))./((xf - posx(3,21)).^2 + (yf - posy(3,21)).^2));
E_yr322 = (rho_r10/(2*pi*e_0)).*((yf - posy(3,22))./((xf - posx(3,22)).^2 + (yf - posy(3,22)).^2));
E_yr323 = (rho_r11/(2*pi*e_0)).*((yf - posy(3,23))./((xf - posx(3,23)).^2 + (yf - posy(3,23)).^2));
E_yr324 = (rho_r12/(2*pi*e_0)).*((yf - posy(3,24))./((xf - posx(3,24)).^2 + (yf - posy(3,24)).^2));

E_yr41 = (-rho_r1/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yr42 = (-rho_r2/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yr43 = (-rho_r3/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yr45 = (-rho_r5/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yr46 = (-rho_r6/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_yr47 = (-rho_r7/(2*pi*e_0)).*((yf - posy(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_yr48 = (-rho_r8/(2*pi*e_0)).*((yf - posy(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_yr49 = (-rho_r9/(2*pi*e_0)).*((yf - posy(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_yr410 = (-rho_r10/(2*pi*e_0)).*((yf - posy(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_yr411 = (-rho_r11/(2*pi*e_0)).*((yf - posy(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_yr412 = (-rho_r12/(2*pi*e_0)).*((yf - posy(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_yr413 = (rho_r1/(2*pi*e_0)).*((yf - posy(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_yr414 = (rho_r2/(2*pi*e_0)).*((yf - posy(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_yr415 = (rho_r3/(2*pi*e_0)).*((yf - posy(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_yr416 = (rho_r4/(2*pi*e_0)).*((yf - posy(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_yr417 = (rho_r5/(2*pi*e_0)).*((yf - posy(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_yr418 = (rho_r6/(2*pi*e_0)).*((yf - posy(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));
E_yr419 = (rho_r7/(2*pi*e_0)).*((yf - posy(4,19))./((xf - posx(4,19)).^2 + (yf - posy(4,19)).^2));
E_yr420 = (rho_r8/(2*pi*e_0)).*((yf - posy(4,20))./((xf - posx(4,20)).^2 + (yf - posy(4,20)).^2));
E_yr421 = (rho_r9/(2*pi*e_0)).*((yf - posy(4,21))./((xf - posx(4,21)).^2 + (yf - posy(4,21)).^2));
E_yr422 = (rho_r10/(2*pi*e_0)).*((yf - posy(4,22))./((xf - posx(4,22)).^2 + (yf - posy(4,22)).^2));
E_yr423 = (rho_r11/(2*pi*e_0)).*((yf - posy(4,23))./((xf - posx(4,23)).^2 + (yf - posy(4,23)).^2));
E_yr424 = (rho_r12/(2*pi*e_0)).*((yf - posy(4,24))./((xf - posx(4,24)).^2 + (yf - posy(4,24)).^2));

E_yr51 = (-rho_r1/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yr52 = (-rho_r2/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yr53 = (-rho_r3/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yr54 = (-rho_r4/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yr56 = (-rho_r6/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_yr57 = (-rho_r7/(2*pi*e_0)).*((yf - posy(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_yr58 = (-rho_r8/(2*pi*e_0)).*((yf - posy(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_yr59 = (-rho_r9/(2*pi*e_0)).*((yf - posy(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_yr510 = (-rho_r10/(2*pi*e_0)).*((yf - posy(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_yr511 = (-rho_r11/(2*pi*e_0)).*((yf - posy(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_yr512 = (-rho_r12/(2*pi*e_0)).*((yf - posy(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_yr513 = (rho_r1/(2*pi*e_0)).*((yf - posy(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_yr514 = (rho_r2/(2*pi*e_0)).*((yf - posy(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_yr515 = (rho_r3/(2*pi*e_0)).*((yf - posy(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_yr516 = (rho_r4/(2*pi*e_0)).*((yf - posy(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_yr517 = (rho_r5/(2*pi*e_0)).*((yf - posy(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_yr518 = (rho_r6/(2*pi*e_0)).*((yf - posy(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));
E_yr519 = (rho_r7/(2*pi*e_0)).*((yf - posy(5,19))./((xf - posx(5,19)).^2 + (yf - posy(5,19)).^2));
E_yr520 = (rho_r8/(2*pi*e_0)).*((yf - posy(5,20))./((xf - posx(5,20)).^2 + (yf - posy(5,20)).^2));
E_yr521 = (rho_r9/(2*pi*e_0)).*((yf - posy(5,21))./((xf - posx(5,21)).^2 + (yf - posy(5,21)).^2));
E_yr522 = (rho_r10/(2*pi*e_0)).*((yf - posy(5,22))./((xf - posx(5,22)).^2 + (yf - posy(5,22)).^2));
E_yr523 = (rho_r11/(2*pi*e_0)).*((yf - posy(5,23))./((xf - posx(5,23)).^2 + (yf - posy(5,23)).^2));
E_yr524 = (rho_r12/(2*pi*e_0)).*((yf - posy(5,24))./((xf - posx(5,24)).^2 + (yf - posy(5,24)).^2));

E_yr61 = (-rho_r1/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yr62 = (-rho_r2/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yr63 = (-rho_r3/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yr64 = (-rho_r4/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yr65 = (-rho_r5/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_yr67 = (-rho_r7/(2*pi*e_0)).*((yf - posy(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_yr68 = (-rho_r8/(2*pi*e_0)).*((yf - posy(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_yr69 = (-rho_r9/(2*pi*e_0)).*((yf - posy(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_yr610 = (-rho_r10/(2*pi*e_0)).*((yf - posy(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_yr611 = (-rho_r11/(2*pi*e_0)).*((yf - posy(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_yr612 = (-rho_r12/(2*pi*e_0)).*((yf - posy(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_yr613 = (rho_r1/(2*pi*e_0)).*((yf - posy(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_yr614 = (rho_r2/(2*pi*e_0)).*((yf - posy(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_yr615 = (rho_r3/(2*pi*e_0)).*((yf - posy(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_yr616 = (rho_r4/(2*pi*e_0)).*((yf - posy(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_yr617 = (rho_r5/(2*pi*e_0)).*((yf - posy(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_yr618 = (rho_r6/(2*pi*e_0)).*((yf - posy(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));
E_yr619 = (rho_r7/(2*pi*e_0)).*((yf - posy(6,19))./((xf - posx(6,19)).^2 + (yf - posy(6,19)).^2));
E_yr620 = (rho_r8/(2*pi*e_0)).*((yf - posy(6,20))./((xf - posx(6,20)).^2 + (yf - posy(6,20)).^2));
E_yr621 = (rho_r9/(2*pi*e_0)).*((yf - posy(6,21))./((xf - posx(6,21)).^2 + (yf - posy(6,21)).^2));
E_yr622 = (rho_r10/(2*pi*e_0)).*((yf - posy(6,22))./((xf - posx(6,22)).^2 + (yf - posy(6,22)).^2));
E_yr623 = (rho_r11/(2*pi*e_0)).*((yf - posy(6,23))./((xf - posx(6,23)).^2 + (yf - posy(6,23)).^2));
E_yr624 = (rho_r12/(2*pi*e_0)).*((yf - posy(6,24))./((xf - posx(6,24)).^2 + (yf - posy(6,24)).^2));

E_yr71 = (-rho_r1/(2*pi*e_0)).*((yf - posy(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_yr72 = (-rho_r2/(2*pi*e_0)).*((yf - posy(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_yr73 = (-rho_r3/(2*pi*e_0)).*((yf - posy(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_yr74 = (-rho_r4/(2*pi*e_0)).*((yf - posy(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_yr75 = (-rho_r5/(2*pi*e_0)).*((yf - posy(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_yr76 = (-rho_r6/(2*pi*e_0)).*((yf - posy(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_yr78 = (-rho_r8/(2*pi*e_0)).*((yf - posy(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_yr79 = (-rho_r9/(2*pi*e_0)).*((yf - posy(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_yr710 = (-rho_r10/(2*pi*e_0)).*((yf - posy(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_yr711 = (-rho_r11/(2*pi*e_0)).*((yf - posy(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_yr712 = (-rho_r12/(2*pi*e_0)).*((yf - posy(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_yr713 = (rho_r1/(2*pi*e_0)).*((yf - posy(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_yr714 = (rho_r2/(2*pi*e_0)).*((yf - posy(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_yr715 = (rho_r3/(2*pi*e_0)).*((yf - posy(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_yr716 = (rho_r4/(2*pi*e_0)).*((yf - posy(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_yr717 = (rho_r5/(2*pi*e_0)).*((yf - posy(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_yr718 = (rho_r6/(2*pi*e_0)).*((yf - posy(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));
E_yr719 = (rho_r7/(2*pi*e_0)).*((yf - posy(7,19))./((xf - posx(7,19)).^2 + (yf - posy(7,19)).^2));
E_yr720 = (rho_r8/(2*pi*e_0)).*((yf - posy(7,20))./((xf - posx(7,20)).^2 + (yf - posy(7,20)).^2));
E_yr721 = (rho_r9/(2*pi*e_0)).*((yf - posy(7,21))./((xf - posx(7,21)).^2 + (yf - posy(7,21)).^2));
E_yr722 = (rho_r10/(2*pi*e_0)).*((yf - posy(7,22))./((xf - posx(7,22)).^2 + (yf - posy(7,22)).^2));
E_yr723 = (rho_r11/(2*pi*e_0)).*((yf - posy(7,23))./((xf - posx(7,23)).^2 + (yf - posy(7,23)).^2));
E_yr724 = (rho_r12/(2*pi*e_0)).*((yf - posy(7,24))./((xf - posx(7,24)).^2 + (yf - posy(7,24)).^2));

E_yr81 = (-rho_r1/(2*pi*e_0)).*((yf - posy(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_yr82 = (-rho_r2/(2*pi*e_0)).*((yf - posy(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_yr83 = (-rho_r3/(2*pi*e_0)).*((yf - posy(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_yr84 = (-rho_r4/(2*pi*e_0)).*((yf - posy(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_yr85 = (-rho_r5/(2*pi*e_0)).*((yf - posy(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_yr86 = (-rho_r6/(2*pi*e_0)).*((yf - posy(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_yr87 = (-rho_r7/(2*pi*e_0)).*((yf - posy(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_yr89 = (-rho_r9/(2*pi*e_0)).*((yf - posy(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_yr810 = (-rho_r10/(2*pi*e_0)).*((yf - posy(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_yr811 = (-rho_r11/(2*pi*e_0)).*((yf - posy(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_yr812 = (-rho_r12/(2*pi*e_0)).*((yf - posy(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_yr813 = (rho_r1/(2*pi*e_0)).*((yf - posy(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_yr814 = (rho_r2/(2*pi*e_0)).*((yf - posy(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_yr815 = (rho_r3/(2*pi*e_0)).*((yf - posy(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_yr816 = (rho_r4/(2*pi*e_0)).*((yf - posy(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_yr817 = (rho_r5/(2*pi*e_0)).*((yf - posy(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_yr818 = (rho_r6/(2*pi*e_0)).*((yf - posy(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));
E_yr819 = (rho_r7/(2*pi*e_0)).*((yf - posy(8,19))./((xf - posx(8,19)).^2 + (yf - posy(8,19)).^2));
E_yr820 = (rho_r8/(2*pi*e_0)).*((yf - posy(8,20))./((xf - posx(8,20)).^2 + (yf - posy(8,20)).^2));
E_yr821 = (rho_r9/(2*pi*e_0)).*((yf - posy(8,21))./((xf - posx(8,21)).^2 + (yf - posy(8,21)).^2));
E_yr822 = (rho_r10/(2*pi*e_0)).*((yf - posy(8,22))./((xf - posx(8,22)).^2 + (yf - posy(8,22)).^2));
E_yr823 = (rho_r11/(2*pi*e_0)).*((yf - posy(8,23))./((xf - posx(8,23)).^2 + (yf - posy(8,23)).^2));
E_yr824 = (rho_r12/(2*pi*e_0)).*((yf - posy(8,24))./((xf - posx(8,24)).^2 + (yf - posy(8,24)).^2));

E_yr91 = (-rho_r1/(2*pi*e_0)).*((yf - posy(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_yr92 = (-rho_r2/(2*pi*e_0)).*((yf - posy(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_yr93 = (-rho_r3/(2*pi*e_0)).*((yf - posy(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_yr94 = (-rho_r4/(2*pi*e_0)).*((yf - posy(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_yr95 = (-rho_r5/(2*pi*e_0)).*((yf - posy(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_yr96 = (-rho_r6/(2*pi*e_0)).*((yf - posy(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_yr97 = (-rho_r7/(2*pi*e_0)).*((yf - posy(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_yr98 = (-rho_r8/(2*pi*e_0)).*((yf - posy(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_yr910 = (-rho_r10/(2*pi*e_0)).*((yf - posy(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_yr911 = (-rho_r11/(2*pi*e_0)).*((yf - posy(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_yr912 = (-rho_r12/(2*pi*e_0)).*((yf - posy(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_yr913 = (rho_r1/(2*pi*e_0)).*((yf - posy(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_yr914 = (rho_r2/(2*pi*e_0)).*((yf - posy(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_yr915 = (rho_r3/(2*pi*e_0)).*((yf - posy(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_yr916 = (rho_r4/(2*pi*e_0)).*((yf - posy(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_yr917 = (rho_r5/(2*pi*e_0)).*((yf - posy(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_yr918 = (rho_r6/(2*pi*e_0)).*((yf - posy(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));
E_yr919 = (rho_r7/(2*pi*e_0)).*((yf - posy(9,19))./((xf - posx(9,19)).^2 + (yf - posy(9,19)).^2));
E_yr920 = (rho_r8/(2*pi*e_0)).*((yf - posy(9,20))./((xf - posx(9,20)).^2 + (yf - posy(9,20)).^2));
E_yr921 = (rho_r9/(2*pi*e_0)).*((yf - posy(9,21))./((xf - posx(9,21)).^2 + (yf - posy(9,21)).^2));
E_yr922 = (rho_r10/(2*pi*e_0)).*((yf - posy(9,22))./((xf - posx(9,22)).^2 + (yf - posy(9,22)).^2));
E_yr923 = (rho_r11/(2*pi*e_0)).*((yf - posy(9,23))./((xf - posx(9,23)).^2 + (yf - posy(9,23)).^2));
E_yr924 = (rho_r12/(2*pi*e_0)).*((yf - posy(9,24))./((xf - posx(9,24)).^2 + (yf - posy(9,24)).^2));

E_yr101 = (-rho_r1/(2*pi*e_0)).*((yf - posy(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_yr102 = (-rho_r2/(2*pi*e_0)).*((yf - posy(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_yr103 = (-rho_r3/(2*pi*e_0)).*((yf - posy(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_yr104 = (-rho_r4/(2*pi*e_0)).*((yf - posy(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_yr105 = (-rho_r5/(2*pi*e_0)).*((yf - posy(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_yr106 = (-rho_r6/(2*pi*e_0)).*((yf - posy(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_yr107 = (-rho_r7/(2*pi*e_0)).*((yf - posy(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_yr108 = (-rho_r8/(2*pi*e_0)).*((yf - posy(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_yr109 = (-rho_r9/(2*pi*e_0)).*((yf - posy(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_yr1011 = (-rho_r11/(2*pi*e_0)).*((yf - posy(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_yr1012 = (-rho_r12/(2*pi*e_0)).*((yf - posy(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_yr1013 = (rho_r1/(2*pi*e_0)).*((yf - posy(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_yr1014 = (rho_r2/(2*pi*e_0)).*((yf - posy(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_yr1015 = (rho_r3/(2*pi*e_0)).*((yf - posy(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_yr1016 = (rho_r4/(2*pi*e_0)).*((yf - posy(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_yr1017 = (rho_r5/(2*pi*e_0)).*((yf - posy(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_yr1018 = (rho_r6/(2*pi*e_0)).*((yf - posy(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));
E_yr10_19 = (rho_r7/(2*pi*e_0)).*((yf - posy(10,19))./((xf - posx(10,19)).^2 + (yf - posy(10,19)).^2));
E_yr10_20 = (rho_r8/(2*pi*e_0)).*((yf - posy(10,20))./((xf - posx(10,20)).^2 + (yf - posy(10,20)).^2));
E_yr10_21 = (rho_r9/(2*pi*e_0)).*((yf - posy(10,21))./((xf - posx(10,21)).^2 + (yf - posy(10,21)).^2));
E_yr10_22 = (rho_r10/(2*pi*e_0)).*((yf - posy(10,22))./((xf - posx(10,22)).^2 + (yf - posy(10,22)).^2));
E_yr10_23 = (rho_r11/(2*pi*e_0)).*((yf - posy(10,23))./((xf - posx(10,23)).^2 + (yf - posy(10,23)).^2));
E_yr10_24 = (rho_r12/(2*pi*e_0)).*((yf - posy(10,24))./((xf - posx(10,24)).^2 + (yf - posy(10,24)).^2));

E_yr11_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_yr11_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_yr113 = (-rho_r3/(2*pi*e_0)).*((yf - posy(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_yr114 = (-rho_r4/(2*pi*e_0)).*((yf - posy(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_yr115 = (-rho_r5/(2*pi*e_0)).*((yf - posy(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_yr116 = (-rho_r6/(2*pi*e_0)).*((yf - posy(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_yr117 = (-rho_r7/(2*pi*e_0)).*((yf - posy(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_yr118 = (-rho_r8/(2*pi*e_0)).*((yf - posy(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_yr119 = (-rho_r9/(2*pi*e_0)).*((yf - posy(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_yr1110 = (-rho_r10/(2*pi*e_0)).*((yf - posy(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_yr1112 = (-rho_r12/(2*pi*e_0)).*((yf - posy(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_yr1113 = (rho_r1/(2*pi*e_0)).*((yf - posy(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_yr1114 = (rho_r2/(2*pi*e_0)).*((yf - posy(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_yr1115 = (rho_r3/(2*pi*e_0)).*((yf - posy(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_yr1116 = (rho_r4/(2*pi*e_0)).*((yf - posy(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_yr1117 = (rho_r5/(2*pi*e_0)).*((yf - posy(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_yr1118 = (rho_r6/(2*pi*e_0)).*((yf - posy(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));
E_yr1019 = (rho_r7/(2*pi*e_0)).*((yf - posy(11,19))./((xf - posx(11,19)).^2 + (yf - posy(11,19)).^2));
E_yr1020 = (rho_r8/(2*pi*e_0)).*((yf - posy(11,20))./((xf - posx(11,20)).^2 + (yf - posy(11,20)).^2));
E_yr1021 = (rho_r9/(2*pi*e_0)).*((yf - posy(11,21))./((xf - posx(11,21)).^2 + (yf - posy(11,21)).^2));
E_yr1022 = (rho_r10/(2*pi*e_0)).*((yf - posy(11,22))./((xf - posx(11,22)).^2 + (yf - posy(11,22)).^2));
E_yr1023 = (rho_r11/(2*pi*e_0)).*((yf - posy(11,23))./((xf - posx(11,23)).^2 + (yf - posy(11,23)).^2));
E_yr1024 = (rho_r12/(2*pi*e_0)).*((yf - posy(11,24))./((xf - posx(11,24)).^2 + (yf - posy(11,24)).^2));

E_yr12_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_yr12_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_yr123 = (-rho_r3/(2*pi*e_0)).*((yf - posy(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_yr124 = (-rho_r4/(2*pi*e_0)).*((yf - posy(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_yr125 = (-rho_r5/(2*pi*e_0)).*((yf - posy(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_yr126 = (-rho_r6/(2*pi*e_0)).*((yf - posy(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_yr127 = (-rho_r7/(2*pi*e_0)).*((yf - posy(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_yr128 = (-rho_r8/(2*pi*e_0)).*((yf - posy(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_yr129 = (-rho_r9/(2*pi*e_0)).*((yf - posy(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_yr1210 = (-rho_r10/(2*pi*e_0)).*((yf - posy(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_yr1211 = (-rho_r11/(2*pi*e_0)).*((yf - posy(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_yr1213 = (rho_r1/(2*pi*e_0)).*((yf - posy(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_yr1214 = (rho_r2/(2*pi*e_0)).*((yf - posy(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_yr1215 = (rho_r3/(2*pi*e_0)).*((yf - posy(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_yr1216 = (rho_r4/(2*pi*e_0)).*((yf - posy(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_yr1217 = (rho_r5/(2*pi*e_0)).*((yf - posy(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_yr1218 = (rho_r6/(2*pi*e_0)).*((yf - posy(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));
E_yr1219 = (rho_r7/(2*pi*e_0)).*((yf - posy(12,19))./((xf - posx(12,19)).^2 + (yf - posy(12,19)).^2));
E_yr1220 = (rho_r8/(2*pi*e_0)).*((yf - posy(12,20))./((xf - posx(12,20)).^2 + (yf - posy(12,20)).^2));
E_yr1221 = (rho_r9/(2*pi*e_0)).*((yf - posy(12,21))./((xf - posx(12,21)).^2 + (yf - posy(12,21)).^2));
E_yr1222 = (rho_r10/(2*pi*e_0)).*((yf - posy(12,22))./((xf - posx(12,22)).^2 + (yf - posy(12,22)).^2));
E_yr1223 = (rho_r11/(2*pi*e_0)).*((yf - posy(12,23))./((xf - posx(12,23)).^2 + (yf - posy(12,23)).^2));
E_yr1224 = (rho_r12/(2*pi*e_0)).*((yf - posy(12,24))./((xf - posx(12,24)).^2 + (yf - posy(12,24)).^2));

E_yr131 = (-rho_r1/(2*pi*e_0)).*((yf - posy(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_yr132 = (-rho_r2/(2*pi*e_0)).*((yf - posy(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_yr133 = (-rho_r3/(2*pi*e_0)).*((yf - posy(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_yr134 = (-rho_r4/(2*pi*e_0)).*((yf - posy(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_yr135 = (-rho_r5/(2*pi*e_0)).*((yf - posy(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_yr136 = (-rho_r6/(2*pi*e_0)).*((yf - posy(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_yr137 = (-rho_r7/(2*pi*e_0)).*((yf - posy(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_yr138 = (-rho_r8/(2*pi*e_0)).*((yf - posy(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_yr139 = (-rho_r9/(2*pi*e_0)).*((yf - posy(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_yr1310 = (-rho_r10/(2*pi*e_0)).*((yf - posy(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_yr1311 = (-rho_r11/(2*pi*e_0)).*((yf - posy(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_yr1312 = (-rho_r12/(2*pi*e_0)).*((yf - posy(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_yr1314 = (rho_r2/(2*pi*e_0)).*((yf - posy(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_yr1315 = (rho_r3/(2*pi*e_0)).*((yf - posy(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_yr1316 = (rho_r4/(2*pi*e_0)).*((yf - posy(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_yr1317 = (rho_r5/(2*pi*e_0)).*((yf - posy(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_yr1318 = (rho_r6/(2*pi*e_0)).*((yf - posy(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));
E_yr1319 = (rho_r7/(2*pi*e_0)).*((yf - posy(13,19))./((xf - posx(13,19)).^2 + (yf - posy(13,19)).^2));
E_yr1320 = (rho_r8/(2*pi*e_0)).*((yf - posy(13,20))./((xf - posx(13,20)).^2 + (yf - posy(13,20)).^2));
E_yr1321 = (rho_r9/(2*pi*e_0)).*((yf - posy(13,21))./((xf - posx(13,21)).^2 + (yf - posy(13,21)).^2));
E_yr1322 = (rho_r10/(2*pi*e_0)).*((yf - posy(13,22))./((xf - posx(13,22)).^2 + (yf - posy(13,22)).^2));
E_yr1323 = (rho_r11/(2*pi*e_0)).*((yf - posy(13,23))./((xf - posx(13,23)).^2 + (yf - posy(13,23)).^2));
E_yr1324 = (rho_r12/(2*pi*e_0)).*((yf - posy(13,24))./((xf - posx(13,24)).^2 + (yf - posy(13,24)).^2));

E_yr141 = (-rho_r1/(2*pi*e_0)).*((yf - posy(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_yr142 = (-rho_r2/(2*pi*e_0)).*((yf - posy(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_yr143 = (-rho_r3/(2*pi*e_0)).*((yf - posy(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_yr144 = (-rho_r4/(2*pi*e_0)).*((yf - posy(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_yr145 = (-rho_r5/(2*pi*e_0)).*((yf - posy(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_yr146 = (-rho_r6/(2*pi*e_0)).*((yf - posy(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_yr147 = (-rho_r7/(2*pi*e_0)).*((yf - posy(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_yr148 = (-rho_r8/(2*pi*e_0)).*((yf - posy(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_yr149 = (-rho_r9/(2*pi*e_0)).*((yf - posy(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_yr1410 = (-rho_r10/(2*pi*e_0)).*((yf - posy(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_yr1411 = (-rho_r11/(2*pi*e_0)).*((yf - posy(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_yr1412 = (-rho_r12/(2*pi*e_0)).*((yf - posy(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_yr1413 = (rho_r1/(2*pi*e_0)).*((yf - posy(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_yr1415 = (rho_r3/(2*pi*e_0)).*((yf - posy(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_yr1416 = (rho_r4/(2*pi*e_0)).*((yf - posy(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_yr1417 = (rho_r5/(2*pi*e_0)).*((yf - posy(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_yr1418 = (rho_r6/(2*pi*e_0)).*((yf - posy(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));
E_yr1419 = (rho_r7/(2*pi*e_0)).*((yf - posy(14,19))./((xf - posx(14,19)).^2 + (yf - posy(14,19)).^2));
E_yr1420 = (rho_r8/(2*pi*e_0)).*((yf - posy(14,20))./((xf - posx(14,20)).^2 + (yf - posy(14,20)).^2));
E_yr1421 = (rho_r9/(2*pi*e_0)).*((yf - posy(14,21))./((xf - posx(14,21)).^2 + (yf - posy(14,21)).^2));
E_yr1422 = (rho_r10/(2*pi*e_0)).*((yf - posy(14,22))./((xf - posx(14,22)).^2 + (yf - posy(14,22)).^2));
E_yr1423 = (rho_r11/(2*pi*e_0)).*((yf - posy(14,23))./((xf - posx(14,23)).^2 + (yf - posy(14,23)).^2));
E_yr1424 = (rho_r12/(2*pi*e_0)).*((yf - posy(14,24))./((xf - posx(14,24)).^2 + (yf - posy(14,24)).^2));

E_yr151 = (-rho_r1/(2*pi*e_0)).*((yf - posy(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_yr152 = (-rho_r2/(2*pi*e_0)).*((yf - posy(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_yr153 = (-rho_r3/(2*pi*e_0)).*((yf - posy(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_yr154 = (-rho_r4/(2*pi*e_0)).*((yf - posy(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_yr155 = (-rho_r5/(2*pi*e_0)).*((yf - posy(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_yr156 = (-rho_r6/(2*pi*e_0)).*((yf - posy(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_yr157 = (-rho_r7/(2*pi*e_0)).*((yf - posy(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_yr158 = (-rho_r8/(2*pi*e_0)).*((yf - posy(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_yr159 = (-rho_r9/(2*pi*e_0)).*((yf - posy(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_yr1510 = (-rho_r10/(2*pi*e_0)).*((yf - posy(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_yr1511 = (-rho_r11/(2*pi*e_0)).*((yf - posy(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_yr1512 = (-rho_r12/(2*pi*e_0)).*((yf - posy(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_yr1513 = (rho_r1/(2*pi*e_0)).*((yf - posy(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_yr1514 = (rho_r2/(2*pi*e_0)).*((yf - posy(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_yr1516 = (rho_r4/(2*pi*e_0)).*((yf - posy(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_yr1517 = (rho_r5/(2*pi*e_0)).*((yf - posy(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_yr1518 = (rho_r6/(2*pi*e_0)).*((yf - posy(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));
E_yr1519 = (rho_r7/(2*pi*e_0)).*((yf - posy(15,19))./((xf - posx(15,19)).^2 + (yf - posy(15,19)).^2));
E_yr1520 = (rho_r8/(2*pi*e_0)).*((yf - posy(15,20))./((xf - posx(15,20)).^2 + (yf - posy(15,20)).^2));
E_yr1521 = (rho_r9/(2*pi*e_0)).*((yf - posy(15,21))./((xf - posx(15,21)).^2 + (yf - posy(15,21)).^2));
E_yr1522 = (rho_r10/(2*pi*e_0)).*((yf - posy(15,22))./((xf - posx(15,22)).^2 + (yf - posy(15,22)).^2));
E_yr1523 = (rho_r11/(2*pi*e_0)).*((yf - posy(15,23))./((xf - posx(15,23)).^2 + (yf - posy(15,23)).^2));
E_yr1524 = (rho_r12/(2*pi*e_0)).*((yf - posy(15,24))./((xf - posx(15,24)).^2 + (yf - posy(15,24)).^2));

E_yr161 = (-rho_r1/(2*pi*e_0)).*((yf - posy(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_yr162 = (-rho_r2/(2*pi*e_0)).*((yf - posy(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_yr163 = (-rho_r3/(2*pi*e_0)).*((yf - posy(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_yr164 = (-rho_r4/(2*pi*e_0)).*((yf - posy(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_yr165 = (-rho_r5/(2*pi*e_0)).*((yf - posy(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_yr166 = (-rho_r6/(2*pi*e_0)).*((yf - posy(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_yr167 = (-rho_r7/(2*pi*e_0)).*((yf - posy(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_yr168 = (-rho_r8/(2*pi*e_0)).*((yf - posy(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_yr169 = (-rho_r9/(2*pi*e_0)).*((yf - posy(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_yr1610 = (-rho_r10/(2*pi*e_0)).*((yf - posy(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_yr1611 = (-rho_r11/(2*pi*e_0)).*((yf - posy(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_yr1612 = (-rho_r12/(2*pi*e_0)).*((yf - posy(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_yr1613 = (rho_r1/(2*pi*e_0)).*((yf - posy(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_yr1614 = (rho_r2/(2*pi*e_0)).*((yf - posy(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_yr1615 = (rho_r3/(2*pi*e_0)).*((yf - posy(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_yr1617 = (rho_r5/(2*pi*e_0)).*((yf - posy(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_yr1618 = (rho_r6/(2*pi*e_0)).*((yf - posy(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));
E_yr1619 = (rho_r7/(2*pi*e_0)).*((yf - posy(16,19))./((xf - posx(16,19)).^2 + (yf - posy(16,19)).^2));
E_yr1620 = (rho_r8/(2*pi*e_0)).*((yf - posy(16,20))./((xf - posx(16,20)).^2 + (yf - posy(16,20)).^2));
E_yr1621 = (rho_r9/(2*pi*e_0)).*((yf - posy(16,21))./((xf - posx(16,21)).^2 + (yf - posy(16,21)).^2));
E_yr1622 = (rho_r10/(2*pi*e_0)).*((yf - posy(16,22))./((xf - posx(16,22)).^2 + (yf - posy(16,22)).^2));
E_yr1623 = (rho_r11/(2*pi*e_0)).*((yf - posy(16,23))./((xf - posx(16,23)).^2 + (yf - posy(16,23)).^2));
E_yr1624 = (rho_r12/(2*pi*e_0)).*((yf - posy(16,24))./((xf - posx(16,24)).^2 + (yf - posy(16,24)).^2));

E_yr171 = (-rho_r1/(2*pi*e_0)).*((yf - posy(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_yr172 = (-rho_r2/(2*pi*e_0)).*((yf - posy(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_yr173 = (-rho_r3/(2*pi*e_0)).*((yf - posy(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_yr174 = (-rho_r4/(2*pi*e_0)).*((yf - posy(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_yr175 = (-rho_r5/(2*pi*e_0)).*((yf - posy(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_yr176 = (-rho_r6/(2*pi*e_0)).*((yf - posy(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_yr177 = (-rho_r7/(2*pi*e_0)).*((yf - posy(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_yr178 = (-rho_r8/(2*pi*e_0)).*((yf - posy(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_yr179 = (-rho_r9/(2*pi*e_0)).*((yf - posy(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_yr1710 = (-rho_r10/(2*pi*e_0)).*((yf - posy(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_yr1711 = (-rho_r11/(2*pi*e_0)).*((yf - posy(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_yr1712 = (-rho_r12/(2*pi*e_0)).*((yf - posy(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_yr1713 = (rho_r1/(2*pi*e_0)).*((yf - posy(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_yr1714 = (rho_r2/(2*pi*e_0)).*((yf - posy(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_yr1715 = (rho_r3/(2*pi*e_0)).*((yf - posy(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_yr1716 = (rho_r4/(2*pi*e_0)).*((yf - posy(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_yr1718 = (rho_r6/(2*pi*e_0)).*((yf - posy(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));
E_yr1719 = (rho_r7/(2*pi*e_0)).*((yf - posy(17,19))./((xf - posx(17,19)).^2 + (yf - posy(17,19)).^2));
E_yr1720 = (rho_r8/(2*pi*e_0)).*((yf - posy(17,20))./((xf - posx(17,20)).^2 + (yf - posy(17,20)).^2));
E_yr1721 = (rho_r9/(2*pi*e_0)).*((yf - posy(17,21))./((xf - posx(17,21)).^2 + (yf - posy(17,21)).^2));
E_yr1722 = (rho_r10/(2*pi*e_0)).*((yf - posy(17,22))./((xf - posx(17,22)).^2 + (yf - posy(17,22)).^2));
E_yr1723 = (rho_r11/(2*pi*e_0)).*((yf - posy(17,23))./((xf - posx(17,23)).^2 + (yf - posy(17,23)).^2));
E_yr1724 = (rho_r12/(2*pi*e_0)).*((yf - posy(17,24))./((xf - posx(17,24)).^2 + (yf - posy(17,24)).^2));

E_yr181 = (-rho_r1/(2*pi*e_0)).*((yf - posy(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_yr182 = (-rho_r2/(2*pi*e_0)).*((yf - posy(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_yr183 = (-rho_r3/(2*pi*e_0)).*((yf - posy(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_yr184 = (-rho_r4/(2*pi*e_0)).*((yf - posy(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_yr185 = (-rho_r5/(2*pi*e_0)).*((yf - posy(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_yr186 = (-rho_r6/(2*pi*e_0)).*((yf - posy(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_yr187 = (-rho_r7/(2*pi*e_0)).*((yf - posy(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_yr188 = (-rho_r8/(2*pi*e_0)).*((yf - posy(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_yr189 = (-rho_r9/(2*pi*e_0)).*((yf - posy(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_yr1810 = (-rho_r10/(2*pi*e_0)).*((yf - posy(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_yr1811 = (-rho_r11/(2*pi*e_0)).*((yf - posy(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_yr1812 = (-rho_r12/(2*pi*e_0)).*((yf - posy(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_yr1813 = (rho_r1/(2*pi*e_0)).*((yf - posy(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_yr1814 = (rho_r2/(2*pi*e_0)).*((yf - posy(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_yr1815 = (rho_r3/(2*pi*e_0)).*((yf - posy(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_yr1816 = (rho_r4/(2*pi*e_0)).*((yf - posy(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_yr1817 = (rho_r5/(2*pi*e_0)).*((yf - posy(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));
E_yr1819 = (rho_r7/(2*pi*e_0)).*((yf - posy(18,19))./((xf - posx(18,19)).^2 + (yf - posy(18,19)).^2));
E_yr1820 = (rho_r8/(2*pi*e_0)).*((yf - posy(18,20))./((xf - posx(18,20)).^2 + (yf - posy(18,20)).^2));
E_yr1821 = (rho_r9/(2*pi*e_0)).*((yf - posy(18,21))./((xf - posx(18,21)).^2 + (yf - posy(18,21)).^2));
E_yr1822 = (rho_r10/(2*pi*e_0)).*((yf - posy(18,22))./((xf - posx(18,22)).^2 + (yf - posy(18,22)).^2));
E_yr1823 = (rho_r11/(2*pi*e_0)).*((yf - posy(18,23))./((xf - posx(18,23)).^2 + (yf - posy(18,23)).^2));
E_yr1824 = (rho_r12/(2*pi*e_0)).*((yf - posy(18,24))./((xf - posx(18,24)).^2 + (yf - posy(18,24)).^2));

E_yr191 = (-rho_r1/(2*pi*e_0)).*((yf - posy(19,1))./((xf - posx(19,1)).^2 + (yf - posy(19,1)).^2));
E_yr192 = (-rho_r2/(2*pi*e_0)).*((yf - posy(19,2))./((xf - posx(19,2)).^2 + (yf - posy(19,2)).^2));
E_yr193 = (-rho_r3/(2*pi*e_0)).*((yf - posy(19,3))./((xf - posx(19,3)).^2 + (yf - posy(19,3)).^2));
E_yr194 = (-rho_r4/(2*pi*e_0)).*((yf - posy(19,4))./((xf - posx(19,4)).^2 + (yf - posy(19,4)).^2));
E_yr195 = (-rho_r5/(2*pi*e_0)).*((yf - posy(19,5))./((xf - posx(19,5)).^2 + (yf - posy(19,5)).^2));
E_yr196 = (-rho_r6/(2*pi*e_0)).*((yf - posy(19,6))./((xf - posx(19,6)).^2 + (yf - posy(19,6)).^2));
E_yr197 = (-rho_r7/(2*pi*e_0)).*((yf - posy(19,7))./((xf - posx(19,7)).^2 + (yf - posy(19,7)).^2));
E_yr198 = (-rho_r8/(2*pi*e_0)).*((yf - posy(19,8))./((xf - posx(19,8)).^2 + (yf - posy(19,8)).^2));
E_yr199 = (-rho_r9/(2*pi*e_0)).*((yf - posy(19,9))./((xf - posx(19,9)).^2 + (yf - posy(19,9)).^2));
E_yr1910 = (-rho_r10/(2*pi*e_0)).*((yf - posy(19,10))./((xf - posx(19,10)).^2 + (yf - posy(19,10)).^2));
E_yr1911 = (-rho_r11/(2*pi*e_0)).*((yf - posy(19,11))./((xf - posx(19,11)).^2 + (yf - posy(19,11)).^2));
E_yr1912 = (-rho_r12/(2*pi*e_0)).*((yf - posy(19,12))./((xf - posx(19,12)).^2 + (yf - posy(19,12)).^2));
E_yr1913 = (rho_r1/(2*pi*e_0)).*((yf - posy(19,13))./((xf - posx(19,13)).^2 + (yf - posy(19,13)).^2));
E_yr1914 = (rho_r2/(2*pi*e_0)).*((yf - posy(19,14))./((xf - posx(19,14)).^2 + (yf - posy(19,14)).^2));
E_yr1915 = (rho_r3/(2*pi*e_0)).*((yf - posy(19,15))./((xf - posx(19,15)).^2 + (yf - posy(19,15)).^2));
E_yr1916 = (rho_r4/(2*pi*e_0)).*((yf - posy(19,16))./((xf - posx(19,16)).^2 + (yf - posy(19,16)).^2));
E_yr1917 = (rho_r5/(2*pi*e_0)).*((yf - posy(19,17))./((xf - posx(19,17)).^2 + (yf - posy(19,17)).^2));
E_yr1918 = (rho_r6/(2*pi*e_0)).*((yf - posy(19,18))./((xf - posx(19,18)).^2 + (yf - posy(19,18)).^2));
E_yr1920 = (rho_r8/(2*pi*e_0)).*((yf - posy(19,20))./((xf - posx(19,20)).^2 + (yf - posy(19,20)).^2));
E_yr1921 = (rho_r9/(2*pi*e_0)).*((yf - posy(19,21))./((xf - posx(19,21)).^2 + (yf - posy(19,21)).^2));
E_yr1922 = (rho_r10/(2*pi*e_0)).*((yf - posy(19,22))./((xf - posx(19,22)).^2 + (yf - posy(19,22)).^2));
E_yr1923 = (rho_r11/(2*pi*e_0)).*((yf - posy(19,23))./((xf - posx(19,23)).^2 + (yf - posy(19,23)).^2));
E_yr1924 = (rho_r12/(2*pi*e_0)).*((yf - posy(19,24))./((xf - posx(19,24)).^2 + (yf - posy(19,24)).^2));

E_yr201 = (-rho_r1/(2*pi*e_0)).*((yf - posy(20,1))./((xf - posx(20,1)).^2 + (yf - posy(20,1)).^2));
E_yr202 = (-rho_r2/(2*pi*e_0)).*((yf - posy(20,2))./((xf - posx(20,2)).^2 + (yf - posy(20,2)).^2));
E_yr203 = (-rho_r3/(2*pi*e_0)).*((yf - posy(20,3))./((xf - posx(20,3)).^2 + (yf - posy(20,3)).^2));
E_yr204 = (-rho_r4/(2*pi*e_0)).*((yf - posy(20,4))./((xf - posx(20,4)).^2 + (yf - posy(20,4)).^2));
E_yr205 = (-rho_r5/(2*pi*e_0)).*((yf - posy(20,5))./((xf - posx(20,5)).^2 + (yf - posy(20,5)).^2));
E_yr206 = (-rho_r6/(2*pi*e_0)).*((yf - posy(20,6))./((xf - posx(20,6)).^2 + (yf - posy(20,6)).^2));
E_yr207 = (-rho_r7/(2*pi*e_0)).*((yf - posy(20,7))./((xf - posx(20,7)).^2 + (yf - posy(20,7)).^2));
E_yr208 = (-rho_r8/(2*pi*e_0)).*((yf - posy(20,8))./((xf - posx(20,8)).^2 + (yf - posy(20,8)).^2));
E_yr209 = (-rho_r9/(2*pi*e_0)).*((yf - posy(20,9))./((xf - posx(20,9)).^2 + (yf - posy(20,9)).^2));
E_yr2010 = (-rho_r10/(2*pi*e_0)).*((yf - posy(20,10))./((xf - posx(20,10)).^2 + (yf - posy(20,10)).^2));
E_yr2011 = (-rho_r11/(2*pi*e_0)).*((yf - posy(20,11))./((xf - posx(20,11)).^2 + (yf - posy(20,11)).^2));
E_yr2012 = (-rho_r12/(2*pi*e_0)).*((yf - posy(20,12))./((xf - posx(20,12)).^2 + (yf - posy(20,12)).^2));
E_yr2013 = (rho_r1/(2*pi*e_0)).*((yf - posy(20,13))./((xf - posx(20,13)).^2 + (yf - posy(20,13)).^2));
E_yr2014 = (rho_r2/(2*pi*e_0)).*((yf - posy(20,14))./((xf - posx(20,14)).^2 + (yf - posy(20,14)).^2));
E_yr2015 = (rho_r3/(2*pi*e_0)).*((yf - posy(20,15))./((xf - posx(20,15)).^2 + (yf - posy(20,15)).^2));
E_yr2016 = (rho_r4/(2*pi*e_0)).*((yf - posy(20,16))./((xf - posx(20,16)).^2 + (yf - posy(20,16)).^2));
E_yr2017 = (rho_r5/(2*pi*e_0)).*((yf - posy(20,17))./((xf - posx(20,17)).^2 + (yf - posy(20,17)).^2));
E_yr2018 = (rho_r6/(2*pi*e_0)).*((yf - posy(20,18))./((xf - posx(20,18)).^2 + (yf - posy(20,18)).^2));
E_yr2019 = (rho_r7/(2*pi*e_0)).*((yf - posy(20,19))./((xf - posx(20,19)).^2 + (yf - posy(20,19)).^2));
E_yr2021 = (rho_r9/(2*pi*e_0)).*((yf - posy(20,21))./((xf - posx(20,21)).^2 + (yf - posy(20,21)).^2));
E_yr2022 = (rho_r10/(2*pi*e_0)).*((yf - posy(20,22))./((xf - posx(20,22)).^2 + (yf - posy(20,22)).^2));
E_yr2023 = (rho_r11/(2*pi*e_0)).*((yf - posy(20,23))./((xf - posx(20,23)).^2 + (yf - posy(20,23)).^2));
E_yr2024 = (rho_r12/(2*pi*e_0)).*((yf - posy(20,24))./((xf - posx(20,24)).^2 + (yf - posy(20,24)).^2));

E_yr21_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(21,1))./((xf - posx(21,1)).^2 + (yf - posy(21,1)).^2));
E_yr21_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(21,2))./((xf - posx(21,2)).^2 + (yf - posy(21,2)).^2));
E_yr21_3 = (-rho_r3/(2*pi*e_0)).*((yf - posy(21,3))./((xf - posx(21,3)).^2 + (yf - posy(21,3)).^2));
E_yr214 = (-rho_r4/(2*pi*e_0)).*((yf - posy(21,4))./((xf - posx(21,4)).^2 + (yf - posy(21,4)).^2));
E_yr215 = (-rho_r5/(2*pi*e_0)).*((yf - posy(21,5))./((xf - posx(21,5)).^2 + (yf - posy(21,5)).^2));
E_yr216 = (-rho_r6/(2*pi*e_0)).*((yf - posy(21,6))./((xf - posx(21,6)).^2 + (yf - posy(21,6)).^2));
E_yr217 = (-rho_r7/(2*pi*e_0)).*((yf - posy(21,7))./((xf - posx(21,7)).^2 + (yf - posy(21,7)).^2));
E_yr218 = (-rho_r8/(2*pi*e_0)).*((yf - posy(21,8))./((xf - posx(21,8)).^2 + (yf - posy(21,8)).^2));
E_yr219 = (-rho_r9/(2*pi*e_0)).*((yf - posy(21,9))./((xf - posx(21,9)).^2 + (yf - posy(21,9)).^2));
E_yr2110 = (-rho_r10/(2*pi*e_0)).*((yf - posy(21,10))./((xf - posx(21,10)).^2 + (yf - posy(21,10)).^2));
E_yr2111 = (-rho_r11/(2*pi*e_0)).*((yf - posy(21,11))./((xf - posx(21,11)).^2 + (yf - posy(21,11)).^2));
E_yr2112 = (-rho_r12/(2*pi*e_0)).*((yf - posy(21,12))./((xf - posx(21,12)).^2 + (yf - posy(21,12)).^2));
E_yr2113 = (rho_r1/(2*pi*e_0)).*((yf - posy(21,13))./((xf - posx(21,13)).^2 + (yf - posy(21,13)).^2));
E_yr2114 = (rho_r2/(2*pi*e_0)).*((yf - posy(21,14))./((xf - posx(21,14)).^2 + (yf - posy(21,14)).^2));
E_yr2115 = (rho_r3/(2*pi*e_0)).*((yf - posy(21,15))./((xf - posx(21,15)).^2 + (yf - posy(21,15)).^2));
E_yr2116 = (rho_r4/(2*pi*e_0)).*((yf - posy(21,16))./((xf - posx(21,16)).^2 + (yf - posy(21,16)).^2));
E_yr2117 = (rho_r5/(2*pi*e_0)).*((yf - posy(21,17))./((xf - posx(21,17)).^2 + (yf - posy(21,17)).^2));
E_yr2118 = (rho_r6/(2*pi*e_0)).*((yf - posy(21,18))./((xf - posx(21,18)).^2 + (yf - posy(21,18)).^2));
E_yr2119 = (rho_r7/(2*pi*e_0)).*((yf - posy(21,19))./((xf - posx(21,19)).^2 + (yf - posy(21,19)).^2));
E_yr2120 = (rho_r8/(2*pi*e_0)).*((yf - posy(21,20))./((xf - posx(21,20)).^2 + (yf - posy(21,20)).^2));
E_yr2122 = (rho_r10/(2*pi*e_0)).*((yf - posy(21,22))./((xf - posx(21,22)).^2 + (yf - posy(21,22)).^2));
E_yr2123 = (rho_r11/(2*pi*e_0)).*((yf - posy(21,23))./((xf - posx(21,23)).^2 + (yf - posy(21,23)).^2));
E_yr2124 = (rho_r12/(2*pi*e_0)).*((yf - posy(21,24))./((xf - posx(21,24)).^2 + (yf - posy(21,24)).^2));

E_yr22_1 = (-rho_r1/(2*pi*e_0)).*((yf - posy(22,1))./((xf - posx(22,1)).^2 + (yf - posy(22,1)).^2));
E_yr22_2 = (-rho_r2/(2*pi*e_0)).*((yf - posy(22,2))./((xf - posx(22,2)).^2 + (yf - posy(22,2)).^2));
E_yr22_3 = (-rho_r3/(2*pi*e_0)).*((yf - posy(22,3))./((xf - posx(22,3)).^2 + (yf - posy(22,3)).^2));
E_yr22_4 = (-rho_r4/(2*pi*e_0)).*((yf - posy(22,4))./((xf - posx(22,4)).^2 + (yf - posy(22,4)).^2));
E_yr225 = (-rho_r5/(2*pi*e_0)).*((yf - posy(22,5))./((xf - posx(22,5)).^2 + (yf - posy(22,5)).^2));
E_yr226 = (-rho_r6/(2*pi*e_0)).*((yf - posy(22,6))./((xf - posx(22,6)).^2 + (yf - posy(22,6)).^2));
E_yr227 = (-rho_r7/(2*pi*e_0)).*((yf - posy(22,7))./((xf - posx(22,7)).^2 + (yf - posy(22,7)).^2));
E_yr228 = (-rho_r8/(2*pi*e_0)).*((yf - posy(22,8))./((xf - posx(22,8)).^2 + (yf - posy(22,8)).^2));
E_yr229 = (-rho_r9/(2*pi*e_0)).*((yf - posy(22,9))./((xf - posx(22,9)).^2 + (yf - posy(22,9)).^2));
E_yr2210 = (-rho_r10/(2*pi*e_0)).*((yf - posy(22,10))./((xf - posx(22,10)).^2 + (yf - posy(22,10)).^2));
E_yr2211 = (-rho_r11/(2*pi*e_0)).*((yf - posy(22,11))./((xf - posx(22,11)).^2 + (yf - posy(22,11)).^2));
E_yr2212 = (-rho_r12/(2*pi*e_0)).*((yf - posy(22,12))./((xf - posx(22,12)).^2 + (yf - posy(22,12)).^2));
E_yr2213 = (rho_r1/(2*pi*e_0)).*((yf - posy(22,13))./((xf - posx(22,13)).^2 + (yf - posy(22,13)).^2));
E_yr2214 = (rho_r2/(2*pi*e_0)).*((yf - posy(22,14))./((xf - posx(22,14)).^2 + (yf - posy(22,14)).^2));
E_yr2215 = (rho_r3/(2*pi*e_0)).*((yf - posy(22,15))./((xf - posx(22,15)).^2 + (yf - posy(22,15)).^2));
E_yr2216 = (rho_r4/(2*pi*e_0)).*((yf - posy(22,16))./((xf - posx(22,16)).^2 + (yf - posy(22,16)).^2));
E_yr2217 = (rho_r5/(2*pi*e_0)).*((yf - posy(22,17))./((xf - posx(22,17)).^2 + (yf - posy(22,17)).^2));
E_yr2218 = (rho_r6/(2*pi*e_0)).*((yf - posy(22,18))./((xf - posx(22,18)).^2 + (yf - posy(22,18)).^2));
E_yr2219 = (rho_r7/(2*pi*e_0)).*((yf - posy(22,19))./((xf - posx(22,19)).^2 + (yf - posy(22,19)).^2));
E_yr2220 = (rho_r8/(2*pi*e_0)).*((yf - posy(22,20))./((xf - posx(22,20)).^2 + (yf - posy(22,20)).^2));
E_yr2221 = (rho_r9/(2*pi*e_0)).*((yf - posy(22,21))./((xf - posx(22,21)).^2 + (yf - posy(22,21)).^2));
E_yr2223 = (rho_r11/(2*pi*e_0)).*((yf - posy(22,23))./((xf - posx(22,23)).^2 + (yf - posy(22,23)).^2));
E_yr2224 = (rho_r12/(2*pi*e_0)).*((yf - posy(22,24))./((xf - posx(22,24)).^2 + (yf - posy(22,24)).^2));

E_yr231 = (-rho_r1/(2*pi*e_0)).*((yf - posy(23,1))./((xf - posx(23,1)).^2 + (yf - posy(23,1)).^2));
E_yr232 = (-rho_r2/(2*pi*e_0)).*((yf - posy(23,2))./((xf - posx(23,2)).^2 + (yf - posy(23,2)).^2));
E_yr233 = (-rho_r3/(2*pi*e_0)).*((yf - posy(23,3))./((xf - posx(23,3)).^2 + (yf - posy(23,3)).^2));
E_yr234 = (-rho_r4/(2*pi*e_0)).*((yf - posy(23,4))./((xf - posx(23,4)).^2 + (yf - posy(23,4)).^2));
E_yr235 = (-rho_r5/(2*pi*e_0)).*((yf - posy(23,5))./((xf - posx(23,5)).^2 + (yf - posy(23,5)).^2));
E_yr236 = (-rho_r6/(2*pi*e_0)).*((yf - posy(23,6))./((xf - posx(23,6)).^2 + (yf - posy(23,6)).^2));
E_yr237 = (-rho_r7/(2*pi*e_0)).*((yf - posy(23,7))./((xf - posx(23,7)).^2 + (yf - posy(23,7)).^2));
E_yr238 = (-rho_r8/(2*pi*e_0)).*((yf - posy(23,8))./((xf - posx(23,8)).^2 + (yf - posy(23,8)).^2));
E_yr239 = (-rho_r9/(2*pi*e_0)).*((yf - posy(23,9))./((xf - posx(23,9)).^2 + (yf - posy(23,9)).^2));
E_yr2310 = (-rho_r10/(2*pi*e_0)).*((yf - posy(23,10))./((xf - posx(23,10)).^2 + (yf - posy(23,10)).^2));
E_yr2311 = (-rho_r11/(2*pi*e_0)).*((yf - posy(23,11))./((xf - posx(23,11)).^2 + (yf - posy(23,11)).^2));
E_yr2312 = (-rho_r12/(2*pi*e_0)).*((yf - posy(23,12))./((xf - posx(23,12)).^2 + (yf - posy(23,12)).^2));
E_yr2313 = (rho_r1/(2*pi*e_0)).*((yf - posy(23,13))./((xf - posx(23,13)).^2 + (yf - posy(23,13)).^2));
E_yr2314 = (rho_r2/(2*pi*e_0)).*((yf - posy(23,14))./((xf - posx(23,14)).^2 + (yf - posy(23,14)).^2));
E_yr2315 = (rho_r3/(2*pi*e_0)).*((yf - posy(23,15))./((xf - posx(23,15)).^2 + (yf - posy(23,15)).^2));
E_yr2316 = (rho_r4/(2*pi*e_0)).*((yf - posy(23,16))./((xf - posx(23,16)).^2 + (yf - posy(23,16)).^2));
E_yr2317 = (rho_r5/(2*pi*e_0)).*((yf - posy(23,17))./((xf - posx(23,17)).^2 + (yf - posy(23,17)).^2));
E_yr2318 = (rho_r6/(2*pi*e_0)).*((yf - posy(23,18))./((xf - posx(23,18)).^2 + (yf - posy(23,18)).^2));
E_yr2319 = (rho_r7/(2*pi*e_0)).*((yf - posy(23,19))./((xf - posx(23,19)).^2 + (yf - posy(23,19)).^2));
E_yr2320 = (rho_r8/(2*pi*e_0)).*((yf - posy(23,20))./((xf - posx(23,20)).^2 + (yf - posy(23,20)).^2));
E_yr2321 = (rho_r9/(2*pi*e_0)).*((yf - posy(23,21))./((xf - posx(23,21)).^2 + (yf - posy(23,21)).^2));
E_yr2322 = (rho_r10/(2*pi*e_0)).*((yf - posy(23,22))./((xf - posx(23,22)).^2 + (yf - posy(23,22)).^2));
E_yr2324 = (rho_r12/(2*pi*e_0)).*((yf - posy(23,24))./((xf - posx(23,24)).^2 + (yf - posy(23,24)).^2));

E_yr241 = (-rho_r1/(2*pi*e_0)).*((yf - posy(24,1))./((xf - posx(24,1)).^2 + (yf - posy(24,1)).^2));
E_yr242 = (-rho_r2/(2*pi*e_0)).*((yf - posy(24,2))./((xf - posx(24,2)).^2 + (yf - posy(24,2)).^2));
E_yr243 = (-rho_r3/(2*pi*e_0)).*((yf - posy(24,3))./((xf - posx(24,3)).^2 + (yf - posy(24,3)).^2));
E_yr244 = (-rho_r4/(2*pi*e_0)).*((yf - posy(24,4))./((xf - posx(24,4)).^2 + (yf - posy(24,4)).^2));
E_yr245 = (-rho_r5/(2*pi*e_0)).*((yf - posy(24,5))./((xf - posx(24,5)).^2 + (yf - posy(24,5)).^2));
E_yr246 = (-rho_r6/(2*pi*e_0)).*((yf - posy(24,6))./((xf - posx(24,6)).^2 + (yf - posy(24,6)).^2));
E_yr247 = (-rho_r7/(2*pi*e_0)).*((yf - posy(24,7))./((xf - posx(24,7)).^2 + (yf - posy(24,7)).^2));
E_yr248 = (-rho_r8/(2*pi*e_0)).*((yf - posy(24,8))./((xf - posx(24,8)).^2 + (yf - posy(24,8)).^2));
E_yr249 = (-rho_r9/(2*pi*e_0)).*((yf - posy(24,9))./((xf - posx(24,9)).^2 + (yf - posy(24,9)).^2));
E_yr2410 = (-rho_r10/(2*pi*e_0)).*((yf - posy(24,10))./((xf - posx(24,10)).^2 + (yf - posy(24,10)).^2));
E_yr2411 = (-rho_r11/(2*pi*e_0)).*((yf - posy(24,11))./((xf - posx(24,11)).^2 + (yf - posy(24,11)).^2));
E_yr2412 = (-rho_r12/(2*pi*e_0)).*((yf - posy(24,12))./((xf - posx(24,12)).^2 + (yf - posy(24,12)).^2));
E_yr2413 = (rho_r1/(2*pi*e_0)).*((yf - posy(24,13))./((xf - posx(24,13)).^2 + (yf - posy(24,13)).^2));
E_yr2414 = (rho_r2/(2*pi*e_0)).*((yf - posy(24,14))./((xf - posx(24,14)).^2 + (yf - posy(24,14)).^2));
E_yr2415 = (rho_r3/(2*pi*e_0)).*((yf - posy(24,15))./((xf - posx(24,15)).^2 + (yf - posy(24,15)).^2));
E_yr2416 = (rho_r4/(2*pi*e_0)).*((yf - posy(24,16))./((xf - posx(24,16)).^2 + (yf - posy(24,16)).^2));
E_yr2417 = (rho_r5/(2*pi*e_0)).*((yf - posy(24,17))./((xf - posx(24,17)).^2 + (yf - posy(24,17)).^2));
E_yr2418 = (rho_r6/(2*pi*e_0)).*((yf - posy(24,18))./((xf - posx(24,18)).^2 + (yf - posy(24,18)).^2));
E_yr2419 = (rho_r7/(2*pi*e_0)).*((yf - posy(24,19))./((xf - posx(24,19)).^2 + (yf - posy(24,19)).^2));
E_yr2420 = (rho_r8/(2*pi*e_0)).*((yf - posy(24,20))./((xf - posx(24,20)).^2 + (yf - posy(24,20)).^2));
E_yr2421 = (rho_r9/(2*pi*e_0)).*((yf - posy(24,21))./((xf - posx(24,21)).^2 + (yf - posy(24,21)).^2));
E_yr2422 = (rho_r10/(2*pi*e_0)).*((yf - posy(24,22))./((xf - posx(24,22)).^2 + (yf - posy(24,22)).^2));
E_yr2423 = (rho_r11/(2*pi*e_0)).*((yf - posy(24,23))./((xf - posx(24,23)).^2 + (yf - posy(24,23)).^2));

%% E_yi21 componente y imaginario campo elétrico condutor 2 fase b

E_yi12 = (-rho_i2/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yi13 = (-rho_i3/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yi14 = (-rho_i4/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yi15 = (-rho_i5/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yi16 = (-rho_i6/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));
E_yi17 = (-rho_i7/(2*pi*e_0)).*((yf - posy(1,7))./((xf - posx(1,7)).^2 + (yf - posy(1,7)).^2));
E_yi18 = (-rho_i8/(2*pi*e_0)).*((yf - posy(1,8))./((xf - posx(1,8)).^2 + (yf - posy(1,8)).^2));
E_yi19 = (-rho_i9/(2*pi*e_0)).*((yf - posy(1,9))./((xf - posx(1,9)).^2 + (yf - posy(1,9)).^2));
E_yi110 = (-rho_i10/(2*pi*e_0)).*((yf - posy(1,10))./((xf - posx(1,10)).^2 + (yf - posy(1,10)).^2));
E_yi1_11 = (-rho_i11/(2*pi*e_0)).*((yf - posy(1,11))./((xf - posx(1,11)).^2 + (yf - posy(1,11)).^2));
E_yi1_12 = (-rho_i12/(2*pi*e_0)).*((yf - posy(1,12))./((xf - posx(1,12)).^2 + (yf - posy(1,12)).^2));
E_yi1_13 = (rho_i1/(2*pi*e_0)).*((yf - posy(1,13))./((xf - posx(1,13)).^2 + (yf - posy(1,13)).^2));
E_yi1_14 = (rho_i2/(2*pi*e_0)).*((yf - posy(1,14))./((xf - posx(1,14)).^2 + (yf - posy(1,14)).^2));
E_yi1_15 = (rho_i3/(2*pi*e_0)).*((yf - posy(1,15))./((xf - posx(1,15)).^2 + (yf - posy(1,15)).^2));
E_yi1_16 = (rho_i4/(2*pi*e_0)).*((yf - posy(1,16))./((xf - posx(1,16)).^2 + (yf - posy(1,16)).^2));
E_yi1_17 = (rho_i5/(2*pi*e_0)).*((yf - posy(1,17))./((xf - posx(1,17)).^2 + (yf - posy(1,17)).^2));
E_yi1_18 = (rho_i6/(2*pi*e_0)).*((yf - posy(1,18))./((xf - posx(1,18)).^2 + (yf - posy(1,18)).^2));
E_yi1_19 = (rho_i7/(2*pi*e_0)).*((yf - posy(1,19))./((xf - posx(1,19)).^2 + (yf - posy(1,19)).^2));
E_yi1_20 = (rho_i8/(2*pi*e_0)).*((yf - posy(1,20))./((xf - posx(1,20)).^2 + (yf - posy(1,20)).^2));
E_yi1_21 = (rho_i9/(2*pi*e_0)).*((yf - posy(1,21))./((xf - posx(1,21)).^2 + (yf - posy(1,21)).^2));
E_yi1_22 = (rho_i10/(2*pi*e_0)).*((yf - posy(1,22))./((xf - posx(1,22)).^2 + (yf - posy(1,22)).^2));
E_yi1_23 = (rho_i11/(2*pi*e_0)).*((yf - posy(1,23))./((xf - posx(1,23)).^2 + (yf - posy(1,23)).^2));
E_yi1_24 = (rho_i12/(2*pi*e_0)).*((yf - posy(1,24))./((xf - posx(1,24)).^2 + (yf - posy(1,24)).^2));

E_yi21 = (-rho_i1/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yi23 = (-rho_i3/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yi24 = (-rho_i4/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yi25 = (-rho_i5/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yi26 = (-rho_i6/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));
E_yi27 = (-rho_i7/(2*pi*e_0)).*((yf - posy(2,7))./((xf - posx(2,7)).^2 + (yf - posy(2,7)).^2));
E_yi28 = (-rho_i8/(2*pi*e_0)).*((yf - posy(2,8))./((xf - posx(2,8)).^2 + (yf - posy(2,8)).^2));
E_yi29 = (-rho_i9/(2*pi*e_0)).*((yf - posy(2,9))./((xf - posx(2,9)).^2 + (yf - posy(2,9)).^2));
E_yi210 = (-rho_i10/(2*pi*e_0)).*((yf - posy(2,10))./((xf - posx(2,10)).^2 + (yf - posy(2,10)).^2));
E_yi2_11 = (-rho_i11/(2*pi*e_0)).*((yf - posy(2,11))./((xf - posx(2,11)).^2 + (yf - posy(2,11)).^2));
E_yi2_12 = (-rho_i12/(2*pi*e_0)).*((yf - posy(2,12))./((xf - posx(2,12)).^2 + (yf - posy(2,12)).^2));
E_yi2_13 = (rho_i1/(2*pi*e_0)).*((yf - posy(2,13))./((xf - posx(2,13)).^2 + (yf - posy(2,13)).^2));
E_yi2_14 = (rho_i2/(2*pi*e_0)).*((yf - posy(2,14))./((xf - posx(2,14)).^2 + (yf - posy(2,14)).^2));
E_yi2_15 = (rho_i3/(2*pi*e_0)).*((yf - posy(2,15))./((xf - posx(2,15)).^2 + (yf - posy(2,15)).^2));
E_yi2_16 = (rho_i4/(2*pi*e_0)).*((yf - posy(2,16))./((xf - posx(2,16)).^2 + (yf - posy(2,16)).^2));
E_yi2_17 = (rho_i5/(2*pi*e_0)).*((yf - posy(2,17))./((xf - posx(2,17)).^2 + (yf - posy(2,17)).^2));
E_yi2_18 = (rho_i6/(2*pi*e_0)).*((yf - posy(2,18))./((xf - posx(2,18)).^2 + (yf - posy(2,18)).^2));
E_yi2_19 = (rho_i7/(2*pi*e_0)).*((yf - posy(2,19))./((xf - posx(2,19)).^2 + (yf - posy(2,19)).^2));
E_yi2_20 = (rho_i8/(2*pi*e_0)).*((yf - posy(2,20))./((xf - posx(2,20)).^2 + (yf - posy(2,20)).^2));
E_yi2_21 = (rho_i9/(2*pi*e_0)).*((yf - posy(2,21))./((xf - posx(2,21)).^2 + (yf - posy(2,21)).^2));
E_yi2_22 = (rho_i10/(2*pi*e_0)).*((yf - posy(2,22))./((xf - posx(2,22)).^2 + (yf - posy(2,22)).^2));
E_yi2_23 = (rho_i11/(2*pi*e_0)).*((yf - posy(2,23))./((xf - posx(2,23)).^2 + (yf - posy(2,23)).^2));
E_yi2_24 = (rho_i12/(2*pi*e_0)).*((yf - posy(2,24))./((xf - posx(2,24)).^2 + (yf - posy(2,24)).^2));

E_yi31 = (-rho_i1/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yi32 = (-rho_i2/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yi34 = (-rho_i4/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yi35 = (-rho_i5/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yi36 = (-rho_i6/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));
E_yi37 = (-rho_i7/(2*pi*e_0)).*((yf - posy(3,7))./((xf - posx(3,7)).^2 + (yf - posy(3,7)).^2));
E_yi38 = (-rho_i8/(2*pi*e_0)).*((yf - posy(3,8))./((xf - posx(3,8)).^2 + (yf - posy(3,8)).^2));
E_yi39 = (-rho_i9/(2*pi*e_0)).*((yf - posy(3,9))./((xf - posx(3,9)).^2 + (yf - posy(3,9)).^2));
E_yi310 = (-rho_i10/(2*pi*e_0)).*((yf - posy(3,10))./((xf - posx(3,10)).^2 + (yf - posy(3,10)).^2));
E_yi311 = (-rho_i11/(2*pi*e_0)).*((yf - posy(3,11))./((xf - posx(3,11)).^2 + (yf - posy(3,11)).^2));
E_yi312 = (-rho_i12/(2*pi*e_0)).*((yf - posy(3,12))./((xf - posx(3,12)).^2 + (yf - posy(3,12)).^2));
E_yi313 = (rho_i1/(2*pi*e_0)).*((yf - posy(3,13))./((xf - posx(3,13)).^2 + (yf - posy(3,13)).^2));
E_yi314 = (rho_i2/(2*pi*e_0)).*((yf - posy(3,14))./((xf - posx(3,14)).^2 + (yf - posy(3,14)).^2));
E_yi315 = (rho_i3/(2*pi*e_0)).*((yf - posy(3,15))./((xf - posx(3,15)).^2 + (yf - posy(3,15)).^2));
E_yi316 = (rho_i4/(2*pi*e_0)).*((yf - posy(3,16))./((xf - posx(3,16)).^2 + (yf - posy(3,16)).^2));
E_yi317 = (rho_i5/(2*pi*e_0)).*((yf - posy(3,17))./((xf - posx(3,17)).^2 + (yf - posy(3,17)).^2));
E_yi318 = (rho_i6/(2*pi*e_0)).*((yf - posy(3,18))./((xf - posx(3,18)).^2 + (yf - posy(3,18)).^2));
E_yi319 = (rho_i7/(2*pi*e_0)).*((yf - posy(3,19))./((xf - posx(3,19)).^2 + (yf - posy(3,19)).^2));
E_yi320 = (rho_i8/(2*pi*e_0)).*((yf - posy(3,20))./((xf - posx(3,20)).^2 + (yf - posy(3,20)).^2));
E_yi321 = (rho_i9/(2*pi*e_0)).*((yf - posy(3,21))./((xf - posx(3,21)).^2 + (yf - posy(3,21)).^2));
E_yi322 = (rho_i10/(2*pi*e_0)).*((yf - posy(3,22))./((xf - posx(3,22)).^2 + (yf - posy(3,22)).^2));
E_yi323 = (rho_i11/(2*pi*e_0)).*((yf - posy(3,23))./((xf - posx(3,23)).^2 + (yf - posy(3,23)).^2));
E_yi324 = (rho_i12/(2*pi*e_0)).*((yf - posy(3,24))./((xf - posx(3,24)).^2 + (yf - posy(3,24)).^2));

E_yi41 = (-rho_i1/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yi42 = (-rho_i2/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yi43 = (-rho_i3/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yi45 = (-rho_i5/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yi46 = (-rho_i6/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));
E_yi47 = (-rho_i7/(2*pi*e_0)).*((yf - posy(4,7))./((xf - posx(4,7)).^2 + (yf - posy(4,7)).^2));
E_yi48 = (-rho_i8/(2*pi*e_0)).*((yf - posy(4,8))./((xf - posx(4,8)).^2 + (yf - posy(4,8)).^2));
E_yi49 = (-rho_i9/(2*pi*e_0)).*((yf - posy(4,9))./((xf - posx(4,9)).^2 + (yf - posy(4,9)).^2));
E_yi410 = (-rho_i10/(2*pi*e_0)).*((yf - posy(4,10))./((xf - posx(4,10)).^2 + (yf - posy(4,10)).^2));
E_yi411 = (-rho_i11/(2*pi*e_0)).*((yf - posy(4,11))./((xf - posx(4,11)).^2 + (yf - posy(4,11)).^2));
E_yi412 = (-rho_i12/(2*pi*e_0)).*((yf - posy(4,12))./((xf - posx(4,12)).^2 + (yf - posy(4,12)).^2));
E_yi413 = (rho_i1/(2*pi*e_0)).*((yf - posy(4,13))./((xf - posx(4,13)).^2 + (yf - posy(4,13)).^2));
E_yi414 = (rho_i2/(2*pi*e_0)).*((yf - posy(4,14))./((xf - posx(4,14)).^2 + (yf - posy(4,14)).^2));
E_yi415 = (rho_i3/(2*pi*e_0)).*((yf - posy(4,15))./((xf - posx(4,15)).^2 + (yf - posy(4,15)).^2));
E_yi416 = (rho_i4/(2*pi*e_0)).*((yf - posy(4,16))./((xf - posx(4,16)).^2 + (yf - posy(4,16)).^2));
E_yi417 = (rho_i5/(2*pi*e_0)).*((yf - posy(4,17))./((xf - posx(4,17)).^2 + (yf - posy(4,17)).^2));
E_yi418 = (rho_i6/(2*pi*e_0)).*((yf - posy(4,18))./((xf - posx(4,18)).^2 + (yf - posy(4,18)).^2));
E_yi419 = (rho_i7/(2*pi*e_0)).*((yf - posy(4,19))./((xf - posx(4,19)).^2 + (yf - posy(4,19)).^2));
E_yi420 = (rho_i8/(2*pi*e_0)).*((yf - posy(4,20))./((xf - posx(4,20)).^2 + (yf - posy(4,20)).^2));
E_yi421 = (rho_i9/(2*pi*e_0)).*((yf - posy(4,21))./((xf - posx(4,21)).^2 + (yf - posy(4,21)).^2));
E_yi422 = (rho_i10/(2*pi*e_0)).*((yf - posy(4,22))./((xf - posx(4,22)).^2 + (yf - posy(4,22)).^2));
E_yi423 = (rho_i11/(2*pi*e_0)).*((yf - posy(4,23))./((xf - posx(4,23)).^2 + (yf - posy(4,23)).^2));
E_yi424 = (rho_i12/(2*pi*e_0)).*((yf - posy(4,24))./((xf - posx(4,24)).^2 + (yf - posy(4,24)).^2));

E_yi51 = (-rho_i1/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yi52 = (-rho_i2/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yi53 = (-rho_i3/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yi54 = (-rho_i4/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yi56 = (-rho_i6/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));
E_yi57 = (-rho_i7/(2*pi*e_0)).*((yf - posy(5,7))./((xf - posx(5,7)).^2 + (yf - posy(5,7)).^2));
E_yi58 = (-rho_i8/(2*pi*e_0)).*((yf - posy(5,8))./((xf - posx(5,8)).^2 + (yf - posy(5,8)).^2));
E_yi59 = (-rho_i9/(2*pi*e_0)).*((yf - posy(5,9))./((xf - posx(5,9)).^2 + (yf - posy(5,9)).^2));
E_yi510 = (-rho_i10/(2*pi*e_0)).*((yf - posy(5,10))./((xf - posx(5,10)).^2 + (yf - posy(5,10)).^2));
E_yi511 = (-rho_i11/(2*pi*e_0)).*((yf - posy(5,11))./((xf - posx(5,11)).^2 + (yf - posy(5,11)).^2));
E_yi512 = (-rho_i12/(2*pi*e_0)).*((yf - posy(5,12))./((xf - posx(5,12)).^2 + (yf - posy(5,12)).^2));
E_yi513 = (rho_i1/(2*pi*e_0)).*((yf - posy(5,13))./((xf - posx(5,13)).^2 + (yf - posy(5,13)).^2));
E_yi514 = (rho_i2/(2*pi*e_0)).*((yf - posy(5,14))./((xf - posx(5,14)).^2 + (yf - posy(5,14)).^2));
E_yi515 = (rho_i3/(2*pi*e_0)).*((yf - posy(5,15))./((xf - posx(5,15)).^2 + (yf - posy(5,15)).^2));
E_yi516 = (rho_i4/(2*pi*e_0)).*((yf - posy(5,16))./((xf - posx(5,16)).^2 + (yf - posy(5,16)).^2));
E_yi517 = (rho_i5/(2*pi*e_0)).*((yf - posy(5,17))./((xf - posx(5,17)).^2 + (yf - posy(5,17)).^2));
E_yi518 = (rho_i6/(2*pi*e_0)).*((yf - posy(5,18))./((xf - posx(5,18)).^2 + (yf - posy(5,18)).^2));
E_yi519 = (rho_i7/(2*pi*e_0)).*((yf - posy(5,19))./((xf - posx(5,19)).^2 + (yf - posy(5,19)).^2));
E_yi520 = (rho_i8/(2*pi*e_0)).*((yf - posy(5,20))./((xf - posx(5,20)).^2 + (yf - posy(5,20)).^2));
E_yi521 = (rho_i9/(2*pi*e_0)).*((yf - posy(5,21))./((xf - posx(5,21)).^2 + (yf - posy(5,21)).^2));
E_yi522 = (rho_i10/(2*pi*e_0)).*((yf - posy(5,22))./((xf - posx(5,22)).^2 + (yf - posy(5,22)).^2));
E_yi523 = (rho_i11/(2*pi*e_0)).*((yf - posy(5,23))./((xf - posx(5,23)).^2 + (yf - posy(5,23)).^2));
E_yi524 = (rho_i12/(2*pi*e_0)).*((yf - posy(5,24))./((xf - posx(5,24)).^2 + (yf - posy(5,24)).^2));

E_yi61 = (-rho_i1/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yi62 = (-rho_i2/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yi63 = (-rho_i3/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yi64 = (-rho_i4/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yi65 = (-rho_i5/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));
E_yi67 = (-rho_i7/(2*pi*e_0)).*((yf - posy(6,7))./((xf - posx(6,7)).^2 + (yf - posy(6,7)).^2));
E_yi68 = (-rho_i8/(2*pi*e_0)).*((yf - posy(6,8))./((xf - posx(6,8)).^2 + (yf - posy(6,8)).^2));
E_yi69 = (-rho_i9/(2*pi*e_0)).*((yf - posy(6,9))./((xf - posx(6,9)).^2 + (yf - posy(6,9)).^2));
E_yi610 = (-rho_i10/(2*pi*e_0)).*((yf - posy(6,10))./((xf - posx(6,10)).^2 + (yf - posy(6,10)).^2));
E_yi611 = (-rho_i11/(2*pi*e_0)).*((yf - posy(6,11))./((xf - posx(6,11)).^2 + (yf - posy(6,11)).^2));
E_yi612 = (-rho_i12/(2*pi*e_0)).*((yf - posy(6,12))./((xf - posx(6,12)).^2 + (yf - posy(6,12)).^2));
E_yi613 = (rho_i1/(2*pi*e_0)).*((yf - posy(6,13))./((xf - posx(6,13)).^2 + (yf - posy(6,13)).^2));
E_yi614 = (rho_i2/(2*pi*e_0)).*((yf - posy(6,14))./((xf - posx(6,14)).^2 + (yf - posy(6,14)).^2));
E_yi615 = (rho_i3/(2*pi*e_0)).*((yf - posy(6,15))./((xf - posx(6,15)).^2 + (yf - posy(6,15)).^2));
E_yi616 = (rho_i4/(2*pi*e_0)).*((yf - posy(6,16))./((xf - posx(6,16)).^2 + (yf - posy(6,16)).^2));
E_yi617 = (rho_i5/(2*pi*e_0)).*((yf - posy(6,17))./((xf - posx(6,17)).^2 + (yf - posy(6,17)).^2));
E_yi618 = (rho_i6/(2*pi*e_0)).*((yf - posy(6,18))./((xf - posx(6,18)).^2 + (yf - posy(6,18)).^2));
E_yi619 = (rho_i7/(2*pi*e_0)).*((yf - posy(6,19))./((xf - posx(6,19)).^2 + (yf - posy(6,19)).^2));
E_yi620 = (rho_i8/(2*pi*e_0)).*((yf - posy(6,20))./((xf - posx(6,20)).^2 + (yf - posy(6,20)).^2));
E_yi621 = (rho_i9/(2*pi*e_0)).*((yf - posy(6,21))./((xf - posx(6,21)).^2 + (yf - posy(6,21)).^2));
E_yi622 = (rho_i10/(2*pi*e_0)).*((yf - posy(6,22))./((xf - posx(6,22)).^2 + (yf - posy(6,22)).^2));
E_yi623 = (rho_i11/(2*pi*e_0)).*((yf - posy(6,23))./((xf - posx(6,23)).^2 + (yf - posy(6,23)).^2));
E_yi624 = (rho_i12/(2*pi*e_0)).*((yf - posy(6,24))./((xf - posx(6,24)).^2 + (yf - posy(6,24)).^2));

E_yi71 = (-rho_i1/(2*pi*e_0)).*((yf - posy(7,1))./((xf - posx(7,1)).^2 + (yf - posy(7,1)).^2));
E_yi72 = (-rho_i2/(2*pi*e_0)).*((yf - posy(7,2))./((xf - posx(7,2)).^2 + (yf - posy(7,2)).^2));
E_yi73 = (-rho_i3/(2*pi*e_0)).*((yf - posy(7,3))./((xf - posx(7,3)).^2 + (yf - posy(7,3)).^2));
E_yi74 = (-rho_i4/(2*pi*e_0)).*((yf - posy(7,4))./((xf - posx(7,4)).^2 + (yf - posy(7,4)).^2));
E_yi75 = (-rho_i5/(2*pi*e_0)).*((yf - posy(7,5))./((xf - posx(7,5)).^2 + (yf - posy(7,5)).^2));
E_yi76 = (-rho_i6/(2*pi*e_0)).*((yf - posy(7,6))./((xf - posx(7,6)).^2 + (yf - posy(7,6)).^2));
E_yi78 = (-rho_i8/(2*pi*e_0)).*((yf - posy(7,8))./((xf - posx(7,8)).^2 + (yf - posy(7,8)).^2));
E_yi79 = (-rho_i9/(2*pi*e_0)).*((yf - posy(7,9))./((xf - posx(7,9)).^2 + (yf - posy(7,9)).^2));
E_yi710 = (-rho_i10/(2*pi*e_0)).*((yf - posy(7,10))./((xf - posx(7,10)).^2 + (yf - posy(7,10)).^2));
E_yi711 = (-rho_i11/(2*pi*e_0)).*((yf - posy(7,11))./((xf - posx(7,11)).^2 + (yf - posy(7,11)).^2));
E_yi712 = (-rho_i12/(2*pi*e_0)).*((yf - posy(7,12))./((xf - posx(7,12)).^2 + (yf - posy(7,12)).^2));
E_yi713 = (rho_i1/(2*pi*e_0)).*((yf - posy(7,13))./((xf - posx(7,13)).^2 + (yf - posy(7,13)).^2));
E_yi714 = (rho_i2/(2*pi*e_0)).*((yf - posy(7,14))./((xf - posx(7,14)).^2 + (yf - posy(7,14)).^2));
E_yi715 = (rho_i3/(2*pi*e_0)).*((yf - posy(7,15))./((xf - posx(7,15)).^2 + (yf - posy(7,15)).^2));
E_yi716 = (rho_i4/(2*pi*e_0)).*((yf - posy(7,16))./((xf - posx(7,16)).^2 + (yf - posy(7,16)).^2));
E_yi717 = (rho_i5/(2*pi*e_0)).*((yf - posy(7,17))./((xf - posx(7,17)).^2 + (yf - posy(7,17)).^2));
E_yi718 = (rho_i6/(2*pi*e_0)).*((yf - posy(7,18))./((xf - posx(7,18)).^2 + (yf - posy(7,18)).^2));
E_yi719 = (rho_i7/(2*pi*e_0)).*((yf - posy(7,19))./((xf - posx(7,19)).^2 + (yf - posy(7,19)).^2));
E_yi720 = (rho_i8/(2*pi*e_0)).*((yf - posy(7,20))./((xf - posx(7,20)).^2 + (yf - posy(7,20)).^2));
E_yi721 = (rho_i9/(2*pi*e_0)).*((yf - posy(7,21))./((xf - posx(7,21)).^2 + (yf - posy(7,21)).^2));
E_yi722 = (rho_i10/(2*pi*e_0)).*((yf - posy(7,22))./((xf - posx(7,22)).^2 + (yf - posy(7,22)).^2));
E_yi723 = (rho_i11/(2*pi*e_0)).*((yf - posy(7,23))./((xf - posx(7,23)).^2 + (yf - posy(7,23)).^2));
E_yi724 = (rho_i12/(2*pi*e_0)).*((yf - posy(7,24))./((xf - posx(7,24)).^2 + (yf - posy(7,24)).^2));

E_yi81 = (-rho_i1/(2*pi*e_0)).*((yf - posy(8,1))./((xf - posx(8,1)).^2 + (yf - posy(8,1)).^2));
E_yi82 = (-rho_i2/(2*pi*e_0)).*((yf - posy(8,2))./((xf - posx(8,2)).^2 + (yf - posy(8,2)).^2));
E_yi83 = (-rho_i3/(2*pi*e_0)).*((yf - posy(8,3))./((xf - posx(8,3)).^2 + (yf - posy(8,3)).^2));
E_yi84 = (-rho_i4/(2*pi*e_0)).*((yf - posy(8,4))./((xf - posx(8,4)).^2 + (yf - posy(8,4)).^2));
E_yi85 = (-rho_i5/(2*pi*e_0)).*((yf - posy(8,5))./((xf - posx(8,5)).^2 + (yf - posy(8,5)).^2));
E_yi86 = (-rho_i6/(2*pi*e_0)).*((yf - posy(8,6))./((xf - posx(8,6)).^2 + (yf - posy(8,6)).^2));
E_yi87 = (-rho_i7/(2*pi*e_0)).*((yf - posy(8,7))./((xf - posx(8,7)).^2 + (yf - posy(8,7)).^2));
E_yi89 = (-rho_i9/(2*pi*e_0)).*((yf - posy(8,9))./((xf - posx(8,9)).^2 + (yf - posy(8,9)).^2));
E_yi810 = (-rho_i10/(2*pi*e_0)).*((yf - posy(8,10))./((xf - posx(8,10)).^2 + (yf - posy(8,10)).^2));
E_yi811 = (-rho_i11/(2*pi*e_0)).*((yf - posy(8,11))./((xf - posx(8,11)).^2 + (yf - posy(8,11)).^2));
E_yi812 = (-rho_i12/(2*pi*e_0)).*((yf - posy(8,12))./((xf - posx(8,12)).^2 + (yf - posy(8,12)).^2));
E_yi813 = (rho_i1/(2*pi*e_0)).*((yf - posy(8,13))./((xf - posx(8,13)).^2 + (yf - posy(8,13)).^2));
E_yi814 = (rho_i2/(2*pi*e_0)).*((yf - posy(8,14))./((xf - posx(8,14)).^2 + (yf - posy(8,14)).^2));
E_yi815 = (rho_i3/(2*pi*e_0)).*((yf - posy(8,15))./((xf - posx(8,15)).^2 + (yf - posy(8,15)).^2));
E_yi816 = (rho_i4/(2*pi*e_0)).*((yf - posy(8,16))./((xf - posx(8,16)).^2 + (yf - posy(8,16)).^2));
E_yi817 = (rho_i5/(2*pi*e_0)).*((yf - posy(8,17))./((xf - posx(8,17)).^2 + (yf - posy(8,17)).^2));
E_yi818 = (rho_i6/(2*pi*e_0)).*((yf - posy(8,18))./((xf - posx(8,18)).^2 + (yf - posy(8,18)).^2));
E_yi819 = (rho_i7/(2*pi*e_0)).*((yf - posy(8,19))./((xf - posx(8,19)).^2 + (yf - posy(8,19)).^2));
E_yi820 = (rho_i8/(2*pi*e_0)).*((yf - posy(8,20))./((xf - posx(8,20)).^2 + (yf - posy(8,20)).^2));
E_yi821 = (rho_i9/(2*pi*e_0)).*((yf - posy(8,21))./((xf - posx(8,21)).^2 + (yf - posy(8,21)).^2));
E_yi822 = (rho_i10/(2*pi*e_0)).*((yf - posy(8,22))./((xf - posx(8,22)).^2 + (yf - posy(8,22)).^2));
E_yi823 = (rho_i11/(2*pi*e_0)).*((yf - posy(8,23))./((xf - posx(8,23)).^2 + (yf - posy(8,23)).^2));
E_yi824 = (rho_i12/(2*pi*e_0)).*((yf - posy(8,24))./((xf - posx(8,24)).^2 + (yf - posy(8,24)).^2));

E_yi91 = (-rho_i1/(2*pi*e_0)).*((yf - posy(9,1))./((xf - posx(9,1)).^2 + (yf - posy(9,1)).^2));
E_yi92 = (-rho_i2/(2*pi*e_0)).*((yf - posy(9,2))./((xf - posx(9,2)).^2 + (yf - posy(9,2)).^2));
E_yi93 = (-rho_i3/(2*pi*e_0)).*((yf - posy(9,3))./((xf - posx(9,3)).^2 + (yf - posy(9,3)).^2));
E_yi94 = (-rho_i4/(2*pi*e_0)).*((yf - posy(9,4))./((xf - posx(9,4)).^2 + (yf - posy(9,4)).^2));
E_yi95 = (-rho_i5/(2*pi*e_0)).*((yf - posy(9,5))./((xf - posx(9,5)).^2 + (yf - posy(9,5)).^2));
E_yi96 = (-rho_i6/(2*pi*e_0)).*((yf - posy(9,6))./((xf - posx(9,6)).^2 + (yf - posy(9,6)).^2));
E_yi97 = (-rho_i7/(2*pi*e_0)).*((yf - posy(9,7))./((xf - posx(9,7)).^2 + (yf - posy(9,7)).^2));
E_yi98 = (-rho_i8/(2*pi*e_0)).*((yf - posy(9,8))./((xf - posx(9,8)).^2 + (yf - posy(9,8)).^2));
E_yi910 = (-rho_i10/(2*pi*e_0)).*((yf - posy(9,10))./((xf - posx(9,10)).^2 + (yf - posy(9,10)).^2));
E_yi911 = (-rho_i11/(2*pi*e_0)).*((yf - posy(9,11))./((xf - posx(9,11)).^2 + (yf - posy(9,11)).^2));
E_yi912 = (-rho_i12/(2*pi*e_0)).*((yf - posy(9,12))./((xf - posx(9,12)).^2 + (yf - posy(9,12)).^2));
E_yi913 = (rho_i1/(2*pi*e_0)).*((yf - posy(9,13))./((xf - posx(9,13)).^2 + (yf - posy(9,13)).^2));
E_yi914 = (rho_i2/(2*pi*e_0)).*((yf - posy(9,14))./((xf - posx(9,14)).^2 + (yf - posy(9,14)).^2));
E_yi915 = (rho_i3/(2*pi*e_0)).*((yf - posy(9,15))./((xf - posx(9,15)).^2 + (yf - posy(9,15)).^2));
E_yi916 = (rho_i4/(2*pi*e_0)).*((yf - posy(9,16))./((xf - posx(9,16)).^2 + (yf - posy(9,16)).^2));
E_yi917 = (rho_i5/(2*pi*e_0)).*((yf - posy(9,17))./((xf - posx(9,17)).^2 + (yf - posy(9,17)).^2));
E_yi918 = (rho_i6/(2*pi*e_0)).*((yf - posy(9,18))./((xf - posx(9,18)).^2 + (yf - posy(9,18)).^2));
E_yi919 = (rho_i7/(2*pi*e_0)).*((yf - posy(9,19))./((xf - posx(9,19)).^2 + (yf - posy(9,19)).^2));
E_yi920 = (rho_i8/(2*pi*e_0)).*((yf - posy(9,20))./((xf - posx(9,20)).^2 + (yf - posy(9,20)).^2));
E_yi921 = (rho_i9/(2*pi*e_0)).*((yf - posy(9,21))./((xf - posx(9,21)).^2 + (yf - posy(9,21)).^2));
E_yi922 = (rho_i10/(2*pi*e_0)).*((yf - posy(9,22))./((xf - posx(9,22)).^2 + (yf - posy(9,22)).^2));
E_yi923 = (rho_i11/(2*pi*e_0)).*((yf - posy(9,23))./((xf - posx(9,23)).^2 + (yf - posy(9,23)).^2));
E_yi924 = (rho_i12/(2*pi*e_0)).*((yf - posy(9,24))./((xf - posx(9,24)).^2 + (yf - posy(9,24)).^2));

E_yi101 = (-rho_i1/(2*pi*e_0)).*((yf - posy(10,1))./((xf - posx(10,1)).^2 + (yf - posy(10,1)).^2));
E_yi102 = (-rho_i2/(2*pi*e_0)).*((yf - posy(10,2))./((xf - posx(10,2)).^2 + (yf - posy(10,2)).^2));
E_yi103 = (-rho_i3/(2*pi*e_0)).*((yf - posy(10,3))./((xf - posx(10,3)).^2 + (yf - posy(10,3)).^2));
E_yi104 = (-rho_i4/(2*pi*e_0)).*((yf - posy(10,4))./((xf - posx(10,4)).^2 + (yf - posy(10,4)).^2));
E_yi105 = (-rho_i5/(2*pi*e_0)).*((yf - posy(10,5))./((xf - posx(10,5)).^2 + (yf - posy(10,5)).^2));
E_yi106 = (-rho_i6/(2*pi*e_0)).*((yf - posy(10,6))./((xf - posx(10,6)).^2 + (yf - posy(10,6)).^2));
E_yi107 = (-rho_i7/(2*pi*e_0)).*((yf - posy(10,7))./((xf - posx(10,7)).^2 + (yf - posy(10,7)).^2));
E_yi108 = (-rho_i8/(2*pi*e_0)).*((yf - posy(10,8))./((xf - posx(10,8)).^2 + (yf - posy(10,8)).^2));
E_yi109 = (-rho_i9/(2*pi*e_0)).*((yf - posy(10,9))./((xf - posx(10,9)).^2 + (yf - posy(10,9)).^2));
E_yi1011 = (-rho_i11/(2*pi*e_0)).*((yf - posy(10,11))./((xf - posx(10,11)).^2 + (yf - posy(10,11)).^2));
E_yi1012 = (-rho_i12/(2*pi*e_0)).*((yf - posy(10,12))./((xf - posx(10,12)).^2 + (yf - posy(10,12)).^2));
E_yi1013 = (rho_i1/(2*pi*e_0)).*((yf - posy(10,13))./((xf - posx(10,13)).^2 + (yf - posy(10,13)).^2));
E_yi1014 = (rho_i2/(2*pi*e_0)).*((yf - posy(10,14))./((xf - posx(10,14)).^2 + (yf - posy(10,14)).^2));
E_yi1015 = (rho_i3/(2*pi*e_0)).*((yf - posy(10,15))./((xf - posx(10,15)).^2 + (yf - posy(10,15)).^2));
E_yi1016 = (rho_i4/(2*pi*e_0)).*((yf - posy(10,16))./((xf - posx(10,16)).^2 + (yf - posy(10,16)).^2));
E_yi1017 = (rho_i5/(2*pi*e_0)).*((yf - posy(10,17))./((xf - posx(10,17)).^2 + (yf - posy(10,17)).^2));
E_yi1018 = (rho_i6/(2*pi*e_0)).*((yf - posy(10,18))./((xf - posx(10,18)).^2 + (yf - posy(10,18)).^2));
E_yi10_19 = (rho_i7/(2*pi*e_0)).*((yf - posy(10,19))./((xf - posx(10,19)).^2 + (yf - posy(10,19)).^2));
E_yi10_20 = (rho_i8/(2*pi*e_0)).*((yf - posy(10,20))./((xf - posx(10,20)).^2 + (yf - posy(10,20)).^2));
E_yi10_21 = (rho_i9/(2*pi*e_0)).*((yf - posy(10,21))./((xf - posx(10,21)).^2 + (yf - posy(10,21)).^2));
E_yi10_22 = (rho_i10/(2*pi*e_0)).*((yf - posy(10,22))./((xf - posx(10,22)).^2 + (yf - posy(10,22)).^2));
E_yi10_23 = (rho_i11/(2*pi*e_0)).*((yf - posy(10,23))./((xf - posx(10,23)).^2 + (yf - posy(10,23)).^2));
E_yi10_24 = (rho_i12/(2*pi*e_0)).*((yf - posy(10,24))./((xf - posx(10,24)).^2 + (yf - posy(10,24)).^2));

E_yi11_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(11,1))./((xf - posx(11,1)).^2 + (yf - posy(11,1)).^2));
E_yi11_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(11,2))./((xf - posx(11,2)).^2 + (yf - posy(11,2)).^2));
E_yi113 = (-rho_i3/(2*pi*e_0)).*((yf - posy(11,3))./((xf - posx(11,3)).^2 + (yf - posy(11,3)).^2));
E_yi114 = (-rho_i4/(2*pi*e_0)).*((yf - posy(11,4))./((xf - posx(11,4)).^2 + (yf - posy(11,4)).^2));
E_yi115 = (-rho_i5/(2*pi*e_0)).*((yf - posy(11,5))./((xf - posx(11,5)).^2 + (yf - posy(11,5)).^2));
E_yi116 = (-rho_i6/(2*pi*e_0)).*((yf - posy(11,6))./((xf - posx(11,6)).^2 + (yf - posy(11,6)).^2));
E_yi117 = (-rho_i7/(2*pi*e_0)).*((yf - posy(11,7))./((xf - posx(11,7)).^2 + (yf - posy(11,7)).^2));
E_yi118 = (-rho_i8/(2*pi*e_0)).*((yf - posy(11,8))./((xf - posx(11,8)).^2 + (yf - posy(11,8)).^2));
E_yi119 = (-rho_i9/(2*pi*e_0)).*((yf - posy(11,9))./((xf - posx(11,9)).^2 + (yf - posy(11,9)).^2));
E_yi1110 = (-rho_i10/(2*pi*e_0)).*((yf - posy(11,10))./((xf - posx(11,10)).^2 + (yf - posy(11,10)).^2));
E_yi1112 = (-rho_i12/(2*pi*e_0)).*((yf - posy(11,12))./((xf - posx(11,12)).^2 + (yf - posy(11,12)).^2));
E_yi1113 = (rho_i1/(2*pi*e_0)).*((yf - posy(11,13))./((xf - posx(11,13)).^2 + (yf - posy(11,13)).^2));
E_yi1114 = (rho_i2/(2*pi*e_0)).*((yf - posy(11,14))./((xf - posx(11,14)).^2 + (yf - posy(11,14)).^2));
E_yi1115 = (rho_i3/(2*pi*e_0)).*((yf - posy(11,15))./((xf - posx(11,15)).^2 + (yf - posy(11,15)).^2));
E_yi1116 = (rho_i4/(2*pi*e_0)).*((yf - posy(11,16))./((xf - posx(11,16)).^2 + (yf - posy(11,16)).^2));
E_yi1117 = (rho_i5/(2*pi*e_0)).*((yf - posy(11,17))./((xf - posx(11,17)).^2 + (yf - posy(11,17)).^2));
E_yi1118 = (rho_i6/(2*pi*e_0)).*((yf - posy(11,18))./((xf - posx(11,18)).^2 + (yf - posy(11,18)).^2));
E_yi1019 = (rho_i7/(2*pi*e_0)).*((yf - posy(11,19))./((xf - posx(11,19)).^2 + (yf - posy(11,19)).^2));
E_yi1020 = (rho_i8/(2*pi*e_0)).*((yf - posy(11,20))./((xf - posx(11,20)).^2 + (yf - posy(11,20)).^2));
E_yi1021 = (rho_i9/(2*pi*e_0)).*((yf - posy(11,21))./((xf - posx(11,21)).^2 + (yf - posy(11,21)).^2));
E_yi1022 = (rho_i10/(2*pi*e_0)).*((yf - posy(11,22))./((xf - posx(11,22)).^2 + (yf - posy(11,22)).^2));
E_yi1023 = (rho_i11/(2*pi*e_0)).*((yf - posy(11,23))./((xf - posx(11,23)).^2 + (yf - posy(11,23)).^2));
E_yi1024 = (rho_i12/(2*pi*e_0)).*((yf - posy(11,24))./((xf - posx(11,24)).^2 + (yf - posy(11,24)).^2));

E_yi12_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(12,1))./((xf - posx(12,1)).^2 + (yf - posy(12,1)).^2));
E_yi12_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(12,2))./((xf - posx(12,2)).^2 + (yf - posy(12,2)).^2));
E_yi123 = (-rho_i3/(2*pi*e_0)).*((yf - posy(12,3))./((xf - posx(12,3)).^2 + (yf - posy(12,3)).^2));
E_yi124 = (-rho_i4/(2*pi*e_0)).*((yf - posy(12,4))./((xf - posx(12,4)).^2 + (yf - posy(12,4)).^2));
E_yi125 = (-rho_i5/(2*pi*e_0)).*((yf - posy(12,5))./((xf - posx(12,5)).^2 + (yf - posy(12,5)).^2));
E_yi126 = (-rho_i6/(2*pi*e_0)).*((yf - posy(12,6))./((xf - posx(12,6)).^2 + (yf - posy(12,6)).^2));
E_yi127 = (-rho_i7/(2*pi*e_0)).*((yf - posy(12,7))./((xf - posx(12,7)).^2 + (yf - posy(12,7)).^2));
E_yi128 = (-rho_i8/(2*pi*e_0)).*((yf - posy(12,8))./((xf - posx(12,8)).^2 + (yf - posy(12,8)).^2));
E_yi129 = (-rho_i9/(2*pi*e_0)).*((yf - posy(12,9))./((xf - posx(12,9)).^2 + (yf - posy(12,9)).^2));
E_yi1210 = (-rho_i10/(2*pi*e_0)).*((yf - posy(12,10))./((xf - posx(12,10)).^2 + (yf - posy(12,10)).^2));
E_yi1211 = (-rho_i11/(2*pi*e_0)).*((yf - posy(12,11))./((xf - posx(12,11)).^2 + (yf - posy(12,11)).^2));
E_yi1213 = (rho_i1/(2*pi*e_0)).*((yf - posy(12,13))./((xf - posx(12,13)).^2 + (yf - posy(12,13)).^2));
E_yi1214 = (rho_i2/(2*pi*e_0)).*((yf - posy(12,14))./((xf - posx(12,14)).^2 + (yf - posy(12,14)).^2));
E_yi1215 = (rho_i3/(2*pi*e_0)).*((yf - posy(12,15))./((xf - posx(12,15)).^2 + (yf - posy(12,15)).^2));
E_yi1216 = (rho_i4/(2*pi*e_0)).*((yf - posy(12,16))./((xf - posx(12,16)).^2 + (yf - posy(12,16)).^2));
E_yi1217 = (rho_i5/(2*pi*e_0)).*((yf - posy(12,17))./((xf - posx(12,17)).^2 + (yf - posy(12,17)).^2));
E_yi1218 = (rho_i6/(2*pi*e_0)).*((yf - posy(12,18))./((xf - posx(12,18)).^2 + (yf - posy(12,18)).^2));
E_yi1219 = (rho_i7/(2*pi*e_0)).*((yf - posy(12,19))./((xf - posx(12,19)).^2 + (yf - posy(12,19)).^2));
E_yi1220 = (rho_i8/(2*pi*e_0)).*((yf - posy(12,20))./((xf - posx(12,20)).^2 + (yf - posy(12,20)).^2));
E_yi1221 = (rho_i9/(2*pi*e_0)).*((yf - posy(12,21))./((xf - posx(12,21)).^2 + (yf - posy(12,21)).^2));
E_yi1222 = (rho_i10/(2*pi*e_0)).*((yf - posy(12,22))./((xf - posx(12,22)).^2 + (yf - posy(12,22)).^2));
E_yi1223 = (rho_i11/(2*pi*e_0)).*((yf - posy(12,23))./((xf - posx(12,23)).^2 + (yf - posy(12,23)).^2));
E_yi1224 = (rho_i12/(2*pi*e_0)).*((yf - posy(12,24))./((xf - posx(12,24)).^2 + (yf - posy(12,24)).^2));

E_yi131 = (-rho_i1/(2*pi*e_0)).*((yf - posy(13,1))./((xf - posx(13,1)).^2 + (yf - posy(13,1)).^2));
E_yi132 = (-rho_i2/(2*pi*e_0)).*((yf - posy(13,2))./((xf - posx(13,2)).^2 + (yf - posy(13,2)).^2));
E_yi133 = (-rho_i3/(2*pi*e_0)).*((yf - posy(13,3))./((xf - posx(13,3)).^2 + (yf - posy(13,3)).^2));
E_yi134 = (-rho_i4/(2*pi*e_0)).*((yf - posy(13,4))./((xf - posx(13,4)).^2 + (yf - posy(13,4)).^2));
E_yi135 = (-rho_i5/(2*pi*e_0)).*((yf - posy(13,5))./((xf - posx(13,5)).^2 + (yf - posy(13,5)).^2));
E_yi136 = (-rho_i6/(2*pi*e_0)).*((yf - posy(13,6))./((xf - posx(13,6)).^2 + (yf - posy(13,6)).^2));
E_yi137 = (-rho_i7/(2*pi*e_0)).*((yf - posy(13,7))./((xf - posx(13,7)).^2 + (yf - posy(13,7)).^2));
E_yi138 = (-rho_i8/(2*pi*e_0)).*((yf - posy(13,8))./((xf - posx(13,8)).^2 + (yf - posy(13,8)).^2));
E_yi139 = (-rho_i9/(2*pi*e_0)).*((yf - posy(13,9))./((xf - posx(13,9)).^2 + (yf - posy(13,9)).^2));
E_yi1310 = (-rho_i10/(2*pi*e_0)).*((yf - posy(13,10))./((xf - posx(13,10)).^2 + (yf - posy(13,10)).^2));
E_yi1311 = (-rho_i11/(2*pi*e_0)).*((yf - posy(13,11))./((xf - posx(13,11)).^2 + (yf - posy(13,11)).^2));
E_yi1312 = (-rho_i12/(2*pi*e_0)).*((yf - posy(13,12))./((xf - posx(13,12)).^2 + (yf - posy(13,12)).^2));
E_yi1314 = (rho_i2/(2*pi*e_0)).*((yf - posy(13,14))./((xf - posx(13,14)).^2 + (yf - posy(13,14)).^2));
E_yi1315 = (rho_i3/(2*pi*e_0)).*((yf - posy(13,15))./((xf - posx(13,15)).^2 + (yf - posy(13,15)).^2));
E_yi1316 = (rho_i4/(2*pi*e_0)).*((yf - posy(13,16))./((xf - posx(13,16)).^2 + (yf - posy(13,16)).^2));
E_yi1317 = (rho_i5/(2*pi*e_0)).*((yf - posy(13,17))./((xf - posx(13,17)).^2 + (yf - posy(13,17)).^2));
E_yi1318 = (rho_i6/(2*pi*e_0)).*((yf - posy(13,18))./((xf - posx(13,18)).^2 + (yf - posy(13,18)).^2));
E_yi1319 = (rho_i7/(2*pi*e_0)).*((yf - posy(13,19))./((xf - posx(13,19)).^2 + (yf - posy(13,19)).^2));
E_yi1320 = (rho_i8/(2*pi*e_0)).*((yf - posy(13,20))./((xf - posx(13,20)).^2 + (yf - posy(13,20)).^2));
E_yi1321 = (rho_i9/(2*pi*e_0)).*((yf - posy(13,21))./((xf - posx(13,21)).^2 + (yf - posy(13,21)).^2));
E_yi1322 = (rho_i10/(2*pi*e_0)).*((yf - posy(13,22))./((xf - posx(13,22)).^2 + (yf - posy(13,22)).^2));
E_yi1323 = (rho_i11/(2*pi*e_0)).*((yf - posy(13,23))./((xf - posx(13,23)).^2 + (yf - posy(13,23)).^2));
E_yi1324 = (rho_i12/(2*pi*e_0)).*((yf - posy(13,24))./((xf - posx(13,24)).^2 + (yf - posy(13,24)).^2));

E_yi141 = (-rho_i1/(2*pi*e_0)).*((yf - posy(14,1))./((xf - posx(14,1)).^2 + (yf - posy(14,1)).^2));
E_yi142 = (-rho_i2/(2*pi*e_0)).*((yf - posy(14,2))./((xf - posx(14,2)).^2 + (yf - posy(14,2)).^2));
E_yi143 = (-rho_i3/(2*pi*e_0)).*((yf - posy(14,3))./((xf - posx(14,3)).^2 + (yf - posy(14,3)).^2));
E_yi144 = (-rho_i4/(2*pi*e_0)).*((yf - posy(14,4))./((xf - posx(14,4)).^2 + (yf - posy(14,4)).^2));
E_yi145 = (-rho_i5/(2*pi*e_0)).*((yf - posy(14,5))./((xf - posx(14,5)).^2 + (yf - posy(14,5)).^2));
E_yi146 = (-rho_i6/(2*pi*e_0)).*((yf - posy(14,6))./((xf - posx(14,6)).^2 + (yf - posy(14,6)).^2));
E_yi147 = (-rho_i7/(2*pi*e_0)).*((yf - posy(14,7))./((xf - posx(14,7)).^2 + (yf - posy(14,7)).^2));
E_yi148 = (-rho_i8/(2*pi*e_0)).*((yf - posy(14,8))./((xf - posx(14,8)).^2 + (yf - posy(14,8)).^2));
E_yi149 = (-rho_i9/(2*pi*e_0)).*((yf - posy(14,9))./((xf - posx(14,9)).^2 + (yf - posy(14,9)).^2));
E_yi1410 = (-rho_i10/(2*pi*e_0)).*((yf - posy(14,10))./((xf - posx(14,10)).^2 + (yf - posy(14,10)).^2));
E_yi1411 = (-rho_i11/(2*pi*e_0)).*((yf - posy(14,11))./((xf - posx(14,11)).^2 + (yf - posy(14,11)).^2));
E_yi1412 = (-rho_i12/(2*pi*e_0)).*((yf - posy(14,12))./((xf - posx(14,12)).^2 + (yf - posy(14,12)).^2));
E_yi1413 = (rho_i1/(2*pi*e_0)).*((yf - posy(14,13))./((xf - posx(14,13)).^2 + (yf - posy(14,13)).^2));
E_yi1415 = (rho_i3/(2*pi*e_0)).*((yf - posy(14,15))./((xf - posx(14,15)).^2 + (yf - posy(14,15)).^2));
E_yi1416 = (rho_i4/(2*pi*e_0)).*((yf - posy(14,16))./((xf - posx(14,16)).^2 + (yf - posy(14,16)).^2));
E_yi1417 = (rho_i5/(2*pi*e_0)).*((yf - posy(14,17))./((xf - posx(14,17)).^2 + (yf - posy(14,17)).^2));
E_yi1418 = (rho_i6/(2*pi*e_0)).*((yf - posy(14,18))./((xf - posx(14,18)).^2 + (yf - posy(14,18)).^2));
E_yi1419 = (rho_i7/(2*pi*e_0)).*((yf - posy(14,19))./((xf - posx(14,19)).^2 + (yf - posy(14,19)).^2));
E_yi1420 = (rho_i8/(2*pi*e_0)).*((yf - posy(14,20))./((xf - posx(14,20)).^2 + (yf - posy(14,20)).^2));
E_yi1421 = (rho_i9/(2*pi*e_0)).*((yf - posy(14,21))./((xf - posx(14,21)).^2 + (yf - posy(14,21)).^2));
E_yi1422 = (rho_i10/(2*pi*e_0)).*((yf - posy(14,22))./((xf - posx(14,22)).^2 + (yf - posy(14,22)).^2));
E_yi1423 = (rho_i11/(2*pi*e_0)).*((yf - posy(14,23))./((xf - posx(14,23)).^2 + (yf - posy(14,23)).^2));
E_yi1424 = (rho_i12/(2*pi*e_0)).*((yf - posy(14,24))./((xf - posx(14,24)).^2 + (yf - posy(14,24)).^2));

E_yi151 = (-rho_i1/(2*pi*e_0)).*((yf - posy(15,1))./((xf - posx(15,1)).^2 + (yf - posy(15,1)).^2));
E_yi152 = (-rho_i2/(2*pi*e_0)).*((yf - posy(15,2))./((xf - posx(15,2)).^2 + (yf - posy(15,2)).^2));
E_yi153 = (-rho_i3/(2*pi*e_0)).*((yf - posy(15,3))./((xf - posx(15,3)).^2 + (yf - posy(15,3)).^2));
E_yi154 = (-rho_i4/(2*pi*e_0)).*((yf - posy(15,4))./((xf - posx(15,4)).^2 + (yf - posy(15,4)).^2));
E_yi155 = (-rho_i5/(2*pi*e_0)).*((yf - posy(15,5))./((xf - posx(15,5)).^2 + (yf - posy(15,5)).^2));
E_yi156 = (-rho_i6/(2*pi*e_0)).*((yf - posy(15,6))./((xf - posx(15,6)).^2 + (yf - posy(15,6)).^2));
E_yi157 = (-rho_i7/(2*pi*e_0)).*((yf - posy(15,7))./((xf - posx(15,7)).^2 + (yf - posy(15,7)).^2));
E_yi158 = (-rho_i8/(2*pi*e_0)).*((yf - posy(15,8))./((xf - posx(15,8)).^2 + (yf - posy(15,8)).^2));
E_yi159 = (-rho_i9/(2*pi*e_0)).*((yf - posy(15,9))./((xf - posx(15,9)).^2 + (yf - posy(15,9)).^2));
E_yi1510 = (-rho_i10/(2*pi*e_0)).*((yf - posy(15,10))./((xf - posx(15,10)).^2 + (yf - posy(15,10)).^2));
E_yi1511 = (-rho_i11/(2*pi*e_0)).*((yf - posy(15,11))./((xf - posx(15,11)).^2 + (yf - posy(15,11)).^2));
E_yi1512 = (-rho_i12/(2*pi*e_0)).*((yf - posy(15,12))./((xf - posx(15,12)).^2 + (yf - posy(15,12)).^2));
E_yi1513 = (rho_i1/(2*pi*e_0)).*((yf - posy(15,13))./((xf - posx(15,13)).^2 + (yf - posy(15,13)).^2));
E_yi1514 = (rho_i2/(2*pi*e_0)).*((yf - posy(15,14))./((xf - posx(15,14)).^2 + (yf - posy(15,14)).^2));
E_yi1516 = (rho_i4/(2*pi*e_0)).*((yf - posy(15,16))./((xf - posx(15,16)).^2 + (yf - posy(15,16)).^2));
E_yi1517 = (rho_i5/(2*pi*e_0)).*((yf - posy(15,17))./((xf - posx(15,17)).^2 + (yf - posy(15,17)).^2));
E_yi1518 = (rho_i6/(2*pi*e_0)).*((yf - posy(15,18))./((xf - posx(15,18)).^2 + (yf - posy(15,18)).^2));
E_yi1519 = (rho_i7/(2*pi*e_0)).*((yf - posy(15,19))./((xf - posx(15,19)).^2 + (yf - posy(15,19)).^2));
E_yi1520 = (rho_i8/(2*pi*e_0)).*((yf - posy(15,20))./((xf - posx(15,20)).^2 + (yf - posy(15,20)).^2));
E_yi1521 = (rho_i9/(2*pi*e_0)).*((yf - posy(15,21))./((xf - posx(15,21)).^2 + (yf - posy(15,21)).^2));
E_yi1522 = (rho_i10/(2*pi*e_0)).*((yf - posy(15,22))./((xf - posx(15,22)).^2 + (yf - posy(15,22)).^2));
E_yi1523 = (rho_i11/(2*pi*e_0)).*((yf - posy(15,23))./((xf - posx(15,23)).^2 + (yf - posy(15,23)).^2));
E_yi1524 = (rho_i12/(2*pi*e_0)).*((yf - posy(15,24))./((xf - posx(15,24)).^2 + (yf - posy(15,24)).^2));

E_yi161 = (-rho_i1/(2*pi*e_0)).*((yf - posy(16,1))./((xf - posx(16,1)).^2 + (yf - posy(16,1)).^2));
E_yi162 = (-rho_i2/(2*pi*e_0)).*((yf - posy(16,2))./((xf - posx(16,2)).^2 + (yf - posy(16,2)).^2));
E_yi163 = (-rho_i3/(2*pi*e_0)).*((yf - posy(16,3))./((xf - posx(16,3)).^2 + (yf - posy(16,3)).^2));
E_yi164 = (-rho_i4/(2*pi*e_0)).*((yf - posy(16,4))./((xf - posx(16,4)).^2 + (yf - posy(16,4)).^2));
E_yi165 = (-rho_i5/(2*pi*e_0)).*((yf - posy(16,5))./((xf - posx(16,5)).^2 + (yf - posy(16,5)).^2));
E_yi166 = (-rho_i6/(2*pi*e_0)).*((yf - posy(16,6))./((xf - posx(16,6)).^2 + (yf - posy(16,6)).^2));
E_yi167 = (-rho_i7/(2*pi*e_0)).*((yf - posy(16,7))./((xf - posx(16,7)).^2 + (yf - posy(16,7)).^2));
E_yi168 = (-rho_i8/(2*pi*e_0)).*((yf - posy(16,8))./((xf - posx(16,8)).^2 + (yf - posy(16,8)).^2));
E_yi169 = (-rho_i9/(2*pi*e_0)).*((yf - posy(16,9))./((xf - posx(16,9)).^2 + (yf - posy(16,9)).^2));
E_yi1610 = (-rho_i10/(2*pi*e_0)).*((yf - posy(16,10))./((xf - posx(16,10)).^2 + (yf - posy(16,10)).^2));
E_yi1611 = (-rho_i11/(2*pi*e_0)).*((yf - posy(16,11))./((xf - posx(16,11)).^2 + (yf - posy(16,11)).^2));
E_yi1612 = (-rho_i12/(2*pi*e_0)).*((yf - posy(16,12))./((xf - posx(16,12)).^2 + (yf - posy(16,12)).^2));
E_yi1613 = (rho_i1/(2*pi*e_0)).*((yf - posy(16,13))./((xf - posx(16,13)).^2 + (yf - posy(16,13)).^2));
E_yi1614 = (rho_i2/(2*pi*e_0)).*((yf - posy(16,14))./((xf - posx(16,14)).^2 + (yf - posy(16,14)).^2));
E_yi1615 = (rho_i3/(2*pi*e_0)).*((yf - posy(16,15))./((xf - posx(16,15)).^2 + (yf - posy(16,15)).^2));
E_yi1617 = (rho_i5/(2*pi*e_0)).*((yf - posy(16,17))./((xf - posx(16,17)).^2 + (yf - posy(16,17)).^2));
E_yi1618 = (rho_i6/(2*pi*e_0)).*((yf - posy(16,18))./((xf - posx(16,18)).^2 + (yf - posy(16,18)).^2));
E_yi1619 = (rho_i7/(2*pi*e_0)).*((yf - posy(16,19))./((xf - posx(16,19)).^2 + (yf - posy(16,19)).^2));
E_yi1620 = (rho_i8/(2*pi*e_0)).*((yf - posy(16,20))./((xf - posx(16,20)).^2 + (yf - posy(16,20)).^2));
E_yi1621 = (rho_i9/(2*pi*e_0)).*((yf - posy(16,21))./((xf - posx(16,21)).^2 + (yf - posy(16,21)).^2));
E_yi1622 = (rho_i10/(2*pi*e_0)).*((yf - posy(16,22))./((xf - posx(16,22)).^2 + (yf - posy(16,22)).^2));
E_yi1623 = (rho_i11/(2*pi*e_0)).*((yf - posy(16,23))./((xf - posx(16,23)).^2 + (yf - posy(16,23)).^2));
E_yi1624 = (rho_i12/(2*pi*e_0)).*((yf - posy(16,24))./((xf - posx(16,24)).^2 + (yf - posy(16,24)).^2));

E_yi171 = (-rho_i1/(2*pi*e_0)).*((yf - posy(17,1))./((xf - posx(17,1)).^2 + (yf - posy(17,1)).^2));
E_yi172 = (-rho_i2/(2*pi*e_0)).*((yf - posy(17,2))./((xf - posx(17,2)).^2 + (yf - posy(17,2)).^2));
E_yi173 = (-rho_i3/(2*pi*e_0)).*((yf - posy(17,3))./((xf - posx(17,3)).^2 + (yf - posy(17,3)).^2));
E_yi174 = (-rho_i4/(2*pi*e_0)).*((yf - posy(17,4))./((xf - posx(17,4)).^2 + (yf - posy(17,4)).^2));
E_yi175 = (-rho_i5/(2*pi*e_0)).*((yf - posy(17,5))./((xf - posx(17,5)).^2 + (yf - posy(17,5)).^2));
E_yi176 = (-rho_i6/(2*pi*e_0)).*((yf - posy(17,6))./((xf - posx(17,6)).^2 + (yf - posy(17,6)).^2));
E_yi177 = (-rho_i7/(2*pi*e_0)).*((yf - posy(17,7))./((xf - posx(17,7)).^2 + (yf - posy(17,7)).^2));
E_yi178 = (-rho_i8/(2*pi*e_0)).*((yf - posy(17,8))./((xf - posx(17,8)).^2 + (yf - posy(17,8)).^2));
E_yi179 = (-rho_i9/(2*pi*e_0)).*((yf - posy(17,9))./((xf - posx(17,9)).^2 + (yf - posy(17,9)).^2));
E_yi1710 = (-rho_i10/(2*pi*e_0)).*((yf - posy(17,10))./((xf - posx(17,10)).^2 + (yf - posy(17,10)).^2));
E_yi1711 = (-rho_i11/(2*pi*e_0)).*((yf - posy(17,11))./((xf - posx(17,11)).^2 + (yf - posy(17,11)).^2));
E_yi1712 = (-rho_i12/(2*pi*e_0)).*((yf - posy(17,12))./((xf - posx(17,12)).^2 + (yf - posy(17,12)).^2));
E_yi1713 = (rho_i1/(2*pi*e_0)).*((yf - posy(17,13))./((xf - posx(17,13)).^2 + (yf - posy(17,13)).^2));
E_yi1714 = (rho_i2/(2*pi*e_0)).*((yf - posy(17,14))./((xf - posx(17,14)).^2 + (yf - posy(17,14)).^2));
E_yi1715 = (rho_i3/(2*pi*e_0)).*((yf - posy(17,15))./((xf - posx(17,15)).^2 + (yf - posy(17,15)).^2));
E_yi1716 = (rho_i4/(2*pi*e_0)).*((yf - posy(17,16))./((xf - posx(17,16)).^2 + (yf - posy(17,16)).^2));
E_yi1718 = (rho_i6/(2*pi*e_0)).*((yf - posy(17,18))./((xf - posx(17,18)).^2 + (yf - posy(17,18)).^2));
E_yi1719 = (rho_i7/(2*pi*e_0)).*((yf - posy(17,19))./((xf - posx(17,19)).^2 + (yf - posy(17,19)).^2));
E_yi1720 = (rho_i8/(2*pi*e_0)).*((yf - posy(17,20))./((xf - posx(17,20)).^2 + (yf - posy(17,20)).^2));
E_yi1721 = (rho_i9/(2*pi*e_0)).*((yf - posy(17,21))./((xf - posx(17,21)).^2 + (yf - posy(17,21)).^2));
E_yi1722 = (rho_i10/(2*pi*e_0)).*((yf - posy(17,22))./((xf - posx(17,22)).^2 + (yf - posy(17,22)).^2));
E_yi1723 = (rho_i11/(2*pi*e_0)).*((yf - posy(17,23))./((xf - posx(17,23)).^2 + (yf - posy(17,23)).^2));
E_yi1724 = (rho_i12/(2*pi*e_0)).*((yf - posy(17,24))./((xf - posx(17,24)).^2 + (yf - posy(17,24)).^2));

E_yi181 = (-rho_i1/(2*pi*e_0)).*((yf - posy(18,1))./((xf - posx(18,1)).^2 + (yf - posy(18,1)).^2));
E_yi182 = (-rho_i2/(2*pi*e_0)).*((yf - posy(18,2))./((xf - posx(18,2)).^2 + (yf - posy(18,2)).^2));
E_yi183 = (-rho_i3/(2*pi*e_0)).*((yf - posy(18,3))./((xf - posx(18,3)).^2 + (yf - posy(18,3)).^2));
E_yi184 = (-rho_i4/(2*pi*e_0)).*((yf - posy(18,4))./((xf - posx(18,4)).^2 + (yf - posy(18,4)).^2));
E_yi185 = (-rho_i5/(2*pi*e_0)).*((yf - posy(18,5))./((xf - posx(18,5)).^2 + (yf - posy(18,5)).^2));
E_yi186 = (-rho_i6/(2*pi*e_0)).*((yf - posy(18,6))./((xf - posx(18,6)).^2 + (yf - posy(18,6)).^2));
E_yi187 = (-rho_i7/(2*pi*e_0)).*((yf - posy(18,7))./((xf - posx(18,7)).^2 + (yf - posy(18,7)).^2));
E_yi188 = (-rho_i8/(2*pi*e_0)).*((yf - posy(18,8))./((xf - posx(18,8)).^2 + (yf - posy(18,8)).^2));
E_yi189 = (-rho_i9/(2*pi*e_0)).*((yf - posy(18,9))./((xf - posx(18,9)).^2 + (yf - posy(18,9)).^2));
E_yi1810 = (-rho_i10/(2*pi*e_0)).*((yf - posy(18,10))./((xf - posx(18,10)).^2 + (yf - posy(18,10)).^2));
E_yi1811 = (-rho_i11/(2*pi*e_0)).*((yf - posy(18,11))./((xf - posx(18,11)).^2 + (yf - posy(18,11)).^2));
E_yi1812 = (-rho_i12/(2*pi*e_0)).*((yf - posy(18,12))./((xf - posx(18,12)).^2 + (yf - posy(18,12)).^2));
E_yi1813 = (rho_i1/(2*pi*e_0)).*((yf - posy(18,13))./((xf - posx(18,13)).^2 + (yf - posy(18,13)).^2));
E_yi1814 = (rho_i2/(2*pi*e_0)).*((yf - posy(18,14))./((xf - posx(18,14)).^2 + (yf - posy(18,14)).^2));
E_yi1815 = (rho_i3/(2*pi*e_0)).*((yf - posy(18,15))./((xf - posx(18,15)).^2 + (yf - posy(18,15)).^2));
E_yi1816 = (rho_i4/(2*pi*e_0)).*((yf - posy(18,16))./((xf - posx(18,16)).^2 + (yf - posy(18,16)).^2));
E_yi1817 = (rho_i5/(2*pi*e_0)).*((yf - posy(18,17))./((xf - posx(18,17)).^2 + (yf - posy(18,17)).^2));
E_yi1819 = (rho_i7/(2*pi*e_0)).*((yf - posy(18,19))./((xf - posx(18,19)).^2 + (yf - posy(18,19)).^2));
E_yi1820 = (rho_i8/(2*pi*e_0)).*((yf - posy(18,20))./((xf - posx(18,20)).^2 + (yf - posy(18,20)).^2));
E_yi1821 = (rho_i9/(2*pi*e_0)).*((yf - posy(18,21))./((xf - posx(18,21)).^2 + (yf - posy(18,21)).^2));
E_yi1822 = (rho_i10/(2*pi*e_0)).*((yf - posy(18,22))./((xf - posx(18,22)).^2 + (yf - posy(18,22)).^2));
E_yi1823 = (rho_i11/(2*pi*e_0)).*((yf - posy(18,23))./((xf - posx(18,23)).^2 + (yf - posy(18,23)).^2));
E_yi1824 = (rho_i12/(2*pi*e_0)).*((yf - posy(18,24))./((xf - posx(18,24)).^2 + (yf - posy(18,24)).^2));

E_yi191 = (-rho_i1/(2*pi*e_0)).*((yf - posy(19,1))./((xf - posx(19,1)).^2 + (yf - posy(19,1)).^2));
E_yi192 = (-rho_i2/(2*pi*e_0)).*((yf - posy(19,2))./((xf - posx(19,2)).^2 + (yf - posy(19,2)).^2));
E_yi193 = (-rho_i3/(2*pi*e_0)).*((yf - posy(19,3))./((xf - posx(19,3)).^2 + (yf - posy(19,3)).^2));
E_yi194 = (-rho_i4/(2*pi*e_0)).*((yf - posy(19,4))./((xf - posx(19,4)).^2 + (yf - posy(19,4)).^2));
E_yi195 = (-rho_i5/(2*pi*e_0)).*((yf - posy(19,5))./((xf - posx(19,5)).^2 + (yf - posy(19,5)).^2));
E_yi196 = (-rho_i6/(2*pi*e_0)).*((yf - posy(19,6))./((xf - posx(19,6)).^2 + (yf - posy(19,6)).^2));
E_yi197 = (-rho_i7/(2*pi*e_0)).*((yf - posy(19,7))./((xf - posx(19,7)).^2 + (yf - posy(19,7)).^2));
E_yi198 = (-rho_i8/(2*pi*e_0)).*((yf - posy(19,8))./((xf - posx(19,8)).^2 + (yf - posy(19,8)).^2));
E_yi199 = (-rho_i9/(2*pi*e_0)).*((yf - posy(19,9))./((xf - posx(19,9)).^2 + (yf - posy(19,9)).^2));
E_yi1910 = (-rho_i10/(2*pi*e_0)).*((yf - posy(19,10))./((xf - posx(19,10)).^2 + (yf - posy(19,10)).^2));
E_yi1911 = (-rho_i11/(2*pi*e_0)).*((yf - posy(19,11))./((xf - posx(19,11)).^2 + (yf - posy(19,11)).^2));
E_yi1912 = (-rho_i12/(2*pi*e_0)).*((yf - posy(19,12))./((xf - posx(19,12)).^2 + (yf - posy(19,12)).^2));
E_yi1913 = (rho_i1/(2*pi*e_0)).*((yf - posy(19,13))./((xf - posx(19,13)).^2 + (yf - posy(19,13)).^2));
E_yi1914 = (rho_i2/(2*pi*e_0)).*((yf - posy(19,14))./((xf - posx(19,14)).^2 + (yf - posy(19,14)).^2));
E_yi1915 = (rho_i3/(2*pi*e_0)).*((yf - posy(19,15))./((xf - posx(19,15)).^2 + (yf - posy(19,15)).^2));
E_yi1916 = (rho_i4/(2*pi*e_0)).*((yf - posy(19,16))./((xf - posx(19,16)).^2 + (yf - posy(19,16)).^2));
E_yi1917 = (rho_i5/(2*pi*e_0)).*((yf - posy(19,17))./((xf - posx(19,17)).^2 + (yf - posy(19,17)).^2));
E_yi1918 = (rho_i6/(2*pi*e_0)).*((yf - posy(19,18))./((xf - posx(19,18)).^2 + (yf - posy(19,18)).^2));
E_yi1920 = (rho_i8/(2*pi*e_0)).*((yf - posy(19,20))./((xf - posx(19,20)).^2 + (yf - posy(19,20)).^2));
E_yi1921 = (rho_i9/(2*pi*e_0)).*((yf - posy(19,21))./((xf - posx(19,21)).^2 + (yf - posy(19,21)).^2));
E_yi1922 = (rho_i10/(2*pi*e_0)).*((yf - posy(19,22))./((xf - posx(19,22)).^2 + (yf - posy(19,22)).^2));
E_yi1923 = (rho_i11/(2*pi*e_0)).*((yf - posy(19,23))./((xf - posx(19,23)).^2 + (yf - posy(19,23)).^2));
E_yi1924 = (rho_i12/(2*pi*e_0)).*((yf - posy(19,24))./((xf - posx(19,24)).^2 + (yf - posy(19,24)).^2));

E_yi201 = (-rho_i1/(2*pi*e_0)).*((yf - posy(20,1))./((xf - posx(20,1)).^2 + (yf - posy(20,1)).^2));
E_yi202 = (-rho_i2/(2*pi*e_0)).*((yf - posy(20,2))./((xf - posx(20,2)).^2 + (yf - posy(20,2)).^2));
E_yi203 = (-rho_i3/(2*pi*e_0)).*((yf - posy(20,3))./((xf - posx(20,3)).^2 + (yf - posy(20,3)).^2));
E_yi204 = (-rho_i4/(2*pi*e_0)).*((yf - posy(20,4))./((xf - posx(20,4)).^2 + (yf - posy(20,4)).^2));
E_yi205 = (-rho_i5/(2*pi*e_0)).*((yf - posy(20,5))./((xf - posx(20,5)).^2 + (yf - posy(20,5)).^2));
E_yi206 = (-rho_i6/(2*pi*e_0)).*((yf - posy(20,6))./((xf - posx(20,6)).^2 + (yf - posy(20,6)).^2));
E_yi207 = (-rho_i7/(2*pi*e_0)).*((yf - posy(20,7))./((xf - posx(20,7)).^2 + (yf - posy(20,7)).^2));
E_yi208 = (-rho_i8/(2*pi*e_0)).*((yf - posy(20,8))./((xf - posx(20,8)).^2 + (yf - posy(20,8)).^2));
E_yi209 = (-rho_i9/(2*pi*e_0)).*((yf - posy(20,9))./((xf - posx(20,9)).^2 + (yf - posy(20,9)).^2));
E_yi2010 = (-rho_i10/(2*pi*e_0)).*((yf - posy(20,10))./((xf - posx(20,10)).^2 + (yf - posy(20,10)).^2));
E_yi2011 = (-rho_i11/(2*pi*e_0)).*((yf - posy(20,11))./((xf - posx(20,11)).^2 + (yf - posy(20,11)).^2));
E_yi2012 = (-rho_i12/(2*pi*e_0)).*((yf - posy(20,12))./((xf - posx(20,12)).^2 + (yf - posy(20,12)).^2));
E_yi2013 = (rho_i1/(2*pi*e_0)).*((yf - posy(20,13))./((xf - posx(20,13)).^2 + (yf - posy(20,13)).^2));
E_yi2014 = (rho_i2/(2*pi*e_0)).*((yf - posy(20,14))./((xf - posx(20,14)).^2 + (yf - posy(20,14)).^2));
E_yi2015 = (rho_i3/(2*pi*e_0)).*((yf - posy(20,15))./((xf - posx(20,15)).^2 + (yf - posy(20,15)).^2));
E_yi2016 = (rho_i4/(2*pi*e_0)).*((yf - posy(20,16))./((xf - posx(20,16)).^2 + (yf - posy(20,16)).^2));
E_yi2017 = (rho_i5/(2*pi*e_0)).*((yf - posy(20,17))./((xf - posx(20,17)).^2 + (yf - posy(20,17)).^2));
E_yi2018 = (rho_i6/(2*pi*e_0)).*((yf - posy(20,18))./((xf - posx(20,18)).^2 + (yf - posy(20,18)).^2));
E_yi2019 = (rho_i7/(2*pi*e_0)).*((yf - posy(20,19))./((xf - posx(20,19)).^2 + (yf - posy(20,19)).^2));
E_yi2021 = (rho_i9/(2*pi*e_0)).*((yf - posy(20,21))./((xf - posx(20,21)).^2 + (yf - posy(20,21)).^2));
E_yi2022 = (rho_i10/(2*pi*e_0)).*((yf - posy(20,22))./((xf - posx(20,22)).^2 + (yf - posy(20,22)).^2));
E_yi2023 = (rho_i11/(2*pi*e_0)).*((yf - posy(20,23))./((xf - posx(20,23)).^2 + (yf - posy(20,23)).^2));
E_yi2024 = (rho_i12/(2*pi*e_0)).*((yf - posy(20,24))./((xf - posx(20,24)).^2 + (yf - posy(20,24)).^2));

E_yi21_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(21,1))./((xf - posx(21,1)).^2 + (yf - posy(21,1)).^2));
E_yi21_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(21,2))./((xf - posx(21,2)).^2 + (yf - posy(21,2)).^2));
E_yi21_3 = (-rho_i3/(2*pi*e_0)).*((yf - posy(21,3))./((xf - posx(21,3)).^2 + (yf - posy(21,3)).^2));
E_yi214 = (-rho_i4/(2*pi*e_0)).*((yf - posy(21,4))./((xf - posx(21,4)).^2 + (yf - posy(21,4)).^2));
E_yi215 = (-rho_i5/(2*pi*e_0)).*((yf - posy(21,5))./((xf - posx(21,5)).^2 + (yf - posy(21,5)).^2));
E_yi216 = (-rho_i6/(2*pi*e_0)).*((yf - posy(21,6))./((xf - posx(21,6)).^2 + (yf - posy(21,6)).^2));
E_yi217 = (-rho_i7/(2*pi*e_0)).*((yf - posy(21,7))./((xf - posx(21,7)).^2 + (yf - posy(21,7)).^2));
E_yi218 = (-rho_i8/(2*pi*e_0)).*((yf - posy(21,8))./((xf - posx(21,8)).^2 + (yf - posy(21,8)).^2));
E_yi219 = (-rho_i9/(2*pi*e_0)).*((yf - posy(21,9))./((xf - posx(21,9)).^2 + (yf - posy(21,9)).^2));
E_yi2110 = (-rho_i10/(2*pi*e_0)).*((yf - posy(21,10))./((xf - posx(21,10)).^2 + (yf - posy(21,10)).^2));
E_yi2111 = (-rho_i11/(2*pi*e_0)).*((yf - posy(21,11))./((xf - posx(21,11)).^2 + (yf - posy(21,11)).^2));
E_yi2112 = (-rho_i12/(2*pi*e_0)).*((yf - posy(21,12))./((xf - posx(21,12)).^2 + (yf - posy(21,12)).^2));
E_yi2113 = (rho_i1/(2*pi*e_0)).*((yf - posy(21,13))./((xf - posx(21,13)).^2 + (yf - posy(21,13)).^2));
E_yi2114 = (rho_i2/(2*pi*e_0)).*((yf - posy(21,14))./((xf - posx(21,14)).^2 + (yf - posy(21,14)).^2));
E_yi2115 = (rho_i3/(2*pi*e_0)).*((yf - posy(21,15))./((xf - posx(21,15)).^2 + (yf - posy(21,15)).^2));
E_yi2116 = (rho_i4/(2*pi*e_0)).*((yf - posy(21,16))./((xf - posx(21,16)).^2 + (yf - posy(21,16)).^2));
E_yi2117 = (rho_i5/(2*pi*e_0)).*((yf - posy(21,17))./((xf - posx(21,17)).^2 + (yf - posy(21,17)).^2));
E_yi2118 = (rho_i6/(2*pi*e_0)).*((yf - posy(21,18))./((xf - posx(21,18)).^2 + (yf - posy(21,18)).^2));
E_yi2119 = (rho_i7/(2*pi*e_0)).*((yf - posy(21,19))./((xf - posx(21,19)).^2 + (yf - posy(21,19)).^2));
E_yi2120 = (rho_i8/(2*pi*e_0)).*((yf - posy(21,20))./((xf - posx(21,20)).^2 + (yf - posy(21,20)).^2));
E_yi2122 = (rho_i10/(2*pi*e_0)).*((yf - posy(21,22))./((xf - posx(21,22)).^2 + (yf - posy(21,22)).^2));
E_yi2123 = (rho_i11/(2*pi*e_0)).*((yf - posy(21,23))./((xf - posx(21,23)).^2 + (yf - posy(21,23)).^2));
E_yi2124 = (rho_i12/(2*pi*e_0)).*((yf - posy(21,24))./((xf - posx(21,24)).^2 + (yf - posy(21,24)).^2));

E_yi22_1 = (-rho_i1/(2*pi*e_0)).*((yf - posy(22,1))./((xf - posx(22,1)).^2 + (yf - posy(22,1)).^2));
E_yi22_2 = (-rho_i2/(2*pi*e_0)).*((yf - posy(22,2))./((xf - posx(22,2)).^2 + (yf - posy(22,2)).^2));
E_yi22_3 = (-rho_i3/(2*pi*e_0)).*((yf - posy(22,3))./((xf - posx(22,3)).^2 + (yf - posy(22,3)).^2));
E_yi22_4 = (-rho_i4/(2*pi*e_0)).*((yf - posy(22,4))./((xf - posx(22,4)).^2 + (yf - posy(22,4)).^2));
E_yi225 = (-rho_i5/(2*pi*e_0)).*((yf - posy(22,5))./((xf - posx(22,5)).^2 + (yf - posy(22,5)).^2));
E_yi226 = (-rho_i6/(2*pi*e_0)).*((yf - posy(22,6))./((xf - posx(22,6)).^2 + (yf - posy(22,6)).^2));
E_yi227 = (-rho_i7/(2*pi*e_0)).*((yf - posy(22,7))./((xf - posx(22,7)).^2 + (yf - posy(22,7)).^2));
E_yi228 = (-rho_i8/(2*pi*e_0)).*((yf - posy(22,8))./((xf - posx(22,8)).^2 + (yf - posy(22,8)).^2));
E_yi229 = (-rho_i9/(2*pi*e_0)).*((yf - posy(22,9))./((xf - posx(22,9)).^2 + (yf - posy(22,9)).^2));
E_yi2210 = (-rho_i10/(2*pi*e_0)).*((yf - posy(22,10))./((xf - posx(22,10)).^2 + (yf - posy(22,10)).^2));
E_yi2211 = (-rho_i11/(2*pi*e_0)).*((yf - posy(22,11))./((xf - posx(22,11)).^2 + (yf - posy(22,11)).^2));
E_yi2212 = (-rho_i12/(2*pi*e_0)).*((yf - posy(22,12))./((xf - posx(22,12)).^2 + (yf - posy(22,12)).^2));
E_yi2213 = (rho_i1/(2*pi*e_0)).*((yf - posy(22,13))./((xf - posx(22,13)).^2 + (yf - posy(22,13)).^2));
E_yi2214 = (rho_i2/(2*pi*e_0)).*((yf - posy(22,14))./((xf - posx(22,14)).^2 + (yf - posy(22,14)).^2));
E_yi2215 = (rho_i3/(2*pi*e_0)).*((yf - posy(22,15))./((xf - posx(22,15)).^2 + (yf - posy(22,15)).^2));
E_yi2216 = (rho_i4/(2*pi*e_0)).*((yf - posy(22,16))./((xf - posx(22,16)).^2 + (yf - posy(22,16)).^2));
E_yi2217 = (rho_i5/(2*pi*e_0)).*((yf - posy(22,17))./((xf - posx(22,17)).^2 + (yf - posy(22,17)).^2));
E_yi2218 = (rho_i6/(2*pi*e_0)).*((yf - posy(22,18))./((xf - posx(22,18)).^2 + (yf - posy(22,18)).^2));
E_yi2219 = (rho_i7/(2*pi*e_0)).*((yf - posy(22,19))./((xf - posx(22,19)).^2 + (yf - posy(22,19)).^2));
E_yi2220 = (rho_i8/(2*pi*e_0)).*((yf - posy(22,20))./((xf - posx(22,20)).^2 + (yf - posy(22,20)).^2));
E_yi2221 = (rho_i9/(2*pi*e_0)).*((yf - posy(22,21))./((xf - posx(22,21)).^2 + (yf - posy(22,21)).^2));
E_yi2223 = (rho_i11/(2*pi*e_0)).*((yf - posy(22,23))./((xf - posx(22,23)).^2 + (yf - posy(22,23)).^2));
E_yi2224 = (rho_i12/(2*pi*e_0)).*((yf - posy(22,24))./((xf - posx(22,24)).^2 + (yf - posy(22,24)).^2));

E_yi231 = (-rho_i1/(2*pi*e_0)).*((yf - posy(23,1))./((xf - posx(23,1)).^2 + (yf - posy(23,1)).^2));
E_yi232 = (-rho_i2/(2*pi*e_0)).*((yf - posy(23,2))./((xf - posx(23,2)).^2 + (yf - posy(23,2)).^2));
E_yi233 = (-rho_i3/(2*pi*e_0)).*((yf - posy(23,3))./((xf - posx(23,3)).^2 + (yf - posy(23,3)).^2));
E_yi234 = (-rho_i4/(2*pi*e_0)).*((yf - posy(23,4))./((xf - posx(23,4)).^2 + (yf - posy(23,4)).^2));
E_yi235 = (-rho_i5/(2*pi*e_0)).*((yf - posy(23,5))./((xf - posx(23,5)).^2 + (yf - posy(23,5)).^2));
E_yi236 = (-rho_i6/(2*pi*e_0)).*((yf - posy(23,6))./((xf - posx(23,6)).^2 + (yf - posy(23,6)).^2));
E_yi237 = (-rho_i7/(2*pi*e_0)).*((yf - posy(23,7))./((xf - posx(23,7)).^2 + (yf - posy(23,7)).^2));
E_yi238 = (-rho_i8/(2*pi*e_0)).*((yf - posy(23,8))./((xf - posx(23,8)).^2 + (yf - posy(23,8)).^2));
E_yi239 = (-rho_i9/(2*pi*e_0)).*((yf - posy(23,9))./((xf - posx(23,9)).^2 + (yf - posy(23,9)).^2));
E_yi2310 = (-rho_i10/(2*pi*e_0)).*((yf - posy(23,10))./((xf - posx(23,10)).^2 + (yf - posy(23,10)).^2));
E_yi2311 = (-rho_i11/(2*pi*e_0)).*((yf - posy(23,11))./((xf - posx(23,11)).^2 + (yf - posy(23,11)).^2));
E_yi2312 = (-rho_i12/(2*pi*e_0)).*((yf - posy(23,12))./((xf - posx(23,12)).^2 + (yf - posy(23,12)).^2));
E_yi2313 = (rho_i1/(2*pi*e_0)).*((yf - posy(23,13))./((xf - posx(23,13)).^2 + (yf - posy(23,13)).^2));
E_yi2314 = (rho_i2/(2*pi*e_0)).*((yf - posy(23,14))./((xf - posx(23,14)).^2 + (yf - posy(23,14)).^2));
E_yi2315 = (rho_i3/(2*pi*e_0)).*((yf - posy(23,15))./((xf - posx(23,15)).^2 + (yf - posy(23,15)).^2));
E_yi2316 = (rho_i4/(2*pi*e_0)).*((yf - posy(23,16))./((xf - posx(23,16)).^2 + (yf - posy(23,16)).^2));
E_yi2317 = (rho_i5/(2*pi*e_0)).*((yf - posy(23,17))./((xf - posx(23,17)).^2 + (yf - posy(23,17)).^2));
E_yi2318 = (rho_i6/(2*pi*e_0)).*((yf - posy(23,18))./((xf - posx(23,18)).^2 + (yf - posy(23,18)).^2));
E_yi2319 = (rho_i7/(2*pi*e_0)).*((yf - posy(23,19))./((xf - posx(23,19)).^2 + (yf - posy(23,19)).^2));
E_yi2320 = (rho_i8/(2*pi*e_0)).*((yf - posy(23,20))./((xf - posx(23,20)).^2 + (yf - posy(23,20)).^2));
E_yi2321 = (rho_i9/(2*pi*e_0)).*((yf - posy(23,21))./((xf - posx(23,21)).^2 + (yf - posy(23,21)).^2));
E_yi2322 = (rho_i10/(2*pi*e_0)).*((yf - posy(23,22))./((xf - posx(23,22)).^2 + (yf - posy(23,22)).^2));
E_yi2324 = (rho_i12/(2*pi*e_0)).*((yf - posy(23,24))./((xf - posx(23,24)).^2 + (yf - posy(23,24)).^2));

E_yi241 = (-rho_i1/(2*pi*e_0)).*((yf - posy(24,1))./((xf - posx(24,1)).^2 + (yf - posy(24,1)).^2));
E_yi242 = (-rho_i2/(2*pi*e_0)).*((yf - posy(24,2))./((xf - posx(24,2)).^2 + (yf - posy(24,2)).^2));
E_yi243 = (-rho_i3/(2*pi*e_0)).*((yf - posy(24,3))./((xf - posx(24,3)).^2 + (yf - posy(24,3)).^2));
E_yi244 = (-rho_i4/(2*pi*e_0)).*((yf - posy(24,4))./((xf - posx(24,4)).^2 + (yf - posy(24,4)).^2));
E_yi245 = (-rho_i5/(2*pi*e_0)).*((yf - posy(24,5))./((xf - posx(24,5)).^2 + (yf - posy(24,5)).^2));
E_yi246 = (-rho_i6/(2*pi*e_0)).*((yf - posy(24,6))./((xf - posx(24,6)).^2 + (yf - posy(24,6)).^2));
E_yi247 = (-rho_i7/(2*pi*e_0)).*((yf - posy(24,7))./((xf - posx(24,7)).^2 + (yf - posy(24,7)).^2));
E_yi248 = (-rho_i8/(2*pi*e_0)).*((yf - posy(24,8))./((xf - posx(24,8)).^2 + (yf - posy(24,8)).^2));
E_yi249 = (-rho_i9/(2*pi*e_0)).*((yf - posy(24,9))./((xf - posx(24,9)).^2 + (yf - posy(24,9)).^2));
E_yi2410 = (-rho_i10/(2*pi*e_0)).*((yf - posy(24,10))./((xf - posx(24,10)).^2 + (yf - posy(24,10)).^2));
E_yi2411 = (-rho_i11/(2*pi*e_0)).*((yf - posy(24,11))./((xf - posx(24,11)).^2 + (yf - posy(24,11)).^2));
E_yi2412 = (-rho_i12/(2*pi*e_0)).*((yf - posy(24,12))./((xf - posx(24,12)).^2 + (yf - posy(24,12)).^2));
E_yi2413 = (rho_i1/(2*pi*e_0)).*((yf - posy(24,13))./((xf - posx(24,13)).^2 + (yf - posy(24,13)).^2));
E_yi2414 = (rho_i2/(2*pi*e_0)).*((yf - posy(24,14))./((xf - posx(24,14)).^2 + (yf - posy(24,14)).^2));
E_yi2415 = (rho_i3/(2*pi*e_0)).*((yf - posy(24,15))./((xf - posx(24,15)).^2 + (yf - posy(24,15)).^2));
E_yi2416 = (rho_i4/(2*pi*e_0)).*((yf - posy(24,16))./((xf - posx(24,16)).^2 + (yf - posy(24,16)).^2));
E_yi2417 = (rho_i5/(2*pi*e_0)).*((yf - posy(24,17))./((xf - posx(24,17)).^2 + (yf - posy(24,17)).^2));
E_yi2418 = (rho_i6/(2*pi*e_0)).*((yf - posy(24,18))./((xf - posx(24,18)).^2 + (yf - posy(24,18)).^2));
E_yi2419 = (rho_i7/(2*pi*e_0)).*((yf - posy(24,19))./((xf - posx(24,19)).^2 + (yf - posy(24,19)).^2));
E_yi2420 = (rho_i8/(2*pi*e_0)).*((yf - posy(24,20))./((xf - posx(24,20)).^2 + (yf - posy(24,20)).^2));
E_yi2421 = (rho_i9/(2*pi*e_0)).*((yf - posy(24,21))./((xf - posx(24,21)).^2 + (yf - posy(24,21)).^2));
E_yi2422 = (rho_i10/(2*pi*e_0)).*((yf - posy(24,22))./((xf - posx(24,22)).^2 + (yf - posy(24,22)).^2));
E_yi2423 = (rho_i11/(2*pi*e_0)).*((yf - posy(24,23))./((xf - posx(24,23)).^2 + (yf - posy(24,23)).^2));

%% aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr12 ;	E_xr13 	;	E_xr14 	;	E_xr15 	;	E_xr16 	;	E_xr17 	;	E_xr18 	;	E_xr19 	;	E_xr110 	;	E_xr1_11 	;	E_xr1_12 	;	E_xr1_13 	;	E_xr1_14 	;	E_xr1_15 	;	E_xr1_16 	;	E_xr1_17 	;	E_xr1_18 	;	E_xr1_19 	;	E_xr1_20 	;	E_xr1_21 	;	E_xr1_22 	;	E_xr1_23 	;	E_xr1_24 	;	E_xr21 	;	E_xr23 	;	E_xr24 	;	E_xr25 	;	E_xr26 	;	E_xr27 	;	E_xr28 	;	E_xr29 	;	E_xr210 	;	E_xr2_11 	;	E_xr2_12 	;	E_xr2_13 	;	E_xr2_14 	;	E_xr2_15 	;	E_xr2_16 	;	E_xr2_17 	;	E_xr2_18 	;	E_xr2_19 	;	E_xr2_20 	;	E_xr2_21 	;	E_xr2_22 	;	E_xr2_23 	;	E_xr2_24 	;	E_xr31 	;	E_xr32 	;	E_xr34 	;	E_xr35 	;	E_xr36 	;	E_xr37 	;	E_xr38 	;	E_xr39 	;	E_xr310 	;	E_xr311 	;	E_xr312 	;	E_xr313 	;	E_xr314 	;	E_xr315 	;	E_xr316 	;	E_xr317 	;	E_xr318 	;	E_xr319 	;	E_xr320 	;	E_xr321 	;	E_xr322 	;	E_xr323 	;	E_xr324 	;	E_xr41 	;	E_xr42 	;	E_xr43 	;	E_xr45 	;	E_xr46 	;	E_xr47 	;	E_xr48 	;	E_xr49 	;	E_xr410 	;	E_xr411 	;	E_xr412 	;	E_xr413 	;	E_xr414 	;	E_xr415 	;	E_xr416 	;	E_xr417 	;	E_xr418 	;	E_xr419 	;	E_xr420 	;	E_xr421 	;	E_xr422 	;	E_xr423 	;	E_xr424 	;	E_xr51 	;	E_xr52 	;	E_xr53 	;	E_xr54 	;	E_xr56 	;	E_xr57 	;	E_xr58 	;	E_xr59 	;	E_xr510 	;	E_xr511 	;	E_xr512 	;	E_xr513 	;	E_xr514 	;	E_xr515 	;	E_xr516 	;	E_xr517 	;	E_xr518 	;	E_xr519 	;	E_xr520 	;	E_xr521 	;	E_xr522 	;	E_xr523 	;	E_xr524 	;	E_xr61 	;	E_xr62 	;	E_xr63 	;	E_xr64 	;	E_xr65 	;	E_xr67 	;	E_xr68 	;	E_xr69 	;	E_xr610 	;	E_xr611 	;	E_xr612 	;	E_xr613 	;	E_xr614 	;	E_xr615 	;	E_xr616 	;	E_xr617 	;	E_xr618 	;	E_xr619 	;	E_xr620 	;	E_xr621 	;	E_xr622 	;	E_xr623 	;	E_xr624 	;	E_xr71 	;	E_xr72 	;	E_xr73 	;	E_xr74 	;	E_xr75 	;	E_xr76 	;	E_xr78 	;	E_xr79 	;	E_xr710 	;	E_xr711 	;	E_xr712 	;	E_xr713 	;	E_xr714 	;	E_xr715 	;	E_xr716 	;	E_xr717 	;	E_xr718 	;	E_xr719 	;	E_xr720 	;	E_xr721 	;	E_xr722 	;	E_xr723 	;	E_xr724 	;	E_xr81 	;	E_xr82 	;	E_xr83 	;	E_xr84 	;	E_xr85 	;	E_xr86 	;	E_xr87 	;	E_xr89 	;	E_xr810 	;	E_xr811 	;	E_xr812 	;	E_xr813 	;	E_xr814 	;	E_xr815 	;	E_xr816 	;	E_xr817 	;	E_xr818 	;	E_xr819 	;	E_xr820 	;	E_xr821 	;	E_xr822 	;	E_xr823 	;	E_xr824 	;	E_xr91 	;	E_xr92 	;	E_xr93 	;	E_xr94 	;	E_xr95 	;	E_xr96 	;	E_xr97 	;	E_xr98 	;	E_xr910 	;	E_xr911 	;	E_xr912 	;	E_xr913 	;	E_xr914 	;	E_xr915 	;	E_xr916 	;	E_xr917 	;	E_xr918 	;	E_xr919 	;	E_xr920 	;	E_xr921 	;	E_xr922 	;	E_xr923 	;	E_xr924 	;	E_xr101 	;	E_xr102 	;	E_xr103 	;	E_xr104 	;	E_xr105 	;	E_xr106 	;	E_xr107 	;	E_xr108 	;	E_xr109 	;	E_xr1011 	;	E_xr1012 	;	E_xr1013 	;	E_xr1014 	;	E_xr1015 	;	E_xr1016 	;	E_xr1017 	;	E_xr1018 	;	E_xr10_19 	;	E_xr10_20 	;	E_xr10_21 	;	E_xr10_22 	;	E_xr10_23 	;	E_xr10_24 	;	E_xr11_1 	;	E_xr11_2 	;	E_xr113 	;	E_xr114 	;	E_xr115 	;	E_xr116 	;	E_xr117 	;	E_xr118 	;	E_xr119 	;	E_xr1110 	;	E_xr1112 	;	E_xr1113 	;	E_xr1114 	;	E_xr1115 	;	E_xr1116 	;	E_xr1117 	;	E_xr1118 	;	E_xr1019 	;	E_xr1020 	;	E_xr1021 	;	E_xr1022 	;	E_xr1023 	;	E_xr1024 	;	E_xr12_1 	;	E_xr12_2 	;	E_xr123 	;	E_xr124 	;	E_xr125 	;	E_xr126 	;	E_xr127 	;	E_xr128 	;	E_xr129 	;	E_xr1210 	;	E_xr1211 	;	E_xr1213 	;	E_xr1214 	;	E_xr1215 	;	E_xr1216 	;	E_xr1217 	;	E_xr1218 	;	E_xr1219 	;	E_xr1220 	;	E_xr1221 	;	E_xr1222 	;	E_xr1223 	;	E_xr1224 	;	E_xr131 	;	E_xr132 	;	E_xr133 	;	E_xr134 	;	E_xr135 	;	E_xr136 	;	E_xr137 	;	E_xr138 	;	E_xr139 	;	E_xr1310 	;	E_xr1311 	;	E_xr1312 	;	E_xr1314 	;	E_xr1315 	;	E_xr1316 	;	E_xr1317 	;	E_xr1318 	;	E_xr1319 	;	E_xr1320 	;	E_xr1321 	;	E_xr1322 	;	E_xr1323 	;	E_xr1324 	;	E_xr141 	;	E_xr142 	;	E_xr143 	;	E_xr144 	;	E_xr145 	;	E_xr146 	;	E_xr147 	;	E_xr148 	;	E_xr149 	;	E_xr1410 	;	E_xr1411 	;	E_xr1412 	;	E_xr1413 	;	E_xr1415 	;	E_xr1416 	;	E_xr1417 	;	E_xr1418 	;	E_xr1419 	;	E_xr1420 	;	E_xr1421 	;	E_xr1422 	;	E_xr1423 	;	E_xr1424 	;	E_xr151 	;	E_xr152 	;	E_xr153 	;	E_xr154 	;	E_xr155 	;	E_xr156 	;	E_xr157 	;	E_xr158 	;	E_xr159 	;	E_xr1510 	;	E_xr1511 	;	E_xr1512 	;	E_xr1513 	;	E_xr1514 	;	E_xr1516 	;	E_xr1517 	;	E_xr1518 	;	E_xr1519 	;	E_xr1520 	;	E_xr1521 	;	E_xr1522 	;	E_xr1523 	;	E_xr1524 	;	E_xr161 	;	E_xr162 	;	E_xr163 	;	E_xr164 	;	E_xr165 	;	E_xr166 	;	E_xr167 	;	E_xr168 	;	E_xr169 	;	E_xr1610 	;	E_xr1611 	;	E_xr1612 	;	E_xr1613 	;	E_xr1614 	;	E_xr1615 	;	E_xr1617 	;	E_xr1618 	;	E_xr1619 	;	E_xr1620 	;	E_xr1621 	;	E_xr1622 	;	E_xr1623 	;	E_xr1624 	;	E_xr171 	;	E_xr172 	;	E_xr173 	;	E_xr174 	;	E_xr175 	;	E_xr176 	;	E_xr177 	;	E_xr178 	;	E_xr179 	;	E_xr1710 	;	E_xr1711 	;	E_xr1712 	;	E_xr1713 	;	E_xr1714 	;	E_xr1715 	;	E_xr1716 	;	E_xr1718 	;	E_xr1719 	;	E_xr1720 	;	E_xr1721 	;	E_xr1722 	;	E_xr1723 	;	E_xr1724 	;	E_xr181 	;	E_xr182 	;	E_xr183 	;	E_xr184 	;	E_xr185 	;	E_xr186 	;	E_xr187 	;	E_xr188 	;	E_xr189 	;	E_xr1810 	;	E_xr1811 	;	E_xr1812 	;	E_xr1813 	;	E_xr1814 	;	E_xr1815 	;	E_xr1816 	;	E_xr1817 	;	E_xr1819 	;	E_xr1820 	;	E_xr1821 	;	E_xr1822 	;	E_xr1823 	;	E_xr1824 	;	E_xr191 	;	E_xr192 	;	E_xr193 	;	E_xr194 	;	E_xr195 	;	E_xr196 	;	E_xr197 	;	E_xr198 	;	E_xr199 	;	E_xr1910 	;	E_xr1911 	;	E_xr1912 	;	E_xr1913 	;	E_xr1914 	;	E_xr1915 	;	E_xr1916 	;	E_xr1917 	;	E_xr1918 	;	E_xr1920 	;	E_xr1921 	;	E_xr1922 	;	E_xr1923 	;	E_xr1924 	;	E_xr201 	;	E_xr202 	;	E_xr203 	;	E_xr204 	;	E_xr205 	;	E_xr206 	;	E_xr207 	;	E_xr208 	;	E_xr209 	;	E_xr2010 	;	E_xr2011 	;	E_xr2012 	;	E_xr2013 	;	E_xr2014 	;	E_xr2015 	;	E_xr2016 	;	E_xr2017 	;	E_xr2018 	;	E_xr2019 	;	E_xr2021 	;	E_xr2022 	;	E_xr2023 	;	E_xr2024 	;	E_xr21_1 	;	E_xr21_2 	;	E_xr21_3 	;	E_xr214 	;	E_xr215 	;	E_xr216 	;	E_xr217 	;	E_xr218 	;	E_xr219 	;	E_xr2110 	;	E_xr2111 	;	E_xr2112 	;	E_xr2113 	;	E_xr2114 	;	E_xr2115 	;	E_xr2116 	;	E_xr2117 	;	E_xr2118 	;	E_xr2119 	;	E_xr2120 	;	E_xr2122 	;	E_xr2123 	;	E_xr2124 	;	E_xr22_1 	;	E_xr22_2 	;	E_xr22_3 	;	E_xr22_4 	;	E_xr225 	;	E_xr226 	;	E_xr227 	;	E_xr228 	;	E_xr229 	;	E_xr2210 	;	E_xr2211 	;	E_xr2212 	;	E_xr2213 	;	E_xr2214 	;	E_xr2215 	;	E_xr2216 	;	E_xr2217 	;	E_xr2218 	;	E_xr2219 	;	E_xr2220 	;	E_xr2221 	;	E_xr2223 	;	E_xr2224 	;	E_xr231 	;	E_xr232 	;	E_xr233 	;	E_xr234 	;	E_xr235 	;	E_xr236 	;	E_xr237 	;	E_xr238 	;	E_xr239 	;	E_xr2310 	;	E_xr2311 	;	E_xr2312 	;	E_xr2313 	;	E_xr2314 	;	E_xr2315 	;	E_xr2316 	;	E_xr2317 	;	E_xr2318 	;	E_xr2319 	;	E_xr2320 	;	E_xr2321 	;	E_xr2322 	;	E_xr2324 	;	E_xr241 	;	E_xr242 	;	E_xr243 	;	E_xr244 	;	E_xr245 	;	E_xr246 	;	E_xr247 	;	E_xr248 	;	E_xr249 	;	E_xr2410 	;	E_xr2411 	;	E_xr2412 	;	E_xr2413 	;	E_xr2414 	;	E_xr2415 	;	E_xr2416 	;	E_xr2417 	;	E_xr2418 	;	E_xr2419 	;	E_xr2420 	;	E_xr2421 	;	E_xr2422 	;	E_xr2423  ];
    
E_xi = [ E_xi12 ;	E_xi13 	;	E_xi14 	;	E_xi15 	;	E_xi16 	;	E_xi17 	;	E_xi18 	;	E_xi19 	;	E_xi110 	;	E_xi1_11 	;	E_xi1_12 	;	E_xi1_13 	;	E_xi1_14 	;	E_xi1_15 	;	E_xi1_16 	;	E_xi1_17 	;	E_xi1_18 	;	E_xi1_19 	;	E_xi1_20 	;	E_xi1_21 	;	E_xi1_22 	;	E_xi1_23 	;	E_xi1_24 	;	E_xi21 	;	E_xi23 	;	E_xi24 	;	E_xi25 	;	E_xi26 	;	E_xi27 	;	E_xi28 	;	E_xi29 	;	E_xi210 	;	E_xi2_11 	;	E_xi2_12 	;	E_xi2_13 	;	E_xi2_14 	;	E_xi2_15 	;	E_xi2_16 	;	E_xi2_17 	;	E_xi2_18 	;	E_xi2_19 	;	E_xi2_20 	;	E_xi2_21 	;	E_xi2_22 	;	E_xi2_23 	;	E_xi2_24 	;	E_xi31 	;	E_xi32 	;	E_xi34 	;	E_xi35 	;	E_xi36 	;	E_xi37 	;	E_xi38 	;	E_xi39 	;	E_xi310 	;	E_xi311 	;	E_xi312 	;	E_xi313 	;	E_xi314 	;	E_xi315 	;	E_xi316 	;	E_xi317 	;	E_xi318 	;	E_xi319 	;	E_xi320 	;	E_xi321 	;	E_xi322 	;	E_xi323 	;	E_xi324 	;	E_xi41 	;	E_xi42 	;	E_xi43 	;	E_xi45 	;	E_xi46 	;	E_xi47 	;	E_xi48 	;	E_xi49 	;	E_xi410 	;	E_xi411 	;	E_xi412 	;	E_xi413 	;	E_xi414 	;	E_xi415 	;	E_xi416 	;	E_xi417 	;	E_xi418 	;	E_xi419 	;	E_xi420 	;	E_xi421 	;	E_xi422 	;	E_xi423 	;	E_xi424 	;	E_xi51 	;	E_xi52 	;	E_xi53 	;	E_xi54 	;	E_xi56 	;	E_xi57 	;	E_xi58 	;	E_xi59 	;	E_xi510 	;	E_xi511 	;	E_xi512 	;	E_xi513 	;	E_xi514 	;	E_xi515 	;	E_xi516 	;	E_xi517 	;	E_xi518 	;	E_xi519 	;	E_xi520 	;	E_xi521 	;	E_xi522 	;	E_xi523 	;	E_xi524 	;	E_xi61 	;	E_xi62 	;	E_xi63 	;	E_xi64 	;	E_xi65 	;	E_xi67 	;	E_xi68 	;	E_xi69 	;	E_xi610 	;	E_xi611 	;	E_xi612 	;	E_xi613 	;	E_xi614 	;	E_xi615 	;	E_xi616 	;	E_xi617 	;	E_xi618 	;	E_xi619 	;	E_xi620 	;	E_xi621 	;	E_xi622 	;	E_xi623 	;	E_xi624 	;	E_xi71 	;	E_xi72 	;	E_xi73 	;	E_xi74 	;	E_xi75 	;	E_xi76 	;	E_xi78 	;	E_xi79 	;	E_xi710 	;	E_xi711 	;	E_xi712 	;	E_xi713 	;	E_xi714 	;	E_xi715 	;	E_xi716 	;	E_xi717 	;	E_xi718 	;	E_xi719 	;	E_xi720 	;	E_xi721 	;	E_xi722 	;	E_xi723 	;	E_xi724 	;	E_xi81 	;	E_xi82 	;	E_xi83 	;	E_xi84 	;	E_xi85 	;	E_xi86 	;	E_xi87 	;	E_xi89 	;	E_xi810 	;	E_xi811 	;	E_xi812 	;	E_xi813 	;	E_xi814 	;	E_xi815 	;	E_xi816 	;	E_xi817 	;	E_xi818 	;	E_xi819 	;	E_xi820 	;	E_xi821 	;	E_xi822 	;	E_xi823 	;	E_xi824 	;	E_xi91 	;	E_xi92 	;	E_xi93 	;	E_xi94 	;	E_xi95 	;	E_xi96 	;	E_xi97 	;	E_xi98 	;	E_xi910 	;	E_xi911 	;	E_xi912 	;	E_xi913 	;	E_xi914 	;	E_xi915 	;	E_xi916 	;	E_xi917 	;	E_xi918 	;	E_xi919 	;	E_xi920 	;	E_xi921 	;	E_xi922 	;	E_xi923 	;	E_xi924 	;	E_xi101 	;	E_xi102 	;	E_xi103 	;	E_xi104 	;	E_xi105 	;	E_xi106 	;	E_xi107 	;	E_xi108 	;	E_xi109 	;	E_xi1011 	;	E_xi1012 	;	E_xi1013 	;	E_xi1014 	;	E_xi1015 	;	E_xi1016 	;	E_xi1017 	;	E_xi1018 	;	E_xi10_19 	;	E_xi10_20 	;	E_xi10_21 	;	E_xi10_22 	;	E_xi10_23 	;	E_xi10_24 	;	E_xi11_1 	;	E_xi11_2 	;	E_xi113 	;	E_xi114 	;	E_xi115 	;	E_xi116 	;	E_xi117 	;	E_xi118 	;	E_xi119 	;	E_xi1110 	;	E_xi1112 	;	E_xi1113 	;	E_xi1114 	;	E_xi1115 	;	E_xi1116 	;	E_xi1117 	;	E_xi1118 	;	E_xi1019 	;	E_xi1020 	;	E_xi1021 	;	E_xi1022 	;	E_xi1023 	;	E_xi1024 	;	E_xi12_1 	;	E_xi12_2 	;	E_xi123 	;	E_xi124 	;	E_xi125 	;	E_xi126 	;	E_xi127 	;	E_xi128 	;	E_xi129 	;	E_xi1210 	;	E_xi1211 	;	E_xi1213 	;	E_xi1214 	;	E_xi1215 	;	E_xi1216 	;	E_xi1217 	;	E_xi1218 	;	E_xi1219 	;	E_xi1220 	;	E_xi1221 	;	E_xi1222 	;	E_xi1223 	;	E_xi1224 	;	E_xi131 	;	E_xi132 	;	E_xi133 	;	E_xi134 	;	E_xi135 	;	E_xi136 	;	E_xi137 	;	E_xi138 	;	E_xi139 	;	E_xi1310 	;	E_xi1311 	;	E_xi1312 	;	E_xi1314 	;	E_xi1315 	;	E_xi1316 	;	E_xi1317 	;	E_xi1318 	;	E_xi1319 	;	E_xi1320 	;	E_xi1321 	;	E_xi1322 	;	E_xi1323 	;	E_xi1324 	;	E_xi141 	;	E_xi142 	;	E_xi143 	;	E_xi144 	;	E_xi145 	;	E_xi146 	;	E_xi147 	;	E_xi148 	;	E_xi149 	;	E_xi1410 	;	E_xi1411 	;	E_xi1412 	;	E_xi1413 	;	E_xi1415 	;	E_xi1416 	;	E_xi1417 	;	E_xi1418 	;	E_xi1419 	;	E_xi1420 	;	E_xi1421 	;	E_xi1422 	;	E_xi1423 	;	E_xi1424 	;	E_xi151 	;	E_xi152 	;	E_xi153 	;	E_xi154 	;	E_xi155 	;	E_xi156 	;	E_xi157 	;	E_xi158 	;	E_xi159 	;	E_xi1510 	;	E_xi1511 	;	E_xi1512 	;	E_xi1513 	;	E_xi1514 	;	E_xi1516 	;	E_xi1517 	;	E_xi1518 	;	E_xi1519 	;	E_xi1520 	;	E_xi1521 	;	E_xi1522 	;	E_xi1523 	;	E_xi1524 	;	E_xi161 	;	E_xi162 	;	E_xi163 	;	E_xi164 	;	E_xi165 	;	E_xi166 	;	E_xi167 	;	E_xi168 	;	E_xi169 	;	E_xi1610 	;	E_xi1611 	;	E_xi1612 	;	E_xi1613 	;	E_xi1614 	;	E_xi1615 	;	E_xi1617 	;	E_xi1618 	;	E_xi1619 	;	E_xi1620 	;	E_xi1621 	;	E_xi1622 	;	E_xi1623 	;	E_xi1624 	;	E_xi171 	;	E_xi172 	;	E_xi173 	;	E_xi174 	;	E_xi175 	;	E_xi176 	;	E_xi177 	;	E_xi178 	;	E_xi179 	;	E_xi1710 	;	E_xi1711 	;	E_xi1712 	;	E_xi1713 	;	E_xi1714 	;	E_xi1715 	;	E_xi1716 	;	E_xi1718 	;	E_xi1719 	;	E_xi1720 	;	E_xi1721 	;	E_xi1722 	;	E_xi1723 	;	E_xi1724 	;	E_xi181 	;	E_xi182 	;	E_xi183 	;	E_xi184 	;	E_xi185 	;	E_xi186 	;	E_xi187 	;	E_xi188 	;	E_xi189 	;	E_xi1810 	;	E_xi1811 	;	E_xi1812 	;	E_xi1813 	;	E_xi1814 	;	E_xi1815 	;	E_xi1816 	;	E_xi1817 	;	E_xi1819 	;	E_xi1820 	;	E_xi1821 	;	E_xi1822 	;	E_xi1823 	;	E_xi1824 	;	E_xi191 	;	E_xi192 	;	E_xi193 	;	E_xi194 	;	E_xi195 	;	E_xi196 	;	E_xi197 	;	E_xi198 	;	E_xi199 	;	E_xi1910 	;	E_xi1911 	;	E_xi1912 	;	E_xi1913 	;	E_xi1914 	;	E_xi1915 	;	E_xi1916 	;	E_xi1917 	;	E_xi1918 	;	E_xi1920 	;	E_xi1921 	;	E_xi1922 	;	E_xi1923 	;	E_xi1924 	;	E_xi201 	;	E_xi202 	;	E_xi203 	;	E_xi204 	;	E_xi205 	;	E_xi206 	;	E_xi207 	;	E_xi208 	;	E_xi209 	;	E_xi2010 	;	E_xi2011 	;	E_xi2012 	;	E_xi2013 	;	E_xi2014 	;	E_xi2015 	;	E_xi2016 	;	E_xi2017 	;	E_xi2018 	;	E_xi2019 	;	E_xi2021 	;	E_xi2022 	;	E_xi2023 	;	E_xi2024 	;	E_xi21_1 	;	E_xi21_2 	;	E_xi21_3 	;	E_xi214 	;	E_xi215 	;	E_xi216 	;	E_xi217 	;	E_xi218 	;	E_xi219 	;	E_xi2110 	;	E_xi2111 	;	E_xi2112 	;	E_xi2113 	;	E_xi2114 	;	E_xi2115 	;	E_xi2116 	;	E_xi2117 	;	E_xi2118 	;	E_xi2119 	;	E_xi2120 	;	E_xi2122 	;	E_xi2123 	;	E_xi2124 	;	E_xi22_1 	;	E_xi22_2 	;	E_xi22_3 	;	E_xi22_4 	;	E_xi225 	;	E_xi226 	;	E_xi227 	;	E_xi228 	;	E_xi229 	;	E_xi2210 	;	E_xi2211 	;	E_xi2212 	;	E_xi2213 	;	E_xi2214 	;	E_xi2215 	;	E_xi2216 	;	E_xi2217 	;	E_xi2218 	;	E_xi2219 	;	E_xi2220 	;	E_xi2221 	;	E_xi2223 	;	E_xi2224 	;	E_xi231 	;	E_xi232 	;	E_xi233 	;	E_xi234 	;	E_xi235 	;	E_xi236 	;	E_xi237 	;	E_xi238 	;	E_xi239 	;	E_xi2310 	;	E_xi2311 	;	E_xi2312 	;	E_xi2313 	;	E_xi2314 	;	E_xi2315 	;	E_xi2316 	;	E_xi2317 	;	E_xi2318 	;	E_xi2319 	;	E_xi2320 	;	E_xi2321 	;	E_xi2322 	;	E_xi2324 	;	E_xi241 	;	E_xi242 	;	E_xi243 	;	E_xi244 	;	E_xi245 	;	E_xi246 	;	E_xi247 	;	E_xi248 	;	E_xi249 	;	E_xi2410 	;	E_xi2411 	;	E_xi2412 	;	E_xi2413 	;	E_xi2414 	;	E_xi2415 	;	E_xi2416 	;	E_xi2417 	;	E_xi2418 	;	E_xi2419 	;	E_xi2420 	;	E_xi2421 	;	E_xi2422 	;	E_xi2423  ];

E_yr = [ E_yr12 ;	E_yr13 	;	E_yr14 	;	E_yr15 	;	E_yr16 	;	E_yr17 	;	E_yr18 	;	E_yr19 	;	E_yr110 	;	E_yr1_11 	;	E_yr1_12 	;	E_yr1_13 	;	E_yr1_14 	;	E_yr1_15 	;	E_yr1_16 	;	E_yr1_17 	;	E_yr1_18 	;	E_yr1_19 	;	E_yr1_20 	;	E_yr1_21 	;	E_yr1_22 	;	E_yr1_23 	;	E_yr1_24 	;	E_yr21 	;	E_yr23 	;	E_yr24 	;	E_yr25 	;	E_yr26 	;	E_yr27 	;	E_yr28 	;	E_yr29 	;	E_yr210 	;	E_yr2_11 	;	E_yr2_12 	;	E_yr2_13 	;	E_yr2_14 	;	E_yr2_15 	;	E_yr2_16 	;	E_yr2_17 	;	E_yr2_18 	;	E_yr2_19 	;	E_yr2_20 	;	E_yr2_21 	;	E_yr2_22 	;	E_yr2_23 	;	E_yr2_24 	;	E_yr31 	;	E_yr32 	;	E_yr34 	;	E_yr35 	;	E_yr36 	;	E_yr37 	;	E_yr38 	;	E_yr39 	;	E_yr310 	;	E_yr311 	;	E_yr312 	;	E_yr313 	;	E_yr314 	;	E_yr315 	;	E_yr316 	;	E_yr317 	;	E_yr318 	;	E_yr319 	;	E_yr320 	;	E_yr321 	;	E_yr322 	;	E_yr323 	;	E_yr324 	;	E_yr41 	;	E_yr42 	;	E_yr43 	;	E_yr45 	;	E_yr46 	;	E_yr47 	;	E_yr48 	;	E_yr49 	;	E_yr410 	;	E_yr411 	;	E_yr412 	;	E_yr413 	;	E_yr414 	;	E_yr415 	;	E_yr416 	;	E_yr417 	;	E_yr418 	;	E_yr419 	;	E_yr420 	;	E_yr421 	;	E_yr422 	;	E_yr423 	;	E_yr424 	;	E_yr51 	;	E_yr52 	;	E_yr53 	;	E_yr54 	;	E_yr56 	;	E_yr57 	;	E_yr58 	;	E_yr59 	;	E_yr510 	;	E_yr511 	;	E_yr512 	;	E_yr513 	;	E_yr514 	;	E_yr515 	;	E_yr516 	;	E_yr517 	;	E_yr518 	;	E_yr519 	;	E_yr520 	;	E_yr521 	;	E_yr522 	;	E_yr523 	;	E_yr524 	;	E_yr61 	;	E_yr62 	;	E_yr63 	;	E_yr64 	;	E_yr65 	;	E_yr67 	;	E_yr68 	;	E_yr69 	;	E_yr610 	;	E_yr611 	;	E_yr612 	;	E_yr613 	;	E_yr614 	;	E_yr615 	;	E_yr616 	;	E_yr617 	;	E_yr618 	;	E_yr619 	;	E_yr620 	;	E_yr621 	;	E_yr622 	;	E_yr623 	;	E_yr624 	;	E_yr71 	;	E_yr72 	;	E_yr73 	;	E_yr74 	;	E_yr75 	;	E_yr76 	;	E_yr78 	;	E_yr79 	;	E_yr710 	;	E_yr711 	;	E_yr712 	;	E_yr713 	;	E_yr714 	;	E_yr715 	;	E_yr716 	;	E_yr717 	;	E_yr718 	;	E_yr719 	;	E_yr720 	;	E_yr721 	;	E_yr722 	;	E_yr723 	;	E_yr724 	;	E_yr81 	;	E_yr82 	;	E_yr83 	;	E_yr84 	;	E_yr85 	;	E_yr86 	;	E_yr87 	;	E_yr89 	;	E_yr810 	;	E_yr811 	;	E_yr812 	;	E_yr813 	;	E_yr814 	;	E_yr815 	;	E_yr816 	;	E_yr817 	;	E_yr818 	;	E_yr819 	;	E_yr820 	;	E_yr821 	;	E_yr822 	;	E_yr823 	;	E_yr824 	;	E_yr91 	;	E_yr92 	;	E_yr93 	;	E_yr94 	;	E_yr95 	;	E_yr96 	;	E_yr97 	;	E_yr98 	;	E_yr910 	;	E_yr911 	;	E_yr912 	;	E_yr913 	;	E_yr914 	;	E_yr915 	;	E_yr916 	;	E_yr917 	;	E_yr918 	;	E_yr919 	;	E_yr920 	;	E_yr921 	;	E_yr922 	;	E_yr923 	;	E_yr924 	;	E_yr101 	;	E_yr102 	;	E_yr103 	;	E_yr104 	;	E_yr105 	;	E_yr106 	;	E_yr107 	;	E_yr108 	;	E_yr109 	;	E_yr1011 	;	E_yr1012 	;	E_yr1013 	;	E_yr1014 	;	E_yr1015 	;	E_yr1016 	;	E_yr1017 	;	E_yr1018 	;	E_yr10_19 	;	E_yr10_20 	;	E_yr10_21 	;	E_yr10_22 	;	E_yr10_23 	;	E_yr10_24 	;	E_yr11_1 	;	E_yr11_2 	;	E_yr113 	;	E_yr114 	;	E_yr115 	;	E_yr116 	;	E_yr117 	;	E_yr118 	;	E_yr119 	;	E_yr1110 	;	E_yr1112 	;	E_yr1113 	;	E_yr1114 	;	E_yr1115 	;	E_yr1116 	;	E_yr1117 	;	E_yr1118 	;	E_yr1019 	;	E_yr1020 	;	E_yr1021 	;	E_yr1022 	;	E_yr1023 	;	E_yr1024 	;	E_yr12_1 	;	E_yr12_2 	;	E_yr123 	;	E_yr124 	;	E_yr125 	;	E_yr126 	;	E_yr127 	;	E_yr128 	;	E_yr129 	;	E_yr1210 	;	E_yr1211 	;	E_yr1213 	;	E_yr1214 	;	E_yr1215 	;	E_yr1216 	;	E_yr1217 	;	E_yr1218 	;	E_yr1219 	;	E_yr1220 	;	E_yr1221 	;	E_yr1222 	;	E_yr1223 	;	E_yr1224 	;	E_yr131 	;	E_yr132 	;	E_yr133 	;	E_yr134 	;	E_yr135 	;	E_yr136 	;	E_yr137 	;	E_yr138 	;	E_yr139 	;	E_yr1310 	;	E_yr1311 	;	E_yr1312 	;	E_yr1314 	;	E_yr1315 	;	E_yr1316 	;	E_yr1317 	;	E_yr1318 	;	E_yr1319 	;	E_yr1320 	;	E_yr1321 	;	E_yr1322 	;	E_yr1323 	;	E_yr1324 	;	E_yr141 	;	E_yr142 	;	E_yr143 	;	E_yr144 	;	E_yr145 	;	E_yr146 	;	E_yr147 	;	E_yr148 	;	E_yr149 	;	E_yr1410 	;	E_yr1411 	;	E_yr1412 	;	E_yr1413 	;	E_yr1415 	;	E_yr1416 	;	E_yr1417 	;	E_yr1418 	;	E_yr1419 	;	E_yr1420 	;	E_yr1421 	;	E_yr1422 	;	E_yr1423 	;	E_yr1424 	;	E_yr151 	;	E_yr152 	;	E_yr153 	;	E_yr154 	;	E_yr155 	;	E_yr156 	;	E_yr157 	;	E_yr158 	;	E_yr159 	;	E_yr1510 	;	E_yr1511 	;	E_yr1512 	;	E_yr1513 	;	E_yr1514 	;	E_yr1516 	;	E_yr1517 	;	E_yr1518 	;	E_yr1519 	;	E_yr1520 	;	E_yr1521 	;	E_yr1522 	;	E_yr1523 	;	E_yr1524 	;	E_yr161 	;	E_yr162 	;	E_yr163 	;	E_yr164 	;	E_yr165 	;	E_yr166 	;	E_yr167 	;	E_yr168 	;	E_yr169 	;	E_yr1610 	;	E_yr1611 	;	E_yr1612 	;	E_yr1613 	;	E_yr1614 	;	E_yr1615 	;	E_yr1617 	;	E_yr1618 	;	E_yr1619 	;	E_yr1620 	;	E_yr1621 	;	E_yr1622 	;	E_yr1623 	;	E_yr1624 	;	E_yr171 	;	E_yr172 	;	E_yr173 	;	E_yr174 	;	E_yr175 	;	E_yr176 	;	E_yr177 	;	E_yr178 	;	E_yr179 	;	E_yr1710 	;	E_yr1711 	;	E_yr1712 	;	E_yr1713 	;	E_yr1714 	;	E_yr1715 	;	E_yr1716 	;	E_yr1718 	;	E_yr1719 	;	E_yr1720 	;	E_yr1721 	;	E_yr1722 	;	E_yr1723 	;	E_yr1724 	;	E_yr181 	;	E_yr182 	;	E_yr183 	;	E_yr184 	;	E_yr185 	;	E_yr186 	;	E_yr187 	;	E_yr188 	;	E_yr189 	;	E_yr1810 	;	E_yr1811 	;	E_yr1812 	;	E_yr1813 	;	E_yr1814 	;	E_yr1815 	;	E_yr1816 	;	E_yr1817 	;	E_yr1819 	;	E_yr1820 	;	E_yr1821 	;	E_yr1822 	;	E_yr1823 	;	E_yr1824 	;	E_yr191 	;	E_yr192 	;	E_yr193 	;	E_yr194 	;	E_yr195 	;	E_yr196 	;	E_yr197 	;	E_yr198 	;	E_yr199 	;	E_yr1910 	;	E_yr1911 	;	E_yr1912 	;	E_yr1913 	;	E_yr1914 	;	E_yr1915 	;	E_yr1916 	;	E_yr1917 	;	E_yr1918 	;	E_yr1920 	;	E_yr1921 	;	E_yr1922 	;	E_yr1923 	;	E_yr1924 	;	E_yr201 	;	E_yr202 	;	E_yr203 	;	E_yr204 	;	E_yr205 	;	E_yr206 	;	E_yr207 	;	E_yr208 	;	E_yr209 	;	E_yr2010 	;	E_yr2011 	;	E_yr2012 	;	E_yr2013 	;	E_yr2014 	;	E_yr2015 	;	E_yr2016 	;	E_yr2017 	;	E_yr2018 	;	E_yr2019 	;	E_yr2021 	;	E_yr2022 	;	E_yr2023 	;	E_yr2024 	;	E_yr21_1 	;	E_yr21_2 	;	E_yr21_3 	;	E_yr214 	;	E_yr215 	;	E_yr216 	;	E_yr217 	;	E_yr218 	;	E_yr219 	;	E_yr2110 	;	E_yr2111 	;	E_yr2112 	;	E_yr2113 	;	E_yr2114 	;	E_yr2115 	;	E_yr2116 	;	E_yr2117 	;	E_yr2118 	;	E_yr2119 	;	E_yr2120 	;	E_yr2122 	;	E_yr2123 	;	E_yr2124 	;	E_yr22_1 	;	E_yr22_2 	;	E_yr22_3 	;	E_yr22_4 	;	E_yr225 	;	E_yr226 	;	E_yr227 	;	E_yr228 	;	E_yr229 	;	E_yr2210 	;	E_yr2211 	;	E_yr2212 	;	E_yr2213 	;	E_yr2214 	;	E_yr2215 	;	E_yr2216 	;	E_yr2217 	;	E_yr2218 	;	E_yr2219 	;	E_yr2220 	;	E_yr2221 	;	E_yr2223 	;	E_yr2224 	;	E_yr231 	;	E_yr232 	;	E_yr233 	;	E_yr234 	;	E_yr235 	;	E_yr236 	;	E_yr237 	;	E_yr238 	;	E_yr239 	;	E_yr2310 	;	E_yr2311 	;	E_yr2312 	;	E_yr2313 	;	E_yr2314 	;	E_yr2315 	;	E_yr2316 	;	E_yr2317 	;	E_yr2318 	;	E_yr2319 	;	E_yr2320 	;	E_yr2321 	;	E_yr2322 	;	E_yr2324 	;	E_yr241 	;	E_yr242 	;	E_yr243 	;	E_yr244 	;	E_yr245 	;	E_yr246 	;	E_yr247 	;	E_yr248 	;	E_yr249 	;	E_yr2410 	;	E_yr2411 	;	E_yr2412 	;	E_yr2413 	;	E_yr2414 	;	E_yr2415 	;	E_yr2416 	;	E_yr2417 	;	E_yr2418 	;	E_yr2419 	;	E_yr2420 	;	E_yr2421 	;	E_yr2422 	;	E_yr2423  ];

E_yi = [ E_yi12 ;	E_yi13 	;	E_yi14 	;	E_yi15 	;	E_yi16 	;	E_yi17 	;	E_yi18 	;	E_yi19 	;	E_yi110 	;	E_yi1_11 	;	E_yi1_12 	;	E_yi1_13 	;	E_yi1_14 	;	E_yi1_15 	;	E_yi1_16 	;	E_yi1_17 	;	E_yi1_18 	;	E_yi1_19 	;	E_yi1_20 	;	E_yi1_21 	;	E_yi1_22 	;	E_yi1_23 	;	E_yi1_24 	;	E_yi21 	;	E_yi23 	;	E_yi24 	;	E_yi25 	;	E_yi26 	;	E_yi27 	;	E_yi28 	;	E_yi29 	;	E_yi210 	;	E_yi2_11 	;	E_yi2_12 	;	E_yi2_13 	;	E_yi2_14 	;	E_yi2_15 	;	E_yi2_16 	;	E_yi2_17 	;	E_yi2_18 	;	E_yi2_19 	;	E_yi2_20 	;	E_yi2_21 	;	E_yi2_22 	;	E_yi2_23 	;	E_yi2_24 	;	E_yi31 	;	E_yi32 	;	E_yi34 	;	E_yi35 	;	E_yi36 	;	E_yi37 	;	E_yi38 	;	E_yi39 	;	E_yi310 	;	E_yi311 	;	E_yi312 	;	E_yi313 	;	E_yi314 	;	E_yi315 	;	E_yi316 	;	E_yi317 	;	E_yi318 	;	E_yi319 	;	E_yi320 	;	E_yi321 	;	E_yi322 	;	E_yi323 	;	E_yi324 	;	E_yi41 	;	E_yi42 	;	E_yi43 	;	E_yi45 	;	E_yi46 	;	E_yi47 	;	E_yi48 	;	E_yi49 	;	E_yi410 	;	E_yi411 	;	E_yi412 	;	E_yi413 	;	E_yi414 	;	E_yi415 	;	E_yi416 	;	E_yi417 	;	E_yi418 	;	E_yi419 	;	E_yi420 	;	E_yi421 	;	E_yi422 	;	E_yi423 	;	E_yi424 	;	E_yi51 	;	E_yi52 	;	E_yi53 	;	E_yi54 	;	E_yi56 	;	E_yi57 	;	E_yi58 	;	E_yi59 	;	E_yi510 	;	E_yi511 	;	E_yi512 	;	E_yi513 	;	E_yi514 	;	E_yi515 	;	E_yi516 	;	E_yi517 	;	E_yi518 	;	E_yi519 	;	E_yi520 	;	E_yi521 	;	E_yi522 	;	E_yi523 	;	E_yi524 	;	E_yi61 	;	E_yi62 	;	E_yi63 	;	E_yi64 	;	E_yi65 	;	E_yi67 	;	E_yi68 	;	E_yi69 	;	E_yi610 	;	E_yi611 	;	E_yi612 	;	E_yi613 	;	E_yi614 	;	E_yi615 	;	E_yi616 	;	E_yi617 	;	E_yi618 	;	E_yi619 	;	E_yi620 	;	E_yi621 	;	E_yi622 	;	E_yi623 	;	E_yi624 	;	E_yi71 	;	E_yi72 	;	E_yi73 	;	E_yi74 	;	E_yi75 	;	E_yi76 	;	E_yi78 	;	E_yi79 	;	E_yi710 	;	E_yi711 	;	E_yi712 	;	E_yi713 	;	E_yi714 	;	E_yi715 	;	E_yi716 	;	E_yi717 	;	E_yi718 	;	E_yi719 	;	E_yi720 	;	E_yi721 	;	E_yi722 	;	E_yi723 	;	E_yi724 	;	E_yi81 	;	E_yi82 	;	E_yi83 	;	E_yi84 	;	E_yi85 	;	E_yi86 	;	E_yi87 	;	E_yi89 	;	E_yi810 	;	E_yi811 	;	E_yi812 	;	E_yi813 	;	E_yi814 	;	E_yi815 	;	E_yi816 	;	E_yi817 	;	E_yi818 	;	E_yi819 	;	E_yi820 	;	E_yi821 	;	E_yi822 	;	E_yi823 	;	E_yi824 	;	E_yi91 	;	E_yi92 	;	E_yi93 	;	E_yi94 	;	E_yi95 	;	E_yi96 	;	E_yi97 	;	E_yi98 	;	E_yi910 	;	E_yi911 	;	E_yi912 	;	E_yi913 	;	E_yi914 	;	E_yi915 	;	E_yi916 	;	E_yi917 	;	E_yi918 	;	E_yi919 	;	E_yi920 	;	E_yi921 	;	E_yi922 	;	E_yi923 	;	E_yi924 	;	E_yi101 	;	E_yi102 	;	E_yi103 	;	E_yi104 	;	E_yi105 	;	E_yi106 	;	E_yi107 	;	E_yi108 	;	E_yi109 	;	E_yi1011 	;	E_yi1012 	;	E_yi1013 	;	E_yi1014 	;	E_yi1015 	;	E_yi1016 	;	E_yi1017 	;	E_yi1018 	;	E_yi10_19 	;	E_yi10_20 	;	E_yi10_21 	;	E_yi10_22 	;	E_yi10_23 	;	E_yi10_24 	;	E_yi11_1 	;	E_yi11_2 	;	E_yi113 	;	E_yi114 	;	E_yi115 	;	E_yi116 	;	E_yi117 	;	E_yi118 	;	E_yi119 	;	E_yi1110 	;	E_yi1112 	;	E_yi1113 	;	E_yi1114 	;	E_yi1115 	;	E_yi1116 	;	E_yi1117 	;	E_yi1118 	;	E_yi1019 	;	E_yi1020 	;	E_yi1021 	;	E_yi1022 	;	E_yi1023 	;	E_yi1024 	;	E_yi12_1 	;	E_yi12_2 	;	E_yi123 	;	E_yi124 	;	E_yi125 	;	E_yi126 	;	E_yi127 	;	E_yi128 	;	E_yi129 	;	E_yi1210 	;	E_yi1211 	;	E_yi1213 	;	E_yi1214 	;	E_yi1215 	;	E_yi1216 	;	E_yi1217 	;	E_yi1218 	;	E_yi1219 	;	E_yi1220 	;	E_yi1221 	;	E_yi1222 	;	E_yi1223 	;	E_yi1224 	;	E_yi131 	;	E_yi132 	;	E_yi133 	;	E_yi134 	;	E_yi135 	;	E_yi136 	;	E_yi137 	;	E_yi138 	;	E_yi139 	;	E_yi1310 	;	E_yi1311 	;	E_yi1312 	;	E_yi1314 	;	E_yi1315 	;	E_yi1316 	;	E_yi1317 	;	E_yi1318 	;	E_yi1319 	;	E_yi1320 	;	E_yi1321 	;	E_yi1322 	;	E_yi1323 	;	E_yi1324 	;	E_yi141 	;	E_yi142 	;	E_yi143 	;	E_yi144 	;	E_yi145 	;	E_yi146 	;	E_yi147 	;	E_yi148 	;	E_yi149 	;	E_yi1410 	;	E_yi1411 	;	E_yi1412 	;	E_yi1413 	;	E_yi1415 	;	E_yi1416 	;	E_yi1417 	;	E_yi1418 	;	E_yi1419 	;	E_yi1420 	;	E_yi1421 	;	E_yi1422 	;	E_yi1423 	;	E_yi1424 	;	E_yi151 	;	E_yi152 	;	E_yi153 	;	E_yi154 	;	E_yi155 	;	E_yi156 	;	E_yi157 	;	E_yi158 	;	E_yi159 	;	E_yi1510 	;	E_yi1511 	;	E_yi1512 	;	E_yi1513 	;	E_yi1514 	;	E_yi1516 	;	E_yi1517 	;	E_yi1518 	;	E_yi1519 	;	E_yi1520 	;	E_yi1521 	;	E_yi1522 	;	E_yi1523 	;	E_yi1524 	;	E_yi161 	;	E_yi162 	;	E_yi163 	;	E_yi164 	;	E_yi165 	;	E_yi166 	;	E_yi167 	;	E_yi168 	;	E_yi169 	;	E_yi1610 	;	E_yi1611 	;	E_yi1612 	;	E_yi1613 	;	E_yi1614 	;	E_yi1615 	;	E_yi1617 	;	E_yi1618 	;	E_yi1619 	;	E_yi1620 	;	E_yi1621 	;	E_yi1622 	;	E_yi1623 	;	E_yi1624 	;	E_yi171 	;	E_yi172 	;	E_yi173 	;	E_yi174 	;	E_yi175 	;	E_yi176 	;	E_yi177 	;	E_yi178 	;	E_yi179 	;	E_yi1710 	;	E_yi1711 	;	E_yi1712 	;	E_yi1713 	;	E_yi1714 	;	E_yi1715 	;	E_yi1716 	;	E_yi1718 	;	E_yi1719 	;	E_yi1720 	;	E_yi1721 	;	E_yi1722 	;	E_yi1723 	;	E_yi1724 	;	E_yi181 	;	E_yi182 	;	E_yi183 	;	E_yi184 	;	E_yi185 	;	E_yi186 	;	E_yi187 	;	E_yi188 	;	E_yi189 	;	E_yi1810 	;	E_yi1811 	;	E_yi1812 	;	E_yi1813 	;	E_yi1814 	;	E_yi1815 	;	E_yi1816 	;	E_yi1817 	;	E_yi1819 	;	E_yi1820 	;	E_yi1821 	;	E_yi1822 	;	E_yi1823 	;	E_yi1824 	;	E_yi191 	;	E_yi192 	;	E_yi193 	;	E_yi194 	;	E_yi195 	;	E_yi196 	;	E_yi197 	;	E_yi198 	;	E_yi199 	;	E_yi1910 	;	E_yi1911 	;	E_yi1912 	;	E_yi1913 	;	E_yi1914 	;	E_yi1915 	;	E_yi1916 	;	E_yi1917 	;	E_yi1918 	;	E_yi1920 	;	E_yi1921 	;	E_yi1922 	;	E_yi1923 	;	E_yi1924 	;	E_yi201 	;	E_yi202 	;	E_yi203 	;	E_yi204 	;	E_yi205 	;	E_yi206 	;	E_yi207 	;	E_yi208 	;	E_yi209 	;	E_yi2010 	;	E_yi2011 	;	E_yi2012 	;	E_yi2013 	;	E_yi2014 	;	E_yi2015 	;	E_yi2016 	;	E_yi2017 	;	E_yi2018 	;	E_yi2019 	;	E_yi2021 	;	E_yi2022 	;	E_yi2023 	;	E_yi2024 	;	E_yi21_1 	;	E_yi21_2 	;	E_yi21_3 	;	E_yi214 	;	E_yi215 	;	E_yi216 	;	E_yi217 	;	E_yi218 	;	E_yi219 	;	E_yi2110 	;	E_yi2111 	;	E_yi2112 	;	E_yi2113 	;	E_yi2114 	;	E_yi2115 	;	E_yi2116 	;	E_yi2117 	;	E_yi2118 	;	E_yi2119 	;	E_yi2120 	;	E_yi2122 	;	E_yi2123 	;	E_yi2124 	;	E_yi22_1 	;	E_yi22_2 	;	E_yi22_3 	;	E_yi22_4 	;	E_yi225 	;	E_yi226 	;	E_yi227 	;	E_yi228 	;	E_yi229 	;	E_yi2210 	;	E_yi2211 	;	E_yi2212 	;	E_yi2213 	;	E_yi2214 	;	E_yi2215 	;	E_yi2216 	;	E_yi2217 	;	E_yi2218 	;	E_yi2219 	;	E_yi2220 	;	E_yi2221 	;	E_yi2223 	;	E_yi2224 	;	E_yi231 	;	E_yi232 	;	E_yi233 	;	E_yi234 	;	E_yi235 	;	E_yi236 	;	E_yi237 	;	E_yi238 	;	E_yi239 	;	E_yi2310 	;	E_yi2311 	;	E_yi2312 	;	E_yi2313 	;	E_yi2314 	;	E_yi2315 	;	E_yi2316 	;	E_yi2317 	;	E_yi2318 	;	E_yi2319 	;	E_yi2320 	;	E_yi2321 	;	E_yi2322 	;	E_yi2324 	;	E_yi241 	;	E_yi242 	;	E_yi243 	;	E_yi244 	;	E_yi245 	;	E_yi246 	;	E_yi247 	;	E_yi248 	;	E_yi249 	;	E_yi2410 	;	E_yi2411 	;	E_yi2412 	;	E_yi2413 	;	E_yi2414 	;	E_yi2415 	;	E_yi2416 	;	E_yi2417 	;	E_yi2418 	;	E_yi2419 	;	E_yi2420 	;	E_yi2421 	;	E_yi2422 	;	E_yi2423  ];

%% aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

Exr = sum(E_xr).^2;

Exi = sum(E_xi).^2;

Eyr = sum(E_yr).^2;

Eyi = sum(E_yi).^2;

%% aqui retira-se o módulo dos quadrados:

% V/m -> kV/cm = 10^-3/10^2
Erms = ((Exr + Exi + Eyr + Eyi).^(1/2))*(10^-5);
thetad = rad2deg(theta); %converte o ângulo de radianos para graus

Ecrit = fcn_supcrit(r*10^(2));
E_crit = linspace(Ecrit,Ecrit,360);

Emed = mean(Erms);
Emax = max(Erms);
Kirreg = Emax/Emed;

E4rms1 = [13.1158493758248 13.1290256929223 13.1423880262588 13.1559321592947 13.1696538185754 13.1835486751771 13.1976123461740 13.2118403961265 13.2262283385849 13.2407716376155 13.2554657093373 13.2703059234840 13.2852876049755 13.3004060355067 13.3156564551517 13.3310340639836 13.3465340237032 13.3621514592850 13.3778814606332 13.3937190842470 13.4096593548990 13.4256972673226 13.4418277879114 13.4580458564166 13.4743463876637 13.4907242732707 13.5071743833709 13.5236915683427 13.5402706605440 13.5569064760462 13.5735938163775 13.5903274702637 13.6071022153709 13.6239128200451 13.6407540450587 13.6576206453546 13.6745073717762 13.6914089728135 13.7083201963316 13.7252357912991 13.7421505095114 13.7590591073136 13.7759563473008 13.7928370000322 13.8096958457272 13.8265276759401 13.8433272952508 13.8600895229283 13.8768091945828 13.8934811638180 13.9101003038537 13.9266615091627 13.9431596970545 13.9595898092926 13.9759468136551 13.9922257055000 14.0084215093227 14.0245292802766 14.0405441056894 14.0564611065614 14.0722754390444 14.0879822959031 14.1035769079507 14.1190545454844 14.1344105196756 14.1496401839547 14.1647389353905 14.1797022160040 14.1945255141117 14.2092043656201 14.2237343553021 14.2381111180420 14.2523303400913 14.2663877602655 14.2802791711326 14.2940004201801 14.3075474109668 14.3209161042219 14.3341025189608 14.3471027335496 14.3599128867476 14.3725291787423 14.3849478721427 14.3971652929583 14.4091778315514 14.4209819435706 14.4325741508441 14.4439510422784 14.4551092746921 14.4660455736674 14.4767567343476 14.4872396222240 14.4974911739063 14.5075083978377 14.5172883750284 14.5268282597352 14.5361252801212 14.5451767389170 14.5539800140157 14.5625325590869 14.5708319041353 14.5788756560652 14.5866614991920 14.5941871957484 14.6014505863771 14.6084495905782 14.6151822071485 14.6216465145931 14.6278406715307 14.6337629170498 14.6394115710625 14.6447850346438 14.6498817903173 14.6547004023586 14.6592395170613 14.6634978629664 14.6674742511004 14.6711675751779 14.6745768117748 14.6777010205123 14.6805393441740 14.6830910088621 14.6853553240749 14.6873316828125 14.6890195616337 14.6904185207061 14.6915282038373 14.6923483384870 14.6928787357539 14.6931192903574 14.6930699805873 14.6927308682422 14.6921020985534 14.6911839000808 14.6899765846028 14.6884805469797 14.6866962650008 14.6846242992181 14.6822652927625 14.6796199711301 14.6766891419794 14.6734736948721 14.6699746010337 14.6661929130754 14.6621297647008 14.6577863704066 14.6531640251555 14.6482641040289 14.6430880618786 14.6376374329531 14.6319138304927 14.6259189463305 14.6196545504638 14.6131224906140 14.6063246917576 14.5992631556566 14.5919399603602 14.5843572596933 14.5765172827264 14.5684223332295 14.5600747891085 14.5514771018235 14.5426317957875 14.5335414677532 14.5242087861760 14.5146364905649 14.5048273908087 14.4947843664948 14.4845103661980 14.4740084067603 14.4632815725495 14.4523330147012 14.4411659503392 14.4297836617833 14.4181894957334 14.4063868624373 14.3943792348448 14.3821701477368 14.3697631968385 14.3571620379149 14.3443703858501 14.3313920137038 14.3182307517489 14.3048904864999 14.2913751597109 14.2776887673632 14.2638353586299 14.2498190348289 14.2356439483496 14.2213143015655 14.2068343457299 14.1922083798509 14.1774407495516 14.1625358459089 14.1474981042781 14.1323320030937 14.1170420626654 14.1016328439412 14.0861089472664 14.0704750111166 14.0547357108220 14.0388957572670 14.0229598955793 14.0069329038035 13.9908195915553 13.9746247986606 13.9583533937827 13.9420102730316 13.9256003585634 13.9091285971536 13.8925999587740 13.8760194351439 13.8593920382732 13.8427227989817 13.8260167654277 13.8092790016116 13.7925145858595 13.7757286093114 13.7589261743973 13.7421123932900 13.7252923863624 13.7084712806362 13.6916542082005 13.6748463046601 13.6580527075373 13.6412785546946 13.6245289827389 13.6078091254222 13.5911241120343 13.5744790658024 13.5578791022711 13.5413293276922 13.5248348374037 13.5084007142086 13.4920320267598 13.4757338279312 13.4595111532008 13.4433690190285 13.4273124212358 13.4113463333847 13.3954757051746 13.3797054608106 13.3640404974240 13.3484856834392 13.3330458570027 13.3177258243737 13.3025303583494 13.2874641966902 13.2725320405383 13.2577385528723 13.2430883569565 13.2285860347853 13.2142361255734 13.2000431242287 13.1860114798426 13.1721455942113 13.1584498203457 13.1449284610084 13.1315857672616 13.1184259370427 13.1054531137343 13.0926713847604 13.0800847802140 13.0676972714762 13.0555127698787 13.0435351253570 13.0317681251604 13.0202154925438 13.0088808855014 12.9977678955267 12.9868800463642 12.9762207928231 12.9657935195841 12.9556015400363 12.9456480951478 12.9359363523363 12.9264694043979 12.9172502684263 12.9082818847821 12.8995671160729 12.8911087461657 12.8829094792202 12.8749719387588 12.8672986667478 12.8598921227186 12.8527546829193 12.8458886394683 12.8392961995724 12.8329794847460 12.8269405300661 12.8211812834583 12.8157036050227 12.8105092663583 12.8055999499446 12.8009772485456 12.7966426646466 12.7925976098956 12.7888434046209 12.7853812773388 12.7822123643046 12.7793377091076 12.7767582622741 12.7744748809177 12.7724883284144 12.7707992741102 12.7694082930572 12.7683158657812 12.7675223780922 12.7670281209071 12.7668332901102 12.7669379864661 12.7673422155230 12.7680458875926 12.7690488177252 12.7703507257408 12.7719512362766 12.7738498788762 12.7760460881110 12.7785392037169 12.7813284707864 12.7844130399821 12.7877919677698 12.7914642167068 12.7954286557431 12.7996840605627 12.8042291139543 12.8090624062187 12.8141824355936 12.8195876087285 12.8252762411790 12.8312465579313 12.8374966939649 12.8440246948430 12.8508285173348 12.8579060300640 12.8652550141937 12.8728731641401 12.8807580883133 12.8889073098945 12.8973182676378 12.9059883166990 12.9149147295050 12.9240946966424 12.9335253277784 12.9432036526088 12.9531266218385 12.9632911081881 12.9736939074244 12.9843317394270 12.9952012492732 13.0062990083557 13.0176215155226 13.0291651982518 13.0409264138377 13.0529014506172 13.0650865292102 13.0774778037909 13.0900713633803 13.1028632331623 13.1158493758248];
E4rms2 = [15.2385193594945 15.2515716570887 15.2643768206039 15.2769310577824 15.2892306535720 15.3012719710806 15.3130514525084 15.3245656200576 15.3358110768122 15.3467845076006 15.3574826798232 15.3679024442674 15.3780407358880 15.3878945745679 15.3974610658522 15.4067374016623 15.4157208609781 15.4244088105028 15.4327987053009 15.4408880894090 15.4486745964264 15.4561559500805 15.4633299647722 15.4701945460848 15.4767476912832 15.4829874897844 15.4889121236008 15.4945198677646 15.4998090907268 15.5047782547298 15.5094259161644 15.5137507258968 15.5177514295757 15.5214268679095 15.5247759769320 15.5277977882416 15.5304914292040 15.5328561231551 15.5348911895631 15.5365960441749 15.5379701991373 15.5390132631044 15.5397249413016 15.5401050355937 15.5401534445173 15.5398701632765 15.5392552837470 15.5383089944348 15.5370315804162 15.5354234232652 15.5334850009402 15.5312168876793 15.5286197538264 15.5256943656921 15.5224415853404 15.5188623703791 15.5149577737349 15.5107289433834 15.5061771220703 15.5013036470103 15.4961099495607 15.4905975548735 15.4847680815179 15.4786232411012 15.4721648378362 15.4653947681063 15.4583150200213 15.4509276729005 15.4432348967919 15.4352389519347 15.4269421882052 15.4183470445285 15.4094560483031 15.4002718147637 15.3907970463345 15.3810345319628 15.3709871464373 15.3606578496513 15.3500496858850 15.3391657830363 15.3280093518290 15.3165836850149 15.3048921565339 15.2929382206609 15.2807254111259 15.2682573402164 15.2555376978402 15.2425702505939 15.2293588407665 15.2159073853628 15.2022198750725 15.1883003732288 15.1741530147527 15.1597820050400 15.1451916188708 15.1303861992616 15.1153701563014 15.1001479659887 15.0847241690015 15.0691033694887 15.0532902338066 15.0372894892592 15.0211059227905 15.0047443796671 14.9882097621525 14.9715070281331 14.9546411897396 14.9376173119394 14.9204405111245 14.9031159536516 14.8856488543799 14.8680444751965 14.8503081234886 14.8324451506359 14.8144609504656 14.7963609576750 14.7781506462623 14.7598355279253 14.7414211504341 14.7229130960147 14.7043169796717 14.6856384475508 14.6668831752279 14.6480568660293 14.6291652493080 14.6102140787178 14.5912091304725 14.5721562015929 14.5530611081337 14.5339296834112 14.5147677762067 14.4955812489651 14.4763759759870 14.4571578416016 14.4379327383416 14.4187065651014 14.3994852252901 14.3802746249795 14.3610806710478 14.3419092693052 14.3227663226408 14.3036577291296 14.2845893801700 14.2655671586002 14.2465969368120 14.2276845748763 14.2088359186573 14.1900567979242 14.1713530244794 14.1527303902800 14.1341946655513 14.1157515969219 14.0974069055546 14.0791662852849 14.0610354007573 14.0430198855812 14.0251253404848 14.0073573314803 13.9897213880386 13.9722230012716 13.9548676221286 13.9376606596007 13.9206074789372 13.9037133998778 13.8869836948945 13.8704235874519 13.8540382502764 13.8378328036503 13.8218123137108 13.8059817907751 13.7903461876796 13.7749103981384 13.7596792551189 13.7446575292408 13.7298499271909 13.7152610901605 13.7008955923089 13.6867579392422 13.6728525665177 13.6591838381714 13.6457560452724 13.6325734044961 13.6196400567235 13.6069600656729 13.5945374165471 13.5823760147142 13.5704796844113 13.5588521674802 13.5474971221252 13.5364181217012 13.5256186535332 13.5151021177595 13.5048718262102 13.4949310013098 13.4852827750136 13.4759301877692 13.4668761875195 13.4581236287226 13.4496752714147 13.4415337802964 13.4337017238577 13.4261815735279 13.4189757028624 13.4120863867629 13.4055158007267 13.3992660201294 13.3933390195435 13.3877366720874 13.3824607488119 13.3775129181099 13.3728947451751 13.3686076914841 13.3646531143162 13.3610322662974 13.3577462949991 13.3547962425577 13.3521830453215 13.3499075335463 13.3479704311239 13.3463723553327 13.3451138166360 13.3441952185155 13.3436168573193 13.3433789221793 13.3434814949274 13.3439245500722 13.3447079547983 13.3458314690040 13.3472947453691 13.3490973294694 13.3512386599103 13.3537180685087 13.3565347804998 13.3596879147819 13.3631764842019 13.3669993958628 13.3711554514778 13.3756433477519 13.3804616767995 13.3856089265926 13.3910834814577 13.3968836225774 13.4030075285698 13.4094532760478 13.4162188402647 13.4233020957518 13.4307008170153 13.4384126792594 13.4464352591280 13.4547660355050 13.4634023903341 13.4723416094544 13.4815808835039 13.4911173088257 13.5009478884117 13.5110695328944 13.5214790615469 13.5321732033249 13.5431485979370 13.5544017969563 13.5659292649408 13.5777273805938 13.5897924379655 13.6021206476597 13.6147081380931 13.6275509567534 13.6406450715238 13.6539863719961 13.6675706708317 13.6813937051546 13.6954511379393 13.7097385594664 13.7242514887708 13.7389853751225 13.7539355995439 13.7690974763205 13.7844662545752 13.8000371198249 13.8158051955873 13.8317655449915 13.8479131724156 13.8642430251404 13.8807499950305 13.8974289202179 13.9142745868170 13.9312817306599 13.9484450390195 13.9657591523873 13.9832186662385 14.0008181328158 14.0185520629310 14.0364149277915 14.0544011608026 14.0725051594136 14.0907212869642 14.1090438745458 14.1274672228406 14.1459856040220 14.1645932636172 14.1832844223869 14.2020532782312 14.2208940080694 14.2398007697451 14.2587677039265 14.2777889360107 14.2968585780261 14.3159707305362 14.3351194845558 14.3542989234454 14.3735031248095 14.3927261624186 14.4119621080772 14.4312050335443 14.4504490124017 14.4696881219488 14.4889164450722 14.5081280721188 14.5273171027633 14.5464776478468 14.5656038312344 14.5846897916512 14.6037296844928 14.6227176836529 14.6416479833163 14.6605147997510 14.6793123730817 14.6980349690594 14.7166768807974 14.7352324305155 14.7536959712556 14.7720618885780 14.7903246022522 14.8084785679255 14.8265182787745 14.8444382671361 14.8622331061234 14.8798974112240 14.8974258418732 14.9148131030166 14.9320539466447 14.9491431733053 14.9660756336073 14.9828462296932 14.9994499166908 15.0158817041452 15.0321366574304 15.0482098991380 15.0640966104372 15.0797920324234 15.0952914674325 15.1105902803391 15.1256838998261 15.1405678196406 15.1552375998100 15.1696888678514 15.1839173199426 15.1979187220795 15.2116889112025 15.2252237963008 15.2385193594945];
E4rms3 = [14.7935150444808 14.7820902751874 14.7704298800819 14.7585372897528 14.7464160074719 14.7340696083080 14.7215017382197 14.7087161131222 14.6957165179263 14.6825068055601 14.6690908959527 14.6554727750101 14.6416564935500 14.6276461662226 14.6134459704009 14.5990601450533 14.5844929895824 14.5697488626486 14.5548321809660 14.5397474180714 14.5244991030741 14.5090918193816 14.4935302034048 14.4778189432294 14.4619627772769 14.4459664929412 14.4298349251962 14.4135729551896 14.3971855088108 14.3806775552379 14.3640541054687 14.3473202108278 14.3304809614550 14.3135414847695 14.2965069439256 14.2793825362470 14.2621734916305 14.2448850709543 14.2275225644525 14.2100912900798 14.1925965918590 14.1750438382211 14.1574384203118 14.1397857503096 14.1220912597177 14.1043603976304 14.0865986290167 14.0688114329716 14.0510043009598 14.0331827350590 14.0153522461776 13.9975183522922 13.9796865766350 13.9618624459264 13.9440514885546 13.9262592327728 13.9084912049004 13.8907529274974 13.8730499175502 13.8553876846545 13.8377717291945 13.8202075405230 13.8027005951343 13.7852563548583 13.7678802650284 13.7505777526701 13.7333542247051 13.7162150661161 13.6991656381675 13.6822112766040 13.6653572898596 13.6486089572653 13.6319715272970 13.6154502157976 13.5990502042214 13.5827766378919 13.5666346242785 13.5506292312563 13.5347654854181 13.5190483703723 13.5034828250584 13.4880737420918 13.4728259661077 13.4577442921311 13.4428334639604 13.4280981725746 13.4135430545437 13.3991726904852 13.3849916035026 13.3710042576831 13.3572150565882 13.3436283417781 13.3302483913658 13.3170794185629 13.3041255702901 13.2913909257808 13.2788794952144 13.2665952183959 13.2545419634218 13.2427235254102 13.2311436252277 13.2198059082679 13.2087139432312 13.1978712209459 13.1872811532260 13.1769470717356 13.1668722268953 13.1570597868113 13.1475128362474 13.1382343756025 13.1292273199341 13.1204944980175 13.1120386514073 13.1038624335624 13.0959684089874 13.0883590523908 13.0810367478996 13.0740037882917 13.0672623742520 13.0608146136866 13.0546625210267 13.0488080166179 13.0432529260867 13.0379989797835 13.0330478122298 13.0284009616095 13.0240598692925 13.0200258793899 13.0163002383366 13.0128840945181 13.0097784979178 13.0069843998035 13.0045026524494 13.0023340088829 13.0004791226745 12.9989385477524 12.9977127382519 12.9968020484020 12.9962067324437 12.9959269445718 12.9959627389332 12.9963140696250 12.9969807907574 12.9979626565314 12.9992593213498 13.0008703399728 13.0027951676955 13.0050331605556 13.0075835755890 13.0104455711097 13.0136182070097 13.0171004451124 13.0208911495468 13.0249890871589 13.0293929279453 13.0341012455320 13.0391125176770 13.0444251268068 13.0500373605827 13.0559474125012 13.0621533825252 13.0686532777437 13.0754450130626 13.0825264119304 13.0898952070877 13.0975490413539 13.1054854684346 13.1137019537712 13.1221958754058 13.1309645248848 13.1400051081880 13.1493147466876 13.1588904781315 13.1687292576595 13.1788279588417 13.1891833747468 13.1997922190395 13.2106511270979 13.2217566571620 13.2331052915033 13.2446934376266 13.2565174294864 13.2685735287337 13.2808579259876 13.2933667421251 13.3060960295956 13.3190417737582 13.3321998942426 13.3455662463244 13.3591366223260 13.3729067530378 13.3868723091556 13.4010289027389 13.4153720886868 13.4298973662307 13.4446001804426 13.4594759237659 13.4745199375522 13.4897275136221 13.5050938958331 13.5206142816674 13.5362838238262 13.5520976318416 13.5680507736987 13.5841382774664 13.6003551329401 13.6166962932945 13.6331566767432 13.6497311682107 13.6664146210018 13.6832018584919 13.7000876758123 13.7170668415471 13.7341340994222 13.7512841700192 13.7685117524795 13.7858115262035 13.8031781525661 13.8206062766299 13.8380905288489 13.8556255267830 13.8732058768116 13.8908261758270 13.9084810129585 13.9261649712580 13.9438726294078 13.9615985634111 13.9793373482819 13.9970835597237 14.0148317758152 14.0325765786717 14.0503125561138 14.0680343033185 14.0857364244635 14.1034135343710 14.1210602601260 14.1386712426996 14.1562411385532 14.1737646212318 14.1912363829453 14.2086511361495 14.2260036150852 14.2432885773503 14.2605008054014 14.2776351080973 14.2946863221837 14.3116493137920 14.3285189799137 14.3452902498462 14.3619580866491 14.3785174885699 14.3949634904376 14.4112911650765 14.4274956246712 14.4435720221224 14.4595155524019 14.4753214538660 14.4909850095627 14.5065015485181 14.5218664470177 14.5370751298453 14.5521230715154 14.5670057974990 14.5817188854067 14.5962579661768 14.6106187252142 14.6247969035519 14.6387882989511 14.6525887670055 14.6661942222319 14.6796006391088 14.6928040531400 14.7058005618656 14.7185863258618 14.7311575697343 14.7435105830615 14.7556417213613 14.7675474069960 14.7792241300883 14.7906684494006 14.8018769932008 14.8128464601060 14.8235736199152 14.8340553144067 14.8442884581299 14.8542700391823 14.8639971199364 14.8734668377903 14.8826764058679 14.8916231137075 14.9003043279344 14.9087174929279 14.9168601314316 14.9247298451794 14.9323243154923 14.9396413038608 14.9466786524780 14.9534342848087 14.9599062060924 14.9660925038405 14.9719913483344 14.9776009930740 14.9829197752299 14.9879461160678 14.9926785213562 14.9971155817535 15.0012559731740 15.0050984571518 15.0086418811605 15.0118851789238 15.0148273707307 15.0174675636826 15.0198049519769 15.0218388171256 15.0235685281859 15.0249935419525 15.0261134031428 15.0269277445637 15.0274362872389 15.0276388405457 15.0275353023193 15.0271256589245 15.0264099853332 15.0253884451631 15.0240612907039 15.0224288629222 15.0204915914520 15.0182499945494 15.0157046790472 15.0128563402772 15.0097057619685 15.0062538161353 15.0025014629389 14.9984497505300 14.9940998148647 14.9894528795062 14.9845102554020 14.9792733406374 14.9737436201731 14.9679226655553 14.9618121346007 14.9554137710716 14.9487294043173 14.9417609488962 14.9345104041742 14.9269798539038 14.9191714657772 14.9110874909545 14.9027302635734 14.8941022002302 14.8852057994401 14.8760436410716 14.8666183857626 14.8569327743019 14.8469896269982 14.8367918430149 14.8263423996880 14.8156443518137 14.8047008309142 14.7935150444808];
E4rms4 = [12.7589173837896 12.7474867033265 12.7362497664670 12.7252101287021 12.7143712838741 12.7037366629893 12.6933096330543 12.6830934959357 12.6730914872389 12.6633067752177 12.6537424596958 12.6444015710274 12.6352870690679 12.6264018421788 12.6177487062520 12.6093304037665 12.6011496028616 12.5932088964443 12.5855108013195 12.5780577573446 12.5708521266145 12.5638961926725 12.5571921597514 12.5507421520327 12.5445482129435 12.5386123044794 12.5329363065502 12.5275220163598 12.5223711478111 12.5174853309388 12.5128661113750 12.5085149498406 12.5044332216663 12.5006222163385 12.4970831370822 12.4938171004718 12.4908251360593 12.4881081860510 12.4856671049984 12.4835026595249 12.4816155280800 12.4800063007298 12.4786754789606 12.4776234755340 12.4768506143600 12.4763571303898 12.4761431695618 12.4762087887599 12.4765539558034 12.4771785494756 12.4780823595667 12.4792650869715 12.4807263437789 12.4824656534370 12.4844824509063 12.4867760828634 12.4893458079398 12.4921907969711 12.4953101332886 12.4987028130375 12.5023677455245 12.5063037535940 12.5105095740274 12.5149838579894 12.5197251714773 12.5247319958155 12.5300027281892 12.5355356821671 12.5413290882986 12.5473810947127 12.5536897677494 12.5602530926102 12.5670689740645 12.5741352371523 12.5814496279259 12.5890098142197 12.5968133864539 12.6048578584392 12.6131406682414 12.6216591790491 12.6304106800675 12.6393923874525 12.6486014452530 12.6580349263879 12.6676898336437 12.6775631007010 12.6876515931714 12.6979521096818 12.7084613829475 12.7191760809030 12.7300928078293 12.7412081055144 12.7525184544414 12.7640202749694 12.7757099285728 12.7875837190721 12.7996378938891 12.8118686453412 12.8242721119177 12.8368443796093 12.8495814832286 12.8624794077703 12.8755340897630 12.8887414186497 12.9020972381933 12.9155973478745 12.9292375043166 12.9430134227183 12.9569207783126 12.9709552078133 12.9851123108897 12.9993876516575 13.0137767601515 13.0282751338409 13.0428782391398 13.0575815129094 13.0723803639953 13.0872701747564 13.1022463025920 13.1173040815018 13.1324388236064 13.1476458207272 13.1629203459131 13.1782576550147 13.1936529882315 13.2091015716734 13.2245986189197 13.2401393325803 13.2557189058475 13.2713325240598 13.2869753662500 13.3026426066972 13.3183294164775 13.3340309650026 13.3497424215638 13.3654589568616 13.3811757445317 13.3968879626678 13.4125907953369 13.4282794340751 13.4439490794020 13.4595949422880 13.4752122456488 13.4907962258087 13.5063421339541 13.5218452375889 13.5373008219663 13.5527041915054 13.5680506712130 13.5833356080816 13.5985543724628 13.6137023594471 13.6287749902194 13.6437677134042 13.6586760063843 13.6734953766211 13.6882213629479 13.7028495368511 13.7173755037331 13.7317949041618 13.7461034151020 13.7602967511288 13.7743706656217 13.7883209519486 13.8021434446230 13.8158340204513 13.8293885996531 13.8428031469761 13.8560736727771 13.8691962340963 13.8821669357079 13.8949819311547 13.9076374237594 13.9201296676220 13.9324549685938 13.9446096852344 13.9565902297529 13.9683930689228 13.9800147249841 13.9914517765206 14.0027008593273 14.0137586672463 14.0246219529918 14.0352875289565 14.0457522679945 14.0560131041866 14.0660670335882 14.0759111149595 14.0855424704718 14.0949582863981 14.1041558137876 14.1131323691167 14.1218853349258 14.1304121604347 14.1387103621416 14.1467775243999 14.1546112999864 14.1622094106386 14.1695696475859 14.1766898720540 14.1835680157591 14.1902020813770 14.1965901430013 14.2027303465818 14.2086209103441 14.2142601251938 14.2196463551033 14.2247780374808 14.2296536835270 14.2342718785628 14.2386312823573 14.2427306294276 14.2465687293262 14.2501444669028 14.2534568025680 14.2565047725308 14.2592874890089 14.2618041404397 14.2640539916727 14.2660363841328 14.2677507359814 14.2691965422614 14.2703733750055 14.2712808833646 14.2719187936825 14.2722869095812 14.2723851120179 14.2722133593292 14.2717716872545 14.2710602089562 14.2700791150065 14.2688286733728 14.2673092293770 14.2655212056404 14.2634651020212 14.2611414955199 14.2585510401827 14.2556944669800 14.2525725836705 14.2491862746452 14.2455365007687 14.2416242991735 14.2374507830841 14.2330171415651 14.2283246393104 14.2233746163684 14.2181684878793 14.2127077437862 14.2069939485129 14.2010287406532 14.1948138326255 14.1883510102935 14.1816421326061 14.1746891311872 14.1674940099141 14.1600588444943 14.1523857819993 14.1444770403908 14.1363349080259 14.1279617431551 14.1193599733775 14.1105320950880 14.1014806729185 14.0922083391329 14.0827177930279 14.0730118002847 14.0630931923391 14.0529648656920 14.0426297812226 14.0320909634855 14.0213514999572 14.0104145403042 13.9992832955983 13.9879610375209 13.9764510975567 13.9647568661400 13.9528817918199 13.9408293803648 13.9286031938755 13.9162068498610 13.9036440203001 13.8909184306787 13.8780338590155 13.8649941348515 13.8518031382334 13.8384647986789 13.8249830940951 13.8113620497153 13.7976057369879 13.7837182724521 13.7697038165972 13.7555665727147 13.7413107856968 13.7269407408476 13.7124607626667 13.6978752136177 13.6831884928520 13.6684050349630 13.6535293086815 13.6385658155649 13.6235190886869 13.6083936912841 13.5931942154050 13.5779252805367 13.5625915322170 13.5471976406308 13.5317482991894 13.5162482231101 13.5007021479612 13.4851148282016 13.4694910357281 13.4538355583688 13.4381531984118 13.4224487710846 13.4067271030491 13.3909930308702 13.3752513994856 13.3595070606683 13.3437648714633 13.3280296926415 13.3123063871363 13.2965998184626 13.2809148491521 13.2652563391683 13.2496291443227 13.2340381146872 13.2184880930086 13.2029839131077 13.1875303982943 13.1721323597716 13.1567945950386 13.1415218863004 13.1263189988779 13.1111906796191 13.0961416553107 13.0811766310978 13.0663002889071 13.0515172858729 13.0368322527723 13.0222497924658 13.0077744783394 12.9934108527673 12.9791634255725 12.9650366725005 12.9510350337013 12.9371629122265 12.9234246725344 12.9098246390056 12.8963670944766 12.8830562787820 12.8698963873115 12.8568915695831 12.8440459278334 12.8313635156163 12.8188483364277 12.8065043423391 12.7943354326546 12.7823454525819 12.7705381919234 12.7589173837896];
E4rms5 = [13.1466881075564 13.1576967924959 13.1689212559759 13.1803579366963 13.1920032087012 13.2038533826221 13.2159047069408 13.2281533692733 13.2405954976707 13.2532271619426 13.2660443749929 13.2790430941809 13.2922192226914 13.3055686109272 13.3190870579146 13.3327703127262 13.3466140759150 13.3606140009669 13.3747656957643 13.3890647240623 13.4035066069766 13.4180868244865 13.4328008169452 13.4476439866015 13.4626116991286 13.4776992851711 13.4929020418862 13.5082152345056 13.5236340978943 13.5391538381226 13.5547696340402 13.5704766388567 13.5862699817268 13.6021447693341 13.6180960874849 13.6341190027030 13.6502085638175 13.6663598035654 13.6825677401825 13.6988273790005 13.7151337140374 13.7314817295973 13.7478664018488 13.7642827004226 13.7807255899927 13.7971900318466 13.8136709854704 13.8301634101117 13.8466622663429 13.8631625176190 13.8796591318193 13.8961470828031 13.9126213519220 13.9290769295658 13.9455088166604 13.9619120261755 13.9782815846285 13.9946125335587 14.0108999310006 14.0271388529479 14.0433243948005 14.0594516728009 14.0755158254526 14.0915120149414 14.1074354285170 14.1232812798795 14.1390448105613 14.1547212912516 14.1703060231580 14.1857943393211 14.2011816059216 14.2164632235602 14.2316346285527 14.2466912941764 14.2616287319093 14.2764424926586 14.2911281679764 14.3056813912323 14.3200978388054 14.3343732312332 14.3485033343454 14.3624839603937 14.3763109691484 14.3899802689855 14.4034878179534 14.4168296248273 14.4300017501288 14.4430003071573 14.4558214629622 14.4684614393385 14.4809165137748 14.4931830203959 14.5052573508955 14.5171359554176 14.5288153434645 14.5402920847536 14.5515628100623 14.5626242120782 14.5734730461893 14.5841061312977 14.5945203505847 14.6047126522850 14.6146800504142 14.6244196254955 14.6339285252783 14.6432039654158 14.6522432301387 14.6610436729076 14.6696027170614 14.6779178564231 14.6859866559059 14.6938067521093 14.7013758538669 14.7086917428140 14.7157522739244 14.7225553760099 14.7290990522379 14.7353813806117 14.7414005144281 14.7471546827490 14.7526421908066 14.7578614204543 14.7628108305377 14.7674889573029 14.7718944147494 14.7760258949893 14.7798821685813 14.7834620848516 14.7867645721929 14.7897886383587 14.7925333707241 14.7949979365453 14.7971815831966 14.7990836383883 14.8007035103770 14.8020406881489 14.8030947415889 14.8038653216415 14.8043521604460 14.8045550714527 14.8044739495402 14.8041087710841 14.8034595940445 14.8025265580098 14.8013098842307 14.7998098756467 14.7980269168809 14.7959614742196 14.7936140955853 14.7909854104838 14.7880761299250 14.7848870463388 14.7814190334697 14.7776730462499 14.7736501206490 14.7693513735161 14.7647780023943 14.7599312853203 14.7548125805984 14.7494233265635 14.7437650413202 14.7378393224578 14.7316478467499 14.7251923698356 14.7184747258734 14.7114968271805 14.7042606638459 14.6967683033289 14.6890218900285 14.6810236448377 14.6727758646720 14.6642809219793 14.6555412642254 14.6465594133566 14.6373379652437 14.6278795890996 14.6181870268764 14.6082630926382 14.5981106719133 14.5877327210208 14.5771322663751 14.5663124037687 14.5552762976275 14.5440271802484 14.5325683510069 14.5209031755486 14.5090350849491 14.4969675748561 14.4847042046060 14.4722485963171 14.4596044339567 14.4467754623884 14.4337654863948 14.4205783696737 14.4072180338152 14.3936884572500 14.3799936741823 14.3661377734911 14.3521248976139 14.3379592414044 14.3236450509724 14.3091866224946 14.2945883010053 14.2798544791695 14.2649895960315 14.2499981357354 14.2348846262355 14.2196536379781 14.2043097825676 14.1888577114039 14.1733021143108 14.1576477181409 14.1418992853582 14.1260616125984 14.1101395292253 14.0941378958616 14.0780616028925 14.0619155689656 14.0457047394765 14.0294340850238 14.0131085998598 13.9967333003331 13.9803132232917 13.9638534245071 13.9473589770556 13.9308349697040 13.9142865052821 13.8977186990355 13.8811366769765 13.8645455742274 13.8479505333449 13.8313567026456 13.8147692345205 13.7981932837377 13.7816340057537 13.7650965550020 13.7485860831899 13.7321077375856 13.7156666593051 13.6992679815894 13.6829168281011 13.6666183111842 13.6503775301743 13.6341995696505 13.6180894977467 13.6020523644204 13.5860931997512 13.5702170122402 13.5544287870934 13.5387334845443 13.5231360381630 13.5076413531612 13.4922543047381 13.4769797364062 13.4618224583346 13.4467872457170 13.4318788371292 13.4171019329090 13.4024611935472 13.3879612381022 13.3736066426083 13.3594019385079 13.3453516111135 13.3314600980597 13.3177317877969 13.3041710180744 13.2907820744784 13.2775691889519 13.2645365383535 13.2516882430433 13.2390283654579 13.2265609087477 13.2142898154057 13.2022189659255 13.1903521774946 13.1786932026822 13.1672457281902 13.1560133735905 13.1449996901135 13.1342081594464 13.1236421925630 13.1133051285762 13.1032002336247 13.0933306997718 13.0836996439442 13.0743101069001 13.0651650522011 13.0562673652460 13.0476198523090 13.0392252396098 13.0310861724182 13.0232052141973 13.0155848457441 13.0082274643902 13.0011353832215 12.9943108303344 12.9877559480936 12.9814727924681 12.9754633323580 12.9697294489613 12.9642729351885 12.9590954950842 12.9541987432949 12.9495842045621 12.9452533132503 12.9412074128990 12.9374477558111 12.9339755026799 12.9307917222315 12.9278973909033 12.9252933925742 12.9229805182875 12.9209594660481 12.9192308406155 12.9177951533507 12.9166528220855 12.9158041710237 12.9152494306820 12.9149887378433 12.9150221355636 12.9153495732016 12.9159709064677 12.9168858975272 12.9180942151140 12.9195954346885 12.9213890386193 12.9234744164037 12.9258508649068 12.9285175886452 12.9314737000934 12.9347182200180 12.9382500778522 12.9420681120930 12.9461710707303 12.9505576117083 12.9552263034110 12.9601756251869 12.9654039678914 12.9709096344727 12.9766908405727 12.9827457151653 12.9890723012216 12.9956685564054 13.0025323537920 13.0096614826196 13.0170536490677 13.0247064770632 13.0326175091102 13.0407842071532 13.0492039534609 13.0578740515417 13.0667917270782 13.0759541288981 13.0853583299581 13.0950013283639 13.1048800484064 13.1149913416276 13.1253319879084 13.1358986965796 13.1466881075564];
E4rms6 = [14.5824601147975 14.5929497002216 14.6032152599678 14.6132538025808 14.6230624070515 14.6326382235012 14.6419784738418 14.6510804524143 14.6599415266025 14.6685591374272 14.6769308001127 14.6850541046381 14.6929267162582 14.7005463760088 14.7079109011857 14.7150181858063 14.7218662010435 14.7284529956450 14.7347766963264 14.7408355081457 14.7466277148536 14.7521516792295 14.7574058433913 14.7623887290865 14.7670989379608 14.7715351518168 14.7756961328372 14.7795807238036 14.7831878482829 14.7865165108055 14.7895657970169 14.7923348738146 14.7948229894667 14.7970294737047 14.7989537378075 14.8005952746662 14.8019536588192 14.8030285464872 14.8038196755761 14.8043268656697 14.8045500179984 14.8044891154046 14.8041442222669 14.8035154844343 14.8026031291288 14.8014074648213 14.7999288811164 14.7981678485983 14.7961249186692 14.7938007233724 14.7911959751867 14.7883114668320 14.7851480710155 14.7817067402124 14.7779885063839 14.7739944807028 14.7697258532689 14.7651838927883 14.7603699462453 14.7552854385640 14.7499318722448 14.7443108269897 14.7384239593032 14.7322730020956 14.7258597642419 14.7191861301438 14.7122540592862 14.7050655857313 14.6976228176519 14.6899279368138 14.6819831980500 14.6737909287068 14.6653535281053 14.6566734669520 14.6477532867422 14.6385955991529 14.6292030854222 14.6195784956852 14.6097246483297 14.5996444293084 14.5893407914362 14.5788167536837 14.5680754004373 14.5571198807502 14.5459534075731 14.5345792569725 14.5230007673149 14.5112213384623 14.4992444309093 14.4870735649425 14.4747123197545 14.4621643325500 14.4494332976438 14.4365229655073 14.4234371418406 14.4101796865944 14.3967545129794 14.3831655864804 14.3694169238121 14.3555125918957 14.3414567067887 14.3272534326227 14.3129069804949 14.2984216073586 14.2838016149043 14.2690513484020 14.2541751955376 14.2391775852294 14.2240629864392 14.2088359069429 14.1935008921003 14.1780625236153 14.1625254182504 14.1468942265588 14.1311736315876 14.1153683475473 14.0994831184968 14.0835227169943 14.0674919427290 14.0513956211649 14.0352386021252 14.0190257584201 14.0027619844025 13.9864521945599 13.9701013220579 13.9537143172925 13.9372961464211 13.9208517898869 13.9043862409256 13.8879045040727 13.8714115936459 13.8549125322302 13.8384123491496 13.8219160789273 13.8054287597460 13.7889554318916 13.7725011361915 13.7560709124552 13.7396697978995 13.7233028255688 13.7069750227644 13.6906914094442 13.6744569966486 13.6582767849019 13.6421557626178 13.6260989045118 13.6101111700008 13.5941975016047 13.5783628233587 13.5626120392186 13.5469500314624 13.5313816591058 13.5159117563174 13.5005451308356 13.4852865623876 13.4701408011214 13.4551125660374 13.4402065434301 13.4254273853322 13.4107797079739 13.3962680902468 13.3818970721742 13.3676711533959 13.3535947916635 13.3396724013440 13.3259083519375 13.3123069666056 13.2988725207168 13.2856092404000 13.2725213011175 13.2596128262497 13.2468878856979 13.2343504945014 13.2220046114715 13.2098541378443 13.1979029159503 13.1861547279022 13.1746132943024 13.1632822729698 13.1521652576866 13.1412657769651 13.1305872928372 13.1201331996620 13.1099068229606 13.0999114182666 13.0901501700070 13.0806261903993 13.0713425183767 13.0623021185362 13.0535078801118 13.0449626159700 13.0366690616337 13.0286298743306 13.0208476320668 13.0133248327274 13.0060638932015 12.9990671485396 12.9923368511315 12.9858751699154 12.9796841896135 12.9737659099982 12.9681222451816 12.9627550229355 12.9576659840434 12.9528567816792 12.9483289808098 12.9440840576357 12.9401233990557 12.9364483021638 12.9330599737707 12.9299595299624 12.9271479956879 12.9246263043721 12.9223952975581 12.9204557245907 12.9188082423245 12.9174534148548 12.9163917132896 12.9156235155548 12.9151491062182 12.9149686763524 12.9150823234368 12.9154900512666 12.9161917699251 12.9171872957593 12.9184763514032 12.9200585658299 12.9219334744282 12.9241005191174 12.9265590484958 12.9293083180124 12.9323474901765 12.9356756347978 12.9392917292522 12.9431946587923 12.9473832168736 12.9518561055243 12.9566119357395 12.9616492279108 12.9669664122791 12.9725618294380 12.9784337308353 12.9845802793507 12.9909995498508 12.9976895298280 13.0046481200267 13.0118731351267 13.0193623044521 13.0271132726918 13.0351236006808 13.0433907661925 13.0519121647505 13.0606851104991 13.0697068370781 13.0789744985314 13.0884851702595 13.0982358499785 13.1082234587194 13.1184448418503 13.1288967701396 13.1395759408239 13.1504789787135 13.1616024373362 13.1729428000841 13.1844964814084 13.1962598280131 13.2082291201087 13.2204005726553 13.2327703366479 13.2453345004322 13.2580890910128 13.2710300754248 13.2841533620962 13.2974548022404 13.3109301912800 13.3245752702647 13.3383857273461 13.3523571992343 13.3664852727007 13.3807654860818 13.3951933308069 13.4097642529391 13.4244736547399 13.4393168962348 13.4542892968054 13.4693861367984 13.4846026591229 13.4999340708945 13.5153755450661 13.5309222220742 13.5465692114999 13.5623115937494 13.5781444217105 13.5940627224477 13.6100614988937 13.6261357315544 13.6422803801874 13.6584903855385 13.6747606710379 13.6910861445099 13.7074616999035 13.7238822189986 13.7403425731291 13.7568376249003 13.7733622299082 13.7899112384507 13.8064794972401 13.8230618511229 13.8396531447745 13.8562482243971 13.8728419394333 13.8894291442299 13.9060046997435 13.9225634751965 13.9391003497519 13.9556102141665 13.9720879724367 13.9885285434412 14.0049268625519 14.0212778832575 14.0375765787654 14.0538179435760 14.0699969950699 14.0861087750570 14.1021483513249 14.1181108191649 14.1339913028899 14.1497849573210 14.1654869692747 14.1810925590216 14.1965969817250 14.2119955288717 14.2272835296755 14.2424563524631 14.2575094060427 14.2724381410466 14.2872380512633 14.3019046749357 14.3164335960550 14.3308204456151 14.3450609028568 14.3591506964873 14.3730856058796 14.3868614622431 14.4004741497781 14.4139196068046 14.4271938268710 14.4402928598326 14.4532128129170 14.4659498517571 14.4785002014073 14.4908601473288 14.5030260363629 14.5149942776664 14.5267613436364 14.5383237708011 14.5496781606944 14.5608211807028 14.5717495648901 14.5824601147975];
E4rms7 = [15.4511712114141 15.4367582802416 15.4221245098728 15.4072742444630 15.3922118952931 15.3769419396059 15.3614689194187 15.3457974403137 15.3299321702046 15.3138778380860 15.2976392327533 15.2812212015116 15.2646286488501 15.2478665351067 15.2309398751028 15.2138537367650 15.1966132397162 15.1792235538562 15.1616898979160 15.1440175379924 15.1262117860626 15.1082779984836 15.0902215744675 15.0720479545388 15.0537626189708 15.0353710862162 15.0168789112991 14.9982916842097 14.9796150282660 14.9608545984700 14.9420160798422 14.9231051857417 14.9041276561736 14.8850892560719 14.8659957735808 14.8468530183181 14.8276668196156 14.8084430247639 14.7891874972307 14.7699061148749 14.7506047681431 14.7312893582706 14.7119657954472 14.6926399970025 14.6733178855673 14.6540053872174 14.6347084296364 14.6154329402501 14.5961848443633 14.5769700632940 14.5577945124904 14.5386640996728 14.5195847229285 14.5005622688567 14.4816026106649 14.4627116062888 14.4438950965194 14.4251589031094 14.4065088268923 14.3879506459084 14.3694901135267 14.3511329565740 14.3328848734633 14.3147515323448 14.2967385692362 14.2788515861791 14.2610961494152 14.2434777875249 14.2260019896311 14.2086742035799 14.1914998341405 14.1744842412045 14.1576327380325 14.1409505894779 14.1244430102362 14.1081151631142 14.0919721573201 14.0760190467422 14.0602608282831 14.0447024401841 14.0293487603718 14.0142046048366 13.9992747260158 13.9845638112060 13.9700764809925 13.9558172877077 13.9417907138943 13.9280011708213 13.9144529969824 13.9011504566602 13.8880977384859 13.8752989540358 13.8627581364628 13.8504792391210 13.8384661342611 13.8267226117206 13.8152523776476 13.8040590532750 13.7931461736841 13.7825171866354 13.7721754514005 13.7621242376483 13.7523667243359 13.7429059986442 13.7337450549538 13.7248867938308 13.7163340210561 13.7080894466838 13.7001556841425 13.6925352493474 13.6852305598606 13.6782439340898 13.6715775904935 13.6652336468511 13.6592141195554 13.6535209229199 13.6481558685513 13.6431206647379 13.6384169158664 13.6340461219018 13.6300096778553 13.6263088733445 13.6229448921305 13.6199188117382 13.6172316030772 13.6148841301176 13.6128771495930 13.6112113107424 13.6098871550818 13.6089051162213 13.6082655197037 13.6079685828899 13.6080144148764 13.6084030164443 13.6091342800529 13.6102079898584 13.6116238217706 13.6133813435534 13.6154800149512 13.6179191878493 13.6206981064870 13.6238159076750 13.6272716210830 13.6310641695367 13.6351923693574 13.6396549307457 13.6444504581871 13.6495774508941 13.6550343032942 13.6608193055451 13.6669306440758 13.6733664021755 13.6801245606142 13.6872029982928 13.6945994929254 13.7023117217634 13.7103372623471 13.7186735932938 13.7273180951129 13.7362680510632 13.7455206480386 13.7550729774812 13.7649220363333 13.7750647280212 13.7854978634665 13.7962181621319 13.8072222530951 13.8185066761597 13.8300678829877 13.8419022382696 13.8540060209192 13.8663754253024 13.8790065624905 13.8918954615435 13.9050380708236 13.9184302593338 13.9320678180857 13.9459464614904 13.9600618287828 13.9744094854631 13.9889849247715 14.0037835691811 14.0188007719205 14.0340318185158 14.0494719283587 14.0651162562960 14.0809598942404 14.0969978728041 14.1132251629517 14.1296366776749 14.1462272736840 14.1629917531206 14.1799248652881 14.1970213083974 14.2142757313324 14.2316827354276 14.2492368762663 14.2669326654874 14.2847645726101 14.3027270268693 14.3208144190659 14.3390211034248 14.3573413994646 14.3757695938808 14.3942999424355 14.4129266718512 14.4316439817193 14.4504460464115 14.4693270169983 14.4882810231674 14.5073021751552 14.5263845656774 14.5455222718603 14.5647093571709 14.5839398733595 14.6032078623956 14.6225073583960 14.6418323895609 14.6611769801107 14.6805351522062 14.6999009278742 14.7192683309352 14.7386313889013 14.7579841349012 14.7773206095683 14.7966348629383 14.8159209563357 14.8351729642439 14.8543849761712 14.8735510985107 14.8926654563769 14.9117221954402 14.9307154837455 14.9496395135114 14.9684885029308 14.9872566979386 15.0059383739780 15.0245278377430 15.0430194289095 15.0614075218397 15.0796865272917 15.0978508940732 15.1158951107279 15.1338137071428 15.1516012561948 15.1692523753303 15.1867617281526 15.2041240259847 15.2213340293885 15.2383865496982 15.2552764505106 15.2719986491402 15.2885481180873 15.3049198864527 15.3211090413360 15.3371107292277 15.3529201573523 15.3685325950015 15.3839433748375 15.3991478941879 15.4141416162898 15.4289200715215 15.4434788586237 15.4578136458686 15.4719201722296 15.4857942484951 15.4994317583992 15.5128286596840 15.5259809851602 15.5388848437482 15.5515364214603 15.5639319824029 15.5760678697184 15.5879405065112 15.5995463967607 15.6108821261769 15.6219443630764 15.6327298591870 15.6432354504609 15.6534580578424 15.6633946880176 15.6730424341365 15.6823984765176 15.6914600833115 15.7002246111538 15.7086895057970 15.7168523026871 15.7247106275581 15.7322621969714 15.7395048188364 15.7464363929124 15.7530549112982 15.7593584588591 15.7653452136665 15.7710134474005 15.7763615257375 15.7813879086782 15.7860911509113 15.7904699021052 15.7945229071873 15.7982490066214 15.8016471366329 15.8047163294291 15.8074557133904 15.8098645132420 15.8119420501950 15.8136877420705 15.8151011034111 15.8161817455486 15.8169293766560 15.8173438018023 15.8174249229366 15.8171727389034 15.8165873453904 15.8156689348838 15.8144177965883 15.8128343163271 15.8109189764292 15.8086723555720 15.8060951286283 15.8031880664821 15.7999520358065 15.7963879988493 15.7924970131703 15.7882802313732 15.7837389008071 15.7788743632556 15.7736880545885 15.7681815044092 15.7623563356705 15.7562142642654 15.7497570986066 15.7429867391768 15.7359051780580 15.7285144984413 15.7208168741092 15.7128145689054 15.7045099361698 15.6959054181676 15.6870035454791 15.6778069363779 15.6683182961847 15.6585404166001 15.6484761750100 15.6381285337736 15.6275005394876 15.6165953222287 15.6054160947697 15.5939661517794 15.5822488689939 15.5702677023704 15.5580261872124 15.5455279372806 15.5327766438717 15.5197760748826 15.5065300738462 15.4930425589483 15.4793175220194 15.4653590275052 15.4511712114141];
E4rms8 = [13.9816690167433 13.9672069009457 13.9529748300957 13.9389772750090 13.9252186335002 13.9117032289252 13.8984353087482 13.8854190431354 13.8726585235725 13.8601577615151 13.8479206870580 13.8359511476434 13.8242529067856 13.8128296428341 13.8016849477594 13.7908223259743 13.7802451931773 13.7699568752344 13.7599606070867 13.7502595316896 13.7408566989820 13.7317550648924 13.7229574903707 13.7144667404544 13.7062854833646 13.6984162896447 13.6908616313153 13.6836238810786 13.6767053115419 13.6701080944859 13.6638343001596 13.6578858966124 13.6522647490604 13.6469726192805 13.6420111650474 13.6373819396041 13.6330863911560 13.6291258624167 13.6255015901741 13.6222147048995 13.6192662303848 13.6166570834281 13.6143880735325 13.6124599026642 13.6108731650334 13.6096283469016 13.6087258264467 13.6081658736447 13.6079486501941 13.6080742094773 13.6085424965468 13.6093533481717 13.6105064928804 13.6120015510854 13.6138380352020 13.6160153498224 13.6185327919333 13.6213895511484 13.6245847099863 13.6281172441879 13.6319860230617 13.6361898098680 13.6407272622323 13.6455969326094 13.6507972687572 13.6563266142635 13.6621832091166 13.6683651902695 13.6748705922884 13.6816973480045 13.6888432892080 13.6963061473638 13.7040835543929 13.7121730434536 13.7205720497665 13.7292779114769 13.7382878705548 13.7475990737017 13.7572085733243 13.7671133285163 13.7773102060726 13.7877959815486 13.7985673403356 13.8096208787740 13.8209531052939 13.8325604415918 13.8444392238198 13.8565857038321 13.8689960504210 13.8816663506204 13.8945926110087 13.9077707590503 13.9211966444725 13.9348660406360 13.9487746459734 13.9629180854218 13.9772919118850 13.9918916077414 14.0067125863367 14.0217501935362 14.0369997092719 14.0524563491352 14.0681152659644 14.0839715514660 14.1000202378662 14.1162562995588 14.1326746547829 14.1492701673149 14.1660376481905 14.1829718574167 14.2000675057172 14.2173192562989 14.2347217266007 14.2522694900924 14.2699570780713 14.2877789814505 14.3057296525923 14.3238035071306 14.3419949257966 14.3602982562840 14.3787078150706 14.3972178893103 14.4158227386663 14.4345165972043 14.4532936752505 14.4721481612759 14.4910742237749 14.5100660131493 14.5291176635873 14.5482232949556 14.5673770146737 14.5865729196020 14.6058050979229 14.6250676310156 14.6443545953377 14.6636600642913 14.6829781100879 14.7023028056152 14.7216282262899 14.7409484519004 14.7602575684614 14.7795496700271 14.7988188605329 14.8180592556006 14.8372649843401 14.8564301911529 14.8755490375070 14.8946157037041 14.9136243906462 14.9325693215800 14.9514447438198 14.9702449304696 14.9889641821293 15.0075968285807 15.0261372304541 15.0445797808911 15.0629189071816 15.0811490723886 15.0992647769484 15.1172605602647 15.1351310022776 15.1528707250114 15.1704743941079 15.1879367203438 15.2052524611214 15.2224164219442 15.2394234578697 15.2562684749479 15.2729464316304 15.2894523401658 15.3057812679709 15.3219283389831 15.3378887349894 15.3536576969329 15.3692305262014 15.3846025858918 15.3997693020526 15.4147261649027 15.4294687300340 15.4439926195834 15.4582935233902 15.4723672001247 15.4862094783996 15.4998162578547 15.5131835102216 15.5263072803651 15.5391836873014 15.5518089251937 15.5641792643262 15.5762910520549 15.5881407137341 15.5997247536227 15.6110397557670 15.6220823848600 15.6328493870794 15.6433375909002 15.6535439078910 15.6634653334797 15.6730989477030 15.6824419159294 15.6914914895647 15.7002450067294 15.7086998929148 15.7168536616228 15.7247039149787 15.7322483443177 15.7394847307585 15.7464109457483 15.7530249515899 15.7593248019390 15.7653086422898 15.7709747104357 15.7763213369037 15.7813469453642 15.7860500530329 15.7904292710432 15.7944833047871 15.7982109542469 15.8016111143074 15.8046827750323 15.8074250219277 15.8098370361929 15.8119180949214 15.8136675713175 15.8150849348579 15.8161697514506 15.8169216835691 15.8173404903564 15.8174260277140 15.8171782483740 15.8165972019353 15.8156830348897 15.8144359906208 15.8128564093762 15.8109447282321 15.8087014810172 15.8061272982278 15.8032229069135 15.7999891305437 15.7964268888423 15.7925371976240 15.7883211685651 15.7837800090099 15.7789150216853 15.7737276044590 15.7682192500208 15.7623915455734 15.7562461724930 15.7497849059447 15.7430096145102 15.7359222597707 15.7285248958506 15.7208196689764 15.7128088169772 15.7044946687698 15.6958796438368 15.6869662516530 15.6777570911024 15.6682548498648 15.6584623037935 15.6483823162409 15.6380178373717 15.6273719034667 15.6164476361715 15.6052482417488 15.5937770102717 15.5820373148362 15.5700326107038 15.5577664344422 15.5452424030483 15.5324642130097 15.5194356393880 15.5061605348414 15.4926428286310 15.4788865256143 15.4648957051820 15.4506745202137 15.4362271959644 15.4215580289571 15.4066713858330 15.3915717021813 15.3762634813426 15.3607512931956 15.3450397729021 15.3291336196429 15.3130375953309 15.2967565232731 15.2802952868470 15.2636588281265 15.2468521464873 15.2298802971960 15.2127483899870 15.1954615875801 15.1780251042106 15.1604442041252 15.1427242000637 15.1248704516866 15.1068883640363 15.0887833859332 15.0705610083614 15.0522267628550 15.0337862198374 15.0152449869597 14.9966087074140 14.9778830582331 14.9590737485653 14.9401865179372 14.9212271345104 14.9022013933026 14.8831151144005 14.8639741411804 14.8447843384716 14.8255515907552 14.8062818003095 14.7869808853702 14.7676547782675 14.7483094235556 14.7289507761389 14.7095847993711 14.6902174631686 14.6708547421075 14.6515026135015 14.6321670554979 14.6128540451468 14.5935695564794 14.5743195585759 14.5551100136382 14.5359468750476 14.5168360854366 14.4977835747512 14.4787952583118 14.4598770348844 14.4410347847473 14.4222743677616 14.4036016214479 14.3850223590627 14.3665423676882 14.3481674063183 14.3299032039634 14.3117554577514 14.2937298310439 14.2758319515612 14.2580674095185 14.2404417557683 14.2229604999602 14.2056291087099 14.1884530037848 14.1714375602979 14.1545881049242 14.1379099141258 14.1214082123979 14.1050881705288 14.0889549038810 14.0730134706865 14.0572688703651 14.0417260418594 14.0263898619926 14.0112651438457 13.9963566351569 13.9816690167433];
E4rms9 = [13.6524395594891 13.6659883518785 13.6797759433918 13.6937980062960 13.7080501406124 13.7225278755739 13.7372266711096 13.7521419193531 13.7672689461698 13.7826030127117 13.7981393169849 13.8138729954479 13.8297991246217 13.8459127227230 13.8622087513143 13.8786821169765 13.8953276729917 13.9121402210493 13.9291145129652 13.9462452524144 13.9635270966814 13.9809546584232 13.9985225074490 14.0162251725015 14.0340571430622 14.0520128711626 14.0700867732034 14.0882732317864 14.1065665975539 14.1249611910328 14.1434513044924 14.1620312038042 14.1806951303069 14.1994373026741 14.2182519187926 14.2371331576425 14.2560751811683 14.2750721361699 14.2941181561804 14.3132073633493 14.3323338703233 14.3514917821324 14.3706751980569 14.3898782135141 14.4090949219294 14.4283194165913 14.4475457925277 14.4667681483573 14.4859805881353 14.5051772232035 14.5243521740143 14.5434995719740 14.5626135612360 14.5816883005331 14.6007179649561 14.6196967477413 14.6386188620551 14.6574785427445 14.6762700480877 14.6949876615309 14.7136256934090 14.7321784826525 14.7506403984725 14.7690058420501 14.7872692481813 14.8054250869235 14.8234678652376 14.8413921285683 14.8591924624576 14.8768634941128 14.8943998939584 14.9117963771603 14.9290477051642 14.9461486871773 14.9630941816438 14.9798790977021 14.9964983966302 15.0129470932391 15.0292202572865 15.0453130148404 15.0612205496242 15.0769381043551 15.0924609820422 15.1077845472746 15.1229042274811 15.1378155141752 15.1525139641603 15.1669952007413 15.1812549148735 15.1952888663287 15.2090928848102 15.2226628710560 15.2359947979281 15.2490847114446 15.2619287318340 15.2745230545321 15.2868639511622 15.2989477705153 15.3107709394644 15.3223299638969 15.3336214295931 15.3446420031075 15.3553884325988 15.3658575486517 15.3760462650867 15.3859515797212 15.3955705751230 15.4049004193334 15.4139383665846 15.4226817579640 15.4311280220758 15.4392746756864 15.4471193243134 15.4546596628306 15.4618934760340 15.4688186391670 15.4754331184544 15.4817349715957 15.4877223482281 15.4933934903972 15.4987467329554 15.5037805040010 15.5084933252302 15.5128838123179 15.5169506752436 15.5206927186085 15.5241088419293 15.5271980399091 15.5299594026789 15.5323921160312 15.5344954616150 15.5362688171179 15.5377116564270 15.5388235497591 15.5396041637794 15.5400532616891 15.5401707032918 15.5399564450425 15.5394105400728 15.5385331381823 15.5373244858355 15.5357849260957 15.5339148985760 15.5317149393445 15.5291856808079 15.5263278515863 15.5231422763525 15.5196298756455 15.5157916656787 15.5116287581152 15.5071423598071 15.5023337725342 15.4972043927076 15.4917557110555 15.4859893122738 15.4799068746705 15.4735101697746 15.4668010619263 15.4597815078417 15.4524535561552 15.4448193469389 15.4368811111967 15.4286411703325 15.4201019356008 15.4112659075257 15.4021356753034 15.3927139161712 15.3830033947649 15.3730069624378 15.3627275565655 15.3521681998222 15.3413319994347 15.3302221464077 15.3188419147306 15.3071946605529 15.2952838213388 15.2831129149996 15.2706855389967 15.2580053694204 15.2450761600460 15.2319017413671 15.2184860195988 15.2048329756572 15.1909466641225 15.1768312121658 15.1624908184590 15.1479297520565 15.1331523512581 15.1181630224401 15.1029662388664 15.0875665394773 15.0719685276498 15.0561768699399 15.0401962947962 15.0240315912538 15.0076876075987 14.9911692500231 14.9744814812403 14.9576293190917 14.9406178351224 14.9234521531417 14.9061374477553 14.8886789428807 14.8710819102428 14.8533516678452 14.8354935784213 14.8175130478692 14.7994155236635 14.7812064932534 14.7628914824300 14.7444760536924 14.7259658045850 14.7073663660218 14.6886834005823 14.6699226008145 14.6510896875082 14.6321904079438 14.6132305341433 14.5942158611057 14.5751522050154 14.5560454014511 14.5369013035851 14.5177257803467 14.4985247146154 14.4793040013621 14.4600695458105 14.4408272615766 14.4215830688020 14.4023428922778 14.3831126595728 14.3638982991396 14.3447057384289 14.3255409019919 14.3064097095789 14.2873180742451 14.2682719004381 14.2492770820989 14.2303395007542 14.2114650236107 14.1926595016470 14.1739287677235 14.1552786346610 14.1367148933738 14.1182433109453 14.0998696287685 14.0815995606448 14.0634387909197 14.0453929726177 14.0274677255685 14.0096686345728 13.9920012475609 13.9744710737462 13.9570835818272 13.9398441981702 13.9227583050117 13.9058312386923 13.8890682878795 13.8724746918183 13.8560556385945 13.8398162634272 13.8237616469550 13.8078968135512 13.7922267296717 13.7767563021951 13.7614903768085 13.7464337363848 13.7315910994191 13.7169671184486 13.7025663785155 13.6883933956595 13.6744526154017 13.6607484112956 13.6472850834701 13.6340668572067 13.6210978815527 13.6083822279348 13.5959238888354 13.5837267764602 13.5717947214589 13.5601314716595 13.5487406908356 13.5376259574999 13.5267907637347 13.5162385140371 13.5059725242068 13.4959960202673 13.4863121373908 13.4769239188912 13.4678343152214 13.4590461830073 13.4505622841177 13.4423852847765 13.4345177546744 13.4269621661435 13.4197208933550 13.4127962115543 13.4061902963013 13.3999052227926 13.3939429651753 13.3883053959060 13.3829942851596 13.3780113002463 13.3733580050801 13.3690358596757 13.3650462196793 13.3613903359322 13.3580693540686 13.3550843141589 13.3524361503683 13.3501256906569 13.3481536565346 13.3465206628083 13.3452272174119 13.3442737212306 13.3436604679871 13.3433876441460 13.3434553288624 13.3438634939676 13.3446120039706 13.3457006161200 13.3471289804891 13.3488966400845 13.3510030310118 13.3534474826582 13.3562292179179 13.3593473534493 13.3628008999729 13.3665887625886 13.3707097411470 13.3751625306425 13.3799457216395 13.3850578007405 13.3904971510844 13.3962620528796 13.4023506839670 13.4087611204213 13.4154913371853 13.4225392087340 13.4299025097780 13.4375789159940 13.4455660047869 13.4538612560941 13.4624620532135 13.4713656836650 13.4805693400845 13.4900701211521 13.4998650325498 13.5099509879470 13.5203248100240 13.5309832315191 13.5419228963092 13.5531403605183 13.5646320936625 13.5763944798103 13.5884238187875 13.6007163273973 13.6132681406769 13.6260753131775 13.6391338202714 13.6524395594891];
E4rms10 = [14.3507897890617 14.3634866625743 14.3759911781728 14.3882996504154 14.4004084534429 14.4123140218908 14.4240128517840 14.4355015014166 14.4467765922080 14.4578348095474 14.4686729036110 14.4792876901718 14.4896760513828 14.4998349365448 14.5097613628564 14.5194524161489 14.5289052515964 14.5381170944148 14.5470852405406 14.5558070572891 14.5642799839984 14.5725015326550 14.5804692885055 14.5881809106386 14.5956341325637 14.6028267627666 14.6097566852451 14.6164218600311 14.6228203236944 14.6289501898263 14.6348096495121 14.6403969717825 14.6457105040482 14.6507486725135 14.6555099825820 14.6599930192422 14.6641964474242 14.6681190123599 14.6717595399098 14.6751169368788 14.6781901913132 14.6809783727883 14.6834806326596 14.6856962043222 14.6876244034390 14.6892646281406 14.6906163592345 14.6916791603774 14.6924526782315 14.6929366426142 14.6931308666117 14.6930352467046 14.6926497628314 14.6919744784879 14.6910095407596 14.6897551803609 14.6882117116625 14.6863795326785 14.6842591250490 14.6818510540028 14.6791559682995 14.6761746001537 14.6729077651326 14.6693563620544 14.6655213728392 14.6614038623581 14.6570049782752 14.6523259508244 14.6473680926188 14.6421327984105 14.6366215448344 14.6308358901206 14.6247774738212 14.6184480164781 14.6118493192828 14.6049832637192 14.5978518111902 14.5904570025945 14.5828009579206 14.5748858757916 14.5667140329896 14.5582877839734 14.5496095603560 14.5406818703694 14.5315072983025 14.5220885039233 14.5124282218596 14.5025292609877 14.4923945037557 14.4820269055249 14.4714294938588 14.4606053678004 14.4495576971340 14.4382897215910 14.4268047500748 14.4151061598303 14.4031973955942 14.3910819687455 14.3787634563895 14.3662455004591 14.3535318067620 14.3406261440290 14.3275323429124 14.3142542949745 14.3007959516631 14.2871613232391 14.2733544776969 14.2593795396521 14.2452406892241 14.2309421608671 14.2164882421977 14.2018832728050 14.1871316430084 14.1722377926288 14.1572062097220 14.1420414292761 14.1267480319138 14.1113306425591 14.0957939290753 14.0801426009104 14.0643814076780 14.0485151377742 14.0325486169186 14.0164867067213 14.0003343032036 13.9840963353117 13.9677777634128 13.9513835777737 13.9349187970168 13.9183884665715 13.9017976570987 13.8851514629057 13.8684550003499 13.8517134062206 13.8349318361168 13.8181154628062 13.8012694745725 13.7843990735554 13.7675094740805 13.7506059009676 13.7336935878557 13.7167777754839 13.6998637099981 13.6829566412322 13.6660618209832 13.6491845012919 13.6323299327095 13.6155033625566 13.5987100331942 13.5819551802834 13.5652440310333 13.5485818024665 13.5319736996734 13.5154249140737 13.4989406216679 13.4825259813056 13.4661861329474 13.4499261959338 13.4337512672584 13.4176664198456 13.4016767008374 13.3857871298843 13.3700026974436 13.3543283630896 13.3387690538294 13.3233296624326 13.3080150457663 13.2928300231497 13.2777793747119 13.2628678397689 13.2481001152129 13.2334808539157 13.2190146631461 13.2047061030058 13.1905596848782 13.1765798698960 13.1627710674292 13.1491376335879 13.1356838697436 13.1224140210722 13.1093322751192 13.0964427603804 13.0837495449049 13.0712566349267 13.0589679735090 13.0468874392176 13.0350188448146 13.0233659359801 13.0119323900521 13.0007218147954 12.9897377471948 12.9789836522730 12.9684629219366 12.9581788738458 12.9481347503126 12.9383337172221 12.9287788629913 12.9194731975423 12.9104196513144 12.9016210742974 12.8930802350995 12.8847998200375 12.8767824322613 12.8690305909071 12.8615467302782 12.8543331990555 12.8473922595405 12.8407260869260 12.8343367686016 12.8282263034786 12.8223966013614 12.8168494823389 12.8115866762107 12.8066098219356 12.8019204671292 12.7975200675821 12.7934099868008 12.7895914955953 12.7860657716950 12.7828338993859 12.7798968691911 12.7772555775831 12.7749108267096 12.7728633241831 12.7711136828685 12.7696624207271 12.7685099606803 12.7676566305090 12.7671026627800 12.7668481948159 12.7668932686817 12.7672378312152 12.7678817340817 12.7688247338615 12.7700664921762 12.7716065758339 12.7734444570176 12.7755795134986 12.7780110288822 12.7807381928824 12.7837601016404 12.7870757580465 12.7906840721353 12.7945838614597 12.7987738515463 12.8032526763411 12.8080188787112 12.8130709109693 12.8184071354147 12.8240258249298 12.8299251635917 12.8361032472998 12.8425580844652 12.8492875967019 12.8562896195522 12.8635619032572 12.8711021135328 12.8789078323851 12.8869765589515 12.8953057103798 12.9038926227139 12.9127345518171 12.9218286743335 12.9311720886521 12.9407618159209 12.9505948010584 12.9606679138269 12.9709779498939 12.9815216319389 12.9922956107883 13.0032964665460 13.0145207097857 13.0259647827368 13.0376250605003 13.0494978522967 13.0615794027099 13.0738658929884 13.0863534423300 13.0990381092130 13.1119158927318 13.1249827339553 13.1382345173023 13.1516670719413 13.1652761731921 13.1790575439583 13.1930068561767 13.2071197322566 13.2213917465703 13.2358184269306 13.2503952560872 13.2651176732401 13.2799810755748 13.2949808197770 13.3101122235874 13.3253705673579 13.3407510956218 13.3562490186435 13.3718595140270 13.3875777282880 13.4033987784452 13.4193177536336 13.4353297166973 13.4514297058045 13.4676127360604 13.4838738011238 13.5002078748234 13.5166099127770 13.5330748540225 13.5495976226302 13.5661731293204 13.5827962731036 13.5994619428732 13.6161650190466 13.6329003751613 13.6496628794961 13.6664473966705 13.6832487892496 13.7000619193460 13.7168816501950 13.7337028477515 13.7505203822677 13.7673291298492 13.7841239740293 13.8008998073140 13.8176515327269 13.8343740653412 13.8510623338089 13.8677112818618 13.8843158698234 13.9008710760951 13.9173718986296 13.9338133563996 13.9501904908505 13.9664983673392 13.9827320765569 13.9988867359415 14.0149574910759 14.0309395170671 14.0468280199177 14.0626182378759 14.0783054427672 14.0938849413220 14.1093520764775 14.1247022286646 14.1399308170785 14.1550333009359 14.1700051807132 14.1848419993629 14.1995393435225 14.2140928446963 14.2284981804251 14.2427510754344 14.2568473027734 14.2707826849217 14.2845530948928 14.2981544573079 14.3115827494594 14.3248340023517 14.3379043017249 14.3507897890617];
E4rms11 = [13.9971000220208 13.9857565379317 13.9742267600114 13.9625140811767 13.9506219511249 13.9385538754603 13.9263134148054 13.9139041838925 13.9013298506326 13.8885941351730 13.8757008089234 13.8626536935786 13.8494566601081 13.8361136277347 13.8226285628905 13.8090054781605 13.7952484311972 13.7813615236268 13.7673488999326 13.7532147463176 13.7389632895534 13.7245987958103 13.7101255694723 13.6955479519248 13.6808703203357 13.6660970864189 13.6512326951733 13.6362816236135 13.6212483794799 13.6061374999328 13.5909535502356 13.5757011224194 13.5603848339333 13.5450093262765 13.5295792636246 13.5140993314414 13.4985742350637 13.4830086982948 13.4674074619697 13.4517752825143 13.4361169304912 13.4204371891435 13.4047408529081 13.3890327259448 13.3733176206426 13.3576003561063 13.3418857566612 13.3261786503309 13.3104838673105 13.2948062384424 13.2791505936708 13.2635217605174 13.2479245625123 13.2323638176705 13.2168443369199 13.2013709225501 13.1859483666652 13.1705814496146 13.1552749384371 13.1400335853034 13.1248621259576 13.1097652781615 13.0947477401340 13.0798141890121 13.0649692792888 13.0502176412724 13.0355638795618 13.0210125714862 13.0065682655992 12.9922354801526 12.9780187015825 12.9639223829938 12.9499509426850 12.9361087626491 12.9224001870959 12.9088295209895 12.8954010286022 12.8821189320553 12.8689874099072 12.8560105957321 12.8431925767129 12.8305373922630 12.8180490326488 12.8057314376349 12.7935884951446 12.7816240399408 12.7698418523116 12.7582456567984 12.7468391209070 12.7356258538750 12.7246094054319 12.7137932645909 12.7031808584682 12.6927755510945 12.6825806422869 12.6725993665137 12.6628348917870 12.6532903185979 12.6439686788388 12.6348729347844 12.6260059780709 12.6173706287209 12.6089696341703 12.6008056683320 12.5928813306951 12.5851991454269 12.5777615605163 12.5705709469347 12.5636295978390 12.5569397277767 12.5505034719340 12.5443228854169 12.5383999425321 12.5327365361278 12.5273344769478 12.5221954930005 12.5173212289802 12.5127132457031 12.5083730195648 12.5043019420513 12.5005013192384 12.4969723713711 12.4937162324194 12.4907339497049 12.4880264835301 12.4855947068483 12.4834394049611 12.4815612752463 12.4799609269083 12.4786388807698 12.4775955690809 12.4768313353643 12.4763464342915 12.4761410315817 12.4762152039395 12.4765689390152 12.4772021353950 12.4781146026258 12.4793060612686 12.4807761429717 12.4825243905969 12.4845502583424 12.4868531119275 12.4894322287877 12.4922867983010 12.4954159220547 12.4988186141305 12.5024938014185 12.5064403239726 12.5106569353892 12.5151423032025 12.5198950093291 12.5249135505304 12.5301963389099 12.5357417024273 12.5415478854576 12.5476130493681 12.5539352731269 12.5605125539393 12.5673428079129 12.5744238707515 12.5817534984751 12.5893293681659 12.5971490787477 12.6052101517853 12.6135100323172 12.6220460897078 12.6308156185364 12.6398158394999 12.6490439003508 12.6584968768570 12.6681717737887 12.6780655259272 12.6881749991015 12.6984969912469 12.7090282334876 12.7197653912467 12.7307050653736 12.7418437932974 12.7531780502007 12.7647042502208 12.7764187476627 12.7883178382402 12.8003977603376 12.8126546962873 12.8250847736683 12.8376840666251 12.8504485972062 12.8633743367151 12.8764572070822 12.8896930822554 12.9030777896022 12.9166071113312 12.9302767859266 12.9440825095975 12.9580199377383 12.9720846864097 12.9862723338220 13.0005784218389 13.0149984574859 13.0295279144745 13.0441622347303 13.0588968299359 13.0737270830793 13.0886483500092 13.1036559609989 13.1187452223155 13.1339114177946 13.1491498104232 13.1644556439166 13.1798241443132 13.1952505215615 13.2107299711139 13.2262576755121 13.2418288059905 13.2574385240697 13.2730819831406 13.2887543300606 13.3044507067461 13.3201662517510 13.3358961018528 13.3516353936349 13.3673792650458 13.3831228569845 13.3988613148460 13.4145897900836 13.4303034417528 13.4459974380499 13.4616669578374 13.4773071921718 13.4929133458053 13.5084806386912 13.5240043074679 13.5394796069347 13.5549018115215 13.5702662167330 13.5855681405919 13.6008029250620 13.6159659374571 13.6310525718343 13.6460582503841 13.6609784247769 13.6758085775382 13.6905442233515 13.7051809103997 13.7197142216454 13.7341397761224 13.7484532302002 13.7626502788173 13.7767266567223 13.7906781396822 13.8045005456560 13.8181897359855 13.8317416165365 13.8451521388264 13.8584173011505 13.8715331496634 13.8844957794528 13.8973013355901 13.9099460141736 13.9224260633277 13.9347377841938 13.9468775319127 13.9588417165620 13.9706268040972 13.9822293172433 13.9936458364043 14.0048730005130 14.0159075078820 14.0267461170393 14.0373856475128 14.0478229806358 14.0580550603001 14.0680788936984 14.0778915520557 14.0874901713140 14.0968719528362 14.1060341640492 14.1149741390978 14.1236892794597 14.1321770545486 14.1404350022932 14.1484607297061 14.1562519134162 14.1638063001972 14.1711217074780 14.1781960238076 14.1850272093399 14.1916132962707 14.1979523892625 14.2040426658546 14.2098823768671 14.2154698467495 14.2208034739465 14.2258817312338 14.2307031660413 14.2352664007291 14.2395701328988 14.2436131356406 14.2473942577776 14.2509124241090 14.2541666356092 14.2571559696288 14.2598795800718 14.2623366975572 14.2645266295615 14.2664487605434 14.2681025520662 14.2694875428811 14.2706033490021 14.2714496637870 14.2720262579559 14.2723329796449 14.2723697544008 14.2721365851888 14.2716335523654 14.2708608136468 14.2698186040625 14.2685072358723 14.2669270984960 14.2650786584138 14.2629624590383 14.2605791205947 14.2579293399670 14.2550138905343 14.2518336219904 14.2483894601527 14.2446824067390 14.2407135391482 14.2364840102131 14.2319950479340 14.2272479552044 14.2222441095170 14.2169849626532 14.2114720403516 14.2057069419667 14.1996913401069 14.1934269802551 14.1869156803763 14.1801593305039 14.1731598923060 14.1659193986454 14.1584399531122 14.1507237295418 14.1427729715139 14.1345899918388 14.1261771720210 14.1175369617047 14.1086718781058 14.0995845054206 14.0902774942186 14.0807535608164 14.0710154866369 14.0610661175415 14.0509083631536 14.0405451961549 14.0299796515695 14.0192148260243 14.0082538769926 13.9971000220208];
E4rms12 = [13.2388696530111 13.2274008090908 13.2161732824708 13.2051906284136 13.1944563268119 13.1839737810015 13.1737463166027 13.1637771803880 13.1540695391711 13.1446264787310 13.1354510027503 13.1265460317958 13.1179144023129 13.1095588656561 13.1014820871440 13.0936866451475 13.0861750301992 13.0789496441398 13.0720127992905 13.0653667176521 13.0590135301389 13.0529552758410 13.0471939013201 13.0417312599261 13.0365691111553 13.0317091200379 13.0271528565504 13.0229017950665 13.0189573138356 13.0153206944925 13.0119931216039 13.0089756822429 13.0062693655963 13.0038750626012 13.0017935656209 13.0000255681524 12.9985716645543 12.9974323498264 12.9966080194078 12.9960989690125 12.9959053944973 12.9960273917684 12.9964649567047 12.9972179851361 12.9982862728434 12.9996695155799 13.0013673091514 13.0033791495108 13.0057044328875 13.0083424559596 13.0112924160432 13.0145534113399 13.0181244411796 13.0220044063426 13.0261921093725 13.0306862549428 13.0354854502631 13.0405882054980 13.0459929342317 13.0516979539633 13.0577014866337 13.0640016591856 13.0705965041486 13.0774839602756 13.0846618731848 13.0921279960499 13.0998799903340 13.1079154265134 13.1162317848809 13.1248264563537 13.1336967433156 13.1428398604819 13.1522529358244 13.1619330114940 13.1718770447854 13.1820819091334 13.1925443951465 13.2032612116411 13.2142289867429 13.2254442689904 13.2369035284687 13.2486031579852 13.2605394742572 13.2727087191357 13.2851070608519 13.2977305952958 13.3105753473060 13.3236372720115 13.3369122561606 13.3503961195126 13.3640846162270 13.3779734362867 13.3920582069520 13.4063344942075 13.4207978042703 13.4354435850903 13.4502672278760 13.4652640686619 13.4804293898572 13.4957584218504 13.5112463446050 13.5268882892958 13.5426793399374 13.5586145350424 13.5746888693050 13.5908972952769 13.6072347250709 13.6236960320726 13.6402760526783 13.6569695880198 13.6737714057178 13.6906762416519 13.7076788017097 13.7247737635807 13.7419557785428 13.7592194732404 13.7765594514969 13.7939702961172 13.8114465706876 13.8289828214076 13.8465735788798 13.8642133599594 13.8818966695462 13.8996180024259 13.9173718450787 13.9351526775028 13.9529549750330 13.9707732101576 13.9886018543252 14.0064353797609 14.0242682612626 14.0420949780003 14.0599100153107 14.0777078664759 14.0954830345054 14.1132300339006 14.1309433924110 14.1486176527855 14.1662473745118 14.1838271355326 14.2013515339771 14.2188151898407 14.2362127466900 14.2535388733296 14.2707882654571 14.2879556473156 14.3050357733173 14.3220234296484 14.3389134358738 14.3557006465154 14.3723799525978 14.3889462831989 14.4053946069698 14.4217199336396 14.4379173154889 14.4539818488199 14.4699086753944 14.4856929838545 14.5013300111200 14.5168150437659 14.5321434193790 14.5473105278881 14.5623118128719 14.5771427728497 14.5917989625408 14.6062759941077 14.6205695383673 14.6346753259906 14.6485891486629 14.6623068602321 14.6758243778266 14.6891376829523 14.7022428225606 14.7151359100956 14.7278131265137 14.7402707212788 14.7525050133369 14.7645123920589 14.7762893181631 14.7878323246091 14.7991380174746 14.8102030767954 14.8210242573885 14.8315983896517 14.8419223803323 14.8519932132735 14.8618079501376 14.8713637311062 14.8806577755492 14.8896873826743 14.8984499321542 14.9069428847235 14.9151637827579 14.9231102508249 14.9307799962134 14.9381708094351 14.9452805647123 14.9521072204300 14.9586488195748 14.9649034901442 14.9708694455383 14.9765449849228 14.9819284935751 14.9870184432060 14.9918133922582 14.9963119861827 15.0005129576949 15.0044151270076 15.0080174020465 15.0113187786321 15.0143183406583 15.0170152602370 15.0194087978277 15.0214983023342 15.0232832112026 15.0247630504864 15.0259374348811 15.0268060677554 15.0273687411600 15.0276253358035 15.0275758210221 15.0272202547302 15.0265587833292 15.0255916416352 15.0243191527473 15.0227417279248 15.0208598664321 15.0186741553678 15.0161852694703 15.0133939709168 15.0103011090866 15.0069076203195 15.0032145276459 14.9992229405001 14.9949340544237 14.9903491507337 14.9854695961874 14.9802968426203 14.9748324265661 14.9690779688555 14.9630351742119 14.9567058307949 14.9500918097756 14.9431950648310 14.9360176316812 14.9285616275583 14.9208292506877 14.9128227797422 14.9045445732607 14.8959970690776 14.8871827837179 14.8781043117551 14.8687643251900 14.8591655727777 14.8493108793420 14.8392031450898 14.8288453448795 14.8182405274839 14.8073918148303 14.7963024012346 14.7849755525915 14.7734146055574 14.7616229667279 14.7496041117685 14.7373615845523 14.7248989962479 14.7122200244280 14.6993284121180 14.6862279668472 14.6729225596870 14.6594161242367 14.6457126556341 14.6318162095137 14.6177309009535 14.6034609034154 14.5890104476366 14.5743838205430 14.5595853641026 14.5446194741896 14.5294905994121 14.5142032399278 14.4987619462366 14.4831713179640 14.4674360026099 14.4515606942930 14.4355501324797 14.4194091006666 14.4031424250864 14.3867549733646 14.3702516531682 14.3536374108381 14.3369172300182 14.3200961302326 14.3031791654771 14.2861714227877 14.2690780207947 14.2519041082350 14.2346548624985 14.2173354881169 14.1999512152496 14.1825072981740 14.1650090137329 14.1474616597900 14.1298705536637 14.1122410305509 14.0945784419365 14.0768881539898 14.0591755459639 14.0414460085638 14.0237049423110 14.0059577559205 13.9882098646257 13.9704666885425 13.9527336509865 13.9350161768096 13.9173196907131 13.8996496155650 13.8820113707137 13.8644103702783 13.8468520214601 13.8293417228397 13.8118848626594 13.7944868171259 13.7771529486954 13.7598886043628 13.7426991139517 13.7255897884081 13.7085659180832 13.6916327710348 13.6747955913215 13.6580595972980 13.6414299799222 13.6249119010636 13.6085104918188 13.5922308508276 13.5760780426038 13.5600570958695 13.5441730018977 13.5284307128669 13.5128351402234 13.4973911530497 13.4821035764553 13.4669771899699 13.4520167259522 13.4372268680104 13.4226122494411 13.4081774516788 13.3939270027604 13.3798653758095 13.3659969875328 13.3523261967355 13.3388573028547 13.3255945445141 13.3125420980889 13.2997040763018 13.2870845268293 13.2746874309358 13.2625167021233 13.2505761848065 13.2388696530111];

h = figure();
plot(thetad,E4rms12,'-','LineWidth',3.5);
ylabel('Campo elétrico (kV/cm)')
xlabel('Ângulo (°)')
grid
hold on
plot(thetad,Erms,':','LineWidth',3.5);
hold on
plot(thetad,E_crit,'--','LineWidth',3.5);
xlim([0 360])
% ylim([13.5 15.5])
legend('Original','Raio menor', 'Campo elétrico crítico')
% title('Campo elétrico superficial condutor um')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'solo_caso1_dist','-dpdf','-r0')