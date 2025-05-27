clear all; close all; clc

% Cálculo do campo elétrico superficial do caso 4 do sistema trifásico

e_0 = 8.854*(10^(-12));
r = 15.980*(10^-3); % raio do condutor
n = 12; % número de condutores do sistema
nc = n*2; % número de condutores total (com imagens)
ci = (2*n-1); % número de cargas imagens por condutor


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

x = [-7.975 -7.025 -7.025 -7.975 -0.475 0.475 0.475 -0.475 7.025 7.975 7.975 7.025 -7.975 -7.025 -7.025 -7.975 -0.475 0.475 0.475 -0.475 7.025 7.975 7.975 7.025];
y = [18.45 18.45 17.5 17.5 25.95 25.95 25 25 18.45 18.45 17.5 17.5 -18.45 -18.45 -17.5 -17.5 -25.95 -25.95 -25 -25 -18.45 -18.45 -17.5 -17.5];

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

theta = linspace(0,2*pi,154); %gera a superfície do condutor em 360 pontos

xcj = x(3); % posição x do centro condutor 1 fase a (para trocar o condutor avaliado deve-se mudar este valor)
ycj = y(3); % posição y do centro condutor 1 fase a (para trocar o condutor avaliado deve-se mudar este valor)

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

%vetores para plot da curva referência condutor 1
xp1 = [ 0	,	1.368821293	,	2.737642586	,	3.650190114	,	5.019011407	,	6.3878327	,	7.756653992	,	8.669201521	,	10.03802281	,	11.40684411	,	12.7756654	,	14.14448669	,	15.51330798	,	16.88212928	,	18.25095057	,	19.61977186	,	20.53231939	,	21.90114068	,	22.81368821	,	24.18250951	,	25.5513308	,	26.92015209	,	28.28897338	,	29.65779468	,	31.02661597	,	32.39543726	,	33.76425856	,	35.13307985	,	36.50190114	,	37.87072243	,	39.23954373	,	40.15209125	,	41.06463878	,	42.43346008	,	44.25855513	,	45.62737643	,	47.45247148	,	48.82129278	,	50.19011407	,	52.01520913	,	53.84030418	,	55.20912548	,	56.57794677	,	58.40304183	,	59.77186312	,	61.14068441	,	62.96577947	,	64.79087452	,	66.15969582	,	67.98479087	,	69.80988593	,	71.63498099	,	73.46007605	,	75.74144487	,	77.56653992	,	79.39163498	,	81.6730038	,	83.95437262	,	86.23574144	,	88.97338403	,	91.71102662	,	94.4486692	,	97.18631179	,	100.3802281	,	103.1178707	,	105.3992395	,	108.1368821	,	110.8745247	,	113.1558935	,	116.3498099	,	118.6311787	,	120.9125475	,	123.6501901	,	126.3878327	,	129.581749	,	132.3193916	,	135.0570342	,	138.2509506	,	140.9885932	,	143.7262357	,	146.0076046	,	148.2889734	,	150.5703422	,	152.851711	,	154.6768061	,	156.9581749	,	158.3269962	,	160.1520913	,	161.9771863	,	163.8022814	,	165.6273764	,	167.9087452	,	169.2775665	,	171.1026616	,	172.9277567	,	174.7528517	,	176.5779468	,	177.9467681	,	179.7718631	,	181.1406844	,	182.5095057	,	183.878327	,	185.7034221	,	187.0722433	,	188.4410646	,	189.8098859	,	191.1787072	,	193.0038023	,	194.3726236	,	195.2851711	,	197.1102662	,	198.4790875	,	199.8479087	,	200.7604563	,	202.1292776	,	203.4980989	,	204.8669202	,	206.2357414	,	207.6045627	,	208.973384	,	210.3422053	,	211.7110266	,	213.0798479	,	214.4486692	,	216.2737643	,	218.0988593	,	219.0114068	,	220.3802281	,	222.2053232	,	224.0304183	,	225.3992395	,	226.7680608	,	228.5931559	,	229.9619772	,	231.3307985	,	233.1558935	,	234.5247148	,	235.8935361	,	237.2623574	,	238.6311787	,	240.4562738	,	242.2813688	,	243.6501901	,	245.0190114	,	246.3878327	,	248.2129278	,	249.581749	,	250.9505703	,	252.7756654	,	254.6007605	,	255.9695817	,	257.7946768	,	259.6197719	,	261.4448669	,	263.269962	,	266.0076046	,	268.2889734	,	270.5703422	,	272.851711	,	275.5893536	,	278.3269962	,	281.0646388	,	284.2585551	,	286.9961977	,	290.6463878	,	293.3840304	,	297.4904943	,	300.6844106	,	304.3346008	,	308.4410646	,	311.634981	,	314.8288973	,	317.1102662	,	320.3041825	,	322.5855513	,	324.8669202	,	327.6045627	,	330.3422053	,	332.6235741	,	334.904943	,	337.1863118	,	339.0114068	,	341.2927757	,	343.5741445	,	345.3992395	,	347.6806084	,	349.9619772	,	351.7870722	,	353.6121673	,	355.4372624	,	356.8060837	,	358.6311787	,	360 ];
yp1 = [ 13.31155116	,	13.33135314	,	13.34785479	,	13.36105611	,	13.38085809	,	13.39735974	,	13.41386139	,	13.42706271	,	13.44356436	,	13.46006601	,	13.47986799	,	13.49636964	,	13.51617162	,	13.53267327	,	13.55247525	,	13.57227723	,	13.58877888	,	13.60858086	,	13.62508251	,	13.64158416	,	13.66468647	,	13.68448845	,	13.70429043	,	13.72409241	,	13.74719472	,	13.76369637	,	13.78349835	,	13.8	,	13.81980198	,	13.83960396	,	13.85610561	,	13.87260726	,	13.88910891	,	13.90561056	,	13.92871287	,	13.95181518	,	13.97161716	,	13.99141914	,	14.01452145	,	14.03762376	,	14.05742574	,	14.08052805	,	14.10033003	,	14.12013201	,	14.13993399	,	14.15973597	,	14.17953795	,	14.19933993	,	14.21914191	,	14.24224422	,	14.2620462	,	14.28184818	,	14.30165017	,	14.32475248	,	14.34455446	,	14.36435644	,	14.38415842	,	14.4039604	,	14.42376238	,	14.44356436	,	14.46336634	,	14.47986799	,	14.49966997	,	14.51617162	,	14.52607261	,	14.5359736	,	14.54587459	,	14.55247525	,	14.55907591	,	14.56237624	,	14.56567657	,	14.56567657	,	14.56567657	,	14.56237624	,	14.55577558	,	14.54917492	,	14.5359736	,	14.52607261	,	14.51617162	,	14.49966997	,	14.48316832	,	14.46666667	,	14.45346535	,	14.4369637	,	14.42046205	,	14.4039604	,	14.39075908	,	14.37425743	,	14.35775578	,	14.34125413	,	14.32475248	,	14.3049505	,	14.28844884	,	14.26864686	,	14.24554455	,	14.22574257	,	14.20924092	,	14.18613861	,	14.16633663	,	14.14653465	,	14.12673267	,	14.10693069	,	14.08712871	,	14.06732673	,	14.05082508	,	14.0310231	,	14.01452145	,	13.99471947	,	13.97821782	,	13.95841584	,	13.93861386	,	13.91881188	,	13.90231023	,	13.88580858	,	13.8660066	,	13.84620462	,	13.82970297	,	13.80990099	,	13.78679868	,	13.7669967	,	13.74719472	,	13.72739274	,	13.70429043	,	13.68118812	,	13.65808581	,	13.63828383	,	13.62178218	,	13.59867987	,	13.57887789	,	13.55907591	,	13.5359736	,	13.51287129	,	13.49306931	,	13.47656766	,	13.45346535	,	13.43366337	,	13.41386139	,	13.39405941	,	13.37425743	,	13.35115512	,	13.33135314	,	13.31155116	,	13.2950495	,	13.27854785	,	13.2620462	,	13.24224422	,	13.22244224	,	13.20594059	,	13.19273927	,	13.17623762	,	13.15643564	,	13.13663366	,	13.11683168	,	13.0970297	,	13.08052805	,	13.05742574	,	13.04422442	,	13.02772277	,	13.01122112	,	12.99141914	,	12.97821782	,	12.96171617	,	12.95181518	,	12.94191419	,	12.9320132	,	12.92541254	,	12.91881188	,	12.91551155	,	12.91551155	,	12.92211221	,	12.92871287	,	12.93861386	,	12.94851485	,	12.96171617	,	12.97161716	,	12.98481848	,	13.00132013	,	13.01782178	,	13.03432343	,	13.05082508	,	13.07062706	,	13.09042904	,	13.11023102	,	13.130033	,	13.14653465	,	13.16963696	,	13.18943894	,	13.21254125	,	13.23234323	,	13.25214521	,	13.26864686	,	13.28844884	,	13.30165017 ];

%vetores para plot da curva referência condutor 2
xp2 = [ 0	,	1.797752809	,	4.04494382	,	5.842696629	,	7.640449438	,	9.438202247	,	11.23595506	,	13.48314607	,	15.73033708	,	17.97752809	,	20.2247191	,	22.47191011	,	25.16853933	,	27.86516854	,	30.56179775	,	33.25842697	,	35.95505618	,	38.20224719	,	41.34831461	,	44.49438202	,	47.19101124	,	49.88764045	,	53.03370787	,	56.17977528	,	59.3258427	,	61.57303371	,	63.82022472	,	66.51685393	,	68.76404494	,	71.01123596	,	73.25842697	,	75.50561798	,	77.30337079	,	79.5505618	,	81.34831461	,	83.14606742	,	84.94382022	,	86.74157303	,	88.53932584	,	90.78651685	,	92.58426966	,	94.38202247	,	96.62921348	,	97.97752809	,	99.7752809	,	101.1235955	,	102.9213483	,	104.7191011	,	106.5168539	,	107.8651685	,	109.2134831	,	111.011236	,	112.3595506	,	113.7078652	,	115.505618	,	117.3033708	,	118.6516854	,	120	,	121.3483146	,	122.6966292	,	124.494382	,	126.2921348	,	127.6404494	,	128.988764	,	129.8876404	,	131.2359551	,	133.0337079	,	134.3820225	,	135.7303371	,	137.5280899	,	139.3258427	,	141.1235955	,	142.4719101	,	144.2696629	,	146.0674157	,	147.8651685	,	149.6629213	,	151.011236	,	152.8089888	,	154.1573034	,	155.9550562	,	157.3033708	,	158.6516854	,	160.4494382	,	161.7977528	,	163.5955056	,	165.8426966	,	167.6404494	,	169.4382022	,	171.2359551	,	173.0337079	,	174.8314607	,	176.6292135	,	179.3258427	,	181.1235955	,	183.3707865	,	185.6179775	,	187.8651685	,	190.1123596	,	192.3595506	,	194.6067416	,	196.8539326	,	199.1011236	,	201.3483146	,	203.5955056	,	206.2921348	,	208.988764	,	211.6853933	,	214.3820225	,	217.5280899	,	220.2247191	,	223.3707865	,	226.0674157	,	229.2134831	,	232.8089888	,	235.9550562	,	239.1011236	,	242.247191	,	244.9438202	,	248.0898876	,	250.7865169	,	253.0337079	,	255.2808989	,	258.4269663	,	260.6741573	,	262.9213483	,	265.1685393	,	267.4157303	,	270.1123596	,	271.9101124	,	273.7078652	,	275.505618	,	277.3033708	,	279.1011236	,	281.3483146	,	283.1460674	,	284.9438202	,	287.1910112	,	289.4382022	,	291.6853933	,	293.4831461	,	295.2808989	,	297.5280899	,	299.7752809	,	302.0224719	,	303.8202247	,	305.6179775	,	307.4157303	,	309.2134831	,	311.011236	,	312.8089888	,	315.0561798	,	316.8539326	,	318.6516854	,	320.4494382	,	322.247191	,	324.494382	,	326.2921348	,	328.0898876	,	329.8876404	,	331.6853933	,	333.9325843	,	336.1797753	,	338.4269663	,	340.6741573	,	342.4719101	,	344.2696629	,	346.5168539	,	348.3146067	,	350.5617978	,	353.258427	,	355.505618	,	357.3033708	,	359.5505618 ];
yp2 = [ 15.20243902	,	15.22195122	,	15.24634146	,	15.27073171	,	15.2902439	,	15.3097561	,	15.32926829	,	15.34878049	,	15.36829268	,	15.38780488	,	15.40243902	,	15.41707317	,	15.43658537	,	15.45121951	,	15.46585366	,	15.47560976	,	15.4804878	,	15.48536585	,	15.4902439	,	15.4902439	,	15.48536585	,	15.4804878	,	15.47560976	,	15.46097561	,	15.44634146	,	15.43658537	,	15.42195122	,	15.40731707	,	15.38780488	,	15.36829268	,	15.34878049	,	15.32439024	,	15.30487805	,	15.28536585	,	15.26585366	,	15.24634146	,	15.22195122	,	15.20243902	,	15.17317073	,	15.14878049	,	15.12439024	,	15.1	,	15.07560976	,	15.05121951	,	15.02195122	,	15.00243902	,	14.97804878	,	14.95365854	,	14.92439024	,	14.9	,	14.87560976	,	14.85121951	,	14.82195122	,	14.79756098	,	14.77317073	,	14.74390244	,	14.7195122	,	14.6902439	,	14.66097561	,	14.64146341	,	14.61219512	,	14.58292683	,	14.55365854	,	14.52926829	,	14.5097561	,	14.48536585	,	14.45609756	,	14.42682927	,	14.40243902	,	14.37317073	,	14.34390244	,	14.3097561	,	14.2804878	,	14.25121951	,	14.22195122	,	14.19756098	,	14.17317073	,	14.13902439	,	14.11463415	,	14.0902439	,	14.06097561	,	14.03658537	,	14.00731707	,	13.98292683	,	13.95853659	,	13.92926829	,	13.9	,	13.87073171	,	13.84146341	,	13.81219512	,	13.78780488	,	13.75853659	,	13.73902439	,	13.71463415	,	13.68536585	,	13.65609756	,	13.63170732	,	13.60731707	,	13.58780488	,	13.56341463	,	13.54390244	,	13.52439024	,	13.50487805	,	13.4902439	,	13.47073171	,	13.45609756	,	13.44634146	,	13.43170732	,	13.42195122	,	13.41219512	,	13.40731707	,	13.40243902	,	13.40243902	,	13.41219512	,	13.42195122	,	13.43658537	,	13.44634146	,	13.46097561	,	13.47560976	,	13.49512195	,	13.5195122	,	13.54390244	,	13.56829268	,	13.59756098	,	13.62195122	,	13.64634146	,	13.67073171	,	13.7	,	13.72926829	,	13.75365854	,	13.77804878	,	13.80243902	,	13.83170732	,	13.86097561	,	13.89512195	,	13.92439024	,	13.95365854	,	13.98780488	,	14.01707317	,	14.05121951	,	14.08536585	,	14.11463415	,	14.15365854	,	14.19268293	,	14.22682927	,	14.26097561	,	14.29512195	,	14.32439024	,	14.36341463	,	14.39268293	,	14.42682927	,	14.46097561	,	14.5	,	14.53414634	,	14.56341463	,	14.59756098	,	14.63170732	,	14.67073171	,	14.7	,	14.73414634	,	14.76829268	,	14.80243902	,	14.83658537	,	14.87560976	,	14.9097561	,	14.94390244	,	14.97804878	,	15.00731707	,	15.03658537	,	15.06585366	,	15.1	,	15.12926829	,	15.15853659	,	15.18292683 ];

%vetores para plot da curva referência condutor 3
xp3 = [ 0	,	1.846153846	,	3.230769231	,	5.076923077	,	6.923076923	,	8.769230769	,	10.61538462	,	12.46153846	,	14.30769231	,	16.61538462	,	18.46153846	,	20.30769231	,	22.15384615	,	24	,	25.84615385	,	27.69230769	,	29.53846154	,	31.38461538	,	33.23076923	,	34.61538462	,	36.46153846	,	37.84615385	,	39.69230769	,	41.53846154	,	43.38461538	,	45.23076923	,	47.07692308	,	48.92307692	,	50.76923077	,	52.61538462	,	54.92307692	,	56.76923077	,	58.15384615	,	60	,	61.84615385	,	63.69230769	,	66	,	68.30769231	,	70.15384615	,	72	,	73.84615385	,	75.69230769	,	77.53846154	,	79.84615385	,	81.23076923	,	83.07692308	,	84.92307692	,	86.76923077	,	89.07692308	,	90.92307692	,	93.23076923	,	95.53846154	,	97.84615385	,	100.1538462	,	102.9230769	,	105.6923077	,	108	,	110.7692308	,	113.5384615	,	116.3076923	,	119.0769231	,	121.8461538	,	125.5384615	,	128.7692308	,	132.4615385	,	134.7692308	,	138	,	141.2307692	,	145.3846154	,	148.6153846	,	151.3846154	,	154.1538462	,	157.3846154	,	160.6153846	,	163.8461538	,	166.6153846	,	169.3846154	,	172.1538462	,	174.9230769	,	177.6923077	,	180	,	182.7692308	,	185.0769231	,	186.9230769	,	189.2307692	,	191.5384615	,	193.8461538	,	196.1538462	,	198.4615385	,	200.7692308	,	202.6153846	,	204.9230769	,	206.7692308	,	208.6153846	,	210.9230769	,	213.6923077	,	216	,	217.8461538	,	220.1538462	,	222.4615385	,	224.3076923	,	226.1538462	,	228	,	230.3076923	,	232.1538462	,	234.4615385	,	236.3076923	,	238.6153846	,	240.9230769	,	243.2307692	,	245.0769231	,	246.9230769	,	248.7692308	,	251.0769231	,	252.9230769	,	255.6923077	,	258	,	260.3076923	,	263.0769231	,	265.3846154	,	268.1538462	,	270.4615385	,	273.2307692	,	275.5384615	,	278.3076923	,	280.6153846	,	282.9230769	,	286.1538462	,	288.9230769	,	291.6923077	,	294.9230769	,	297.6923077	,	300.4615385	,	303.2307692	,	306.4615385	,	310.1538462	,	313.3846154	,	316.6153846	,	319.8461538	,	323.0769231	,	325.8461538	,	329.0769231	,	332.3076923	,	335.5384615	,	338.3076923	,	341.0769231	,	343.8461538	,	346.6153846	,	348.9230769	,	351.6923077	,	353.5384615	,	355.8461538	,	358.1538462	,	360 ];
yp3 = [ 14.76126878	,	14.74373957	,	14.72621035	,	14.70283806	,	14.67946578	,	14.65609349	,	14.62687813	,	14.59766277	,	14.56844741	,	14.53923205	,	14.51001669	,	14.48664441	,	14.45742905	,	14.42821369	,	14.39899833	,	14.36978297	,	14.34056761	,	14.31135225	,	14.28213689	,	14.25292154	,	14.22370618	,	14.19449082	,	14.15943239	,	14.12437396	,	14.0951586	,	14.06594324	,	14.03088481	,	14.00166945	,	13.96661102	,	13.93155259	,	13.89065109	,	13.86143573	,	13.83222037	,	13.80300501	,	13.77378965	,	13.73873122	,	13.70367279	,	13.66277129	,	13.63355593	,	13.60434057	,	13.57512521	,	13.54590985	,	13.51085142	,	13.48163606	,	13.4524207	,	13.42904841	,	13.39983306	,	13.37646077	,	13.34724541	,	13.31803005	,	13.29465776	,	13.2654424	,	13.23622705	,	13.21285476	,	13.18948247	,	13.16026711	,	13.13689482	,	13.11352254	,	13.09599332	,	13.07846411	,	13.06093489	,	13.04340568	,	13.02587646	,	13.00834725	,	13.00250417	,	13.00250417	,	13.00250417	,	13.00250417	,	13.00834725	,	13.02003339	,	13.03171953	,	13.04924875	,	13.06677796	,	13.08430718	,	13.10767947	,	13.13689482	,	13.16026711	,	13.1836394	,	13.21285476	,	13.24207012	,	13.27712855	,	13.30634391	,	13.33555927	,	13.36477462	,	13.39983306	,	13.42320534	,	13.45826377	,	13.4933222	,	13.52838063	,	13.56343907	,	13.5984975	,	13.63355593	,	13.66861436	,	13.70367279	,	13.74457429	,	13.77963272	,	13.81469115	,	13.84974958	,	13.88480801	,	13.92570952	,	13.96076795	,	14.00166945	,	14.03088481	,	14.06594324	,	14.10100167	,	14.14190317	,	14.18280467	,	14.21786311	,	14.25292154	,	14.28797997	,	14.3230384	,	14.35809683	,	14.39315526	,	14.42821369	,	14.46327212	,	14.49833055	,	14.53338898	,	14.56844741	,	14.60350584	,	14.63856427	,	14.6736227	,	14.70868114	,	14.73789649	,	14.76711185	,	14.79632721	,	14.8196995	,	14.84891486	,	14.87813022	,	14.90734558	,	14.92487479	,	14.94824708	,	14.96577629	,	14.98914858	,	15.0066778	,	15.01836394	,	15.02420701	,	15.03005008	,	15.03005008	,	15.03005008	,	15.02420701	,	15.02420701	,	15.01252087	,	15.00083472	,	14.98914858	,	14.96577629	,	14.94824708	,	14.93071786	,	14.91318865	,	14.88981636	,	14.86644407	,	14.84307179	,	14.82554257	,	14.80217028	,	14.78464107 ];

%vetores para plot da curva referência condutor 4
xp4 = [ 0	,	1.86770428	,	3.73540856	,	5.60311284	,	7.470817121	,	9.338521401	,	11.67315175	,	13.54085603	,	15.87548638	,	18.21011673	,	21.01167315	,	23.3463035	,	26.14785992	,	28.94941634	,	31.75097276	,	35.01945525	,	37.82101167	,	40.62256809	,	42.95719844	,	46.22568093	,	49.49416342	,	51.82879377	,	54.63035019	,	57.43190661	,	60.70038911	,	63.50194553	,	66.30350195	,	69.10505837	,	71.90661479	,	75.17509728	,	78.44357977	,	81.71206226	,	84.51361868	,	87.78210117	,	91.05058366	,	93.85214008	,	96.18677043	,	98.98832685	,	101.7898833	,	104.1245136	,	106.92607	,	109.7276265	,	112.0622568	,	114.3968872	,	117.1984436	,	119.5330739	,	121.8677043	,	124.2023346	,	126.536965	,	128.4046693	,	130.7392996	,	132.6070039	,	134.4747082	,	136.3424125	,	138.6770428	,	141.0116732	,	143.8132296	,	146.1478599	,	148.4824903	,	150.3501946	,	152.6848249	,	155.0194553	,	157.3540856	,	159.688716	,	162.0233463	,	164.3579767	,	166.692607	,	169.0272374	,	171.8287938	,	174.6303502	,	177.4319066	,	180.233463	,	183.0350195	,	185.8365759	,	188.6381323	,	191.4396887	,	193.7743191	,	196.5758755	,	199.3774319	,	202.1789883	,	205.4474708	,	208.7159533	,	211.0505837	,	214.3190661	,	217.1206226	,	220.3891051	,	224.1245136	,	227.3929961	,	230.6614786	,	233.463035	,	237.1984436	,	240.9338521	,	244.6692607	,	248.4046693	,	252.1400778	,	255.8754864	,	259.1439689	,	262.4124514	,	265.6809339	,	269.4163424	,	272.6848249	,	275.4863813	,	278.2879377	,	280.6225681	,	282.4902724	,	285.2918288	,	287.6264591	,	290.4280156	,	293.229572	,	295.5642023	,	298.3657588	,	300.7003891	,	303.0350195	,	305.3696498	,	307.7042802	,	310.0389105	,	312.3735409	,	314.7081712	,	317.5097276	,	319.3774319	,	321.2451362	,	323.5797665	,	325.9143969	,	328.2490272	,	330.5836576	,	332.4513619	,	335.2529183	,	337.5875486	,	339.4552529	,	341.3229572	,	343.6575875	,	345.5252918	,	347.3929961	,	349.2607004	,	351.5953307	,	353.9299611	,	355.7976654	,	358.1322957	,	360 ];
yp4 = [ 12.97969543	,	12.96192893	,	12.93824027	,	12.91455161	,	12.89086294	,	12.86717428	,	12.84940778	,	12.83164129	,	12.80795262	,	12.78426396	,	12.7605753	,	12.7428088	,	12.71912014	,	12.69543147	,	12.67766497	,	12.65989848	,	12.64213198	,	12.63028765	,	12.61844332	,	12.60659898	,	12.60067682	,	12.59475465	,	12.58883249	,	12.58883249	,	12.58883249	,	12.58883249	,	12.59475465	,	12.60067682	,	12.61252115	,	12.63028765	,	12.64213198	,	12.65397631	,	12.67174281	,	12.69543147	,	12.71912014	,	12.73688663	,	12.75465313	,	12.77834179	,	12.80795262	,	12.83164129	,	12.85532995	,	12.88494078	,	12.91455161	,	12.93824027	,	12.9678511	,	12.99746193	,	13.02115059	,	13.05076142	,	13.08037225	,	13.10998308	,	13.13367174	,	13.15736041	,	13.18104907	,	13.2106599	,	13.24027073	,	13.27580372	,	13.31133672	,	13.34094755	,	13.37055838	,	13.40609137	,	13.4357022	,	13.47123519	,	13.50084602	,	13.53045685	,	13.56006768	,	13.59560068	,	13.62521151	,	13.6607445	,	13.6962775	,	13.73181049	,	13.76734349	,	13.79695431	,	13.83248731	,	13.86209814	,	13.89170897	,	13.9213198	,	13.94500846	,	13.97461929	,	14.00423012	,	14.02791878	,	14.05160745	,	14.06937394	,	14.08121827	,	14.09898477	,	14.11675127	,	14.13451777	,	14.1463621	,	14.1641286	,	14.17005076	,	14.18189509	,	14.18189509	,	14.18189509	,	14.17597293	,	14.1641286	,	14.15228426	,	14.14043993	,	14.1285956	,	14.11675127	,	14.09306261	,	14.06937394	,	14.04568528	,	14.02199662	,	13.99830795	,	13.98054146	,	13.95685279	,	13.92724196	,	13.9035533	,	13.87986464	,	13.85025381	,	13.82064298	,	13.79103215	,	13.76142132	,	13.73181049	,	13.70812183	,	13.678511	,	13.64890017	,	13.61928934	,	13.58967851	,	13.55414552	,	13.52453469	,	13.49492386	,	13.45939086	,	13.42385787	,	13.39424704	,	13.35871404	,	13.32910321	,	13.29949239	,	13.26395939	,	13.24027073	,	13.21658206	,	13.18697124	,	13.15736041	,	13.13367174	,	13.10406091	,	13.08037225	,	13.05076142	,	13.02707276	,	13.00338409	,	12.99153976 ];

%vetores para plot da curva referência condutor 5
xp5 = [ 0	,	1.375796178	,	3.210191083	,	5.044585987	,	6.878980892	,	8.713375796	,	10.5477707	,	12.38216561	,	14.21656051	,	16.50955414	,	18.80254777	,	21.0955414	,	23.38853503	,	25.68152866	,	27.51592357	,	29.35031847	,	31.18471338	,	33.47770701	,	35.31210191	,	37.14649682	,	39.43949045	,	41.27388535	,	43.56687898	,	45.85987261	,	48.15286624	,	50.44585987	,	52.7388535	,	55.49044586	,	57.78343949	,	60.07643312	,	62.36942675	,	64.66242038	,	66.95541401	,	69.24840764	,	71.54140127	,	73.8343949	,	76.58598726	,	78.87898089	,	81.17197452	,	83.92356688	,	86.67515924	,	89.42675159	,	92.17834395	,	94.47133758	,	97.22292994	,	100.433121	,	103.6433121	,	106.3949045	,	109.6050955	,	112.8152866	,	116.0254777	,	118.7770701	,	122.4458599	,	126.1146497	,	129.7834395	,	133.4522293	,	136.6624204	,	139.4140127	,	142.1656051	,	144.9171975	,	147.6687898	,	150.4203822	,	153.1719745	,	155.4649682	,	158.2165605	,	160.5095541	,	162.8025478	,	165.5541401	,	167.8471338	,	170.1401274	,	172.8917197	,	175.1847134	,	177.0191083	,	179.3121019	,	181.6050955	,	182.9808917	,	185.2738854	,	187.566879	,	189.4012739	,	191.2356688	,	193.0700637	,	194.9044586	,	196.2802548	,	198.1146497	,	199.9490446	,	202.2420382	,	204.0764331	,	205.4522293	,	207.2866242	,	209.1210191	,	210.955414	,	212.7898089	,	214.6242038	,	216.9171975	,	219.2101911	,	221.044586	,	222.4203822	,	224.7133758	,	227.0063694	,	228.8407643	,	230.6751592	,	232.5095541	,	234.343949	,	236.1783439	,	238.0127389	,	239.8471338	,	241.6815287	,	243.5159236	,	245.8089172	,	247.6433121	,	249.9363057	,	251.7707006	,	254.0636943	,	256.3566879	,	258.6496815	,	260.9426752	,	263.2356688	,	265.5286624	,	267.8216561	,	270.1146497	,	272.4076433	,	275.1592357	,	277.910828	,	280.6624204	,	283.4140127	,	285.7070064	,	288.9171975	,	292.1273885	,	294.8789809	,	298.089172	,	301.2993631	,	304.9681529	,	308.1783439	,	311.388535	,	315.0573248	,	318.7261146	,	322.3949045	,	325.6050955	,	328.8152866	,	332.0254777	,	335.2356688	,	338.4458599	,	342.1146497	,	344.866242	,	348.0764331	,	350.3694268	,	352.6624204	,	354.955414	,	357.2484076	,	359.0828025	,	360 ];
yp5 = [ 13.37810945	,	13.40298507	,	13.4278607	,	13.44776119	,	13.47263682	,	13.49253731	,	13.52238806	,	13.54726368	,	13.5721393	,	13.60199005	,	13.62686567	,	13.65671642	,	13.68656716	,	13.71641791	,	13.74626866	,	13.7761194	,	13.80099502	,	13.83084577	,	13.86567164	,	13.89054726	,	13.92039801	,	13.95024876	,	13.9800995	,	14.00995025	,	14.04477612	,	14.07960199	,	14.11442786	,	14.14427861	,	14.17910448	,	14.21393035	,	14.24378109	,	14.27363184	,	14.30348259	,	14.33333333	,	14.35820896	,	14.38308458	,	14.41293532	,	14.43781095	,	14.46268657	,	14.49253731	,	14.51741294	,	14.54228856	,	14.56716418	,	14.58706468	,	14.60696517	,	14.62686567	,	14.64179104	,	14.66666667	,	14.67661692	,	14.69154229	,	14.70149254	,	14.71144279	,	14.71144279	,	14.71144279	,	14.71144279	,	14.70646766	,	14.70149254	,	14.69154229	,	14.68159204	,	14.66169154	,	14.65174129	,	14.63681592	,	14.62189055	,	14.60199005	,	14.58208955	,	14.56218905	,	14.54228856	,	14.52238806	,	14.49751244	,	14.47263682	,	14.45273632	,	14.42288557	,	14.39303483	,	14.37313433	,	14.34825871	,	14.32338308	,	14.30348259	,	14.27860697	,	14.25373134	,	14.2238806	,	14.2039801	,	14.17910448	,	14.14925373	,	14.12437811	,	14.09950249	,	14.07462687	,	14.04975124	,	14.02487562	,	14	,	13.97512438	,	13.94527363	,	13.92039801	,	13.89054726	,	13.86567164	,	13.8358209	,	13.80597015	,	13.7761194	,	13.75124378	,	13.72139303	,	13.69651741	,	13.66666667	,	13.64179104	,	13.6119403	,	13.58706468	,	13.56716418	,	13.54228856	,	13.51243781	,	13.48756219	,	13.45771144	,	13.43283582	,	13.4079602	,	13.38308458	,	13.3681592	,	13.33830846	,	13.30845771	,	13.28855721	,	13.26368159	,	13.23880597	,	13.21393035	,	13.19402985	,	13.17412935	,	13.15422886	,	13.13432836	,	13.11940299	,	13.10447761	,	13.08457711	,	13.06965174	,	13.05472637	,	13.04477612	,	13.03482587	,	13.02985075	,	13.02985075	,	13.02985075	,	13.039801	,	13.04477612	,	13.04975124	,	13.06467662	,	13.07960199	,	13.09452736	,	13.11442786	,	13.13930348	,	13.16915423	,	13.19402985	,	13.21890547	,	13.24378109	,	13.26865672	,	13.29353234	,	13.31343284	,	13.34328358	,	13.36318408	,	13.37810945 ];

%vetores para plot da curva referência condutor 6
xp6 = [ 0	,	1.38996139	,	2.77992278	,	5.096525097	,	7.413127413	,	9.72972973	,	11.58301158	,	13.8996139	,	16.21621622	,	18.06949807	,	20.38610039	,	22.7027027	,	25.48262548	,	28.26254826	,	31.04247104	,	33.82239382	,	37.06563707	,	40.30888031	,	43.55212355	,	46.33204633	,	49.57528958	,	52.35521236	,	54.67181467	,	57.45173745	,	59.76833977	,	62.08494208	,	64.86486486	,	67.64478764	,	69.4980695	,	72.27799228	,	74.59459459	,	76.44787645	,	78.3011583	,	80.61776062	,	82.00772201	,	83.86100386	,	86.17760618	,	88.03088803	,	89.88416988	,	91.73745174	,	93.59073359	,	94.98069498	,	96.83397683	,	98.68725869	,	100.5405405	,	101.9305019	,	103.7837838	,	105.6370656	,	107.4903475	,	109.3436293	,	110.7335907	,	112.5868726	,	113.976834	,	115.3667954	,	116.7567568	,	118.1467181	,	119.0733591	,	120.9266409	,	122.3166023	,	123.7065637	,	125.5598456	,	126.9498069	,	128.3397683	,	129.2664093	,	130.6563707	,	132.046332	,	133.4362934	,	134.8262548	,	136.6795367	,	138.0694981	,	139.4594595	,	141.3127413	,	142.7027027	,	144.5559846	,	145.9459459	,	147.3359073	,	148.7258687	,	150.5791506	,	151.969112	,	153.3590734	,	154.7490347	,	156.1389961	,	157.992278	,	159.8455598	,	161.6988417	,	163.5521236	,	165.4054054	,	166.7953668	,	168.1853282	,	170.5019305	,	171.8918919	,	173.2818533	,	175.1351351	,	176.988417	,	178.8416988	,	180.6949807	,	183.011583	,	185.3281853	,	187.1814672	,	189.4980695	,	191.8146718	,	194.1312741	,	196.4478764	,	198.7644788	,	201.0810811	,	203.3976834	,	206.1776062	,	208.957529	,	212.6640927	,	215.4440154	,	218.6872587	,	221.4671815	,	224.2471042	,	227.953668	,	230.7335907	,	233.976834	,	237.2200772	,	240	,	242.3166023	,	245.0965251	,	247.4131274	,	250.6563707	,	253.4362934	,	256.2162162	,	258.0694981	,	260.3861004	,	262.7027027	,	265.019305	,	267.3359073	,	269.6525097	,	271.5057915	,	273.3590734	,	275.2123552	,	277.0656371	,	278.4555985	,	280.3088803	,	282.1621622	,	284.4787645	,	285.8687259	,	287.7220077	,	289.5752896	,	291.4285714	,	293.2818533	,	295.1351351	,	296.988417	,	298.8416988	,	300.6949807	,	302.5482625	,	304.4015444	,	305.7915058	,	307.6447876	,	309.034749	,	310.4247104	,	312.2779923	,	314.1312741	,	315.984556	,	317.8378378	,	319.6911197	,	321.5444015	,	323.3976834	,	324.7876448	,	326.6409266	,	328.4942085	,	330.3474903	,	332.2007722	,	334.0540541	,	335.9073359	,	337.7606178	,	339.6138996	,	341.4671815	,	343.3204633	,	345.1737452	,	347.027027	,	348.8803089	,	351.1969112	,	353.0501931	,	355.3667954	,	357.2200772	,	358.6100386	,	360 ];
yp6 = [ 14.56302521	,	14.58403361	,	14.60084034	,	14.62184874	,	14.64705882	,	14.66806723	,	14.68907563	,	14.71008403	,	14.72689076	,	14.74369748	,	14.7605042	,	14.77310924	,	14.78991597	,	14.80252101	,	14.81512605	,	14.82352941	,	14.83193277	,	14.83613445	,	14.83613445	,	14.83613445	,	14.83613445	,	14.82773109	,	14.81932773	,	14.80672269	,	14.79831933	,	14.78571429	,	14.76890756	,	14.75210084	,	14.73529412	,	14.71848739	,	14.70168067	,	14.68067227	,	14.66386555	,	14.64285714	,	14.6302521	,	14.6092437	,	14.58403361	,	14.55882353	,	14.53361345	,	14.51680672	,	14.49159664	,	14.46638655	,	14.44537815	,	14.42016807	,	14.39915966	,	14.37394958	,	14.35294118	,	14.32773109	,	14.29831933	,	14.26890756	,	14.24369748	,	14.21848739	,	14.19747899	,	14.17226891	,	14.1512605	,	14.1302521	,	14.10504202	,	14.07983193	,	14.06302521	,	14.03781513	,	14.01260504	,	13.99159664	,	13.96638655	,	13.94537815	,	13.92016807	,	13.89495798	,	13.8697479	,	13.84453782	,	13.81932773	,	13.79411765	,	13.76470588	,	13.74369748	,	13.71428571	,	13.68907563	,	13.66386555	,	13.64285714	,	13.61764706	,	13.59243697	,	13.57142857	,	13.54621849	,	13.52521008	,	13.5	,	13.47478992	,	13.44957983	,	13.42016807	,	13.3907563	,	13.36554622	,	13.34453782	,	13.32352941	,	13.29831933	,	13.27310924	,	13.25210084	,	13.23109244	,	13.20588235	,	13.18487395	,	13.16386555	,	13.13865546	,	13.11344538	,	13.09243697	,	13.07142857	,	13.04621849	,	13.02941176	,	13.00840336	,	12.99159664	,	12.97478992	,	12.95798319	,	12.94117647	,	12.92436975	,	12.91176471	,	12.90336134	,	12.89495798	,	12.8907563	,	12.8907563	,	12.89495798	,	12.90336134	,	12.91176471	,	12.92016807	,	12.93697479	,	12.94957983	,	12.96638655	,	12.98319328	,	13.00420168	,	13.02521008	,	13.04621849	,	13.06722689	,	13.08823529	,	13.1092437	,	13.13445378	,	13.15546218	,	13.18067227	,	13.20588235	,	13.23109244	,	13.25630252	,	13.28151261	,	13.30672269	,	13.33193277	,	13.35714286	,	13.38655462	,	13.41176471	,	13.43697479	,	13.45798319	,	13.48739496	,	13.51680672	,	13.54621849	,	13.57563025	,	13.60504202	,	13.63445378	,	13.66386555	,	13.69327731	,	13.72268908	,	13.75210084	,	13.77731092	,	13.80672269	,	13.83613445	,	13.86554622	,	13.89495798	,	13.92436975	,	13.95378151	,	13.98739496	,	14.01680672	,	14.04621849	,	14.07563025	,	14.1092437	,	14.13445378	,	14.16386555	,	14.19327731	,	14.22268908	,	14.25210084	,	14.28151261	,	14.30672269	,	14.33613445	,	14.36554622	,	14.39495798	,	14.42016807	,	14.44537815	,	14.47058824	,	14.5	,	14.5210084	,	14.54201681	,	14.55882353 ];

%vetores para plot da curva referência condutor 7
xp7 = [ 0	,	1.358490566	,	2.716981132	,	4.075471698	,	5.433962264	,	6.79245283	,	8.150943396	,	9.962264151	,	11.32075472	,	12.67924528	,	14.03773585	,	15.39622642	,	16.75471698	,	18.11320755	,	19.47169811	,	20.37735849	,	21.73584906	,	23.54716981	,	24.90566038	,	26.26415094	,	27.16981132	,	28.52830189	,	29.88679245	,	31.69811321	,	33.05660377	,	34.41509434	,	35.32075472	,	36.67924528	,	38.49056604	,	39.8490566	,	41.20754717	,	42.56603774	,	43.9245283	,	45.73584906	,	47.09433962	,	48.90566038	,	50.26415094	,	52.0754717	,	53.88679245	,	55.69811321	,	56.60377358	,	57.96226415	,	59.32075472	,	60.67924528	,	62.49056604	,	63.8490566	,	65.20754717	,	67.47169811	,	69.28301887	,	71.09433962	,	72.45283019	,	73.81132075	,	75.16981132	,	76.98113208	,	78.79245283	,	80.1509434	,	81.50943396	,	83.32075472	,	84.67924528	,	86.03773585	,	87.8490566	,	89.66037736	,	91.9245283	,	94.18867925	,	96	,	98.26415094	,	100.5283019	,	102.3396226	,	104.1509434	,	105.9622642	,	108.2264151	,	110.9433962	,	113.6603774	,	116.3773585	,	119.0943396	,	122.2641509	,	125.8867925	,	129.509434	,	132.6792453	,	135.8490566	,	139.0188679	,	142.6415094	,	145.8113208	,	148.9811321	,	151.6981132	,	154.4150943	,	156.6792453	,	159.3962264	,	162.1132075	,	164.8301887	,	167.0943396	,	169.3584906	,	171.6226415	,	173.4339623	,	175.245283	,	177.509434	,	179.3207547	,	181.1320755	,	183.3962264	,	185.2075472	,	186.5660377	,	188.3773585	,	190.6415094	,	192.4528302	,	194.2641509	,	196.0754717	,	197.8867925	,	199.6981132	,	201.509434	,	202.8679245	,	204.6792453	,	206.490566	,	208.3018868	,	210.1132075	,	211.4716981	,	213.2830189	,	214.6415094	,	216.4528302	,	217.8113208	,	219.6226415	,	220.9811321	,	222.7924528	,	224.6037736	,	225.9622642	,	227.7735849	,	229.5849057	,	231.3962264	,	232.754717	,	234.5660377	,	236.3773585	,	238.1886792	,	240	,	241.8113208	,	243.6226415	,	245.4339623	,	247.245283	,	249.0566038	,	250.8679245	,	252.6792453	,	254.490566	,	256.754717	,	258.5660377	,	260.3773585	,	262.6415094	,	264.9056604	,	267.1698113	,	268.9811321	,	271.245283	,	273.0566038	,	275.3207547	,	277.5849057	,	280.3018868	,	283.0188679	,	286.1886792	,	289.8113208	,	292.9811321	,	296.1509434	,	299.7735849	,	303.3962264	,	307.0188679	,	310.1886792	,	313.8113208	,	317.8867925	,	321.509434	,	324.2264151	,	327.3962264	,	330.5660377	,	333.2830189	,	336	,	338.7169811	,	341.4339623	,	344.1509434	,	346.8679245	,	349.1320755	,	351.3962264	,	352.754717	,	354.5660377	,	355.9245283	,	358.1886792	,	360 ];
yp7 = [ 15.41707718	,	15.39655172	,	15.38013136	,	15.35960591	,	15.33908046	,	15.31855501	,	15.29802956	,	15.27339901	,	15.24876847	,	15.23234811	,	15.21182266	,	15.18719212	,	15.16256158	,	15.14203612	,	15.12151067	,	15.10098522	,	15.08045977	,	15.05582923	,	15.03530378	,	15.01067323	,	14.99014778	,	14.96962233	,	14.9408867	,	14.91625616	,	14.89162562	,	14.86699507	,	14.84236453	,	14.8136289	,	14.78078818	,	14.76026273	,	14.73152709	,	14.70689655	,	14.68226601	,	14.65353038	,	14.62479475	,	14.59605911	,	14.57142857	,	14.54269294	,	14.51395731	,	14.48932677	,	14.46469622	,	14.44006568	,	14.41543514	,	14.3908046	,	14.36617406	,	14.34154351	,	14.31691297	,	14.28407225	,	14.25944171	,	14.23070608	,	14.20197044	,	14.1773399	,	14.15681445	,	14.136289	,	14.11165846	,	14.08702791	,	14.06239737	,	14.04187192	,	14.01724138	,	13.99671593	,	13.97208539	,	13.95155993	,	13.93103448	,	13.90640394	,	13.8817734	,	13.86124795	,	13.8407225	,	13.82430213	,	13.80788177	,	13.78735632	,	13.77093596	,	13.7545156	,	13.73399015	,	13.71756979	,	13.70525452	,	13.69293924	,	13.68062397	,	13.67651888	,	13.67651888	,	13.67651888	,	13.68883415	,	13.70114943	,	13.7134647	,	13.72577997	,	13.74220033	,	13.75862069	,	13.77504105	,	13.7955665	,	13.82019704	,	13.84482759	,	13.86535304	,	13.88587849	,	13.91050903	,	13.93513957	,	13.96387521	,	13.98850575	,	14.01313629	,	14.03776683	,	14.06239737	,	14.08702791	,	14.11576355	,	14.14449918	,	14.16912972	,	14.19376026	,	14.2183908	,	14.25123153	,	14.28407225	,	14.31280788	,	14.3456486	,	14.37027915	,	14.39901478	,	14.4318555	,	14.46469622	,	14.49343186	,	14.52216749	,	14.55090312	,	14.57553366	,	14.60426929	,	14.63300493	,	14.66584565	,	14.69458128	,	14.727422	,	14.75615764	,	14.78489327	,	14.8136289	,	14.84236453	,	14.87520525	,	14.90394089	,	14.93267652	,	14.96551724	,	14.99425287	,	15.02298851	,	15.05582923	,	15.08866995	,	15.12151067	,	15.14614122	,	15.17487685	,	15.20361248	,	15.2364532	,	15.26929392	,	15.29802956	,	15.33087028	,	15.35960591	,	15.39244663	,	15.41707718	,	15.44581281	,	15.46633826	,	15.4909688	,	15.51559934	,	15.54433498	,	15.56896552	,	15.59770115	,	15.62233169	,	15.64696223	,	15.67159278	,	15.69211823	,	15.70853859	,	15.72495895	,	15.74137931	,	15.7454844	,	15.74958949	,	15.74958949	,	15.7454844	,	15.73727422	,	15.72906404	,	15.71674877	,	15.7044335	,	15.68801314	,	15.66338259	,	15.64285714	,	15.6182266	,	15.59770115	,	15.5771757	,	15.55254516	,	15.5320197	,	15.51149425	,	15.49507389	,	15.47454844	,	15.45402299	,	15.43349754 ];

%vetores para plot da curva referência condutor 8
xp8 = [ 0	,	1.858064516	,	3.716129032	,	5.574193548	,	7.432258065	,	9.75483871	,	12.07741935	,	13.93548387	,	16.25806452	,	18.11612903	,	20.43870968	,	22.76129032	,	25.08387097	,	27.87096774	,	30.19354839	,	32.98064516	,	35.76774194	,	38.55483871	,	41.80645161	,	44.59354839	,	47.38064516	,	50.16774194	,	53.88387097	,	57.13548387	,	60.38709677	,	64.10322581	,	67.35483871	,	71.07096774	,	74.32258065	,	77.57419355	,	80.82580645	,	83.61290323	,	85.93548387	,	89.18709677	,	92.43870968	,	95.22580645	,	97.08387097	,	98.94193548	,	101.7290323	,	104.0516129	,	106.3741935	,	108.6967742	,	111.0193548	,	113.3419355	,	115.2	,	117.5225806	,	119.8451613	,	122.1677419	,	124.4903226	,	127.2774194	,	129.6	,	131.9225806	,	134.2451613	,	137.0322581	,	139.3548387	,	141.6774194	,	144	,	146.3225806	,	148.1806452	,	150.5032258	,	152.8258065	,	155.1483871	,	157.9354839	,	160.2580645	,	162.116129	,	164.4387097	,	167.2258065	,	170.0129032	,	172.8	,	175.1225806	,	177.9096774	,	180.6967742	,	183.9483871	,	187.2	,	189.9870968	,	193.2387097	,	196.4903226	,	199.2774194	,	201.6	,	204.3870968	,	207.6387097	,	210.8903226	,	214.6064516	,	217.3935484	,	220.6451613	,	224.8258065	,	228.5419355	,	232.2580645	,	235.5096774	,	239.2258065	,	243.4064516	,	247.1225806	,	250.8387097	,	254.5548387	,	257.8064516	,	261.5225806	,	264.7741935	,	268.0258065	,	271.2774194	,	274.5290323	,	277.316129	,	280.5677419	,	283.8193548	,	286.6064516	,	289.3935484	,	291.716129	,	294.0387097	,	295.8967742	,	298.683871	,	301.0064516	,	302.8645161	,	305.1870968	,	307.0451613	,	309.8322581	,	312.6193548	,	314.4774194	,	316.8	,	319.1225806	,	321.4451613	,	323.3032258	,	325.6258065	,	327.9483871	,	330.2709677	,	332.5935484	,	334.4516129	,	336.3096774	,	339.0967742	,	341.4193548	,	343.7419355	,	345.6	,	347.4580645	,	349.7806452	,	352.1032258	,	354.4258065	,	356.7483871	,	359.0709677	,	360 ];
yp8 = [ 14.17905405	,	14.16131757	,	14.13766892	,	14.11402027	,	14.08445946	,	14.06081081	,	14.03716216	,	14.01351351	,	13.98986486	,	13.96621622	,	13.94847973	,	13.93074324	,	13.90709459	,	13.88344595	,	13.86570946	,	13.84797297	,	13.83614865	,	13.81841216	,	13.80658784	,	13.79476351	,	13.78885135	,	13.78293919	,	13.78293919	,	13.78293919	,	13.78885135	,	13.79476351	,	13.8125	,	13.83023649	,	13.84206081	,	13.86570946	,	13.88935811	,	13.90709459	,	13.93074324	,	13.95439189	,	13.9839527	,	14.00760135	,	14.03125	,	14.06081081	,	14.09037162	,	14.11993243	,	14.14949324	,	14.18496622	,	14.21452703	,	14.24408784	,	14.27364865	,	14.30320946	,	14.33277027	,	14.36824324	,	14.40962838	,	14.44510135	,	14.48648649	,	14.52195946	,	14.55743243	,	14.59290541	,	14.63429054	,	14.66976351	,	14.71114865	,	14.75253378	,	14.78800676	,	14.82347973	,	14.86486486	,	14.90033784	,	14.93581081	,	14.97128378	,	15.00675676	,	15.04222973	,	15.08361486	,	15.125	,	15.16638514	,	15.20185811	,	15.23733108	,	15.27871622	,	15.32010135	,	15.36148649	,	15.3910473	,	15.42060811	,	15.45016892	,	15.47972973	,	15.50337838	,	15.52702703	,	15.55658784	,	15.58023649	,	15.59797297	,	15.62162162	,	15.63344595	,	15.64527027	,	15.65709459	,	15.66300676	,	15.66300676	,	15.65709459	,	15.65118243	,	15.63344595	,	15.61570946	,	15.60388514	,	15.58614865	,	15.55658784	,	15.52702703	,	15.49746622	,	15.46790541	,	15.43243243	,	15.39695946	,	15.36739865	,	15.33192568	,	15.30236486	,	15.27280405	,	15.23141892	,	15.19594595	,	15.16047297	,	15.125	,	15.08952703	,	15.05996622	,	15.03040541	,	14.99493243	,	14.95945946	,	14.92398649	,	14.88851351	,	14.85304054	,	14.81756757	,	14.78209459	,	14.74662162	,	14.71114865	,	14.66976351	,	14.63429054	,	14.59881757	,	14.56334459	,	14.52787162	,	14.49239865	,	14.45692568	,	14.4214527	,	14.39189189	,	14.35641892	,	14.32094595	,	14.28547297	,	14.25591216	,	14.22635135	,	14.2027027	,	14.19087838 ];

%vetores para plot da curva referência condutor 9
xp9 = [ 0	,	1.355081556	,	3.161856964	,	4.968632371	,	6.775407779	,	8.582183187	,	10.84065245	,	12.64742785	,	14.45420326	,	16.26097867	,	18.06775408	,	19.87452949	,	21.68130489	,	23.4880803	,	24.84316186	,	26.19824341	,	28.00501882	,	29.36010038	,	31.16687578	,	32.52195734	,	33.8770389	,	35.6838143	,	37.03889586	,	38.84567127	,	40.20075282	,	42.00752823	,	43.81430364	,	45.62107905	,	47.42785445	,	48.78293601	,	50.58971142	,	52.39648683	,	54.20326223	,	56.01003764	,	57.81681305	,	59.62358846	,	61.43036386	,	63.23713927	,	65.04391468	,	66.85069009	,	68.20577164	,	70.01254705	,	71.81932246	,	73.62609787	,	75.43287327	,	77.23964868	,	79.49811794	,	81.7565872	,	84.01505646	,	86.27352572	,	88.98368883	,	91.24215809	,	93.9523212	,	96.21079046	,	99.37264743	,	102.5345044	,	105.2446675	,	108.4065245	,	112.0200753	,	115.6336261	,	119.2471769	,	122.8607277	,	126.4742785	,	129.6361355	,	133.2496863	,	136.4115433	,	140.0250941	,	144.0903388	,	147.7038896	,	151.3174404	,	154.0276035	,	156.7377666	,	158.9962359	,	161.2547051	,	163.9648683	,	166.2233375	,	168.0301129	,	170.2885822	,	172.9987453	,	175.7089084	,	177.9673777	,	180.2258469	,	182.4843162	,	184.2910916	,	186.097867	,	188.3563363	,	190.1631117	,	192.4215809	,	194.2283563	,	196.4868256	,	198.293601	,	200.1003764	,	201.9071518	,	203.7139272	,	205.9723965	,	207.7791719	,	209.5859473	,	211.3927227	,	213.1994981	,	215.0062735	,	216.8130489	,	218.6198243	,	220.4265997	,	222.2333752	,	224.0401506	,	226.2986198	,	228.5570891	,	230.8155583	,	232.6223338	,	233.9774153	,	235.3324969	,	237.1392723	,	238.9460477	,	241.2045169	,	243.0112923	,	244.8180678	,	247.076537	,	248.8833124	,	251.1417817	,	252.9485571	,	254.7553325	,	257.0138018	,	258.8205772	,	261.0790464	,	263.3375157	,	265.5959849	,	267.8544542	,	270.5646173	,	272.8230866	,	275.0815558	,	277.3400251	,	280.0501882	,	281.8569636	,	284.5671267	,	287.2772898	,	290.8908407	,	294.0526976	,	297.6662484	,	301.7314931	,	305.3450439	,	308.9585947	,	312.5721455	,	316.1856964	,	319.7992472	,	322.9611041	,	325.6712673	,	328.8331242	,	331.9949812	,	334.7051443	,	337.4153074	,	340.1254705	,	342.3839398	,	345.0941029	,	347.3525721	,	349.6110414	,	352.3212045	,	354.5796738	,	356.838143	,	358.6449184	,	360 ];
yp9 = [ 13.76229508	,	13.78688525	,	13.81639344	,	13.84590164	,	13.88032787	,	13.90983607	,	13.9442623	,	13.97377049	,	14.00327869	,	14.03770492	,	14.07213115	,	14.10163934	,	14.13606557	,	14.1704918	,	14.2	,	14.22459016	,	14.25409836	,	14.27868852	,	14.30819672	,	14.33770492	,	14.36721311	,	14.39672131	,	14.42622951	,	14.46065574	,	14.49016393	,	14.52459016	,	14.55409836	,	14.58360656	,	14.61311475	,	14.64754098	,	14.68196721	,	14.71147541	,	14.74098361	,	14.77540984	,	14.80491803	,	14.83934426	,	14.87377049	,	14.90327869	,	14.93770492	,	14.96721311	,	14.99672131	,	15.02622951	,	15.0557377	,	15.0852459	,	15.1147541	,	15.1442623	,	15.17377049	,	15.20327869	,	15.23770492	,	15.27213115	,	15.30163934	,	15.33114754	,	15.3557377	,	15.38032787	,	15.40983607	,	15.43934426	,	15.46393443	,	15.48852459	,	15.50819672	,	15.52786885	,	15.54262295	,	15.55245902	,	15.55737705	,	15.55737705	,	15.55737705	,	15.54754098	,	15.53278689	,	15.51803279	,	15.49836066	,	15.47377049	,	15.44918033	,	15.4295082	,	15.40491803	,	15.37540984	,	15.3557377	,	15.33114754	,	15.30163934	,	15.27213115	,	15.23770492	,	15.20327869	,	15.16885246	,	15.13442623	,	15.1	,	15.07540984	,	15.04098361	,	15.00655738	,	14.96721311	,	14.93770492	,	14.90327869	,	14.86885246	,	14.83442623	,	14.8	,	14.7704918	,	14.74098361	,	14.70655738	,	14.67213115	,	14.63278689	,	14.59344262	,	14.56393443	,	14.5295082	,	14.49016393	,	14.45081967	,	14.41147541	,	14.37704918	,	14.33770492	,	14.29836066	,	14.26393443	,	14.21967213	,	14.19016393	,	14.16065574	,	14.13114754	,	14.09672131	,	14.06229508	,	14.02786885	,	13.99344262	,	13.95409836	,	13.92459016	,	13.89508197	,	13.8557377	,	13.82622951	,	13.79180328	,	13.75737705	,	13.72786885	,	13.69344262	,	13.65901639	,	13.6295082	,	13.6	,	13.57540984	,	13.54590164	,	13.51639344	,	13.49672131	,	13.47213115	,	13.44754098	,	13.42786885	,	13.40327869	,	13.38360656	,	13.36885246	,	13.35409836	,	13.3442623	,	13.33934426	,	13.33442623	,	13.33934426	,	13.34918033	,	13.35901639	,	13.37377049	,	13.38852459	,	13.40327869	,	13.43278689	,	13.45737705	,	13.47704918	,	13.50163934	,	13.52622951	,	13.5557377	,	13.5852459	,	13.6147541	,	13.6442623	,	13.67868852	,	13.70327869	,	13.73770492	,	13.76229508 ];

%vetores para plot da curva referência condutor 10
xp10 = [ 0	,	1.811320755	,	3.622641509	,	5.433962264	,	7.245283019	,	9.962264151	,	12.67924528	,	14.94339623	,	17.66037736	,	20.83018868	,	23.54716981	,	26.71698113	,	29.88679245	,	33.50943396	,	37.58490566	,	40.30188679	,	43.9245283	,	48.45283019	,	52.0754717	,	55.69811321	,	59.77358491	,	63.39622642	,	66.56603774	,	69.28301887	,	72	,	74.71698113	,	77.43396226	,	80.60377358	,	83.32075472	,	86.03773585	,	88.75471698	,	91.01886792	,	93.28301887	,	96	,	98.71698113	,	100.9811321	,	103.245283	,	105.9622642	,	108.2264151	,	110.490566	,	112.3018868	,	114.1132075	,	115.9245283	,	117.7358491	,	119.5471698	,	121.3584906	,	123.1698113	,	124.9811321	,	126.7924528	,	128.6037736	,	130.4150943	,	132.2264151	,	134.490566	,	136.3018868	,	138.5660377	,	140.3773585	,	142.1886792	,	144.9056604	,	146.7169811	,	148.9811321	,	150.7924528	,	153.0566038	,	155.3207547	,	157.1320755	,	159.8490566	,	161.6603774	,	163.9245283	,	166.1886792	,	168	,	170.2641509	,	172.9811321	,	175.245283	,	177.509434	,	180.2264151	,	182.9433962	,	185.2075472	,	187.4716981	,	190.6415094	,	193.8113208	,	196.5283019	,	199.245283	,	201.9622642	,	204.6792453	,	207.8490566	,	211.9245283	,	216	,	220.0754717	,	223.6981132	,	227.3207547	,	231.3962264	,	235.0188679	,	239.0943396	,	242.7169811	,	246.3396226	,	249.9622642	,	253.5849057	,	256.3018868	,	259.0188679	,	261.7358491	,	264.9056604	,	268.0754717	,	270.7924528	,	273.509434	,	275.7735849	,	278.0377358	,	280.3018868	,	282.5660377	,	285.2830189	,	287.5471698	,	289.8113208	,	292.5283019	,	294.7924528	,	297.0566038	,	299.3207547	,	301.5849057	,	303.3962264	,	305.6603774	,	307.9245283	,	310.1886792	,	312	,	313.8113208	,	316.0754717	,	318.7924528	,	321.0566038	,	322.8679245	,	325.1320755	,	327.3962264	,	329.2075472	,	331.4716981	,	333.2830189	,	336	,	338.7169811	,	340.5283019	,	342.3396226	,	344.1509434	,	346.4150943	,	348.2264151	,	350.490566	,	353.2075472	,	355.9245283	,	358.1886792	,	360 ];
yp10 = [ 14.35526316	,	14.375	,	14.39473684	,	14.41940789	,	14.44407895	,	14.46875	,	14.49342105	,	14.51809211	,	14.54276316	,	14.5625	,	14.58717105	,	14.60690789	,	14.62171053	,	14.63651316	,	14.65131579	,	14.66118421	,	14.67105263	,	14.67598684	,	14.67105263	,	14.66118421	,	14.65131579	,	14.63651316	,	14.62171053	,	14.60690789	,	14.58717105	,	14.56743421	,	14.54769737	,	14.52796053	,	14.50328947	,	14.47368421	,	14.44901316	,	14.42434211	,	14.39473684	,	14.36513158	,	14.33059211	,	14.30098684	,	14.27631579	,	14.24671053	,	14.21710526	,	14.18256579	,	14.15789474	,	14.12828947	,	14.09868421	,	14.07401316	,	14.04934211	,	14.01973684	,	13.99013158	,	13.96052632	,	13.93092105	,	13.90131579	,	13.87664474	,	13.84210526	,	13.80756579	,	13.77302632	,	13.74342105	,	13.70888158	,	13.67434211	,	13.63980263	,	13.60526316	,	13.57565789	,	13.54111842	,	13.50657895	,	13.47697368	,	13.44736842	,	13.40789474	,	13.37335526	,	13.33881579	,	13.30427632	,	13.26973684	,	13.24013158	,	13.21052632	,	13.17598684	,	13.14638158	,	13.11184211	,	13.08223684	,	13.05263158	,	13.02796053	,	12.99835526	,	12.96875	,	12.94407895	,	12.91940789	,	12.89967105	,	12.875	,	12.85526316	,	12.83552632	,	12.81578947	,	12.80098684	,	12.79605263	,	12.79605263	,	12.79605263	,	12.80098684	,	12.81085526	,	12.82565789	,	12.84539474	,	12.86513158	,	12.88980263	,	12.90953947	,	12.93421053	,	12.95394737	,	12.97861842	,	13.00822368	,	13.03782895	,	13.07236842	,	13.09703947	,	13.12664474	,	13.16118421	,	13.19078947	,	13.22039474	,	13.25	,	13.28453947	,	13.31414474	,	13.34868421	,	13.38322368	,	13.41776316	,	13.45230263	,	13.48684211	,	13.52138158	,	13.56085526	,	13.59046053	,	13.625	,	13.65953947	,	13.69407895	,	13.73355263	,	13.76809211	,	13.79769737	,	13.83223684	,	13.87171053	,	13.90625	,	13.94078947	,	13.97532895	,	14.00986842	,	14.04440789	,	14.07894737	,	14.10855263	,	14.13815789	,	14.16776316	,	14.19736842	,	14.22697368	,	14.25657895	,	14.29111842	,	14.32072368	,	14.34046053 ];

%vetores para plot da curva referência condutor 11
xp11 = [ 0	,	1.829733164	,	4.116899619	,	6.404066074	,	8.691232529	,	10.97839898	,	13.72299873	,	16.01016518	,	17.83989835	,	20.1270648	,	21.95679797	,	23.78653113	,	26.07369759	,	28.36086404	,	30.6480305	,	32.93519695	,	35.22236341	,	37.50952986	,	39.33926302	,	41.62642948	,	43.91359593	,	45.7433291	,	48.03049555	,	49.86022872	,	52.14739517	,	54.43456163	,	56.72172808	,	59.00889454	,	61.29606099	,	63.58322745	,	65.41296061	,	67.70012706	,	69.98729352	,	71.81702668	,	74.56162643	,	77.30622618	,	79.59339263	,	81.88055909	,	84.16772554	,	86.45489199	,	88.74205845	,	91.94409149	,	95.14612452	,	97.89072427	,	100.1778907	,	103.3799238	,	106.5819568	,	109.7839898	,	112.5285896	,	116.1880559	,	119.8475222	,	123.0495553	,	127.1664549	,	130.3684879	,	133.570521	,	136.772554	,	140.4320203	,	144.0914867	,	147.2935197	,	150.952986	,	154.6124524	,	158.2719187	,	161.931385	,	165.133418	,	168.7928844	,	171.9949174	,	175.1969504	,	177.9415502	,	180.6861499	,	183.4307497	,	185.7179161	,	188.0050826	,	190.7496823	,	193.0368488	,	195.7814485	,	198.5260483	,	200.8132147	,	203.5578145	,	205.8449809	,	208.5895807	,	210.4193139	,	212.7064803	,	214.9936468	,	217.7382465	,	220.4828463	,	222.7700127	,	225.0571792	,	227.3443456	,	230.0889454	,	232.3761118	,	234.6632783	,	237.407878	,	240.1524778	,	242.8970775	,	245.6416773	,	247.9288437	,	250.2160102	,	252.9606099	,	255.7052097	,	258.4498094	,	261.6518424	,	263.9390089	,	267.1410419	,	270.343075	,	273.0876747	,	276.2897078	,	279.4917408	,	282.6937738	,	286.8106734	,	290.0127065	,	293.6721728	,	297.3316391	,	300.9911055	,	304.1931385	,	307.3951715	,	310.5972046	,	314.2566709	,	317.4587039	,	320.660737	,	323.4053367	,	327.064803	,	330.7242694	,	333.9263024	,	336.6709022	,	339.8729352	,	343.0749682	,	346.2770013	,	349.4790343	,	352.6810673	,	354.9682338	,	357.2554003	,	360 ];
yp11 = [ 14	,	13.97350993	,	13.94701987	,	13.91390728	,	13.8807947	,	13.84768212	,	13.81456954	,	13.78145695	,	13.74834437	,	13.71523179	,	13.68874172	,	13.66225166	,	13.62913907	,	13.59602649	,	13.56291391	,	13.52980132	,	13.49668874	,	13.45695364	,	13.41721854	,	13.38410596	,	13.35099338	,	13.31788079	,	13.28476821	,	13.25165563	,	13.21854305	,	13.18543046	,	13.15231788	,	13.1192053	,	13.08609272	,	13.05298013	,	13.01986755	,	12.98675497	,	12.95364238	,	12.9205298	,	12.88741722	,	12.86092715	,	12.82781457	,	12.79470199	,	12.76821192	,	12.74172185	,	12.70860927	,	12.68211921	,	12.64900662	,	12.62251656	,	12.59602649	,	12.57615894	,	12.55629139	,	12.53642384	,	12.51655629	,	12.50331126	,	12.49006623	,	12.47682119	,	12.47019868	,	12.47019868	,	12.47019868	,	12.47019868	,	12.47682119	,	12.49006623	,	12.50331126	,	12.52317881	,	12.54304636	,	12.56953642	,	12.59602649	,	12.62251656	,	12.64900662	,	12.68211921	,	12.71523179	,	12.74834437	,	12.78145695	,	12.81456954	,	12.84768212	,	12.87417219	,	12.91390728	,	12.94701987	,	12.98675497	,	13.02649007	,	13.06622517	,	13.09933775	,	13.13907285	,	13.17218543	,	13.20529801	,	13.24503311	,	13.2781457	,	13.31788079	,	13.36423841	,	13.40397351	,	13.43708609	,	13.47682119	,	13.50993377	,	13.54966887	,	13.59602649	,	13.63576159	,	13.67549669	,	13.71523179	,	13.74834437	,	13.78807947	,	13.82119205	,	13.85430464	,	13.88741722	,	13.92715232	,	13.96688742	,	14	,	14.03311258	,	14.06622517	,	14.09933775	,	14.13245033	,	14.1589404	,	14.18543046	,	14.21192053	,	14.23178808	,	14.25165563	,	14.27152318	,	14.28476821	,	14.29139073	,	14.29801325	,	14.29801325	,	14.29801325	,	14.29139073	,	14.29139073	,	14.2781457	,	14.26490066	,	14.25165563	,	14.23178808	,	14.21854305	,	14.1986755	,	14.17218543	,	14.14569536	,	14.1192053	,	14.09271523	,	14.06622517	,	14.0397351	,	14.01324503 ];

%vetores para plot da curva referência condutor 12
xp12 = [ 0	,	1.358490566	,	2.716981132	,	4.075471698	,	5.886792453	,	7.698113208	,	9.509433962	,	11.32075472	,	13.58490566	,	15.8490566	,	18.11320755	,	20.37735849	,	22.64150943	,	24.45283019	,	26.71698113	,	28.98113208	,	31.69811321	,	34.41509434	,	37.58490566	,	40.75471698	,	44.37735849	,	47.54716981	,	50.71698113	,	53.88679245	,	57.50943396	,	60.67924528	,	63.8490566	,	67.47169811	,	70.18867925	,	73.35849057	,	76.0754717	,	78.79245283	,	81.50943396	,	84.22641509	,	86.49056604	,	88.75471698	,	91.47169811	,	93.28301887	,	95.54716981	,	97.35849057	,	99.16981132	,	101.4339623	,	103.245283	,	105.0566038	,	106.4150943	,	108.2264151	,	110.0377358	,	111.8490566	,	113.2075472	,	115.0188679	,	116.3773585	,	118.1886792	,	120	,	121.3584906	,	123.1698113	,	124.9811321	,	126.7924528	,	128.1509434	,	129.9622642	,	131.3207547	,	133.1320755	,	134.9433962	,	136.754717	,	138.5660377	,	140.3773585	,	142.1886792	,	144	,	145.8113208	,	147.1698113	,	148.9811321	,	150.3396226	,	152.1509434	,	153.9622642	,	156.2264151	,	157.5849057	,	159.3962264	,	161.6603774	,	163.9245283	,	165.7358491	,	167.5471698	,	169.3584906	,	171.1698113	,	172.9811321	,	174.7924528	,	177.0566038	,	179.3207547	,	181.5849057	,	183.3962264	,	185.2075472	,	187.0188679	,	189.2830189	,	191.5471698	,	193.8113208	,	196.5283019	,	199.245283	,	201.509434	,	203.7735849	,	206.490566	,	210.1132075	,	213.7358491	,	217.3584906	,	220.5283019	,	223.6981132	,	227.3207547	,	230.490566	,	234.1132075	,	238.1886792	,	240.9056604	,	244.0754717	,	247.6981132	,	250.4150943	,	252.6792453	,	254.9433962	,	257.2075472	,	259.4716981	,	261.7358491	,	264	,	266.2641509	,	268.0754717	,	270.3396226	,	272.6037736	,	274.4150943	,	276.2264151	,	278.0377358	,	279.8490566	,	281.6603774	,	283.4716981	,	285.2830189	,	287.5471698	,	289.8113208	,	291.6226415	,	293.4339623	,	295.245283	,	297.0566038	,	298.8679245	,	300.6792453	,	302.490566	,	303.8490566	,	305.2075472	,	306.5660377	,	308.3773585	,	310.1886792	,	311.5471698	,	312.9056604	,	314.7169811	,	316.5283019	,	317.8867925	,	319.245283	,	321.0566038	,	322.8679245	,	324.2264151	,	325.5849057	,	327.3962264	,	329.2075472	,	330.5660377	,	331.9245283	,	333.7358491	,	335.5471698	,	337.3584906	,	339.1698113	,	340.9811321	,	343.245283	,	345.0566038	,	347.3207547	,	349.1320755	,	350.9433962	,	352.3018868	,	354.1132075	,	355.9245283	,	358.1886792	,	360 ];
yp12 = [ 13.46229508	,	13.44590164	,	13.42622951	,	13.40983607	,	13.39016393	,	13.3704918	,	13.35081967	,	13.33114754	,	13.31147541	,	13.29180328	,	13.27213115	,	13.2557377	,	13.23934426	,	13.22295082	,	13.20983607	,	13.19344262	,	13.18360656	,	13.1704918	,	13.15737705	,	13.15081967	,	13.1442623	,	13.14098361	,	13.13770492	,	13.13770492	,	13.14098361	,	13.15081967	,	13.16393443	,	13.17704918	,	13.19344262	,	13.20655738	,	13.22295082	,	13.24262295	,	13.26229508	,	13.28196721	,	13.30491803	,	13.32786885	,	13.34754098	,	13.3704918	,	13.39344262	,	13.41639344	,	13.43606557	,	13.45901639	,	13.48196721	,	13.50491803	,	13.52459016	,	13.55081967	,	13.57704918	,	13.59672131	,	13.61639344	,	13.63606557	,	13.65901639	,	13.68196721	,	13.70491803	,	13.72786885	,	13.75409836	,	13.78032787	,	13.80655738	,	13.83278689	,	13.85245902	,	13.87540984	,	13.90163934	,	13.93114754	,	13.95737705	,	13.98688525	,	14.01639344	,	14.04262295	,	14.06885246	,	14.09508197	,	14.12131148	,	14.1442623	,	14.1704918	,	14.19672131	,	14.21967213	,	14.24918033	,	14.27540984	,	14.30491803	,	14.33114754	,	14.35737705	,	14.38688525	,	14.41311475	,	14.43934426	,	14.46557377	,	14.49180328	,	14.5147541	,	14.54098361	,	14.56721311	,	14.59016393	,	14.61311475	,	14.63606557	,	14.65901639	,	14.68196721	,	14.70491803	,	14.72459016	,	14.74754098	,	14.7704918	,	14.78688525	,	14.80655738	,	14.82622951	,	14.84262295	,	14.85901639	,	14.87540984	,	14.8852459	,	14.89508197	,	14.90163934	,	14.90491803	,	14.90491803	,	14.89508197	,	14.88852459	,	14.87868852	,	14.86557377	,	14.85245902	,	14.83606557	,	14.82295082	,	14.80983607	,	14.79344262	,	14.77704918	,	14.75737705	,	14.73770492	,	14.71803279	,	14.69836066	,	14.67868852	,	14.65901639	,	14.63606557	,	14.61639344	,	14.59344262	,	14.57377049	,	14.55081967	,	14.52459016	,	14.50163934	,	14.47540984	,	14.44918033	,	14.42295082	,	14.39672131	,	14.3704918	,	14.34754098	,	14.32131148	,	14.29508197	,	14.27540984	,	14.25245902	,	14.2295082	,	14.20983607	,	14.18688525	,	14.16393443	,	14.13770492	,	14.11803279	,	14.09180328	,	14.06885246	,	14.04590164	,	14.02295082	,	14	,	13.97704918	,	13.95081967	,	13.92459016	,	13.89836066	,	13.87213115	,	13.84918033	,	13.82622951	,	13.8	,	13.7704918	,	13.74098361	,	13.7147541	,	13.6852459	,	13.66229508	,	13.63606557	,	13.60983607	,	13.58688525	,	13.56393443	,	13.5442623	,	13.52131148	,	13.49836066	,	13.47868852 ];

figure();
plot(thetad,Erms,'LineWidth',4.5);
ylabel('Campo elétrico (kV/cm)')
xlabel('Ângulo (°)')
%title('Campo elétrico superficial condutor três') %título (trocar ao mudar condutor avaliado)
grid
hold on
% yy1 = smooth(xp1,yp1,'rloess'); %comando suavizar a curva referência condutor 1
% yy2 = smooth(xp2,yp2,'rloess'); %comando suavizar a curva referência condutor 2
% yy3 = smooth(xp3,yp3,'rloess'); %comando suavizar a curva referência condutor 3
% yy4 = smooth(xp4,yp4,'rloess'); %comando suavizar a curva referência condutor 4
% yy5 = smooth(xp5,yp5,'rloess'); %comando suavizar a curva referência condutor 5
% yy6 = smooth(xp6,yp6,'rloess'); %comando suavizar a curva referência condutor 6
% yy7 = smooth(xp7,yp7,'rloess'); %comando suavizar a curva referência condutor 7
% yy8 = smooth(xp8,yp8,'rloess'); %comando suavizar a curva referência condutor 8
% yy9 = smooth(xp9,yp9,'rloess'); %comando suavizar a curva referência condutor 9
% yy10 = smooth(xp10,yp10,'rloess'); %comando suavizar a curva referência condutor 10
% yy11 = smooth(xp11,yp11,'rloess'); %comando suavizar a curva referência condutor 11
% yy12 = smooth(xp12,yp12,'rloess'); %comando suavizar a curva referência condutor 12

% plot(xp1,yy1,'--','LineWidth',2.5); %plot curva referência (trocar variáveis ao mudar condutor avaliado)
plot(xp3,yp3,'--','LineWidth',2.5); %plot curva referência (trocar variáveis ao mudar condutor avaliado)

xlim([0 360])
legend('Simulado', '(PAGANOTTI, 2012)')
% 
% figure();
% plot(xr,yr,'o','color','blue')
% legend('Cabos condutores fase A, B e C');
% xlabel('Faixa de passagem (m)')
% ylabel('Altura do cabo (m)')
% %title('Configuração geométrica dos condutores')
% grid
% xlim([-8.5 8.5])
% ylim([17 27])

err = abs((yp3 - Erms)/yp3)*100;