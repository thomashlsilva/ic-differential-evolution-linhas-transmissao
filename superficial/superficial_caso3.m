clear all; close all; clc

% Cálculo do campo elétrico superficial do caso 3 do sistema trifásico

e_0 = 8.854*(10^(-12));
r = 14.37*(10^-3); % raio do condutor
n = 9; % número de condutores do sistema
nc = n*2; % número de condutores total (com imagens)
ci = (2*n-1); % número de cargas imagens por condutor


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

%tensão de operação máxima do sistema - neste caso SANTOS usou a nominal
V = 500*10^3;

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

x = [-10.4785 -10.25 -10.0215 -0.2285 0 0.2285 10.0215 10.25 10.4785 -10.4785 -10.25 -10.0215 -0.2285 0 0.2285 10.0215 10.25 10.4785]; % posição x dos condutores (com imagens)
y = [16.532 16.728 16.532 16.532 16.728 16.532 16.532 16.728 16.532 -16.532 -16.728 -16.532 -16.532 -16.728 -16.532 -16.532 -16.728 -16.532]; % posição y dos condutores (com imagens)

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

%% Pontos de avaliação (toda superfície do condutor)

theta = linspace(0,2*pi,142); %gera a superfície do condutor em 360 pontos

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

%vetores para plot da curva referência condutor 1
xp1 = [ 0.48582996	,	2.914979757	,	5.829959514	,	8.744939271	,	12.63157895	,	17.00404858	,	21.86234818	,	25.74898785	,	30.60728745	,	35.46558704	,	40.32388664	,	44.69635628	,	48.58299595	,	51.98380567	,	55.38461538	,	58.7854251	,	61.2145749	,	64.12955466	,	67.04453441	,	68.98785425	,	70.93117409	,	73.36032389	,	75.78947368	,	77.73279352	,	80.16194332	,	82.59109312	,	84.53441296	,	86.47773279	,	88.90688259	,	90.85020243	,	92.79352227	,	95.22267206	,	96.68016194	,	98.62348178	,	101.0526316	,	102.9959514	,	104.9392713	,	107.3684211	,	109.7975709	,	111.7408907	,	113.6842105	,	116.1133603	,	118.0566802	,	120	,	121.9433198	,	123.8866397	,	126.8016194	,	129.2307692	,	131.659919	,	134.0890688	,	136.5182186	,	138.9473684	,	140.8906883	,	142.8340081	,	145.7489879	,	148.1781377	,	150.6072874	,	153.0364372	,	155.465587	,	158.3805668	,	161.7813765	,	164.6963563	,	168.097166	,	170.5263158	,	173.9271255	,	176.8421053	,	180.242915	,	184.6153846	,	189.4736842	,	194.3319838	,	198.7044534	,	202.5910931	,	206.9635628	,	211.8218623	,	215.708502	,	220.0809717	,	224.4534413	,	228.340081	,	232.2267206	,	235.1417004	,	237.5708502	,	240.48583	,	243.4008097	,	246.8016194	,	249.2307692	,	251.659919	,	254.5748988	,	257.0040486	,	259.4331984	,	262.8340081	,	264.7773279	,	267.2064777	,	269.6356275	,	272.5506073	,	274.9797571	,	277.4089069	,	280.3238866	,	282.7530364	,	284.6963563	,	286.6396761	,	288.582996	,	291.4979757	,	293.4412955	,	296.3562753	,	298.2995951	,	300.242915	,	302.1862348	,	303.6437247	,	306.0728745	,	308.5020243	,	310.4453441	,	312.8744939	,	314.8178138	,	317.7327935	,	320.1619433	,	321.6194332	,	324.048583	,	325.9919028	,	328.9068826	,	331.3360324	,	333.2793522	,	335.2226721	,	337.6518219	,	339.5951417	,	342.0242915	,	344.4534413	,	346.8825911	,	349.7975709	,	352.7125506	,	355.1417004	,	358.0566802 ];
yp1 = [ 13.00675325	,	12.96983302	,	12.9329128	,	12.89599258	,	12.8664564	,	12.83692022	,	12.81476809	,	12.81476809	,	12.82215213	,	12.84430427	,	12.88860853	,	12.94029685	,	12.99198516	,	13.05105751	,	13.11012987	,	13.16920223	,	13.22827458	,	13.29473098	,	13.37595547	,	13.43502783	,	13.49410019	,	13.56055659	,	13.63439703	,	13.69346939	,	13.76730983	,	13.83376623	,	13.90022263	,	13.96667904	,	14.03313544	,	14.09959184	,	14.17343228	,	14.23988868	,	14.29896104	,	14.37280148	,	14.43925788	,	14.52048237	,	14.59432282	,	14.67554731	,	14.76415584	,	14.83061224	,	14.90445269	,	14.97090909	,	15.03736549	,	15.11858998	,	15.19243043	,	15.26627087	,	15.34749536	,	15.4361039	,	15.52471243	,	15.60593692	,	15.67239332	,	15.73146568	,	15.79053803	,	15.85699443	,	15.93821892	,	16.00467532	,	16.06374768	,	16.12282004	,	16.18189239	,	16.25573284	,	16.31480519	,	16.3812616	,	16.447718	,	16.49202226	,	16.54371058	,	16.59539889	,	16.6470872	,	16.68400742	,	16.72831169	,	16.76523191	,	16.77261596	,	16.78	,	16.78	,	16.76523191	,	16.74307978	,	16.70615955	,	16.66185529	,	16.62493506	,	16.56586271	,	16.52155844	,	16.47725417	,	16.42556586	,	16.35910946	,	16.30003711	,	16.24096475	,	16.19666048	,	16.13020408	,	16.06374768	,	15.98990724	,	15.91606679	,	15.86437848	,	15.79792208	,	15.70931354	,	15.64285714	,	15.56163265	,	15.47302412	,	15.39918367	,	15.31795918	,	15.23673469	,	15.1555102	,	15.08166976	,	14.99306122	,	14.92660482	,	14.84538033	,	14.77153989	,	14.70508349	,	14.64601113	,	14.57217069	,	14.48356215	,	14.40972171	,	14.35064935	,	14.26204082	,	14.19558442	,	14.11435993	,	14.04051948	,	13.96667904	,	13.90022263	,	13.83376623	,	13.76730983	,	13.69346939	,	13.61962894	,	13.56055659	,	13.50886827	,	13.44979592	,	13.39072356	,	13.31688312	,	13.25781076	,	13.1987384	,	13.15443414	,	13.10274583	,	13.05105751 ];

%vetores para plot da curva referência condutor 2
xp2 = [ 0.240639255	,	2.685607568	,	4.640112666	,	6.104154548	,	8.058659645	,	9.032238312	,	11.47628815	,	13.43079325	,	14.89391666	,	16.84842176	,	18.80292686	,	20.75835043	,	22.22239231	,	24.17506047	,	26.12956556	,	28.08498913	,	30.5281205	,	32.48078866	,	34.43621223	,	36.39071733	,	38.34522242	,	40.78927227	,	43.23332211	,	45.67920889	,	47.63463246	,	50.56914552	,	53.50641399	,	55.9513823	,	59.37727704	,	62.80409026	,	65.74044025	,	69.65863515	,	74.0663748	,	77.49686189	,	81.4178122	,	85.33968098	,	89.26246824	,	92.69754768	,	96.13262713	,	100.0600067	,	103.4978416	,	106.934758	,	109.8821296	,	112.3381196	,	115.2882466	,	118.2365368	,	120.6962006	,	124.1358724	,	127.0859994	,	129.5447448	,	132.4957903	,	134.4640725	,	136.9237363	,	138.892937	,	140.8621376	,	142.8313382	,	144.3082387	,	146.2792762	,	149.2312402	,	151.2004409	,	153.6628601	,	155.6320607	,	157.6003429	,	159.0790803	,	160.5559808	,	162.5251814	,	164.0039188	,	166.4635826	,	167.9423201	,	169.9115207	,	171.8825582	,	173.3612957	,	174.8391146	,	177.3006154	,	178.7775158	,	179.7648716	,	181.7340722	,	182.7205094	,	184.6897101	,	186.6589107	,	188.1385666	,	189.6163855	,	191.0942045	,	192.5720234	,	195.0335242	,	197.0027248	,	198.4832991	,	200.9447999	,	202.9121636	,	205.3755013	,	207.3456204	,	209.8062027	,	212.7572483	,	214.7273674	,	216.696568	,	219.1562318	,	221.6158957	,	225.0574044	,	227.0256866	,	228.9948872	,	232.4354775	,	236.3674494	,	238.8252763	,	241.7744849	,	244.722775	,	247.6701466	,	252.5802896	,	256.5067508	,	260.9245936	,	264.8492178	,	268.7729235	,	274.6575636	,	279.5594403	,	283.9690169	,	287.3976671	,	291.3167805	,	295.2349754	,	297.6808621	,	301.1076754	,	304.0449438	,	306.0003674	,	308.9367174	,	311.3826042	,	313.8275725	,	316.7620855	,	319.2070539	,	321.161559	,	323.1169825	,	325.0705691	,	327.0259927	,	328.9804978	,	330.9359214	,	333.3790528	,	334.8412577	,	337.7748523	,	339.2388942	,	341.6820255	,	344.6156201	,	346.0787435	,	348.0323302	,	350.4754615	,	351.9395034	,	353.8940085	,	355.8457582	,	357.8011818	,	358.7756789 ];
yp2 = [ 13.01973783	,	13.07531835	,	13.13089888	,	13.1864794	,	13.24205993	,	13.29764045	,	13.36016854	,	13.41574906	,	13.47827715	,	13.53385768	,	13.5894382	,	13.63807116	,	13.69365169	,	13.76312734	,	13.81870787	,	13.86734082	,	13.93681648	,	14.00629213	,	14.05492509	,	14.11050562	,	14.16608614	,	14.22861423	,	14.29114232	,	14.33977528	,	14.38840824	,	14.45093633	,	14.49262172	,	14.54820225	,	14.60378277	,	14.65241573	,	14.70104869	,	14.74273408	,	14.79136704	,	14.81220974	,	14.83305243	,	14.84694757	,	14.85389513	,	14.84	,	14.82610487	,	14.79831461	,	14.76357678	,	14.73578652	,	14.70104869	,	14.67325843	,	14.6176779	,	14.57599251	,	14.52041199	,	14.47177903	,	14.4161985	,	14.36756554	,	14.30503745	,	14.25640449	,	14.20082397	,	14.14524345	,	14.08966292	,	14.0340824	,	13.992397	,	13.92292135	,	13.85344569	,	13.79786517	,	13.72144195	,	13.66586142	,	13.61722846	,	13.56164794	,	13.51996255	,	13.46438202	,	13.4088015	,	13.35322097	,	13.29764045	,	13.24205993	,	13.17258427	,	13.11700375	,	13.06837079	,	12.99889513	,	12.95720974	,	12.90857678	,	12.85299625	,	12.81131086	,	12.75573034	,	12.70014981	,	12.63762172	,	12.58898876	,	12.54035581	,	12.49172285	,	12.42224719	,	12.36666667	,	12.29719101	,	12.22771536	,	12.18602996	,	12.10265918	,	12.04013109	,	11.977603	,	11.91507491	,	11.85254682	,	11.79696629	,	11.74138577	,	11.68580524	,	11.62327715	,	11.57464419	,	11.51906367	,	11.46348315	,	11.40095506	,	11.35926966	,	11.3106367	,	11.26895131	,	11.23421348	,	11.19252809	,	11.17168539	,	11.14389513	,	11.13694757	,	11.13694757	,	11.14389513	,	11.16473783	,	11.19947566	,	11.23421348	,	11.26895131	,	11.3106367	,	11.35926966	,	11.40790262	,	11.44958801	,	11.49822097	,	11.54685393	,	11.59548689	,	11.65106742	,	11.71359551	,	11.76917603	,	11.82475655	,	11.87338951	,	11.9359176	,	11.98455056	,	12.04013109	,	12.08876404	,	12.1582397	,	12.22771536	,	12.29719101	,	12.35277154	,	12.42224719	,	12.49172285	,	12.55425094	,	12.61677903	,	12.68625468	,	12.74183521	,	12.79741573	,	12.87383895	,	12.92247191	,	12.97110487 ];

%vetores para plot da curva referência condutor 3
xp3 = [ 0	,	1.980742779	,	4.951856946	,	6.932599725	,	9.408528198	,	11.88445667	,	14.85557084	,	17.82668501	,	20.79779917	,	22.77854195	,	25.74965612	,	27.7303989	,	29.71114168	,	32.68225585	,	34.66299862	,	36.6437414	,	39.11966988	,	41.59559835	,	43.08115543	,	45.06189821	,	47.04264099	,	48.52819807	,	50.50894085	,	52.48968363	,	54.47042641	,	56.45116919	,	57.93672627	,	59.91746905	,	61.89821183	,	63.87895461	,	65.36451169	,	67.34525447	,	68.83081155	,	70.31636864	,	72.29711142	,	74.2778542	,	76.25859697	,	78.23933975	,	79.72489684	,	81.70563961	,	83.68638239	,	85.17193948	,	86.65749656	,	88.63823934	,	90.12379642	,	92.1045392	,	94.08528198	,	96.56121045	,	98.04676754	,	100.0275103	,	102.0082531	,	103.9889959	,	106.4649243	,	107.9504814	,	109.9312242	,	111.911967	,	113.8927098	,	115.8734525	,	118.8445667	,	121.3204952	,	123.7964237	,	127.2627235	,	129.2434663	,	132.2145805	,	135.1856946	,	138.1568088	,	141.6231087	,	144.5942228	,	148.0605227	,	153.0123796	,	156.9738652	,	161.4305365	,	165.392022	,	168.8583219	,	171.829436	,	175.7909216	,	179.7524072	,	183.218707	,	186.1898212	,	189.656121	,	192.1320495	,	195.1031637	,	198.0742779	,	200.5502063	,	202.5309491	,	204.5116919	,	206.9876204	,	208.9683631	,	211.4442916	,	213.4250344	,	215.4057772	,	217.8817056	,	219.8624484	,	221.8431912	,	223.823934	,	225.8046768	,	228.2806052	,	230.261348	,	232.7372765	,	234.7180193	,	237.1939477	,	238.6795048	,	240.6602476	,	242.6409904	,	245.1169188	,	246.6024759	,	248.5832187	,	250.5639615	,	253.03989	,	255.0206327	,	257.4965612	,	259.477304	,	261.4580468	,	263.9339752	,	266.4099037	,	268.8858322	,	271.3617607	,	273.3425034	,	276.3136176	,	278.2943604	,	281.7606602	,	284.2365887	,	287.2077029	,	289.6836314	,	293.1499312	,	296.6162311	,	300.0825309	,	303.0536451	,	306.0247593	,	309.4910591	,	312.957359	,	317.4140303	,	321.8707015	,	325.8321871	,	329.7936726	,	333.7551582	,	337.7166437	,	343.6588721	,	347.6203576	,	351.5818432	,	355.5433287	,	359.0096286 ];
yp3 = [ 17.15348485	,	17.12136364	,	17.09727273	,	17.05712121	,	17.0169697	,	16.97681818	,	16.92060606	,	16.86439394	,	16.80015152	,	16.7519697	,	16.68772727	,	16.63954545	,	16.56727273	,	16.5030303	,	16.43075758	,	16.35848485	,	16.30227273	,	16.2380303	,	16.18984848	,	16.12560606	,	16.06136364	,	15.99712121	,	15.94090909	,	15.86863636	,	15.80439394	,	15.73212121	,	15.67590909	,	15.61166667	,	15.54742424	,	15.47515152	,	15.41893939	,	15.35469697	,	15.29045455	,	15.22621212	,	15.15393939	,	15.07363636	,	15.00136364	,	14.93712121	,	14.87287879	,	14.79257576	,	14.72833333	,	14.66409091	,	14.59181818	,	14.51954545	,	14.45530303	,	14.3830303	,	14.32681818	,	14.25454545	,	14.19833333	,	14.13409091	,	14.06984848	,	14.00560606	,	13.92530303	,	13.87712121	,	13.81287879	,	13.74863636	,	13.68439394	,	13.62015152	,	13.56393939	,	13.50772727	,	13.44348485	,	13.37924242	,	13.315	,	13.26681818	,	13.21060606	,	13.17045455	,	13.13833333	,	13.10621212	,	13.07409091	,	13.05	,	13.0419697	,	13.05	,	13.07409091	,	13.10621212	,	13.13833333	,	13.17848485	,	13.24272727	,	13.29090909	,	13.35515152	,	13.42742424	,	13.48363636	,	13.55590909	,	13.62015152	,	13.68439394	,	13.74863636	,	13.81287879	,	13.87712121	,	13.94136364	,	14.01363636	,	14.08590909	,	14.16621212	,	14.23848485	,	14.31075758	,	14.3669697	,	14.43924242	,	14.51954545	,	14.59984848	,	14.67212121	,	14.75242424	,	14.83272727	,	14.905	,	14.97727273	,	15.04954545	,	15.12984848	,	15.19409091	,	15.25833333	,	15.33060606	,	15.41090909	,	15.49121212	,	15.57151515	,	15.65181818	,	15.71606061	,	15.80439394	,	15.88469697	,	15.965	,	16.04530303	,	16.11757576	,	16.18984848	,	16.27015152	,	16.35045455	,	16.43878788	,	16.51106061	,	16.57530303	,	16.65560606	,	16.72787879	,	16.80015152	,	16.88045455	,	16.93666667	,	16.99287879	,	17.05712121	,	17.12136364	,	17.16954545	,	17.21772727	,	17.24984848	,	17.2819697	,	17.2819697	,	17.29	,	17.2819697	,	17.27393939	,	17.24181818	,	17.20969697	,	17.17757576 ];

%vetores para plot da curva referência condutor 4
xp4 = [ 0	,	1.983471074	,	4.958677686	,	7.933884298	,	10.90909091	,	13.38842975	,	16.85950413	,	20.82644628	,	23.30578512	,	27.27272727	,	31.23966942	,	35.20661157	,	37.68595041	,	41.15702479	,	43.63636364	,	46.61157025	,	49.58677686	,	52.56198347	,	55.53719008	,	58.51239669	,	60.49586777	,	62.97520661	,	65.45454545	,	67.9338843	,	69.91735537	,	72.39669421	,	74.38016529	,	76.85950413	,	78.84297521	,	80.82644628	,	82.80991736	,	83.80165289	,	85.78512397	,	87.76859504	,	89.75206612	,	91.23966942	,	93.71900826	,	95.70247934	,	97.68595041	,	100.1652893	,	102.1487603	,	103.1404959	,	104.6280992	,	106.6115702	,	108.0991736	,	110.5785124	,	112.5619835	,	114.0495868	,	116.5289256	,	119.0082645	,	120.4958678	,	122.4793388	,	124.9586777	,	126.446281	,	128.4297521	,	129.9173554	,	131.9008264	,	134.3801653	,	136.3636364	,	138.3471074	,	140.3305785	,	142.8099174	,	144.7933884	,	146.7768595	,	150.2479339	,	152.231405	,	155.2066116	,	158.1818182	,	160.661157	,	164.1322314	,	167.107438	,	170.0826446	,	174.0495868	,	177.0247934	,	180.4958678	,	184.9586777	,	187.9338843	,	191.9008264	,	195.8677686	,	200.8264463	,	205.2892562	,	209.2561983	,	213.2231405	,	217.1900826	,	220.661157	,	224.6280992	,	227.6033058	,	230.5785124	,	233.553719	,	236.5289256	,	239.5041322	,	242.4793388	,	245.9504132	,	248.9256198	,	251.4049587	,	253.3884298	,	255.8677686	,	258.3471074	,	260.8264463	,	262.8099174	,	265.2892562	,	267.768595	,	270.2479339	,	272.231405	,	274.214876	,	276.1983471	,	278.1818182	,	279.6694215	,	281.1570248	,	283.1404959	,	285.1239669	,	286.6115702	,	288.5950413	,	291.0743802	,	293.553719	,	295.0413223	,	297.0247934	,	299.0082645	,	300	,	301.9834711	,	303.4710744	,	305.4545455	,	306.9421488	,	308.9256198	,	311.4049587	,	313.3884298	,	315.3719008	,	317.3553719	,	318.8429752	,	320.8264463	,	322.3140496	,	324.7933884	,	327.2727273	,	329.7520661	,	331.2396694	,	333.2231405	,	335.7024793	,	338.1818182	,	340.661157	,	342.6446281	,	344.6280992	,	347.6033058	,	350.5785124	,	353.553719	,	356.0330579	,	358.5123967 ];
yp4 = [ 14.22070076	,	14.18698864	,	14.15327652	,	14.11113636	,	14.08585227	,	14.05214015	,	14.02685606	,	14.01842803	,	14.01842803	,	14.02685606	,	14.03528409	,	14.06056818	,	14.08585227	,	14.11956439	,	14.17013258	,	14.22070076	,	14.27126894	,	14.32183712	,	14.38926136	,	14.45668561	,	14.52410985	,	14.59153409	,	14.65895833	,	14.71795455	,	14.78537879	,	14.85280303	,	14.9286553	,	14.99607955	,	15.06350379	,	15.13092803	,	15.19835227	,	15.24892045	,	15.3163447	,	15.400625	,	15.46804924	,	15.53547348	,	15.62818182	,	15.70403409	,	15.77988636	,	15.86416667	,	15.93159091	,	15.99058712	,	16.05801136	,	16.12543561	,	16.20971591	,	16.27714015	,	16.36984848	,	16.43727273	,	16.52155303	,	16.61426136	,	16.69011364	,	16.75753788	,	16.83339015	,	16.89238636	,	16.97666667	,	17.03566288	,	17.11151515	,	17.19579545	,	17.25479167	,	17.33064394	,	17.40649621	,	17.47392045	,	17.53291667	,	17.60034091	,	17.68462121	,	17.75204545	,	17.82789773	,	17.90375	,	17.96274621	,	18.03017045	,	18.0975947	,	18.16501894	,	18.22401515	,	18.28301136	,	18.33357955	,	18.3757197	,	18.40943182	,	18.44314394	,	18.46842803	,	18.46842803	,	18.46842803	,	18.46	,	18.43471591	,	18.40100379	,	18.36729167	,	18.30829545	,	18.27458333	,	18.23244318	,	18.16501894	,	18.11445076	,	18.05545455	,	17.97960227	,	17.91217803	,	17.86160985	,	17.79418561	,	17.72676136	,	17.65933712	,	17.59191288	,	17.53291667	,	17.45706439	,	17.38964015	,	17.30535985	,	17.22950758	,	17.1536553	,	17.07780303	,	17.01037879	,	16.94295455	,	16.86710227	,	16.81653409	,	16.74910985	,	16.68168561	,	16.61426136	,	16.54683712	,	16.46255682	,	16.37827652	,	16.31085227	,	16.235	,	16.15914773	,	16.10015152	,	16.03272727	,	15.97373106	,	15.90630682	,	15.83888258	,	15.7630303	,	15.68717803	,	15.61132576	,	15.54390152	,	15.46804924	,	15.39219697	,	15.3163447	,	15.25734848	,	15.18149621	,	15.09721591	,	15.02136364	,	14.95393939	,	14.89494318	,	14.83594697	,	14.7600947	,	14.69267045	,	14.60839015	,	14.55782197	,	14.49039773	,	14.42297348	,	14.35554924	,	14.30498106	,	14.27126894 ];

%vetores para plot da curva referência condutor 5
xp5 = [ 0	,	1.504178273	,	3.509749304	,	5.515320334	,	7.019498607	,	8.52367688	,	10.52924791	,	12.03342618	,	14.03899721	,	15.54317549	,	17.54874652	,	19.55431755	,	21.55988858	,	23.56545961	,	25.06963788	,	27.07520891	,	29.08077994	,	31.08635097	,	33.09192201	,	34.59610028	,	36.60167131	,	38.60724234	,	40.61281337	,	42.6183844	,	44.62395543	,	47.13091922	,	49.63788301	,	52.1448468	,	54.65181058	,	57.66016713	,	60.16713092	,	63.17548747	,	65.1810585	,	68.18941504	,	71.19777159	,	74.70752089	,	77.71587744	,	81.22562674	,	84.23398329	,	88.24512535	,	91.75487465	,	95.76601671	,	98.77437326	,	103.2869081	,	106.2952646	,	109.8050139	,	112.3119777	,	114.8189415	,	118.8300836	,	121.3370474	,	123.3426184	,	125.8495822	,	128.356546	,	130.362117	,	132.367688	,	134.8746518	,	137.3816156	,	139.8885794	,	141.8941504	,	143.8997214	,	145.9052925	,	147.9108635	,	149.9164345	,	152.4233983	,	153.9275766	,	155.9331476	,	157.9387187	,	159.4428969	,	161.448468	,	162.9526462	,	164.9582173	,	167.4651811	,	169.4707521	,	171.9777159	,	173.4818942	,	174.9860724	,	176.4902507	,	177.994429	,	180.5013928	,	182.005571	,	184.0111421	,	186.0167131	,	188.0222841	,	190.5292479	,	192.0334262	,	194.54039	,	196.0445682	,	198.0501393	,	200.0557103	,	202.0612813	,	204.0668524	,	206.5738162	,	208.0779944	,	210.0835655	,	212.5905292	,	214.5961003	,	216.6016713	,	219.6100279	,	222.1169916	,	224.6239554	,	226.6295265	,	229.1364903	,	232.1448468	,	235.1532033	,	238.1615599	,	241.1699164	,	244.178273	,	247.6880223	,	252.2005571	,	256.2116992	,	260.2228412	,	264.2339833	,	268.2451253	,	272.2562674	,	276.7688022	,	280.7799443	,	284.7910864	,	288.8022284	,	292.8133705	,	296.3231198	,	299.3314763	,	302.3398329	,	305.3481894	,	308.356546	,	311.3649025	,	314.3732591	,	316.3788301	,	318.8857939	,	320.8913649	,	323.3983287	,	325.9052925	,	327.9108635	,	330.4178273	,	332.4233983	,	334.4289694	,	336.4345404	,	338.9415042	,	341.448468	,	343.454039	,	345.45961	,	347.9665738	,	349.9721448	,	351.9777159	,	353.9832869	,	355.9888579	,	357.4930362	,	359.4986072 ];
yp5 = [ 14.06573614	,	14.10426386	,	14.16590822	,	14.22755258	,	14.28919694	,	14.3508413	,	14.39707457	,	14.45871893	,	14.53577438	,	14.59741874	,	14.6590631	,	14.71300191	,	14.78235182	,	14.84399618	,	14.90564054	,	14.96728489	,	15.02892925	,	15.08286807	,	15.14451243	,	15.19845124	,	15.2600956	,	15.31403442	,	15.37567878	,	15.43732314	,	15.49126195	,	15.55290631	,	15.61455067	,	15.66848948	,	15.71472275	,	15.77636711	,	15.83030593	,	15.8765392	,	15.91506692	,	15.96130019	,	15.99212237	,	16.0306501	,	16.06147228	,	16.07688337	,	16.09229446	,	16.10770554	,	16.10770554	,	16.1	,	16.08458891	,	16.04606119	,	16.01523901	,	15.97671128	,	15.93818356	,	15.89965583	,	15.85342256	,	15.81489484	,	15.76866157	,	15.70701721	,	15.64537285	,	15.60684512	,	15.54520076	,	15.49126195	,	15.42961759	,	15.36797323	,	15.31403442	,	15.24468451	,	15.18304015	,	15.12910134	,	15.07516252	,	15.01351816	,	14.94416826	,	14.8825239	,	14.82858509	,	14.76694073	,	14.68988528	,	14.62053537	,	14.56659656	,	14.49724665	,	14.42789675	,	14.3508413	,	14.28919694	,	14.21984704	,	14.15820268	,	14.08885277	,	14.01179732	,	13.95015296	,	13.87309751	,	13.79604207	,	13.73439771	,	13.67275335	,	13.60340344	,	13.54175908	,	13.48011472	,	13.41076482	,	13.33370937	,	13.27206501	,	13.21042065	,	13.14877629	,	13.08713193	,	13.02548757	,	12.97154876	,	12.90219885	,	12.84055449	,	12.76349904	,	12.70956023	,	12.64791587	,	12.58627151	,	12.54003824	,	12.47839388	,	12.40904398	,	12.34739962	,	12.2934608	,	12.24722753	,	12.20869981	,	12.16246654	,	12.12393881	,	12.10082218	,	12.09311663	,	12.07770554	,	12.07770554	,	12.09311663	,	12.10852772	,	12.1393499	,	12.18558317	,	12.2241109	,	12.27804971	,	12.31657744	,	12.3782218	,	12.43216061	,	12.48609943	,	12.54774379	,	12.60938815	,	12.6710325	,	12.72497132	,	12.77891013	,	12.84055449	,	12.91760994	,	12.9792543	,	13.04860421	,	13.11795411	,	13.18730402	,	13.24894837	,	13.32600382	,	13.39535373	,	13.46470363	,	13.54175908	,	13.61881453	,	13.68045889	,	13.75751434	,	13.82686424	,	13.89621415	,	13.96556405	,	14.01950287 ];

%vetores para plot da curva referência condutor 6
xp6 = [ 0	,	3.554301834	,	5.585331453	,	8.631875882	,	11.17066291	,	14.21720733	,	16.75599436	,	18.78702398	,	20.8180536	,	23.86459803	,	26.40338505	,	28.43441467	,	30.46544429	,	32.49647391	,	35.03526093	,	37.06629055	,	38.58956276	,	40.11283498	,	42.1438646	,	44.68265162	,	46.71368124	,	48.74471086	,	50.77574048	,	52.8067701	,	54.33004231	,	56.36107193	,	58.39210155	,	60.42313117	,	61.43864598	,	62.96191819	,	64.99294781	,	67.02397743	,	68.54724965	,	70.57827927	,	72.60930889	,	74.1325811	,	76.67136812	,	78.70239774	,	80.73342736	,	82.76445698	,	84.7954866	,	86.82651622	,	88.85754584	,	90.38081805	,	92.41184767	,	94.9506347	,	97.48942172	,	99.52045134	,	101.551481	,	103.5825106	,	106.1212976	,	108.6600846	,	110.6911142	,	112.7221439	,	115.7686883	,	117.7997179	,	120.3385049	,	123.3850494	,	125.9238364	,	129.9858956	,	133.5401975	,	137.0944993	,	140.1410437	,	143.6953456	,	147.2496474	,	150.8039492	,	154.8660085	,	158.9280677	,	162.9901269	,	166.5444288	,	170.606488	,	174.1607898	,	177.2073343	,	180.2538787	,	183.3004231	,	186.854725	,	189.9012694	,	192.9478138	,	195.4866008	,	197.5176305	,	200.0564175	,	202.5952045	,	205.1339915	,	207.1650212	,	209.1960508	,	211.2270804	,	213.25811	,	215.2891396	,	217.3201693	,	219.3511989	,	221.3822285	,	222.9055007	,	224.4287729	,	226.4598025	,	228.4908322	,	230.0141044	,	231.5373766	,	234.0761636	,	235.5994358	,	237.6304654	,	239.1537377	,	240.6770099	,	243.2157969	,	245.2468265	,	247.7856135	,	249.8166432	,	252.3554302	,	254.3864598	,	255.909732	,	258.448519	,	259.9717913	,	262.0028209	,	264.5416079	,	266.5726375	,	269.1114245	,	271.1424542	,	273.1734838	,	275.7122708	,	278.2510578	,	280.2820874	,	282.8208745	,	285.3596615	,	288.4062059	,	290.9449929	,	293.48378	,	296.022567	,	299.0691114	,	302.1156559	,	305.1622003	,	307.7009873	,	310.7475317	,	314.809591	,	318.8716502	,	322.9337094	,	326.9957687	,	331.0578279	,	334.1043724	,	337.6586742	,	341.212976	,	345.7827927	,	349.8448519	,	353.3991537	,	355.9379408	,	359.4922426 ];
yp6 = [ 18.30416342	,	18.26953307	,	18.23490272	,	18.1829572	,	18.13101167	,	18.07040856	,	18.01846304	,	17.96651751	,	17.91457198	,	17.85396887	,	17.79336576	,	17.73276265	,	17.66350195	,	17.61155642	,	17.54229572	,	17.48169261	,	17.43840467	,	17.36914397	,	17.29988327	,	17.22196498	,	17.15270428	,	17.08344358	,	17.00552529	,	16.93626459	,	16.86700389	,	16.7890856	,	16.71116732	,	16.64190661	,	16.59861868	,	16.53801556	,	16.46875486	,	16.39949416	,	16.33023346	,	16.24365759	,	16.1657393	,	16.08782101	,	16.00124514	,	15.91466926	,	15.84540856	,	15.77614786	,	15.68957198	,	15.6116537	,	15.52507782	,	15.45581712	,	15.37789883	,	15.29132296	,	15.20474708	,	15.12682879	,	15.04891051	,	14.97099222	,	14.90173152	,	14.82381323	,	14.74589494	,	14.67663424	,	14.60737354	,	14.53811284	,	14.47750973	,	14.40824903	,	14.33033074	,	14.26107004	,	14.20912451	,	14.1485214	,	14.10523346	,	14.06194553	,	14.03597276	,	14.01	,	14.01	,	14.01	,	14.02731518	,	14.05328794	,	14.08791829	,	14.13120623	,	14.17449416	,	14.22643969	,	14.2870428	,	14.3563035	,	14.41690661	,	14.4948249	,	14.5640856	,	14.6333463	,	14.69394942	,	14.7718677	,	14.8411284	,	14.91904669	,	14.97964981	,	15.04891051	,	15.11817121	,	15.19608949	,	15.26535019	,	15.33461089	,	15.39521401	,	15.46447471	,	15.52507782	,	15.60299611	,	15.67225681	,	15.72420233	,	15.80212062	,	15.88003891	,	15.9579572	,	16.03587549	,	16.0964786	,	16.1657393	,	16.26097276	,	16.33023346	,	16.40815175	,	16.49472763	,	16.5813035	,	16.66787938	,	16.74579767	,	16.81505837	,	16.88431907	,	16.96223735	,	17.04015564	,	17.11807393	,	17.19599222	,	17.27391051	,	17.33451362	,	17.40377432	,	17.48169261	,	17.56826848	,	17.63752918	,	17.71544747	,	17.78470817	,	17.87128405	,	17.92322957	,	18.00114786	,	18.06175097	,	18.13101167	,	18.17429961	,	18.22624514	,	18.27819066	,	18.33013619	,	18.37342412	,	18.41671206	,	18.44268482	,	18.46865759	,	18.47731518	,	18.47731518	,	18.46865759	,	18.44268482	,	18.41671206	,	18.3907393	,	18.36476654	,	18.33013619 ];

%vetores para plot da curva referência condutor 7
xp7 = [ 0	,	1.477428181	,	3.447332421	,	5.909712722	,	8.372093023	,	11.32694938	,	14.77428181	,	18.22161423	,	22.16142271	,	26.10123119	,	30.04103967	,	33.98084815	,	37.92065663	,	40.875513	,	43.83036936	,	46.78522572	,	49.24760602	,	51.70998632	,	54.66484268	,	57.12722298	,	59.09712722	,	61.55950752	,	63.52941176	,	65.99179207	,	68.45417237	,	70.42407661	,	72.88645691	,	75.34883721	,	77.31874145	,	79.28864569	,	81.25854993	,	83.22845417	,	85.19835841	,	87.16826265	,	89.13816689	,	91.10807114	,	93.57045144	,	95.54035568	,	98.00273598	,	99.97264022	,	101.9425445	,	103.9124487	,	106.374829	,	108.8372093	,	110.8071135	,	112.7770178	,	115.2393981	,	117.2093023	,	119.1792066	,	121.6415869	,	123.6114911	,	125.0889193	,	127.5512996	,	129.5212038	,	131.9835841	,	133.9534884	,	135.9233926	,	138.3857729	,	140.3556772	,	142.3255814	,	144.7879617	,	147.7428181	,	150.6976744	,	153.1600547	,	156.1149111	,	159.0697674	,	162.0246238	,	164.9794802	,	168.4268126	,	171.874145	,	175.8139535	,	179.753762	,	183.6935705	,	187.6333789	,	191.5731874	,	195.0205198	,	198.9603283	,	203.3926129	,	207.3324213	,	211.2722298	,	215.2120383	,	219.1518468	,	223.0916553	,	227.5239398	,	230.4787962	,	232.9411765	,	235.4035568	,	237.8659371	,	240.8207934	,	243.7756498	,	246.2380301	,	249.1928865	,	251.6552668	,	253.625171	,	256.0875513	,	258.5499316	,	261.0123119	,	262.9822161	,	264.9521204	,	267.4145007	,	270.369357	,	272.3392613	,	274.3091655	,	276.2790698	,	278.248974	,	280.2188782	,	282.1887825	,	284.1586867	,	286.621067	,	289.0834473	,	291.0533516	,	293.0232558	,	294.500684	,	296.9630643	,	298.4404925	,	300.4103967	,	302.380301	,	304.3502052	,	306.3201094	,	307.7975376	,	309.7674419	,	311.7373461	,	313.7072503	,	315.1846785	,	317.6470588	,	319.6169631	,	321.5868673	,	323.5567715	,	325.5266758	,	327.49658	,	329.4664843	,	331.9288646	,	333.8987688	,	336.3611491	,	338.8235294	,	341.2859097	,	343.74829	,	346.7031464	,	348.6730506	,	351.1354309	,	353.5978112	,	356.0601915	,	358.5225718 ];
yp7 = [ 13.24924812	,	13.21736842	,	13.19345865	,	13.16157895	,	13.12969925	,	13.09781955	,	13.07390977	,	13.05796992	,	13.05	,	13.05796992	,	13.06593985	,	13.08984962	,	13.12969925	,	13.16157895	,	13.20142857	,	13.24924812	,	13.30503759	,	13.35285714	,	13.40067669	,	13.46443609	,	13.52022556	,	13.57601504	,	13.62383459	,	13.69556391	,	13.75932331	,	13.82308271	,	13.89481203	,	13.96654135	,	14.03827068	,	14.09406015	,	14.15781955	,	14.22954887	,	14.30924812	,	14.37300752	,	14.44473684	,	14.52443609	,	14.59616541	,	14.68383459	,	14.76353383	,	14.85120301	,	14.93090226	,	15.00263158	,	15.08233083	,	15.17	,	15.24969925	,	15.3293985	,	15.40909774	,	15.48879699	,	15.56052632	,	15.63225564	,	15.71195489	,	15.78368421	,	15.85541353	,	15.93511278	,	16.00684211	,	16.07857143	,	16.15030075	,	16.22203008	,	16.28578947	,	16.3575188	,	16.4212782	,	16.50097744	,	16.58067669	,	16.65240602	,	16.73210526	,	16.79586466	,	16.86759398	,	16.93932331	,	16.99511278	,	17.05887218	,	17.11466165	,	17.17045113	,	17.21030075	,	17.24218045	,	17.26609023	,	17.29	,	17.29796992	,	17.29796992	,	17.29	,	17.27406015	,	17.24218045	,	17.21030075	,	17.1624812	,	17.10669173	,	17.07481203	,	17.02699248	,	16.97120301	,	16.92338346	,	16.86759398	,	16.81180451	,	16.75601504	,	16.68428571	,	16.62849624	,	16.57270677	,	16.50894737	,	16.44518797	,	16.38142857	,	16.30172932	,	16.22203008	,	16.15827068	,	16.07857143	,	16.01481203	,	15.94308271	,	15.87135338	,	15.80759398	,	15.72789474	,	15.66413534	,	15.59240602	,	15.51270677	,	15.44097744	,	15.35330827	,	15.28954887	,	15.22578947	,	15.14609023	,	15.0743609	,	15.0106015	,	14.94684211	,	14.88308271	,	14.81932331	,	14.74759398	,	14.68383459	,	14.61210526	,	14.53240602	,	14.46864662	,	14.39691729	,	14.33315789	,	14.26142857	,	14.2056391	,	14.12593985	,	14.06218045	,	13.9824812	,	13.9187218	,	13.85496241	,	13.79917293	,	13.72744361	,	13.65571429	,	13.59195489	,	13.53616541	,	13.48037594	,	13.41661654	,	13.36879699	,	13.32097744	,	13.28112782 ];

%vetores para plot da curva referência condutor 8
xp8 = [ 0	,	2.036775106	,	4.073550212	,	6.110325318	,	8.147100424	,	10.18387553	,	12.22065064	,	14.25742574	,	15.78500707	,	17.82178218	,	19.34936351	,	20.87694484	,	22.91371994	,	24.95049505	,	26.98727016	,	29.02404526	,	31.06082037	,	33.09759547	,	35.13437058	,	37.17114569	,	39.71711457	,	41.75388967	,	44.29985856	,	46.84582744	,	48.88260255	,	51.93776521	,	54.48373409	,	57.02970297	,	59.57567185	,	62.12164074	,	65.17680339	,	68.23196605	,	72.30551627	,	75.8698727	,	79.43422914	,	82.99858557	,	86.56294201	,	89.61810467	,	92.16407355	,	95.72842999	,	99.29278642	,	102.8571429	,	105.9123055	,	109.476662	,	112.0226308	,	114.5685997	,	118.1329562	,	121.6973126	,	124.2432815	,	126.7892504	,	129.844413	,	132.3903819	,	134.427157	,	136.9731259	,	139.009901	,	141.5558699	,	143.592645	,	145.1202263	,	147.1570014	,	149.1937765	,	151.2305516	,	153.2673267	,	155.3041018	,	156.8316832	,	158.8684583	,	160.9052334	,	162.9420085	,	164.4695898	,	166.5063649	,	168.0339463	,	170.0707214	,	172.1074965	,	174.1442716	,	175.6718529	,	177.1994342	,	178.7270156	,	180.7637907	,	182.8005658	,	184.3281471	,	186.3649222	,	188.4016973	,	190.4384724	,	191.9660537	,	193.4936351	,	195.5304102	,	197.5671853	,	199.6039604	,	201.6407355	,	203.1683168	,	205.2050919	,	207.7510608	,	209.7878359	,	212.3338048	,	214.8797737	,	216.9165488	,	219.4625177	,	221.4992928	,	224.0452617	,	226.0820368	,	229.1371994	,	231.6831683	,	234.2291372	,	237.7934936	,	240.8486563	,	244.4130127	,	246.9589816	,	250.0141443	,	253.5785007	,	257.1428571	,	262.2347949	,	266.3083451	,	270.3818953	,	274.9646393	,	280.0565771	,	284.1301273	,	287.18529	,	290.7496464	,	293.2956153	,	296.8599717	,	300.4243281	,	303.4794908	,	306.5346535	,	309.5898161	,	312.135785	,	315.7001414	,	318.7553041	,	321.8104668	,	324.3564356	,	326.3932107	,	328.9391796	,	331.4851485	,	333.5219236	,	335.5586987	,	338.1046676	,	339.6322489	,	342.1782178	,	344.2149929	,	346.251768	,	348.2885431	,	349.8161245	,	351.8528996	,	353.3804809	,	355.417256	,	357.4540311	,	359.4908062 ];
yp8 = [ 12.91023033	,	12.97431862	,	13.03128599	,	13.09537428	,	13.15234165	,	13.21642994	,	13.28763916	,	13.34460653	,	13.4015739	,	13.46566219	,	13.51550864	,	13.57247601	,	13.62944338	,	13.68641075	,	13.74337812	,	13.81458733	,	13.87867562	,	13.93564299	,	13.99261036	,	14.04957774	,	14.11366603	,	14.1706334	,	14.22760077	,	14.28456814	,	14.34153551	,	14.39850288	,	14.45547025	,	14.51243762	,	14.55516315	,	14.6050096	,	14.65485605	,	14.69758157	,	14.73318618	,	14.76879079	,	14.79727447	,	14.82575816	,	14.84	,	14.84712092	,	14.84712092	,	14.84712092	,	14.82575816	,	14.81151631	,	14.78303263	,	14.74742802	,	14.72606526	,	14.68333973	,	14.64773512	,	14.59788868	,	14.54092131	,	14.49107486	,	14.44122841	,	14.38426104	,	14.33441459	,	14.27744722	,	14.22760077	,	14.17775432	,	14.12790787	,	14.06381958	,	14.01397313	,	13.96412668	,	13.90715931	,	13.85731286	,	13.79322457	,	13.74337812	,	13.68641075	,	13.62944338	,	13.56535509	,	13.50838772	,	13.45142035	,	13.4015739	,	13.34460653	,	13.27339731	,	13.2021881	,	13.14522073	,	13.08825336	,	13.03128599	,	12.97431862	,	12.9031094	,	12.85326296	,	12.78917466	,	12.71796545	,	12.66099808	,	12.60403071	,	12.54706334	,	12.49009597	,	12.42600768	,	12.36904031	,	12.29783109	,	12.25510557	,	12.20525912	,	12.1340499	,	12.06996161	,	12.00587332	,	11.94890595	,	11.89193858	,	11.83497121	,	11.77800384	,	11.72103647	,	11.6640691	,	11.59998081	,	11.53589251	,	11.47892514	,	11.42907869	,	11.37211132	,	11.32226488	,	11.28666027	,	11.24393474	,	11.20833013	,	11.17272553	,	11.15136276	,	11.12287908	,	11.12287908	,	11.12287908	,	11.14424184	,	11.16560461	,	11.17984645	,	11.21545106	,	11.25817658	,	11.29378119	,	11.33650672	,	11.37923225	,	11.43619962	,	11.49316699	,	11.55013436	,	11.60710173	,	11.67119002	,	11.73527831	,	11.80648752	,	11.86345489	,	11.92042226	,	11.98451056	,	12.04147793	,	12.10556622	,	12.16965451	,	12.24086372	,	12.30495202	,	12.36904031	,	12.4331286	,	12.49009597	,	12.54706334	,	12.60403071	,	12.668119	,	12.72508637	,	12.78917466	,	12.85326296 ];

%vetores para plot da curva referência condutor 9
xp9 = [ 0	,	2.479338843	,	5.454545455	,	8.429752066	,	11.40495868	,	13.88429752	,	16.36363636	,	18.34710744	,	20.82644628	,	23.30578512	,	25.78512397	,	28.26446281	,	30.24793388	,	32.72727273	,	34.21487603	,	36.69421488	,	38.18181818	,	40.16528926	,	42.6446281	,	45.12396694	,	47.10743802	,	49.09090909	,	51.07438017	,	53.05785124	,	55.04132231	,	56.52892562	,	58.01652893	,	60	,	61.98347107	,	63.96694215	,	65.95041322	,	67.9338843	,	69.91735537	,	71.90082645	,	73.88429752	,	75.8677686	,	77.85123967	,	79.83471074	,	81.81818182	,	83.30578512	,	85.2892562	,	87.27272727	,	89.75206612	,	91.73553719	,	93.71900826	,	95.70247934	,	97.68595041	,	99.66942149	,	101.1570248	,	103.6363636	,	105.6198347	,	108.0991736	,	111.0743802	,	113.0578512	,	115.0413223	,	117.5206612	,	120.4958678	,	123.4710744	,	126.446281	,	129.4214876	,	132.892562	,	135.8677686	,	139.338843	,	143.3057851	,	148.2644628	,	152.231405	,	156.1983471	,	160.661157	,	165.1239669	,	169.5867769	,	173.0578512	,	176.5289256	,	180.4958678	,	183.9669421	,	186.9421488	,	189.4214876	,	192.892562	,	195.3719008	,	198.3471074	,	200.8264463	,	203.8016529	,	206.2809917	,	208.2644628	,	210.7438017	,	212.7272727	,	214.7107438	,	216.6942149	,	219.1735537	,	221.1570248	,	223.1404959	,	225.6198347	,	228.0991736	,	230.0826446	,	232.0661157	,	234.0495868	,	235.5371901	,	237.5206612	,	240	,	241.4876033	,	243.4710744	,	245.4545455	,	247.4380165	,	249.9173554	,	251.9008264	,	254.3801653	,	256.3636364	,	257.8512397	,	260.3305785	,	262.3140496	,	265.2892562	,	268.2644628	,	270.2479339	,	272.7272727	,	275.2066116	,	277.6859504	,	280.1652893	,	282.6446281	,	285.1239669	,	288.0991736	,	291.5702479	,	294.0495868	,	297.0247934	,	300.4958678	,	303.9669421	,	306.9421488	,	310.4132231	,	313.8842975	,	317.8512397	,	321.8181818	,	326.2809917	,	329.7520661	,	333.7190083	,	338.1818182	,	342.1487603	,	345.6198347	,	349.5867769	,	353.0578512	,	356.0330579	,	359.5041322 ];
yp9 = [ 16.62952741	,	16.59190926	,	16.55429112	,	16.50914934	,	16.47153119	,	16.41886578	,	16.37372401	,	16.32858223	,	16.26839319	,	16.22325142	,	16.16306238	,	16.1179206	,	16.05773157	,	15.9900189	,	15.9373535	,	15.87716446	,	15.81697543	,	15.75678639	,	15.69659735	,	15.63640832	,	15.56869565	,	15.50098299	,	15.44079395	,	15.38812854	,	15.32793951	,	15.26775047	,	15.21508507	,	15.1473724	,	15.07213611	,	15.00442344	,	14.92918715	,	14.86147448	,	14.79376181	,	14.71852552	,	14.65081285	,	14.59062382	,	14.52291115	,	14.46272212	,	14.39500945	,	14.31977316	,	14.25206049	,	14.18434783	,	14.10911153	,	14.04139887	,	13.9736862	,	13.90597353	,	13.8457845	,	13.78559546	,	13.72540643	,	13.65769376	,	13.5899811	,	13.52226843	,	13.45455577	,	13.3868431	,	13.31913043	,	13.2589414	,	13.20627599	,	13.13856333	,	13.08589792	,	13.03323251	,	12.97304348	,	12.9279017	,	12.88275992	,	12.85266541	,	12.81504726	,	12.80752363	,	12.8	,	12.81504726	,	12.84514178	,	12.87523629	,	12.91285444	,	12.95799622	,	13.01818526	,	13.07085066	,	13.1310397	,	13.19875236	,	13.26646503	,	13.33417769	,	13.40189036	,	13.46960302	,	13.53731569	,	13.60502836	,	13.67274102	,	13.74797732	,	13.82321361	,	13.89092628	,	13.95863894	,	14.01882798	,	14.08654064	,	14.15425331	,	14.22196597	,	14.29720227	,	14.37243856	,	14.44767486	,	14.50786389	,	14.56805293	,	14.6357656	,	14.70347826	,	14.77871456	,	14.84642722	,	14.91413989	,	14.97432892	,	15.04204159	,	15.11727788	,	15.19251418	,	15.26022684	,	15.33546314	,	15.41069943	,	15.48593573	,	15.56117202	,	15.63640832	,	15.69659735	,	15.77183365	,	15.84706994	,	15.91478261	,	15.99754253	,	16.05773157	,	16.12544423	,	16.1931569	,	16.26086957	,	16.3210586	,	16.38877127	,	16.4489603	,	16.50914934	,	16.55429112	,	16.59943289	,	16.64457467	,	16.68971645	,	16.71981096	,	16.74990548	,	16.77247637	,	16.78	,	16.78	,	16.77247637	,	16.74990548	,	16.72733459	,	16.70476371	,	16.68219282	,	16.63705104 ];

figure();
plot(thetad,Erms,'LineWidth',4.5);
ylabel('Campo elétrico (kV/cm)')
xlabel('Ângulo (°)')
%title('Campo elétrico superficial condutor um') %título (trocar ao mudar condutor avaliado)
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

% plot(xp6,yy6,'--','LineWidth',2.5); %plot curva referência (trocar variáveis ao mudar condutor avaliado)
plot(xp3,yp3,'--','LineWidth',2.5); %plot curva referência (trocar variáveis ao mudar condutor avaliado)
xlim([0 360])
legend('Simulado', '(SANTOS, 2017)')

figure(); 
plot(x,y,'o','color','blue')
legend('Cabos condutores fase A, B e C');
xlabel('Faixa de passagem (m)')
ylabel('Altura do cabo (m)')
%title('Configuração geométrica dos condutores')
grid
xlim([-15 15])
ylim([16.5 16.85])

%%
% PAGANOTTI, 2012
% %vetores para plot da curva referência condutor 1
% xp1 = [ 0	,	1.806775408	,	3.613550816	,	5.872020075	,	8.130489335	,	10.38895859	,	13.09912171	,	15.80928482	,	18.97114178	,	21.68130489	,	23.93977415	,	26.64993726	,	28.90840652	,	32.07026349	,	35.6838143	,	39.29736512	,	42.45922208	,	46.0727729	,	48.78293601	,	51.94479297	,	54.20326223	,	56.91342535	,	59.1718946	,	60.97867001	,	62.78544542	,	65.04391468	,	67.30238394	,	69.10915935	,	71.36762861	,	73.17440402	,	74.98117942	,	76.78795483	,	78.59473024	,	80.40150565	,	82.20828105	,	84.01505646	,	85.82183187	,	87.62860728	,	89.43538269	,	91.24215809	,	93.0489335	,	94.85570891	,	96.21079046	,	97.56587202	,	99.37264743	,	101.1794228	,	102.5345044	,	104.3412798	,	106.1480552	,	107.9548306	,	109.3099122	,	110.6649937	,	112.4717691	,	113.8268507	,	115.6336261	,	117.4404015	,	119.2471769	,	120.6022585	,	122.4090339	,	124.2158093	,	126.0225847	,	128.281054	,	130.0878294	,	132.3462986	,	134.153074	,	135.9598494	,	138.2183187	,	140.0250941	,	142.2835634	,	144.0903388	,	145.8971142	,	148.1555834	,	150.4140527	,	152.672522	,	154.9309912	,	157.1894605	,	159.4479297	,	162.1580928	,	164.4165621	,	166.6750314	,	169.3851945	,	172.0953576	,	173.902133	,	176.6122961	,	179.3224592	,	182.4843162	,	185.6461731	,	188.3563363	,	191.5181932	,	194.6800502	,	197.8419072	,	201.455458	,	204.1656211	,	207.7791719	,	211.3927227	,	214.5545797	,	218.1681305	,	222.2333752	,	225.3952321	,	228.5570891	,	230.8155583	,	233.9774153	,	236.2358846	,	238.9460477	,	240.7528231	,	243.0112923	,	245.2697616	,	247.9799247	,	249.7867001	,	252.0451694	,	254.3036386	,	256.1104141	,	257.9171895	,	259.7239649	,	261.5307403	,	263.3375157	,	265.1442911	,	266.9510665	,	268.7578419	,	270.1129235	,	271.9196989	,	274.1781681	,	275.9849435	,	277.3400251	,	278.6951066	,	279.5984944	,	281.4052698	,	283.2120452	,	285.0188206	,	286.825596	,	288.6323714	,	289.9874529	,	291.3425345	,	293.1493099	,	294.9560853	,	296.3111669	,	298.1179423	,	299.9247177	,	301.2797992	,	302.6348808	,	304.4416562	,	306.2484316	,	308.055207	,	309.8619824	,	311.217064	,	313.4755332	,	315.2823087	,	317.0890841	,	319.3475533	,	321.1543287	,	322.9611041	,	324.7678795	,	326.1229611	,	327.9297365	,	329.7365119	,	331.0915935	,	332.8983689	,	334.7051443	,	336.5119197	,	338.770389	,	340.1254705	,	341.9322459	,	344.1907152	,	345.5457967	,	347.804266	,	349.6110414	,	351.8695107	,	353.2245922	,	355.0313676	,	357.2898369	,	360 ];
% yp1 = [ 13.71311475	,	13.68852459	,	13.6557377	,	13.6147541	,	13.59016393	,	13.55737705	,	13.52459016	,	13.50819672	,	13.49180328	,	13.47540984	,	13.46721311	,	13.45901639	,	13.46721311	,	13.48360656	,	13.50819672	,	13.54098361	,	13.57377049	,	13.6147541	,	13.66393443	,	13.70491803	,	13.75409836	,	13.78688525	,	13.83606557	,	13.87704918	,	13.93442623	,	13.99180328	,	14.04098361	,	14.09016393	,	14.13934426	,	14.18852459	,	14.23770492	,	14.28688525	,	14.3442623	,	14.39344262	,	14.45081967	,	14.5	,	14.56557377	,	14.62295082	,	14.68032787	,	14.73770492	,	14.79508197	,	14.86065574	,	14.90983607	,	14.96721311	,	15.02459016	,	15.09016393	,	15.14754098	,	15.21311475	,	15.2704918	,	15.3442623	,	15.40163934	,	15.45901639	,	15.51639344	,	15.57377049	,	15.62295082	,	15.68852459	,	15.75409836	,	15.81147541	,	15.86885246	,	15.93442623	,	16	,	16.06557377	,	16.13114754	,	16.19672131	,	16.26229508	,	16.32786885	,	16.3852459	,	16.45081967	,	16.50819672	,	16.57377049	,	16.63114754	,	16.68032787	,	16.74590164	,	16.80327869	,	16.86065574	,	16.91803279	,	16.97540984	,	17.03278689	,	17.09016393	,	17.13934426	,	17.18032787	,	17.2295082	,	17.2704918	,	17.31147541	,	17.36065574	,	17.39344262	,	17.42622951	,	17.45901639	,	17.48360656	,	17.50819672	,	17.53278689	,	17.54098361	,	17.54098361	,	17.54098361	,	17.52459016	,	17.50819672	,	17.49180328	,	17.45901639	,	17.43442623	,	17.39344262	,	17.35245902	,	17.31147541	,	17.26229508	,	17.22131148	,	17.18032787	,	17.13114754	,	17.09016393	,	17.04918033	,	17.00819672	,	16.95901639	,	16.90163934	,	16.85245902	,	16.80327869	,	16.75409836	,	16.70491803	,	16.64754098	,	16.60655738	,	16.55737705	,	16.5	,	16.45081967	,	16.40163934	,	16.3442623	,	16.28688525	,	16.24590164	,	16.19672131	,	16.14754098	,	16.09016393	,	16.02459016	,	15.96721311	,	15.91803279	,	15.86885246	,	15.81147541	,	15.75409836	,	15.69672131	,	15.63934426	,	15.59016393	,	15.53278689	,	15.47540984	,	15.42622951	,	15.37704918	,	15.31967213	,	15.24590164	,	15.18032787	,	15.12295082	,	15.07377049	,	14.99180328	,	14.93442623	,	14.86885246	,	14.81147541	,	14.75409836	,	14.68852459	,	14.62295082	,	14.57377049	,	14.51639344	,	14.46721311	,	14.42622951	,	14.36065574	,	14.30327869	,	14.25409836	,	14.19672131	,	14.14754098	,	14.09836066	,	14.04918033	,	14	,	13.95081967	,	13.90983607	,	13.86885246	,	13.83606557	,	13.80327869	,	13.76229508	,	13.7295082 ];
% 
% %vetores para plot da curva referência condutor 2
% xp2 = [ 0.456273764	,	2.281368821	,	4.562737643	,	6.844106464	,	8.669201521	,	10.95057034	,	13.23193916	,	15.51330798	,	17.79467681	,	20.07604563	,	22.35741445	,	24.63878327	,	26.92015209	,	29.20152091	,	31.48288973	,	33.76425856	,	36.04562738	,	38.3269962	,	40.60836502	,	43.3460076	,	46.08365019	,	48.82129278	,	52.01520913	,	55.20912548	,	58.40304183	,	61.14068441	,	64.33460076	,	67.98479087	,	71.17870722	,	74.37262357	,	78.02281369	,	81.21673004	,	85.32319392	,	88.51711027	,	91.71102662	,	94.90494297	,	98.09885932	,	101.2927757	,	104.9429658	,	107.6806084	,	110.8745247	,	113.6121673	,	116.3498099	,	119.0874525	,	121.8250951	,	124.5627376	,	127.3003802	,	129.581749	,	131.8631179	,	134.1444867	,	136.8821293	,	139.1634981	,	141.4448669	,	143.7262357	,	146.0076046	,	148.2889734	,	150.5703422	,	152.3954373	,	154.2205323	,	156.5019011	,	158.78327	,	160.608365	,	162.8897338	,	165.1711027	,	166.9961977	,	169.2775665	,	171.1026616	,	172.9277567	,	174.7528517	,	177.0342205	,	178.8593156	,	180.6844106	,	182.5095057	,	184.7908745	,	186.6159696	,	188.4410646	,	190.7224335	,	193.0038023	,	195.2851711	,	197.5665399	,	199.8479087	,	202.1292776	,	204.4106464	,	207.148289	,	208.973384	,	211.7110266	,	213.9923954	,	216.2737643	,	218.5551331	,	221.2927757	,	224.0304183	,	226.3117871	,	229.0494297	,	231.3307985	,	234.0684411	,	236.3498099	,	239.0874525	,	242.2813688	,	244.5627376	,	247.3003802	,	250.0380228	,	253.2319392	,	256.4258555	,	260.5323194	,	264.1825095	,	267.3764259	,	271.026616	,	274.6768061	,	278.3269962	,	281.9771863	,	285.6273764	,	288.365019	,	291.5589354	,	294.7528517	,	297.4904943	,	300.2281369	,	303.4220532	,	307.0722433	,	309.3536122	,	312.0912548	,	314.8288973	,	317.5665399	,	319.8479087	,	322.1292776	,	324.8669202	,	327.6045627	,	329.4296578	,	332.1673004	,	334.4486692	,	336.730038	,	339.0114068	,	341.2927757	,	343.5741445	,	346.3117871	,	348.5931559	,	350.8745247	,	353.1558935	,	354.9809886	,	357.2623574	,	359.0874525 ];
% yp2 = [ 14.00496689	,	14.07450331	,	14.15562914	,	14.22516556	,	14.29470199	,	14.37582781	,	14.44536424	,	14.52649007	,	14.60761589	,	14.67715232	,	14.75827815	,	14.82781457	,	14.88576159	,	14.95529801	,	15.02483444	,	15.09437086	,	15.15231788	,	15.2218543	,	15.29139073	,	15.34933775	,	15.41887417	,	15.4884106	,	15.54635762	,	15.61589404	,	15.67384106	,	15.72019868	,	15.7781457	,	15.82450331	,	15.87086093	,	15.90562914	,	15.94039735	,	15.96357616	,	15.97516556	,	15.97516556	,	15.97516556	,	15.97516556	,	15.96357616	,	15.92880795	,	15.89403974	,	15.87086093	,	15.83609272	,	15.7897351	,	15.74337748	,	15.69701987	,	15.65066225	,	15.60430464	,	15.55794702	,	15.5115894	,	15.46523179	,	15.39569536	,	15.33774834	,	15.27980132	,	15.2218543	,	15.15231788	,	15.09437086	,	15.02483444	,	14.95529801	,	14.89735099	,	14.83940397	,	14.76986755	,	14.70033113	,	14.6307947	,	14.56125828	,	14.48013245	,	14.41059603	,	14.3410596	,	14.27152318	,	14.21357616	,	14.15562914	,	14.07450331	,	14.01655629	,	13.94701987	,	13.88907285	,	13.81953642	,	13.75	,	13.68046358	,	13.61092715	,	13.54139073	,	13.4718543	,	13.37913907	,	13.29801325	,	13.21688742	,	13.14735099	,	13.07781457	,	13.00827815	,	12.93874172	,	12.8692053	,	12.81125828	,	12.74172185	,	12.68377483	,	12.61423841	,	12.55629139	,	12.49834437	,	12.44039735	,	12.38245033	,	12.34768212	,	12.2897351	,	12.24337748	,	12.19701987	,	12.17384106	,	12.13907285	,	12.11589404	,	12.08112583	,	12.05794702	,	12.03476821	,	12.02317881	,	12.02317881	,	12.03476821	,	12.05794702	,	12.08112583	,	12.11589404	,	12.13907285	,	12.17384106	,	12.23178808	,	12.2781457	,	12.32450331	,	12.38245033	,	12.44039735	,	12.50993377	,	12.56788079	,	12.62582781	,	12.69536424	,	12.76490066	,	12.82284768	,	12.89238411	,	12.96192053	,	13.03145695	,	13.10099338	,	13.18211921	,	13.25165563	,	13.32119205	,	13.40231788	,	13.48344371	,	13.56456954	,	13.63410596	,	13.70364238	,	13.78476821	,	13.84271523	,	13.91225166	,	13.98178808 ];
% 
% %vetores para plot da curva referência condutor 3
% xp3 = [ 0	,	2.241594022	,	4.03486924	,	6.276463263	,	8.518057285	,	10.3113325	,	12.55292653	,	14.79452055	,	17.03611457	,	19.27770859	,	21.51930262	,	23.76089664	,	25.55417186	,	27.79576588	,	29.5890411	,	31.38231631	,	33.17559153	,	34.96886675	,	36.76214197	,	38.55541719	,	39.9003736	,	41.69364882	,	43.48692403	,	45.28019925	,	47.07347447	,	48.41843088	,	49.7633873	,	51.55666252	,	53.34993773	,	55.14321295	,	56.93648817	,	58.72976339	,	60.0747198	,	61.41967621	,	62.76463263	,	64.10958904	,	65.45454545	,	67.24782067	,	69.04109589	,	70.83437111	,	72.62764633	,	73.97260274	,	75.31755915	,	76.66251557	,	78.45579078	,	79.8007472	,	81.14570361	,	82.49066002	,	83.83561644	,	85.18057285	,	86.52552927	,	87.87048568	,	89.21544209	,	91.00871731	,	92.80199253	,	94.14694894	,	95.94022416	,	97.73349938	,	99.5267746	,	100.871731	,	102.6650062	,	104.4582814	,	106.6998755	,	108.4931507	,	110.2864259	,	111.6313823	,	113.4246575	,	115.6662516	,	117.9078456	,	119.7011208	,	121.494396	,	123.2876712	,	125.0809465	,	127.7708593	,	130.4607721	,	133.1506849	,	136.2889166	,	139.4271482	,	143.0136986	,	146.6002491	,	149.7384807	,	152.4283935	,	155.1183064	,	159.1531756	,	162.739726	,	165.8779577	,	169.4645081	,	172.6027397	,	175.2926526	,	177.9825654	,	180.6724782	,	183.362391	,	186.0523039	,	187.8455791	,	190.5354919	,	192.7770859	,	194.5703611	,	196.8119552	,	199.0535492	,	201.2951432	,	203.5367372	,	205.7783313	,	208.0199253	,	210.2615193	,	212.0547945	,	213.8480697	,	216.0896638	,	217.882939	,	220.124533	,	221.9178082	,	224.1594022	,	225.5043587	,	227.2976339	,	229.0909091	,	230.8841843	,	232.6774595	,	234.4707347	,	236.26401	,	238.0572852	,	239.8505604	,	241.6438356	,	243.4371108	,	245.2303861	,	247.0236613	,	248.8169365	,	250.6102117	,	252.4034869	,	254.1967621	,	255.9900374	,	257.7833126	,	259.5765878	,	261.369863	,	263.611457	,	265.4047323	,	267.1980075	,	268.9912827	,	270.7845579	,	272.5778331	,	274.8194271	,	277.0610212	,	278.8542964	,	281.0958904	,	282.8891656	,	285.1307597	,	287.3723537	,	290.0622665	,	292.7521793	,	294.9937733	,	297.6836862	,	299.9252802	,	302.615193	,	305.7534247	,	307.9950187	,	310.6849315	,	313.3748443	,	316.513076	,	319.6513076	,	322.7895392	,	326.8244085	,	330.4109589	,	334.4458281	,	338.4806974	,	342.0672478	,	345.6537983	,	349.2403487	,	352.3785803	,	355.516812	,	358.2067248	,	360 ];
% yp3 = [ 17.8495935	,	17.82520325	,	17.79268293	,	17.7601626	,	17.7195122	,	17.67886179	,	17.63821138	,	17.58943089	,	17.53252033	,	17.48373984	,	17.42682927	,	17.37804878	,	17.32926829	,	17.2804878	,	17.23170732	,	17.17479675	,	17.11788618	,	17.06097561	,	17.00406504	,	16.95528455	,	16.89837398	,	16.84146341	,	16.79268293	,	16.73577236	,	16.68699187	,	16.6300813	,	16.58130081	,	16.52439024	,	16.46747967	,	16.40243902	,	16.33739837	,	16.27235772	,	16.21544715	,	16.15853659	,	16.10162602	,	16.05284553	,	15.99593496	,	15.93902439	,	15.88211382	,	15.82520325	,	15.76829268	,	15.71138211	,	15.65447154	,	15.59756098	,	15.54065041	,	15.48373984	,	15.43495935	,	15.38617886	,	15.32926829	,	15.2804878	,	15.23170732	,	15.18292683	,	15.12601626	,	15.08536585	,	15.0203252	,	14.97154472	,	14.90650407	,	14.84146341	,	14.79268293	,	14.72764228	,	14.67073171	,	14.62195122	,	14.56504065	,	14.51626016	,	14.46747967	,	14.42682927	,	14.3699187	,	14.32113821	,	14.27235772	,	14.21544715	,	14.17479675	,	14.13414634	,	14.09349593	,	14.04471545	,	13.99593496	,	13.94715447	,	13.91463415	,	13.88211382	,	13.84146341	,	13.81707317	,	13.80081301	,	13.78455285	,	13.78455285	,	13.79268293	,	13.81707317	,	13.8495935	,	13.88211382	,	13.91463415	,	13.96341463	,	14.00406504	,	14.04471545	,	14.09349593	,	14.14227642	,	14.19105691	,	14.24796748	,	14.30487805	,	14.36178862	,	14.41869919	,	14.48373984	,	14.54065041	,	14.59756098	,	14.66260163	,	14.72764228	,	14.79268293	,	14.85772358	,	14.92276423	,	14.98780488	,	15.05284553	,	15.11788618	,	15.19105691	,	15.25609756	,	15.31300813	,	15.37804878	,	15.44308943	,	15.50813008	,	15.57317073	,	15.63821138	,	15.70325203	,	15.77642276	,	15.84146341	,	15.90650407	,	15.97154472	,	16.02845528	,	16.09349593	,	16.16666667	,	16.23170732	,	16.28861789	,	16.35365854	,	16.41869919	,	16.48373984	,	16.54878049	,	16.61382114	,	16.67886179	,	16.74390244	,	16.80894309	,	16.86585366	,	16.93089431	,	16.98780488	,	17.04471545	,	17.1097561	,	17.16666667	,	17.22357724	,	17.2804878	,	17.34552846	,	17.41056911	,	17.45934959	,	17.50813008	,	17.57317073	,	17.6300813	,	17.67886179	,	17.72764228	,	17.77642276	,	17.81707317	,	17.85772358	,	17.8902439	,	17.93089431	,	17.95528455	,	17.98780488	,	18.01219512	,	18.03658537	,	18.03658537	,	18.03658537	,	18.02845528	,	18.01219512	,	17.9796748	,	17.94715447	,	17.91463415	,	17.88211382	,	17.86585366 ];
% 
% %vetores para plot da curva referência condutor 4
% xp4 = [ 0	,	1.8	,	4.05	,	6.75	,	9.45	,	12.15	,	14.85	,	17.55	,	20.7	,	23.4	,	27	,	30.15	,	32.85	,	36.9	,	40.5	,	44.55	,	48.15	,	51.75	,	54.45	,	57.6	,	61.2	,	63.9	,	67.05	,	69.75	,	72.9	,	75.15	,	78.3	,	81	,	84.15	,	87.3	,	90.45	,	93.15	,	95.85	,	98.55	,	101.25	,	103.95	,	106.65	,	109.35	,	111.6	,	114.3	,	117	,	119.7	,	121.95	,	124.65	,	127.35	,	129.6	,	132.3	,	135	,	137.7	,	139.95	,	142.65	,	145.35	,	148.05	,	150.3	,	153	,	155.7	,	158.4	,	161.55	,	164.25	,	166.95	,	169.65	,	172.8	,	175.95	,	179.1	,	182.25	,	185.85	,	188.55	,	191.7	,	194.85	,	198	,	201.15	,	204.75	,	208.35	,	211.95	,	215.1	,	218.7	,	222.3	,	225.45	,	228.15	,	232.2	,	235.35	,	238.95	,	242.55	,	246.15	,	249.3	,	252.9	,	255.6	,	258.75	,	261.9	,	264.15	,	267.3	,	270.45	,	273.15	,	275.85	,	278.1	,	280.8	,	283.05	,	285.75	,	288.45	,	290.7	,	292.95	,	295.2	,	297.9	,	300.15	,	302.85	,	305.55	,	308.25	,	310.5	,	312.75	,	315.45	,	318.15	,	320.4	,	323.1	,	325.35	,	328.05	,	330.75	,	333	,	335.7	,	338.4	,	340.65	,	342.9	,	345.15	,	347.85	,	350.55	,	352.35	,	355.05	,	357.75	,	360 ];
% yp4 = [ 15.98039216	,	15.94607843	,	15.91176471	,	15.87745098	,	15.83169935	,	15.79738562	,	15.7745098	,	15.74019608	,	15.71732026	,	15.69444444	,	15.67156863	,	15.66013072	,	15.64869281	,	15.64869281	,	15.66013072	,	15.67156863	,	15.68300654	,	15.70588235	,	15.72875817	,	15.7630719	,	15.79738562	,	15.82026144	,	15.85457516	,	15.9003268	,	15.93464052	,	15.98039216	,	16.02614379	,	16.07189542	,	16.11764706	,	16.1748366	,	16.23202614	,	16.28921569	,	16.34640523	,	16.39215686	,	16.44934641	,	16.51797386	,	16.5751634	,	16.63235294	,	16.68954248	,	16.74673203	,	16.80392157	,	16.86111111	,	16.91830065	,	16.9754902	,	17.03267974	,	17.08986928	,	17.15849673	,	17.21568627	,	17.28431373	,	17.34150327	,	17.39869281	,	17.45588235	,	17.5130719	,	17.57026144	,	17.62745098	,	17.69607843	,	17.75326797	,	17.82189542	,	17.87908497	,	17.93627451	,	18.00490196	,	18.0620915	,	18.10784314	,	18.17647059	,	18.22222222	,	18.27941176	,	18.31372549	,	18.35947712	,	18.39379085	,	18.42810458	,	18.47385621	,	18.49673203	,	18.51960784	,	18.54248366	,	18.56535948	,	18.57679739	,	18.58823529	,	18.58823529	,	18.57679739	,	18.56535948	,	18.55392157	,	18.53104575	,	18.50816993	,	18.47385621	,	18.42810458	,	18.38235294	,	18.33660131	,	18.29084967	,	18.25653595	,	18.21078431	,	18.15359477	,	18.10784314	,	18.05065359	,	18.00490196	,	17.94771242	,	17.89052288	,	17.83333333	,	17.77614379	,	17.71895425	,	17.6503268	,	17.59313725	,	17.53594771	,	17.46732026	,	17.39869281	,	17.33006536	,	17.26143791	,	17.19281046	,	17.13562092	,	17.06699346	,	16.99836601	,	16.91830065	,	16.86111111	,	16.78104575	,	16.72385621	,	16.66666667	,	16.59803922	,	16.52941176	,	16.47222222	,	16.40359477	,	16.35784314	,	16.3120915	,	16.26633987	,	16.20915033	,	16.16339869	,	16.11764706	,	16.06045752	,	16.01470588	,	15.99183007 ];
% 
% %vetores para plot da curva referência condutor 5
% xp5 = [ 0	,	1.806775408	,	3.613550816	,	5.420326223	,	7.227101631	,	9.033877039	,	10.84065245	,	12.64742785	,	14.45420326	,	16.26097867	,	18.06775408	,	20.32622334	,	22.13299875	,	24.39146801	,	26.64993726	,	28.45671267	,	30.71518193	,	32.97365119	,	35.6838143	,	38.39397742	,	40.65244668	,	43.36260979	,	46.0727729	,	49.23462986	,	51.94479297	,	54.65495609	,	57.81681305	,	60.97867001	,	64.14052698	,	67.30238394	,	70.91593476	,	74.07779172	,	76.78795483	,	79.94981179	,	82.20828105	,	84.91844417	,	87.62860728	,	90.33877039	,	93.0489335	,	95.30740276	,	97.56587202	,	99.82434128	,	102.0828105	,	104.3412798	,	106.1480552	,	107.9548306	,	109.761606	,	112.0200753	,	114.2785445	,	116.0853199	,	117.8920954	,	119.6988708	,	121.0539523	,	122.4090339	,	124.6675031	,	126.0225847	,	127.8293601	,	129.1844417	,	130.5395232	,	131.8946048	,	133.2496863	,	135.0564617	,	136.4115433	,	138.2183187	,	140.0250941	,	141.3801757	,	143.1869511	,	144.9937265	,	146.8005019	,	148.1555834	,	149.510665	,	151.3174404	,	152.672522	,	154.0276035	,	155.3826851	,	156.7377666	,	158.544542	,	159.8996236	,	161.706399	,	163.5131744	,	164.868256	,	166.6750314	,	168.0301129	,	169.3851945	,	171.1919699	,	172.5470514	,	173.902133	,	175.2572146	,	177.06399	,	178.4190715	,	180.2258469	,	181.5809285	,	183.3877039	,	184.7427854	,	186.5495609	,	188.8080301	,	190.6148055	,	192.4215809	,	194.2283563	,	196.4868256	,	197.8419072	,	199.6486826	,	201.455458	,	203.7139272	,	206.4240903	,	208.6825596	,	210.489335	,	212.7478043	,	215.0062735	,	217.7164366	,	220.4265997	,	223.1367629	,	225.846926	,	228.5570891	,	231.2672522	,	234.4291092	,	237.5909661	,	241.2045169	,	245.2697616	,	248.8833124	,	252.4968632	,	255.6587202	,	259.272271	,	262.434128	,	265.1442911	,	268.3061481	,	270.5646173	,	273.2747804	,	275.9849435	,	278.2434128	,	280.5018821	,	282.7603513	,	285.0188206	,	287.2772898	,	289.0840652	,	291.3425345	,	293.6010038	,	294.9560853	,	296.7628607	,	299.02133	,	300.8281054	,	302.6348808	,	304.4416562	,	305.7967378	,	307.6035132	,	309.4102886	,	311.217064	,	312.5721455	,	314.378921	,	315.7340025	,	317.5407779	,	319.3475533	,	320.7026349	,	322.5094103	,	322.5094103	,	322.5094103	,	324.3161857	,	326.1229611	,	328.3814304	,	329.7365119	,	331.5432873	,	332.8983689	,	334.7051443	,	336.5119197	,	338.3186951	,	340.1254705	,	341.9322459	,	343.7390213	,	345.5457967	,	347.804266	,	349.6110414	,	351.4178168	,	352.7728984	,	354.5796738	,	356.3864492	,	358.1932246	,	360 ];
% yp5 = [ 16.2295082	,	16.28852459	,	16.34754098	,	16.41639344	,	16.4852459	,	16.5442623	,	16.60327869	,	16.66229508	,	16.72131148	,	16.7704918	,	16.81967213	,	16.87868852	,	16.93770492	,	16.98688525	,	17.04590164	,	17.09508197	,	17.1442623	,	17.19344262	,	17.24262295	,	17.29180328	,	17.33114754	,	17.3704918	,	17.40983607	,	17.44918033	,	17.47868852	,	17.50819672	,	17.52786885	,	17.54754098	,	17.55737705	,	17.55737705	,	17.54754098	,	17.52786885	,	17.50819672	,	17.48852459	,	17.45901639	,	17.4295082	,	17.39016393	,	17.35081967	,	17.31147541	,	17.27213115	,	17.21311475	,	17.16393443	,	17.1147541	,	17.0557377	,	17.00655738	,	16.95737705	,	16.90819672	,	16.83934426	,	16.79016393	,	16.72131148	,	16.67213115	,	16.62295082	,	16.56393443	,	16.5147541	,	16.44590164	,	16.38688525	,	16.32786885	,	16.26885246	,	16.20983607	,	16.16065574	,	16.09180328	,	16.03278689	,	15.97377049	,	15.9147541	,	15.84590164	,	15.77704918	,	15.69836066	,	15.6295082	,	15.5704918	,	15.50163934	,	15.45245902	,	15.38360656	,	15.3147541	,	15.2557377	,	15.18688525	,	15.11803279	,	15.05901639	,	14.99016393	,	14.92131148	,	14.85245902	,	14.78360656	,	14.7147541	,	14.64590164	,	14.57704918	,	14.50819672	,	14.43934426	,	14.3704918	,	14.32131148	,	14.25245902	,	14.20327869	,	14.1442623	,	14.0852459	,	14.00655738	,	13.93770492	,	13.86885246	,	13.8	,	13.73114754	,	13.66229508	,	13.60327869	,	13.55409836	,	13.49508197	,	13.43606557	,	13.37704918	,	13.32786885	,	13.26885246	,	13.20983607	,	13.15081967	,	13.10163934	,	13.05245902	,	13.00327869	,	12.95409836	,	12.9147541	,	12.8852459	,	12.84590164	,	12.82622951	,	12.81639344	,	12.79672131	,	12.78688525	,	12.79672131	,	12.81639344	,	12.84590164	,	12.8852459	,	12.93442623	,	12.97377049	,	13.02295082	,	13.07213115	,	13.13114754	,	13.18032787	,	13.23934426	,	13.29836066	,	13.34754098	,	13.41639344	,	13.47540984	,	13.53442623	,	13.60327869	,	13.66229508	,	13.72131148	,	13.78032787	,	13.84918033	,	13.90819672	,	13.96721311	,	14.02622951	,	14.0852459	,	14.15409836	,	14.21311475	,	14.28196721	,	14.34098361	,	14.4	,	14.45901639	,	14.52786885	,	14.59672131	,	14.66557377	,	14.72459016	,	14.79344262	,	14.79344262	,	14.79344262	,	14.87213115	,	14.94098361	,	15.00983607	,	15.06885246	,	15.13770492	,	15.19672131	,	15.27540984	,	15.3442623	,	15.41311475	,	15.47213115	,	15.54098361	,	15.61967213	,	15.69836066	,	15.77704918	,	15.84590164	,	15.9147541	,	15.97377049	,	16.03278689	,	16.09180328	,	16.15081967	,	16.20983607 ];
% 
% %vetores para plot da curva referência condutor 6
% xp6 = [ 0	,	2.683229814	,	4.919254658	,	7.155279503	,	9.391304348	,	12.07453416	,	14.31055901	,	16.99378882	,	19.67701863	,	22.80745342	,	25.04347826	,	27.72670807	,	29.51552795	,	31.7515528	,	33.98757764	,	36.22360248	,	38.45962733	,	40.2484472	,	42.48447205	,	44.72049689	,	46.95652174	,	48.74534161	,	50.53416149	,	52.32298137	,	54.11180124	,	55.90062112	,	57.68944099	,	59.47826087	,	61.71428571	,	63.95031056	,	65.73913043	,	67.52795031	,	69.31677019	,	71.10559006	,	72.89440994	,	74.68322981	,	76.91925466	,	78.70807453	,	80.49689441	,	81.83850932	,	83.18012422	,	84.9689441	,	86.75776398	,	88.54658385	,	90.33540373	,	92.1242236	,	94.36024845	,	96.14906832	,	97.9378882	,	99.72670807	,	101.515528	,	103.7515528	,	105.9875776	,	107.7763975	,	110.0124224	,	112.6956522	,	114.931677	,	117.6149068	,	119.8509317	,	122.5341615	,	125.2173913	,	127.9006211	,	130.5838509	,	133.2670807	,	135.9503106	,	139.5279503	,	143.1055901	,	146.2360248	,	149.3664596	,	152.9440994	,	156.0745342	,	159.2049689	,	162.3354037	,	165.9130435	,	168.5962733	,	171.7267081	,	174.8571429	,	177.5403727	,	180.6708075	,	183.3540373	,	186.0372671	,	188.7204969	,	191.4037267	,	193.6397516	,	196.3229814	,	199.0062112	,	201.689441	,	203.9254658	,	206.1614907	,	208.8447205	,	211.0807453	,	213.3167702	,	215.552795	,	218.2360248	,	220.9192547	,	223.1552795	,	225.8385093	,	228.0745342	,	230.310559	,	232.5465839	,	234.3354037	,	236.5714286	,	238.8074534	,	240.5962733	,	242.8322981	,	245.068323	,	247.3043478	,	249.5403727	,	251.3291925	,	253.5652174	,	255.8012422	,	258.0372671	,	260.2732919	,	262.5093168	,	264.7453416	,	266.9813665	,	269.2173913	,	271.4534161	,	273.689441	,	275.9254658	,	278.1614907	,	280.8447205	,	283.9751553	,	286.6583851	,	289.7888199	,	292.9192547	,	296.0496894	,	299.6273292	,	302.757764	,	305.8881988	,	309.4658385	,	313.0434783	,	316.621118	,	319.3043478	,	322.4347826	,	325.5652174	,	329.1428571	,	332.2732919	,	336.2981366	,	339.8757764	,	343.4534161	,	346.5838509	,	350.1614907	,	353.7391304	,	356.8695652	,	360 ];
% yp6 = [ 19.0097561	,	18.97560976	,	18.9300813	,	18.89593496	,	18.83902439	,	18.79349593	,	18.74796748	,	18.6796748	,	18.62276423	,	18.55447154	,	18.48617886	,	18.41788618	,	18.36097561	,	18.29268293	,	18.23577236	,	18.15609756	,	18.08780488	,	18.0195122	,	17.95121951	,	17.87154472	,	17.81463415	,	17.73495935	,	17.67804878	,	17.6097561	,	17.54146341	,	17.47317073	,	17.40487805	,	17.33658537	,	17.26829268	,	17.18861789	,	17.1203252	,	17.04065041	,	16.97235772	,	16.90406504	,	16.83577236	,	16.76747967	,	16.71056911	,	16.64227642	,	16.57398374	,	16.51707317	,	16.44878049	,	16.3804878	,	16.31219512	,	16.25528455	,	16.18699187	,	16.1300813	,	16.06178862	,	15.99349593	,	15.92520325	,	15.85691057	,	15.78861789	,	15.7203252	,	15.65203252	,	15.59512195	,	15.52682927	,	15.4699187	,	15.41300813	,	15.34471545	,	15.29918699	,	15.23089431	,	15.16260163	,	15.11707317	,	15.0601626	,	15.02601626	,	14.99186992	,	14.95772358	,	14.93495935	,	14.90081301	,	14.87804878	,	14.86666667	,	14.87804878	,	14.88943089	,	14.90081301	,	14.93495935	,	14.96910569	,	15.00325203	,	15.04878049	,	15.10569106	,	15.16260163	,	15.2195122	,	15.26504065	,	15.33333333	,	15.40162602	,	15.45853659	,	15.52682927	,	15.60650407	,	15.67479675	,	15.74308943	,	15.82276423	,	15.90243902	,	15.99349593	,	16.07317073	,	16.15284553	,	16.23252033	,	16.32357724	,	16.40325203	,	16.49430894	,	16.57398374	,	16.66504065	,	16.75609756	,	16.83577236	,	16.90406504	,	16.99512195	,	17.06341463	,	17.13170732	,	17.22276423	,	17.29105691	,	17.37073171	,	17.4504065	,	17.5300813	,	17.6097561	,	17.68943089	,	17.76910569	,	17.83739837	,	17.90569106	,	17.98536585	,	18.05365854	,	18.12195122	,	18.20162602	,	18.2699187	,	18.3495935	,	18.42926829	,	18.50894309	,	18.57723577	,	18.65691057	,	18.72520325	,	18.79349593	,	18.86178862	,	18.91869919	,	18.97560976	,	19.02113821	,	19.06666667	,	19.10081301	,	19.13495935	,	19.15772358	,	19.1804878	,	19.21463415	,	19.22601626	,	19.22601626	,	19.21463415	,	19.1804878	,	19.14634146	,	19.12357724	,	19.07804878	,	19.04390244	,	19.0097561 ];
% 
% %vetores para plot da curva referência condutor 7
% xp7 = [ 0	,	2.241594022	,	4.03486924	,	5.828144458	,	8.069738481	,	10.3113325	,	12.55292653	,	15.24283935	,	17.48443337	,	20.62266501	,	23.76089664	,	27.34744707	,	30.93399751	,	34.07222914	,	37.21046077	,	40.79701121	,	43.93524284	,	47.52179328	,	51.10834371	,	54.24657534	,	57.38480697	,	60.0747198	,	63.21295143	,	65.90286426	,	68.59277709	,	70.83437111	,	73.07596513	,	75.76587796	,	78.00747198	,	80.249066	,	82.04234122	,	83.83561644	,	86.07721046	,	88.31880448	,	90.56039851	,	92.80199253	,	95.04358655	,	96.83686177	,	99.07845579	,	101.3200498	,	103.113325	,	104.9066002	,	106.6998755	,	108.4931507	,	110.2864259	,	112.5280199	,	114.3212951	,	116.5628892	,	118.3561644	,	120.1494396	,	121.9427148	,	123.73599	,	125.5292653	,	127.3225405	,	129.5641345	,	131.3574097	,	133.1506849	,	134.9439601	,	136.7372354	,	138.5305106	,	140.7721046	,	142.5653798	,	144.358655	,	146.6002491	,	148.3935243	,	150.6351183	,	152.4283935	,	154.2216687	,	156.014944	,	158.256538	,	160.498132	,	162.739726	,	164.98132	,	167.2229141	,	169.4645081	,	171.7061021	,	173.4993773	,	175.7409714	,	177.9825654	,	179.7758406	,	182.0174346	,	184.7073474	,	186.9489415	,	188.7422167	,	190.9838107	,	193.2254047	,	195.9153176	,	198.6052304	,	201.743462	,	204.4333748	,	207.5716065	,	211.1581569	,	214.2963885	,	217.882939	,	221.4694894	,	224.607721	,	228.1942715	,	231.7808219	,	235.3673724	,	238.0572852	,	241.1955168	,	243.4371108	,	245.2303861	,	248.3686177	,	250.1618929	,	251.9551681	,	254.1967621	,	256.4383562	,	258.6799502	,	260.9215442	,	262.7148194	,	264.5080946	,	266.3013699	,	268.0946451	,	269.8879203	,	271.6811955	,	273.4744707	,	275.7160648	,	277.0610212	,	278.8542964	,	280.6475716	,	282.4408468	,	284.234122	,	286.0273973	,	287.8206725	,	289.1656289	,	290.5105853	,	292.3038605	,	293.6488169	,	295.4420922	,	296.7870486	,	298.5803238	,	299.9252802	,	301.7185554	,	303.5118306	,	304.856787	,	306.2017435	,	307.5466999	,	309.3399751	,	310.6849315	,	312.4782067	,	314.2714819	,	316.0647572	,	317.4097136	,	319.2029888	,	320.5479452	,	321.8929016	,	323.6861768	,	325.4794521	,	326.8244085	,	328.6176837	,	330.4109589	,	332.2042341	,	333.9975093	,	335.3424658	,	337.135741	,	338.9290162	,	340.7222914	,	342.5155666	,	344.3088418	,	346.1021171	,	347.8953923	,	349.6886675	,	351.9302615	,	353.7235367	,	355.9651308	,	358.2067248	,	360 ];
% yp7 = [ 14.87928222	,	14.85073409	,	14.8278956	,	14.79934747	,	14.76508972	,	14.74225122	,	14.7137031	,	14.6908646	,	14.66231648	,	14.63947798	,	14.6223491	,	14.60522023	,	14.59380098	,	14.58809135	,	14.58238173	,	14.58809135	,	14.5995106	,	14.61663948	,	14.63947798	,	14.65660685	,	14.67944535	,	14.70228385	,	14.7365416	,	14.77079935	,	14.8050571	,	14.83360522	,	14.86786297	,	14.89641109	,	14.93066884	,	14.97063622	,	14.99918434	,	15.03344209	,	15.07340946	,	15.11337684	,	15.15334421	,	15.19331158	,	15.23327896	,	15.27324633	,	15.31892333	,	15.3588907	,	15.39885808	,	15.43882545	,	15.47879282	,	15.51305057	,	15.55872757	,	15.60440457	,	15.64437194	,	15.68433931	,	15.72430669	,	15.76427406	,	15.79853181	,	15.84420881	,	15.88417618	,	15.92414356	,	15.96982055	,	16.01549755	,	16.06117455	,	16.10114192	,	16.14681892	,	16.19249592	,	16.2324633	,	16.27243067	,	16.31810767	,	16.35807504	,	16.39804241	,	16.44371941	,	16.48368679	,	16.52365416	,	16.56362153	,	16.60929853	,	16.65497553	,	16.70065253	,	16.7406199	,	16.78058728	,	16.82626427	,	16.86623165	,	16.90619902	,	16.94045677	,	16.97471452	,	17.00897227	,	17.04323002	,	17.07748777	,	17.10603589	,	17.14029364	,	17.16884176	,	17.19738989	,	17.22593801	,	17.24877651	,	17.27732463	,	17.30016313	,	17.32300163	,	17.34013051	,	17.35154976	,	17.362969	,	17.36867863	,	17.36867863	,	17.362969	,	17.35154976	,	17.33442088	,	17.31729201	,	17.29445351	,	17.27732463	,	17.25448613	,	17.23164763	,	17.20880914	,	17.18597064	,	17.16313214	,	17.13458401	,	17.10032626	,	17.06606852	,	17.03752039	,	17.00326264	,	16.96900489	,	16.93474715	,	16.9004894	,	16.86623165	,	16.82626427	,	16.79200653	,	16.75774878	,	16.72349103	,	16.6949429	,	16.66068515	,	16.61500816	,	16.57504078	,	16.54078303	,	16.50652529	,	16.46655791	,	16.42659054	,	16.38662316	,	16.35236542	,	16.31239804	,	16.27243067	,	16.23817292	,	16.19820555	,	16.15823817	,	16.12398042	,	16.08401305	,	16.04404568	,	16.0040783	,	15.95840131	,	15.91843393	,	15.87846656	,	15.83278956	,	15.79282219	,	15.74714519	,	15.70717781	,	15.66721044	,	15.62724307	,	15.59298532	,	15.55872757	,	15.5187602	,	15.47879282	,	15.43882545	,	15.39885808	,	15.3588907	,	15.3132137	,	15.27324633	,	15.23898858	,	15.19902121	,	15.16476346	,	15.12479608	,	15.09053834	,	15.05628059	,	15.02202284	,	14.98776509	,	14.94779772	,	14.91924959	,	14.89070147 ];
% 
% %vetores para plot da curva referência condutor 8
% xp8 = [ 0	,	1.793275218	,	4.03486924	,	6.276463263	,	8.96637609	,	11.65628892	,	13.89788294	,	16.13947696	,	18.82938979	,	21.96762142	,	24.65753425	,	27.79576588	,	30.4856787	,	33.17559153	,	35.86550436	,	39.00373599	,	42.59028643	,	46.17683686	,	49.31506849	,	52.90161893	,	56.48816936	,	60.0747198	,	63.66127024	,	67.24782067	,	70.3860523	,	73.97260274	,	78.00747198	,	81.59402242	,	85.18057285	,	88.76712329	,	91.90535492	,	94.59526775	,	97.73349938	,	100.871731	,	104.0099626	,	107.1481943	,	110.2864259	,	112.9763387	,	115.2179328	,	117.9078456	,	120.5977584	,	122.8393524	,	125.0809465	,	127.3225405	,	129.1158157	,	131.3574097	,	134.0473225	,	136.2889166	,	138.9788294	,	141.2204234	,	143.4620174	,	145.7036115	,	147.9452055	,	150.1867995	,	152.4283935	,	154.6699875	,	156.4632628	,	158.256538	,	160.498132	,	162.739726	,	164.98132	,	166.7745953	,	169.0161893	,	171.2577833	,	173.0510585	,	174.8443337	,	177.0859278	,	179.3275218	,	181.5691158	,	183.8107098	,	185.6039851	,	188.2938979	,	190.0871731	,	192.3287671	,	194.5703611	,	197.260274	,	199.501868	,	201.743462	,	204.4333748	,	207.1232877	,	210.2615193	,	213.3997509	,	216.0896638	,	219.2278954	,	222.366127	,	226.4009963	,	229.9875467	,	233.5740971	,	237.1606476	,	240.747198	,	244.7820672	,	247.9202989	,	251.9551681	,	255.0933998	,	257.7833126	,	261.8181818	,	264.5080946	,	268.0946451	,	270.7845579	,	273.4744707	,	277.0610212	,	279.750934	,	282.8891656	,	285.5790785	,	288.7173101	,	291.4072229	,	294.0971357	,	296.3387298	,	298.5803238	,	301.2702366	,	303.9601494	,	306.2017435	,	308.4433375	,	310.6849315	,	313.3748443	,	315.6164384	,	318.3063512	,	320.996264	,	323.6861768	,	325.9277709	,	328.6176837	,	330.8592777	,	333.1008717	,	335.3424658	,	337.5840598	,	339.8256538	,	342.5155666	,	345.2054795	,	347.4470735	,	350.1369863	,	352.3785803	,	354.6201743	,	357.3100872	,	360 ];
% yp8 = [ 14.91544715	,	14.9804878	,	15.04552846	,	15.12357724	,	15.20162602	,	15.2796748	,	15.34471545	,	15.42276423	,	15.48780488	,	15.5398374	,	15.61788618	,	15.68292683	,	15.74796748	,	15.81300813	,	15.87804878	,	15.94308943	,	15.99512195	,	16.03414634	,	16.0601626	,	16.09918699	,	16.12520325	,	16.13821138	,	16.15121951	,	16.15121951	,	16.15121951	,	16.13821138	,	16.11219512	,	16.09918699	,	16.07317073	,	16.02113821	,	15.96910569	,	15.9300813	,	15.86504065	,	15.8	,	15.73495935	,	15.65691057	,	15.59186992	,	15.51382114	,	15.43577236	,	15.37073171	,	15.30569106	,	15.22764228	,	15.16260163	,	15.08455285	,	15.00650407	,	14.92845528	,	14.8504065	,	14.77235772	,	14.68130081	,	14.60325203	,	14.49918699	,	14.42113821	,	14.34308943	,	14.25203252	,	14.17398374	,	14.08292683	,	13.99186992	,	13.92682927	,	13.83577236	,	13.75772358	,	13.66666667	,	13.58861789	,	13.49756098	,	13.4195122	,	13.32845528	,	13.26341463	,	13.19837398	,	13.10731707	,	13.01626016	,	12.93821138	,	12.87317073	,	12.79512195	,	12.7300813	,	12.65203252	,	12.57398374	,	12.49593496	,	12.41788618	,	12.35284553	,	12.28780488	,	12.22276423	,	12.14471545	,	12.09268293	,	12.04065041	,	11.97560976	,	11.92357724	,	11.88455285	,	11.85853659	,	11.8195122	,	11.80650407	,	11.80650407	,	11.80650407	,	11.8195122	,	11.83252033	,	11.85853659	,	11.88455285	,	11.93658537	,	11.98861789	,	12.05365854	,	12.10569106	,	12.15772358	,	12.22276423	,	12.28780488	,	12.36585366	,	12.43089431	,	12.49593496	,	12.57398374	,	12.63902439	,	12.71707317	,	12.80813008	,	12.87317073	,	12.93821138	,	13.02926829	,	13.10731707	,	13.18536585	,	13.27642276	,	13.36747967	,	13.45853659	,	13.5495935	,	13.64065041	,	13.73170732	,	13.8097561	,	13.88780488	,	13.97886179	,	14.05691057	,	14.14796748	,	14.23902439	,	14.31707317	,	14.40813008	,	14.49918699	,	14.57723577	,	14.65528455	,	14.73333333	,	14.81138211	,	14.88943089 ];
% 
% %vetores para plot da curva referência condutor 9
% xp9 = [ 0	,	1.860465116	,	3.720930233	,	5.581395349	,	7.441860465	,	9.302325581	,	11.1627907	,	13.02325581	,	14.41860465	,	15.81395349	,	17.6744186	,	19.53488372	,	21.39534884	,	23.25581395	,	24.65116279	,	26.51162791	,	28.37209302	,	29.76744186	,	31.62790698	,	33.02325581	,	34.41860465	,	35.81395349	,	37.20930233	,	38.60465116	,	40	,	41.39534884	,	42.79069767	,	44.18604651	,	45.58139535	,	46.97674419	,	48.37209302	,	49.76744186	,	51.1627907	,	52.55813953	,	53.48837209	,	54.88372093	,	56.27906977	,	57.6744186	,	59.06976744	,	60.46511628	,	61.86046512	,	63.25581395	,	64.18604651	,	65.58139535	,	66.97674419	,	67.90697674	,	69.30232558	,	70.69767442	,	72.09302326	,	73.48837209	,	74.88372093	,	76.27906977	,	77.6744186	,	79.06976744	,	80.46511628	,	81.39534884	,	82.79069767	,	84.18604651	,	85.11627907	,	86.51162791	,	87.90697674	,	89.30232558	,	90.69767442	,	92.09302326	,	93.48837209	,	94.88372093	,	96.27906977	,	98.13953488	,	99.53488372	,	101.3953488	,	103.255814	,	105.1162791	,	106.5116279	,	108.372093	,	109.7674419	,	111.627907	,	113.0232558	,	114.8837209	,	116.744186	,	118.1395349	,	120.4651163	,	122.3255814	,	124.6511628	,	126.5116279	,	128.8372093	,	131.1627907	,	133.9534884	,	136.744186	,	140	,	143.255814	,	146.0465116	,	149.3023256	,	152.5581395	,	155.8139535	,	159.0697674	,	162.3255814	,	165.5813953	,	168.8372093	,	171.627907	,	174.8837209	,	177.6744186	,	180	,	182.3255814	,	184.6511628	,	186.9767442	,	188.8372093	,	190.6976744	,	192.5581395	,	194.4186047	,	196.2790698	,	198.6046512	,	200.9302326	,	202.7906977	,	204.6511628	,	206.5116279	,	208.372093	,	210.2325581	,	212.0930233	,	213.9534884	,	215.8139535	,	217.6744186	,	219.0697674	,	220.9302326	,	222.3255814	,	224.1860465	,	225.5813953	,	226.9767442	,	228.372093	,	230.2325581	,	232.0930233	,	233.4883721	,	234.8837209	,	236.744186	,	238.1395349	,	240	,	241.3953488	,	242.7906977	,	244.1860465	,	246.0465116	,	247.4418605	,	249.3023256	,	251.1627907	,	252.5581395	,	254.4186047	,	255.8139535	,	257.6744186	,	259.0697674	,	260.4651163	,	262.3255814	,	263.7209302	,	265.5813953	,	267.4418605	,	269.3023256	,	271.1627907	,	273.0232558	,	274.8837209	,	276.744186	,	279.0697674	,	280.9302326	,	283.255814	,	285.5813953	,	287.9069767	,	290.2325581	,	292.5581395	,	294.8837209	,	297.6744186	,	300	,	302.3255814	,	305.1162791	,	307.9069767	,	311.1627907	,	313.9534884	,	317.2093023	,	320.4651163	,	323.7209302	,	326.9767442	,	330.2325581	,	333.9534884	,	338.6046512	,	342.3255814	,	346.0465116	,	349.3023256	,	352.0930233	,	355.8139535	,	358.1395349	,	360 ];
% yp9 = [ 17.32236842	,	17.28947368	,	17.26315789	,	17.23026316	,	17.19736842	,	17.16447368	,	17.13157895	,	17.09868421	,	17.06578947	,	17.03289474	,	16.99342105	,	16.95394737	,	16.91447368	,	16.86842105	,	16.82894737	,	16.78947368	,	16.75	,	16.70394737	,	16.65789474	,	16.61184211	,	16.56578947	,	16.52631579	,	16.48026316	,	16.43421053	,	16.39473684	,	16.34868421	,	16.30921053	,	16.26315789	,	16.21710526	,	16.17105263	,	16.13157895	,	16.07894737	,	16.03947368	,	16	,	15.96710526	,	15.92105263	,	15.875	,	15.82236842	,	15.77631579	,	15.73026316	,	15.68421053	,	15.63815789	,	15.59868421	,	15.55921053	,	15.51973684	,	15.47368421	,	15.42763158	,	15.38157895	,	15.33552632	,	15.28289474	,	15.23684211	,	15.18421053	,	15.13815789	,	15.08552632	,	15.03947368	,	14.99342105	,	14.94078947	,	14.90131579	,	14.86184211	,	14.81578947	,	14.76973684	,	14.72368421	,	14.67763158	,	14.63815789	,	14.59210526	,	14.53947368	,	14.49342105	,	14.44078947	,	14.39473684	,	14.34868421	,	14.30263158	,	14.25657895	,	14.21052632	,	14.16447368	,	14.11842105	,	14.07236842	,	14.02631579	,	13.97368421	,	13.93421053	,	13.89473684	,	13.85526316	,	13.81578947	,	13.76973684	,	13.73684211	,	13.69736842	,	13.65789474	,	13.61842105	,	13.59210526	,	13.55921053	,	13.53289474	,	13.51315789	,	13.5	,	13.5	,	13.5	,	13.50657895	,	13.52631579	,	13.55263158	,	13.58552632	,	13.61842105	,	13.65131579	,	13.69736842	,	13.73684211	,	13.77631579	,	13.81578947	,	13.86184211	,	13.90789474	,	13.95394737	,	14	,	14.04605263	,	14.09868421	,	14.15789474	,	14.21052632	,	14.26315789	,	14.32236842	,	14.375	,	14.42763158	,	14.48684211	,	14.54605263	,	14.59868421	,	14.65789474	,	14.71710526	,	14.76973684	,	14.82894737	,	14.875	,	14.93421053	,	14.98684211	,	15.03947368	,	15.08552632	,	15.13815789	,	15.19736842	,	15.25	,	15.30263158	,	15.35526316	,	15.40789474	,	15.46052632	,	15.51973684	,	15.57236842	,	15.61842105	,	15.67105263	,	15.72368421	,	15.78289474	,	15.84210526	,	15.90131579	,	15.94736842	,	16	,	16.05263158	,	16.10526316	,	16.15789474	,	16.21710526	,	16.26973684	,	16.31578947	,	16.38157895	,	16.44078947	,	16.5	,	16.54605263	,	16.60526316	,	16.65789474	,	16.71052632	,	16.76973684	,	16.82894737	,	16.88815789	,	16.94078947	,	16.98684211	,	17.03947368	,	17.08552632	,	17.13157895	,	17.17763158	,	17.22368421	,	17.26973684	,	17.31578947	,	17.34868421	,	17.38815789	,	17.42105263	,	17.44078947	,	17.47368421	,	17.48684211	,	17.5	,	17.50657895	,	17.5	,	17.48684211	,	17.46710526	,	17.44078947	,	17.41447368	,	17.38815789	,	17.35526316	,	17.32894737 ];

err = abs((yp3 - Erms)/yp3)*100;