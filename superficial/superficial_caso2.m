clear all; close all; clc

% Cálculo do campo elétrico superficial do caso 2 do sistema trifásico

e_0 = 8.854*(10^(-12));
r = 14.37*(10^-3); % raio do condutor
n = 6; % número de condutores do sistema
nc = n*2; % número de condutores total (com imagens)
ci = (2*n-1); % número de cargas imagens


%% matriz P

% posição dos condutores reais do sistema
xr = [-10.478 -10.25 -0.228 0 10.022 10.25];
yr = [14.290 14.290 14.290 14.290 14.290 14.290];
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
y = [14.290 14.290 14.290 14.290 14.290 14.290 -14.290 -14.290 -14.290 -14.290 -14.290 -14.290]; % posição y dos condutores (com imagens)

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

theta = linspace(0,2*pi,152); %gera a superfície do condutor em 360 pontos

xcj = x(3); % posição x do centro condutor 1 fase a (para trocar o condutor avaliado deve-se mudar este valor)
ycj = y(3); % posição y do centro condutor 1 fase a

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

%vetores para plot da curva referência condutor 1
xp1 = [ 0	,	2.723833544	,	4.99369483	,	7.717528373	,	10.44136192	,	13.16519546	,	14.98108449	,	17.70491803	,	19.97477932	,	22.24464061	,	24.96847415	,	27.23833544	,	29.50819672	,	31.77805801	,	33.59394704	,	35.40983607	,	36.77175284	,	38.58764187	,	40.4035309	,	41.76544767	,	43.12736444	,	44.94325347	,	46.7591425	,	48.57503153	,	49.9369483	,	51.29886507	,	53.1147541	,	54.47667087	,	55.83858764	,	57.20050441	,	59.01639344	,	60.83228247	,	62.19419924	,	63.55611602	,	64.91803279	,	66.27994956	,	67.64186633	,	69.0037831	,	70.36569987	,	71.72761665	,	73.08953342	,	74.45145019	,	75.81336696	,	77.62925599	,	78.99117276	,	80.35308953	,	81.71500631	,	83.53089533	,	85.34678436	,	86.70870113	,	88.52459016	,	89.88650694	,	91.70239596	,	93.06431274	,	94.42622951	,	95.78814628	,	98.05800757	,	99.41992434	,	101.2358134	,	103.0517024	,	104.8675914	,	106.6834805	,	108.4993695	,	110.3152585	,	112.1311475	,	113.9470366	,	116.2168979	,	118.0327869	,	120.3026482	,	122.1185372	,	124.3883985	,	126.2042875	,	128.4741488	,	130.7440101	,	133.4678436	,	135.7377049	,	137.5535939	,	140.2774275	,	142.5472888	,	145.2711223	,	148.4489281	,	151.6267339	,	154.3505675	,	157.5283733	,	161.1601513	,	164.7919294	,	168.4237074	,	172.0554855	,	175.6872636	,	179.3190416	,	182.9508197	,	186.1286255	,	188.852459	,	192.0302648	,	195.6620429	,	198.8398487	,	202.0176545	,	205.1954603	,	207.9192938	,	210.1891551	,	213.3669609	,	216.0907945	,	219.2686003	,	221.9924338	,	224.7162673	,	227.4401009	,	230.1639344	,	231.9798235	,	234.703657	,	237.4274905	,	239.2433796	,	241.0592686	,	242.8751576	,	244.6910467	,	247.4148802	,	249.2307692	,	251.0466583	,	253.3165195	,	254.6784363	,	257.4022699	,	259.2181589	,	260.5800757	,	262.8499369	,	264.665826	,	266.481715	,	267.3896595	,	269.2055485	,	271.0214376	,	272.3833544	,	274.1992434	,	276.0151324	,	277.8310214	,	279.6469105	,	281.0088272	,	282.8247163	,	284.186633	,	285.5485498	,	287.3644388	,	288.7263556	,	290.5422446	,	292.3581337	,	293.7200504	,	295.9899117	,	297.3518285	,	299.1677175	,	300.9836066	,	303.2534678	,	304.6153846	,	306.4312736	,	308.2471627	,	310.0630517	,	311.8789407	,	313.6948298	,	315.964691	,	318.2345523	,	320.0504414	,	322.7742749	,	324.5901639	,	326.8600252	,	329.1298865	,	331.8537201	,	334.1235813	,	335.9394704	,	338.6633039	,	341.3871375	,	344.5649433	,	348.1967213	,	351.8284994	,	355.4602774	,	359.0920555 ];
yp1 = [ 11.67874794	,	11.684514	,	11.69604613	,	11.70757825	,	11.72487644	,	11.74217463	,	11.76523888	,	11.78830313	,	11.81713344	,	11.85172982	,	11.88632619	,	11.92668863	,	11.96128501	,	11.99588138	,	12.03047776	,	12.06507414	,	12.09967051	,	12.13426689	,	12.16886326	,	12.20345964	,	12.24382208	,	12.28995058	,	12.32454695	,	12.36490939	,	12.41103789	,	12.45140033	,	12.49176277	,	12.53212521	,	12.57248764	,	12.60708402	,	12.65321252	,	12.69934102	,	12.73970346	,	12.7800659	,	12.82042834	,	12.86079077	,	12.90115321	,	12.94728171	,	12.98187809	,	13.01647446	,	13.0568369	,	13.1029654	,	13.14332784	,	13.18369028	,	13.22981878	,	13.27594728	,	13.32207578	,	13.36820428	,	13.40856672	,	13.45469522	,	13.50082372	,	13.54695222	,	13.59308072	,	13.63344316	,	13.67957166	,	13.72570016	,	13.77182867	,	13.8121911	,	13.8583196	,	13.91021417	,	13.95057661	,	13.99670511	,	14.04283361	,	14.09472817	,	14.14085667	,	14.18121911	,	14.23311367	,	14.27924217	,	14.31960461	,	14.36573311	,	14.41186161	,	14.44645799	,	14.49258649	,	14.53871499	,	14.57331137	,	14.61943987	,	14.65980231	,	14.70016474	,	14.73476112	,	14.76359143	,	14.79818781	,	14.83278418	,	14.86738056	,	14.90197694	,	14.93080725	,	14.95963756	,	14.97693575	,	14.98846787	,	14.99423394	,	15	,	15	,	15	,	14.99423394	,	14.97116969	,	14.9538715	,	14.93080725	,	14.907743	,	14.87891269	,	14.85584843	,	14.82701812	,	14.79242175	,	14.75782537	,	14.71746293	,	14.67133443	,	14.63097199	,	14.58484349	,	14.54448105	,	14.50411862	,	14.45799012	,	14.41762768	,	14.37149918	,	14.33113674	,	14.28500824	,	14.2446458	,	14.1985173	,	14.15815486	,	14.1062603	,	14.05436573	,	14.00247117	,	13.95057661	,	13.89868204	,	13.85255354	,	13.80065898	,	13.74876442	,	13.70263591	,	13.65074135	,	13.59884679	,	13.55848435	,	13.51235585	,	13.45469522	,	13.40280066	,	13.35667216	,	13.31054366	,	13.26441516	,	13.21828666	,	13.17215815	,	13.12026359	,	13.07413509	,	13.03377265	,	12.98764415	,	12.94151565	,	12.88385502	,	12.8261944	,	12.7800659	,	12.7339374	,	12.6878089	,	12.6416804	,	12.58401977	,	12.52635914	,	12.48023064	,	12.43410214	,	12.38797364	,	12.33607908	,	12.27841845	,	12.23228995	,	12.18616145	,	12.13426689	,	12.09390445	,	12.04200988	,	11.99588138	,	11.95551895	,	11.92092257	,	11.88056013	,	11.84596376	,	11.79983526	,	11.76523888	,	11.74217463	,	11.71911038	,	11.69604613	,	11.69028007 ];

%vetores para plot da curva referência condutor 2
xp2 = [ 0.459183673	,	4.132653061	,	7.806122449	,	11.02040816	,	14.69387755	,	17.44897959	,	20.20408163	,	23.41836735	,	26.17346939	,	28.46938776	,	30.76530612	,	32.60204082	,	34.43877551	,	36.73469388	,	39.03061224	,	40.86734694	,	43.16326531	,	45.91836735	,	48.67346939	,	50.51020408	,	52.80612245	,	55.10204082	,	57.39795918	,	59.23469388	,	61.07142857	,	62.90816327	,	64.74489796	,	66.58163265	,	69.33673469	,	71.63265306	,	73.46938776	,	75.30612245	,	77.14285714	,	78.97959184	,	80.81632653	,	82.65306122	,	84.48979592	,	86.32653061	,	88.62244898	,	90.91836735	,	92.75510204	,	94.59183673	,	96.42857143	,	98.26530612	,	99.64285714	,	101.0204082	,	102.8571429	,	104.6938776	,	106.5306122	,	108.3673469	,	110.2040816	,	112.0408163	,	113.4183673	,	114.7959184	,	116.6326531	,	118.0102041	,	119.8469388	,	121.2244898	,	122.6020408	,	123.5204082	,	125.3571429	,	127.1938776	,	128.5714286	,	130.4081633	,	132.244898	,	134.0816327	,	135.9183673	,	138.2142857	,	140.9693878	,	142.8061224	,	145.5612245	,	148.3163265	,	150.6122449	,	152.9081633	,	155.6632653	,	157.9591837	,	160.7142857	,	163.4693878	,	167.1428571	,	170.8163265	,	174.0306122	,	177.7040816	,	180.9183673	,	183.6734694	,	186.8877551	,	189.1836735	,	191.9387755	,	194.6938776	,	197.4489796	,	200.2040816	,	202.0408163	,	204.7959184	,	207.5510204	,	209.8469388	,	212.6020408	,	215.3571429	,	217.1938776	,	219.4897959	,	222.244898	,	224.5408163	,	226.8367347	,	229.1326531	,	231.4285714	,	233.2653061	,	235.1020408	,	237.3979592	,	239.2346939	,	241.0714286	,	242.9081633	,	244.744898	,	246.122449	,	247.9591837	,	249.7959184	,	251.6326531	,	253.4693878	,	255.3061224	,	257.1428571	,	259.4387755	,	260.8163265	,	262.6530612	,	264.9489796	,	266.7857143	,	269.0816327	,	271.377551	,	273.2142857	,	275.5102041	,	277.3469388	,	279.6428571	,	281.4795918	,	283.7755102	,	285.6122449	,	287.9081633	,	289.744898	,	292.0408163	,	293.877551	,	295.7142857	,	297.5510204	,	299.3877551	,	301.6836735	,	303.9795918	,	306.2755102	,	308.5714286	,	311.3265306	,	314.0816327	,	316.8367347	,	319.5918367	,	322.3469388	,	325.1020408	,	327.8571429	,	331.5306122	,	334.744898	,	338.877551	,	343.4693878	,	347.1428571	,	350.8163265	,	354.9489796	,	358.622449 ];
yp2 = [ 15.22757475	,	15.21428571	,	15.20099668	,	15.18770764	,	15.1744186	,	15.15448505	,	15.1345515	,	15.10797342	,	15.08803987	,	15.05481728	,	15.0282392	,	15.00166113	,	14.96843854	,	14.93521595	,	14.90199336	,	14.86212625	,	14.82225914	,	14.78239203	,	14.74252492	,	14.70265781	,	14.64950166	,	14.60299003	,	14.55647841	,	14.52325581	,	14.47009967	,	14.42358804	,	14.38372093	,	14.3372093	,	14.28405316	,	14.23089701	,	14.18438538	,	14.14451827	,	14.08471761	,	14.03156146	,	13.9717608	,	13.91860465	,	13.85215947	,	13.79900332	,	13.74584718	,	13.67940199	,	13.61960133	,	13.56644518	,	13.51328904	,	13.45348837	,	13.40697674	,	13.3538206	,	13.30730897	,	13.25415282	,	13.20099668	,	13.14119601	,	13.08139535	,	13.02159468	,	12.98172757	,	12.93521595	,	12.88870432	,	12.84883721	,	12.80232558	,	12.76910299	,	12.71594684	,	12.67607973	,	12.63621262	,	12.589701	,	12.53654485	,	12.49667774	,	12.45016611	,	12.39700997	,	12.34385382	,	12.29069767	,	12.23754153	,	12.18438538	,	12.12458472	,	12.07142857	,	12.02491694	,	11.98504983	,	11.94518272	,	11.90531561	,	11.8654485	,	11.8255814	,	11.79900332	,	11.76578073	,	11.75249169	,	11.73920266	,	11.73920266	,	11.73920266	,	11.74584718	,	11.76578073	,	11.78571429	,	11.81229236	,	11.83887043	,	11.8654485	,	11.8986711	,	11.93853821	,	11.98504983	,	12.02491694	,	12.07807309	,	12.12458472	,	12.17109635	,	12.23089701	,	12.29069767	,	12.35049834	,	12.39700997	,	12.45016611	,	12.5166113	,	12.56976744	,	12.62956811	,	12.68272425	,	12.74252492	,	12.80232558	,	12.85548173	,	12.90863787	,	12.9551495	,	13.00166113	,	13.06146179	,	13.12126246	,	13.18106312	,	13.23421927	,	13.30066445	,	13.3538206	,	13.40697674	,	13.47342193	,	13.53322259	,	13.59302326	,	13.65282392	,	13.7192691	,	13.77906977	,	13.83887043	,	13.8986711	,	13.95847176	,	14.01827243	,	14.07142857	,	14.13787375	,	14.1910299	,	14.25083056	,	14.30398671	,	14.35049834	,	14.39036545	,	14.43023256	,	14.47674419	,	14.53654485	,	14.58305648	,	14.63621262	,	14.69601329	,	14.74252492	,	14.78239203	,	14.83554817	,	14.8820598	,	14.93521595	,	14.98837209	,	15.02159468	,	15.05481728	,	15.08803987	,	15.12790698	,	15.16112957	,	15.18770764	,	15.21428571	,	15.22757475	,	15.23421927 ];

%vetores para plot da curva referência condutor 3
xp3 = [ 0.406777131	,	2.775041051	,	6.412897447	,	10.0522466	,	13.69159576	,	17.3316913	,	20.51798776	,	23.24973877	,	26.43603523	,	28.71324078	,	30.99119272	,	33.26914465	,	35.54784296	,	37.36975668	,	39.64770861	,	41.47111509	,	43.74906702	,	46.02776534	,	47.8519182	,	50.13061651	,	51.95402299	,	54.2327213	,	56.05687416	,	58.33557247	,	59.7044335	,	61.52858636	,	63.35199283	,	65.17614569	,	67.00104493	,	69.28048963	,	71.55918794	,	73.38408718	,	75.66353187	,	77.48843111	,	79.31333035	,	81.13748321	,	82.96088969	,	84.78504254	,	86.6091954	,	88.4318555	,	89.80071652	,	91.624123	,	92.99373041	,	95.2731751	,	97.55261979	,	99.83206449	,	101.6562173	,	103.935662	,	105.7590685	,	107.582475	,	109.8611733	,	112.1391252	,	113.9617853	,	116.2397373	,	118.5184356	,	121.2516794	,	123.0750858	,	124.8977459	,	127.1764442	,	129.9089416	,	132.1868936	,	134.919391	,	137.6518883	,	141.2927303	,	144.0252276	,	146.7562323	,	149.9417824	,	152.672787	,	155.8575907	,	159.0423944	,	162.6817435	,	165.8650545	,	168.5945664	,	172.2316764	,	175.8695328	,	179.5058964	,	183.14226	,	186.7778773	,	189.9589491	,	193.5930736	,	197.6817435	,	201.3151217	,	204.9477534	,	207.6705478	,	210.8478877	,	214.025974	,	217.2040603	,	219.9268548	,	222.6496492	,	225.3716973	,	228.5497835	,	230.3627407	,	232.6302433	,	235.351545	,	237.619794	,	239.4334975	,	241.7010001	,	243.9685028	,	246.235259	,	248.0482162	,	249.8604269	,	252.1279295	,	254.3946858	,	256.661442	,	258.4729064	,	259.8320645	,	261.6450216	,	263.4579788	,	265.724735	,	267.5369458	,	269.803702	,	271.6159128	,	273.8834154	,	275.6956262	,	277.9608897	,	279.3178086	,	281.1300194	,	282.9422302	,	284.7536946	,	286.5666517	,	288.3788625	,	290.6456187	,	292.4570831	,	294.2692939	,	295.6291984	,	296.9868637	,	298.7990745	,	301.0658307	,	302.8780415	,	304.6909987	,	306.0501567	,	307.4085684	,	309.6760711	,	311.4875355	,	313.7542917	,	316.4755934	,	318.2892969	,	320.1022541	,	322.3690103	,	324.6372593	,	326.9040155	,	329.6260636	,	332.348858	,	335.5247052	,	338.2474996	,	340.516495	,	343.2400358	,	345.964323	,	349.1431557	,	352.3219884	,	355.5023138	,	359.1386774 ];
yp3 = [ 12.4729064	,	12.4729064	,	12.48768473	,	12.51724138	,	12.54679803	,	12.58374384	,	12.62807882	,	12.67241379	,	12.71674877	,	12.76108374	,	12.81280788	,	12.86453202	,	12.92364532	,	12.96059113	,	13.01231527	,	13.06403941	,	13.11576355	,	13.17487685	,	13.23399015	,	13.29310345	,	13.34482759	,	13.40394089	,	13.46305419	,	13.52216749	,	13.57389163	,	13.63300493	,	13.68472906	,	13.74384236	,	13.81034483	,	13.87684729	,	13.93596059	,	14.00246305	,	14.06896552	,	14.13546798	,	14.20197044	,	14.26108374	,	14.31280788	,	14.37192118	,	14.43103448	,	14.47536946	,	14.5270936	,	14.57881773	,	14.63793103	,	14.7044335	,	14.77093596	,	14.83743842	,	14.89655172	,	14.96305419	,	15.01477833	,	15.06650246	,	15.12561576	,	15.1773399	,	15.22167488	,	15.27339901	,	15.33251232	,	15.39162562	,	15.44334975	,	15.48768473	,	15.54679803	,	15.59852217	,	15.65024631	,	15.70197044	,	15.75369458	,	15.79802956	,	15.84975369	,	15.88669951	,	15.92364532	,	15.96059113	,	15.99014778	,	16.01970443	,	16.04926108	,	16.06403941	,	16.0862069	,	16.09359606	,	16.10837438	,	16.10837438	,	16.10837438	,	16.10098522	,	16.09359606	,	16.07142857	,	16.04926108	,	16.01970443	,	15.98275862	,	15.93842365	,	15.89408867	,	15.85714286	,	15.82019704	,	15.77586207	,	15.73152709	,	15.67980296	,	15.64285714	,	15.591133	,	15.53940887	,	15.48029557	,	15.43596059	,	15.39162562	,	15.33990148	,	15.28817734	,	15.22906404	,	15.1773399	,	15.1182266	,	15.06650246	,	15.00738916	,	14.94827586	,	14.8817734	,	14.83743842	,	14.78571429	,	14.73399015	,	14.67487685	,	14.61576355	,	14.55665025	,	14.49753695	,	14.44581281	,	14.38669951	,	14.31280788	,	14.24630542	,	14.18719212	,	14.12807882	,	14.06157635	,	14.00985222	,	13.95073892	,	13.89162562	,	13.82512315	,	13.76600985	,	13.72906404	,	13.66995074	,	13.61083744	,	13.55172414	,	13.49261084	,	13.4408867	,	13.39655172	,	13.34482759	,	13.29310345	,	13.22660099	,	13.16748768	,	13.10837438	,	13.06403941	,	13.01231527	,	12.95320197	,	12.908867	,	12.84975369	,	12.79802956	,	12.75369458	,	12.69458128	,	12.65024631	,	12.61330049	,	12.57635468	,	12.54679803	,	12.51724138	,	12.48768473	,	12.4729064	,	12.4729064 ];

%vetores para plot da curva referência condutor 4
xp4 = [ 0	,	3.272727273	,	6.545454545	,	9.818181818	,	13.55844156	,	16.36363636	,	19.63636364	,	22.90909091	,	26.18181818	,	28.98701299	,	31.32467532	,	34.12987013	,	36.46753247	,	38.80519481	,	41.61038961	,	43.94805195	,	46.28571429	,	48.62337662	,	50.96103896	,	53.2987013	,	55.16883117	,	57.03896104	,	58.90909091	,	60.77922078	,	62.64935065	,	64.51948052	,	66.85714286	,	69.19480519	,	70.5974026	,	72.93506494	,	74.80519481	,	76.20779221	,	77.61038961	,	79.48051948	,	81.35064935	,	83.22077922	,	84.62337662	,	86.49350649	,	88.83116883	,	90.7012987	,	92.1038961	,	93.97402597	,	95.37662338	,	97.24675325	,	99.11688312	,	100.987013	,	101.9220779	,	103.3246753	,	104.7272727	,	106.5974026	,	108.4675325	,	109.8701299	,	111.7402597	,	113.6103896	,	115.4805195	,	116.8831169	,	118.7532468	,	120.1558442	,	122.025974	,	124.3636364	,	126.2337662	,	128.1038961	,	129.974026	,	131.8441558	,	133.7142857	,	136.0519481	,	137.9220779	,	139.7922078	,	141.6623377	,	144	,	146.3376623	,	148.6753247	,	151.012987	,	153.8181818	,	157.0909091	,	160.3636364	,	163.6363636	,	166.9090909	,	170.1818182	,	174.3896104	,	177.6623377	,	180.9350649	,	184.2077922	,	187.9480519	,	191.6883117	,	194.961039	,	198.2337662	,	201.038961	,	203.8441558	,	207.1168831	,	210.3896104	,	213.1948052	,	216	,	217.8701299	,	220.2077922	,	222.5454545	,	225.3506494	,	227.2207792	,	229.0909091	,	231.4285714	,	233.2987013	,	236.1038961	,	238.4415584	,	240.3116883	,	242.6493506	,	244.5194805	,	246.3896104	,	248.2597403	,	250.5974026	,	252.4675325	,	254.3376623	,	255.7402597	,	257.6103896	,	259.4805195	,	261.3506494	,	263.2207792	,	265.0909091	,	267.4285714	,	269.7662338	,	271.6363636	,	273.974026	,	275.8441558	,	277.7142857	,	280.0519481	,	282.3896104	,	284.2597403	,	286.5974026	,	288.4675325	,	290.8051948	,	292.6753247	,	295.4805195	,	297.8181818	,	300.1558442	,	302.4935065	,	305.2987013	,	307.6363636	,	310.9090909	,	313.7142857	,	316.987013	,	320.2597403	,	323.5324675	,	326.3376623	,	330.0779221	,	333.3506494	,	336.6233766	,	340.8311688	,	344.5714286	,	348.7792208	,	352.5194805	,	355.7922078	,	360 ];
yp4 = [ 16.10472973	,	16.10472973	,	16.09712838	,	16.08192568	,	16.05912162	,	16.03631757	,	16.01351351	,	15.98310811	,	15.94510135	,	15.90709459	,	15.86908784	,	15.83108108	,	15.80067568	,	15.75506757	,	15.71706081	,	15.6714527	,	15.63344595	,	15.58783784	,	15.54983108	,	15.51182432	,	15.46621622	,	15.42060811	,	15.375	,	15.32939189	,	15.29138514	,	15.23817568	,	15.18496622	,	15.13935811	,	15.08614865	,	15.02533784	,	14.97212838	,	14.92652027	,	14.88091216	,	14.83530405	,	14.77449324	,	14.72128378	,	14.66047297	,	14.60726351	,	14.53885135	,	14.49324324	,	14.44003378	,	14.37922297	,	14.32601351	,	14.28040541	,	14.21959459	,	14.17398649	,	14.12837838	,	14.08277027	,	14.03716216	,	13.9839527	,	13.93074324	,	13.86993243	,	13.80912162	,	13.74831081	,	13.7027027	,	13.65709459	,	13.60388514	,	13.55067568	,	13.49746622	,	13.43665541	,	13.38344595	,	13.32263514	,	13.26182432	,	13.20861486	,	13.15540541	,	13.1097973	,	13.05658784	,	13.00337838	,	12.95016892	,	12.89695946	,	12.84375	,	12.79814189	,	12.74493243	,	12.70692568	,	12.66131757	,	12.60810811	,	12.5625	,	12.53209459	,	12.50168919	,	12.48648649	,	12.47128378	,	12.47128378	,	12.47888514	,	12.50168919	,	12.52449324	,	12.55489865	,	12.58530405	,	12.62331081	,	12.66891892	,	12.71452703	,	12.76773649	,	12.82094595	,	12.87415541	,	12.92736486	,	12.98817568	,	13.04898649	,	13.1097973	,	13.17060811	,	13.23141892	,	13.29222973	,	13.36064189	,	13.42905405	,	13.50506757	,	13.55827703	,	13.62668919	,	13.6875	,	13.74831081	,	13.80912162	,	13.86993243	,	13.93074324	,	13.99155405	,	14.04476351	,	14.09797297	,	14.15878378	,	14.21959459	,	14.27280405	,	14.34881757	,	14.41722973	,	14.47804054	,	14.53885135	,	14.61486486	,	14.68327703	,	14.74408784	,	14.8125	,	14.87331081	,	14.93412162	,	14.98733108	,	15.04814189	,	15.1089527	,	15.16976351	,	15.23057432	,	15.29138514	,	15.34459459	,	15.39780405	,	15.45861486	,	15.52702703	,	15.58783784	,	15.64864865	,	15.71706081	,	15.76266892	,	15.80827703	,	15.86148649	,	15.90709459	,	15.9527027	,	15.99831081	,	16.02871622	,	16.05152027	,	16.07432432	,	16.08952703	,	16.10472973	,	16.11233108 ];

%vetores para plot da curva referência condutor 5
xp5 = [ 0	,	3.146067416	,	6.292134831	,	9.438202247	,	13.03370787	,	16.17977528	,	18.42696629	,	21.12359551	,	23.82022472	,	26.51685393	,	29.21348315	,	31.46067416	,	33.25842697	,	35.05617978	,	37.30337079	,	39.5505618	,	41.34831461	,	43.14606742	,	45.39325843	,	47.64044944	,	49.43820225	,	51.23595506	,	53.48314607	,	55.28089888	,	57.07865169	,	58.87640449	,	60.2247191	,	62.02247191	,	63.82022472	,	65.61797753	,	67.41573034	,	69.21348315	,	71.01123596	,	72.80898876	,	75.05617978	,	77.30337079	,	78.65168539	,	80.4494382	,	82.24719101	,	84.49438202	,	86.29213483	,	88.08988764	,	89.88764045	,	91.68539326	,	93.48314607	,	95.28089888	,	97.52808989	,	99.3258427	,	101.1235955	,	102.9213483	,	104.7191011	,	106.9662921	,	109.2134831	,	111.4606742	,	113.258427	,	115.505618	,	117.752809	,	120.4494382	,	122.6966292	,	124.9438202	,	127.6404494	,	129.4382022	,	132.1348315	,	134.3820225	,	137.5280899	,	140.6741573	,	144.2696629	,	147.4157303	,	150.5617978	,	154.1573034	,	158.2022472	,	161.3483146	,	164.9438202	,	169.4382022	,	173.4831461	,	177.9775281	,	182.0224719	,	185.6179775	,	190.1123596	,	193.7078652	,	197.752809	,	201.3483146	,	204.9438202	,	208.0898876	,	210.7865169	,	213.4831461	,	216.1797753	,	218.8764045	,	221.5730337	,	224.7191011	,	226.9662921	,	229.2134831	,	231.9101124	,	234.1573034	,	236.8539326	,	239.1011236	,	241.3483146	,	243.5955056	,	246.2921348	,	248.5393258	,	250.3370787	,	252.1348315	,	254.3820225	,	256.6292135	,	258.4269663	,	260.6741573	,	262.4719101	,	264.2696629	,	266.5168539	,	268.3146067	,	270.1123596	,	271.9101124	,	273.258427	,	275.0561798	,	276.8539326	,	279.1011236	,	280.8988764	,	282.6966292	,	284.494382	,	285.8426966	,	287.6404494	,	288.988764	,	291.2359551	,	293.0337079	,	294.3820225	,	296.1797753	,	297.9775281	,	299.7752809	,	301.5730337	,	303.8202247	,	305.6179775	,	306.9662921	,	309.2134831	,	311.4606742	,	313.7078652	,	315.9550562	,	317.752809	,	320	,	322.247191	,	324.0449438	,	326.2921348	,	328.5393258	,	331.2359551	,	333.9325843	,	336.6292135	,	340.2247191	,	342.4719101	,	346.0674157	,	350.1123596	,	352.3595506	,	356.4044944	,	360 ];
yp5 = [ 11.73414634	,	11.74065041	,	11.74715447	,	11.77317073	,	11.79918699	,	11.82520325	,	11.85772358	,	11.8902439	,	11.92926829	,	11.97479675	,	12.0203252	,	12.06585366	,	12.10487805	,	12.1504065	,	12.19593496	,	12.24796748	,	12.29349593	,	12.33902439	,	12.39105691	,	12.44308943	,	12.49512195	,	12.54715447	,	12.59918699	,	12.65772358	,	12.7097561	,	12.76178862	,	12.81382114	,	12.86585366	,	12.91788618	,	12.9699187	,	13.03495935	,	13.09349593	,	13.14552846	,	13.20406504	,	13.26260163	,	13.32764228	,	13.3796748	,	13.43821138	,	13.49674797	,	13.54878049	,	13.60731707	,	13.66585366	,	13.71788618	,	13.7699187	,	13.82195122	,	13.87398374	,	13.92601626	,	13.98455285	,	14.0300813	,	14.07560976	,	14.14065041	,	14.19918699	,	14.26422764	,	14.3097561	,	14.36178862	,	14.4203252	,	14.47886179	,	14.53089431	,	14.58292683	,	14.63495935	,	14.6804878	,	14.72601626	,	14.77804878	,	14.82357724	,	14.86910569	,	14.92764228	,	14.97317073	,	15.02520325	,	15.06422764	,	15.10325203	,	15.13577236	,	15.16178862	,	15.18780488	,	15.21382114	,	15.22682927	,	15.2398374	,	15.2398374	,	15.23333333	,	15.2203252	,	15.20081301	,	15.16829268	,	15.14227642	,	15.10325203	,	15.06422764	,	15.02520325	,	14.98617886	,	14.95365854	,	14.91463415	,	14.86260163	,	14.81707317	,	14.77154472	,	14.72601626	,	14.6804878	,	14.63495935	,	14.57642276	,	14.52439024	,	14.47235772	,	14.42682927	,	14.37479675	,	14.32276423	,	14.26422764	,	14.21219512	,	14.1601626	,	14.10162602	,	14.04308943	,	13.98455285	,	13.93902439	,	13.88699187	,	13.82195122	,	13.76341463	,	13.71138211	,	13.65284553	,	13.60081301	,	13.54878049	,	13.49674797	,	13.44471545	,	13.39268293	,	13.33414634	,	13.27560976	,	13.21707317	,	13.17154472	,	13.12601626	,	13.07398374	,	13.01544715	,	12.9699187	,	12.90487805	,	12.85284553	,	12.78780488	,	12.72926829	,	12.67723577	,	12.62520325	,	12.57317073	,	12.52113821	,	12.46910569	,	12.40406504	,	12.35203252	,	12.30650407	,	12.24796748	,	12.20243902	,	12.1504065	,	12.09837398	,	12.04634146	,	12.00081301	,	11.95528455	,	11.91626016	,	11.87073171	,	11.83170732	,	11.79268293	,	11.7601626	,	11.74715447	,	11.74065041	,	11.73414634 ];

%vetores para plot da curva referência condutor 6
xp6 = [ 0	,	3.311432326	,	6.622864652	,	9.934296978	,	13.71879106	,	17.03022339	,	19.86859396	,	22.70696452	,	26.01839685	,	28.85676741	,	32.6412615	,	35.0065703	,	37.37187911	,	40.21024967	,	42.57555848	,	44.46780552	,	47.30617608	,	50.14454665	,	52.03679369	,	53.45597898	,	55.82128778	,	57.71353482	,	60.07884363	,	61.49802891	,	63.39027595	,	65.282523	,	66.70170828	,	68.59395532	,	70.95926413	,	72.85151117	,	74.74375821	,	76.63600526	,	78.5282523	,	80.42049934	,	82.31274639	,	83.73193167	,	85.15111695	,	87.04336399	,	88.46254928	,	89.88173456	,	91.30091984	,	92.72010512	,	94.61235217	,	96.03153745	,	97.45072273	,	98.86990802	,	100.7621551	,	102.6544021	,	104.5466491	,	105.9658344	,	107.3850197	,	108.804205	,	110.2233903	,	112.1156373	,	113.5348226	,	115.4270696	,	117.3193167	,	118.738502	,	120.1576873	,	121.5768725	,	122.9960578	,	124.4152431	,	126.3074901	,	127.7266754	,	129.6189225	,	131.0381078	,	132.9303548	,	134.8226018	,	136.2417871	,	138.6070959	,	140.499343	,	142.39159	,	144.2838371	,	146.1760841	,	148.5413929	,	150.9067017	,	152.7989488	,	155.1642576	,	157.5295664	,	160.8409987	,	164.152431	,	166.9908016	,	170.3022339	,	173.6136662	,	177.3981603	,	180.7095926	,	184.9671485	,	189.2247043	,	193.4822602	,	197.2667543	,	201.0512484	,	203.4165572	,	206.2549277	,	209.0932983	,	211.9316689	,	214.2969777	,	216.1892247	,	218.5545335	,	220.9198423	,	222.8120894	,	225.1773982	,	227.0696452	,	229.434954	,	231.3272011	,	233.2194481	,	235.1116951	,	237.0039422	,	238.8961892	,	240.3153745	,	242.2076216	,	244.0998686	,	245.9921156	,	247.8843627	,	249.7766097	,	251.195795	,	253.088042	,	254.9802891	,	256.8725361	,	258.7647832	,	260.6570302	,	262.5492773	,	264.4415243	,	266.8068331	,	268.6990802	,	270.5913272	,	272.4835742	,	274.3758213	,	275.7950066	,	277.6872536	,	280.0525624	,	281.9448095	,	284.3101183	,	286.2023653	,	288.5676741	,	290.4599212	,	292.3521682	,	294.717477	,	297.0827858	,	299.4480946	,	301.3403417	,	303.7056505	,	306.0709593	,	308.9093298	,	311.7477004	,	314.1130092	,	316.478318	,	318.8436268	,	321.6819974	,	324.5203679	,	327.3587385	,	330.6701708	,	333.9816032	,	338.239159	,	341.5505913	,	345.3350854	,	349.5926413	,	352.9040736	,	356.2155059	,	359.5269382 ];
yp6 = [ 15	,	15	,	14.99398625	,	14.98797251	,	14.96993127	,	14.95189003	,	14.9338488	,	14.90378007	,	14.87371134	,	14.83762887	,	14.80154639	,	14.77147766	,	14.72938144	,	14.68728522	,	14.645189	,	14.60910653	,	14.56701031	,	14.52491409	,	14.48281787	,	14.4467354	,	14.40463918	,	14.36254296	,	14.32044674	,	14.28436426	,	14.24226804	,	14.20618557	,	14.16408935	,	14.12199313	,	14.07388316	,	14.0257732	,	13.98367698	,	13.92955326	,	13.8814433	,	13.83333333	,	13.79123711	,	13.75515464	,	13.71305842	,	13.66494845	,	13.62285223	,	13.58676976	,	13.54467354	,	13.49054983	,	13.44845361	,	13.3943299	,	13.35223368	,	13.31013746	,	13.26202749	,	13.20790378	,	13.15979381	,	13.11168385	,	13.06357388	,	13.02147766	,	12.97938144	,	12.92525773	,	12.87714777	,	12.8290378	,	12.78092784	,	12.73883162	,	12.69072165	,	12.64261168	,	12.60652921	,	12.55841924	,	12.51632302	,	12.46821306	,	12.42611684	,	12.37800687	,	12.3419244	,	12.29381443	,	12.24570447	,	12.1975945	,	12.15549828	,	12.11340206	,	12.07731959	,	12.02920962	,	11.9871134	,	11.95103093	,	11.90893471	,	11.86683849	,	11.83676976	,	11.78865979	,	11.75257732	,	11.72250859	,	11.70446735	,	11.68041237	,	11.66838488	,	11.66237113	,	11.66838488	,	11.69243986	,	11.72250859	,	11.76460481	,	11.80670103	,	11.84278351	,	11.88487973	,	11.92697595	,	11.98109966	,	12.02319588	,	12.0652921	,	12.11340206	,	12.16151203	,	12.20360825	,	12.25773196	,	12.31185567	,	12.36597938	,	12.41408935	,	12.46219931	,	12.51030928	,	12.55841924	,	12.61254296	,	12.65463918	,	12.70876289	,	12.7628866	,	12.81701031	,	12.87714777	,	12.93127148	,	12.98539519	,	13.0395189	,	13.08762887	,	13.14175258	,	13.19587629	,	13.25601375	,	13.31013746	,	13.37027491	,	13.43642612	,	13.49656357	,	13.54467354	,	13.59278351	,	13.64089347	,	13.69501718	,	13.74312715	,	13.8032646	,	13.86340206	,	13.92353952	,	13.98367698	,	14.03178694	,	14.0919244	,	14.14003436	,	14.19415808	,	14.24226804	,	14.29037801	,	14.33848797	,	14.38659794	,	14.44072165	,	14.49484536	,	14.54295533	,	14.59707904	,	14.63316151	,	14.67525773	,	14.71735395	,	14.75945017	,	14.80154639	,	14.83762887	,	14.87972509	,	14.91580756	,	14.93986254	,	14.96993127	,	14.98797251	,	14.99398625	,	15	,	15 ];

figure();
plot(thetad,Erms,'LineWidth',4.5);
ylabel('Campo elétrico (kV/cm)')
xlabel('Ângulo (°)')
% title('Campo elétrico superficial condutor seis') %título (trocar ao mudar condutor avaliado)
grid
hold on
% yy1 = smooth(xp1,yp1,'rloess'); %comando suavizar a curva referência condutor 1
% yy2 = smooth(xp2,yp2,'rloess'); %comando suavizar a curva referência condutor 2
% yy3 = smooth(xp3,yp3,'rloess'); %comando suavizar a curva referência condutor 3
% yy4 = smooth(xp4,yp4,'rloess'); %comando suavizar a curva referência condutor 4
% yy5 = smooth(xp5,yp5,'rloess'); %comando suavizar a curva referência condutor 5
% yy6 = smooth(xp6,yp6,'rloess'); %comando suavizar a curva referência condutor 6

% plot(xp6,yy6,'--','LineWidth',2.5); %plot curva referência (trocar variáveis ao mudar condutor avaliado)
plot(xp3,yp3,'--','LineWidth',2.5); %plot curva referência (trocar variáveis ao mudar condutor avaliado)
xlim([0 360])
legend('Simulado', '(PAGANOTTI, 2012)')

figure();
plot(xr,yr,'o','color','blue')
legend('Cabos condutores fase A, B e C');
xlabel('Faixa de passagem (m)')
ylabel('Altura do cabo (m)')
grid
xlim([-15 15])
ylim([13 15.5])

err = abs((yp3 - Erms)/yp3)*100;