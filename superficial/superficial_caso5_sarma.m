clear all; close all; clc

% C�lculo do campo el�trico superficial do caso 4 do sistema trif�sico

e_0 = 8.854*(10^(-12));
r = 0.035052/2; % raio do condutor
n = 12; % n�mero de condutores do sistema
nc = n*2; % n�mero de condutores total (com imagens)
ci = (2*n-1); % n�mero de cargas imagens por condutor


%% matriz P

% posi��o dos condutores reais do sistema
xr = [-15.4686 -15.4686 -15.0114 -15.0114 -0.2286 -0.2286 0.2286 0.2286 15.0114 15.0114 15.4686 15.4686];
yr = [19.5834 20.0406 19.5834 20.0406 19.5834 20.0406 19.5834 20.0406 19.5834 20.0406 19.5834 20.0406];
yi = -yr;

% c�lculo matriz P
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

%% C�lculo tens�o por fase

%tens�o de opera��o m�xima do sistema
V = 1;

%tens�o condutores fase a
V_ra = V*(cos(2*pi/3));
V_ia = 1i*V*(sin(2*pi/3));

%tens�o condutores fase b
V_rb = V;
V_ib = 0;

%tens�o condutores fase c
V_rc = V*(cos(-2*pi/3));
V_ic = 1i*V*(sin(-2*pi/3));

Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ];

%% C�lculo densidade de carga

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


%% Dist�ncia entre condutores

x = [-15.4686 -15.4686 -15.0114 -15.0114 -0.2286 -0.2286 0.2286 0.2286 15.0114 15.0114 15.4686 15.4686 -15.4686 -15.4686 -15.0114 -15.0114 -0.2286 -0.2286 0.2286 0.2286 15.0114 15.0114 15.4686 15.4686];
y = [19.5834 20.0406 19.5834 20.0406 19.5834 20.0406 19.5834 20.0406 19.5834 20.0406 19.5834 20.0406 -19.5834 -20.0406 -19.5834 -20.0406 -19.5834 -20.0406 -19.5834 -20.0406 -19.5834 -20.0406 -19.5834 -20.0406];

%sendo: i -> o n�mero do primeiro condutor e j -> o n�mero do segundo condutor
%exemplo: 1,2 � a dist�ncia entre condutor 1 e o condutor 2, os pr�ximos
%c�lculos seguem esta mesma l�gica
%sendo 10 a imagem do condutor 1, 12 a imagem do condutor 3 e 11 a imagem do
%condutor 2

distance = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j
            distance(i,j) = 0; %como n�o existe dist�ncia no centro do pr�prio condutor ent�o � zero
        else
            distance(i,j) = sqrt(((x(j)-x(i))^2)+((y(j)-y(i))^2)); %c�lculo dist�ncia entre dois pontos (entre centro de dois condutores)
        end
    end
end

% plot(x,y,'o')
% grid

%% C�lculo de delta

delta = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j
            delta(i,j) = 0; %delta para o centro do pr�prio condutor � zero
        else
            delta(i,j) = (r^2)./(distance(i,j)); %c�lculo de delta utilizando a f�rmula e a dist�ncia calculada acima
        end
    end
end

%% C�lculo das posi��es em x das cargas imagens
 
phix = zeros(nc,nc);
posx = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j || x(j)==x(i) %c�lculo de posi��o para cargas posicionadas no centro dos condutores
            posx(i,j) = x(i);
        elseif y(j)==y(i) %c�lculo para cargas de mesma altura
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

%% C�lculo das posi��es em y das cargas imagens

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

%% Pontos de avalia��o (toda superf�cie do condutor)

theta = linspace(0,2*pi,221); %gera a superf�cie do condutor em 360 pontos

xcj = x(5); % posi��o x do centro condutor 1 fase a (para trocar o condutor avaliado deve-se mudar este valor)
ycj = y(5); % posi��o y do centro condutor 1 fase a (para trocar o condutor avaliado deve-se mudar este valor)

xf = r.*cos(theta) + xcj; % eixo x do ponto de avalia��o
yf = r.*sin(theta) + ycj; % eixo y do ponto de avalia��o

% plot(xf,yf)

%% E_xr componente x real campo el�trico condutor 1 fase a assim segue:

%busca o valor da carga imagem no condutor 1 real, valor do eixo x dos
%pontos de avalia��o e as posi��es geradas pelas matrizes posx e posy das
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

%% E_xi componente x imaginario campo el�trico condutor 2 fase b assim segue:

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

%% E_yr componente y real campo el�trico condutor 2 fase b

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

%% E_yi21 componente y imaginario campo el�trico condutor 2 fase b

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

%% aqui em cada matriz dessa tem as componentes geradas na dire��o x e y, pela parte real e imagin�ria de cada carga el�trica:

E_xr = [ E_xr12 ;	E_xr13 	;	E_xr14 	;	E_xr15 	;	E_xr16 	;	E_xr17 	;	E_xr18 	;	E_xr19 	;	E_xr110 	;	E_xr1_11 	;	E_xr1_12 	;	E_xr1_13 	;	E_xr1_14 	;	E_xr1_15 	;	E_xr1_16 	;	E_xr1_17 	;	E_xr1_18 	;	E_xr1_19 	;	E_xr1_20 	;	E_xr1_21 	;	E_xr1_22 	;	E_xr1_23 	;	E_xr1_24 	;	E_xr21 	;	E_xr23 	;	E_xr24 	;	E_xr25 	;	E_xr26 	;	E_xr27 	;	E_xr28 	;	E_xr29 	;	E_xr210 	;	E_xr2_11 	;	E_xr2_12 	;	E_xr2_13 	;	E_xr2_14 	;	E_xr2_15 	;	E_xr2_16 	;	E_xr2_17 	;	E_xr2_18 	;	E_xr2_19 	;	E_xr2_20 	;	E_xr2_21 	;	E_xr2_22 	;	E_xr2_23 	;	E_xr2_24 	;	E_xr31 	;	E_xr32 	;	E_xr34 	;	E_xr35 	;	E_xr36 	;	E_xr37 	;	E_xr38 	;	E_xr39 	;	E_xr310 	;	E_xr311 	;	E_xr312 	;	E_xr313 	;	E_xr314 	;	E_xr315 	;	E_xr316 	;	E_xr317 	;	E_xr318 	;	E_xr319 	;	E_xr320 	;	E_xr321 	;	E_xr322 	;	E_xr323 	;	E_xr324 	;	E_xr41 	;	E_xr42 	;	E_xr43 	;	E_xr45 	;	E_xr46 	;	E_xr47 	;	E_xr48 	;	E_xr49 	;	E_xr410 	;	E_xr411 	;	E_xr412 	;	E_xr413 	;	E_xr414 	;	E_xr415 	;	E_xr416 	;	E_xr417 	;	E_xr418 	;	E_xr419 	;	E_xr420 	;	E_xr421 	;	E_xr422 	;	E_xr423 	;	E_xr424 	;	E_xr51 	;	E_xr52 	;	E_xr53 	;	E_xr54 	;	E_xr56 	;	E_xr57 	;	E_xr58 	;	E_xr59 	;	E_xr510 	;	E_xr511 	;	E_xr512 	;	E_xr513 	;	E_xr514 	;	E_xr515 	;	E_xr516 	;	E_xr517 	;	E_xr518 	;	E_xr519 	;	E_xr520 	;	E_xr521 	;	E_xr522 	;	E_xr523 	;	E_xr524 	;	E_xr61 	;	E_xr62 	;	E_xr63 	;	E_xr64 	;	E_xr65 	;	E_xr67 	;	E_xr68 	;	E_xr69 	;	E_xr610 	;	E_xr611 	;	E_xr612 	;	E_xr613 	;	E_xr614 	;	E_xr615 	;	E_xr616 	;	E_xr617 	;	E_xr618 	;	E_xr619 	;	E_xr620 	;	E_xr621 	;	E_xr622 	;	E_xr623 	;	E_xr624 	;	E_xr71 	;	E_xr72 	;	E_xr73 	;	E_xr74 	;	E_xr75 	;	E_xr76 	;	E_xr78 	;	E_xr79 	;	E_xr710 	;	E_xr711 	;	E_xr712 	;	E_xr713 	;	E_xr714 	;	E_xr715 	;	E_xr716 	;	E_xr717 	;	E_xr718 	;	E_xr719 	;	E_xr720 	;	E_xr721 	;	E_xr722 	;	E_xr723 	;	E_xr724 	;	E_xr81 	;	E_xr82 	;	E_xr83 	;	E_xr84 	;	E_xr85 	;	E_xr86 	;	E_xr87 	;	E_xr89 	;	E_xr810 	;	E_xr811 	;	E_xr812 	;	E_xr813 	;	E_xr814 	;	E_xr815 	;	E_xr816 	;	E_xr817 	;	E_xr818 	;	E_xr819 	;	E_xr820 	;	E_xr821 	;	E_xr822 	;	E_xr823 	;	E_xr824 	;	E_xr91 	;	E_xr92 	;	E_xr93 	;	E_xr94 	;	E_xr95 	;	E_xr96 	;	E_xr97 	;	E_xr98 	;	E_xr910 	;	E_xr911 	;	E_xr912 	;	E_xr913 	;	E_xr914 	;	E_xr915 	;	E_xr916 	;	E_xr917 	;	E_xr918 	;	E_xr919 	;	E_xr920 	;	E_xr921 	;	E_xr922 	;	E_xr923 	;	E_xr924 	;	E_xr101 	;	E_xr102 	;	E_xr103 	;	E_xr104 	;	E_xr105 	;	E_xr106 	;	E_xr107 	;	E_xr108 	;	E_xr109 	;	E_xr1011 	;	E_xr1012 	;	E_xr1013 	;	E_xr1014 	;	E_xr1015 	;	E_xr1016 	;	E_xr1017 	;	E_xr1018 	;	E_xr10_19 	;	E_xr10_20 	;	E_xr10_21 	;	E_xr10_22 	;	E_xr10_23 	;	E_xr10_24 	;	E_xr11_1 	;	E_xr11_2 	;	E_xr113 	;	E_xr114 	;	E_xr115 	;	E_xr116 	;	E_xr117 	;	E_xr118 	;	E_xr119 	;	E_xr1110 	;	E_xr1112 	;	E_xr1113 	;	E_xr1114 	;	E_xr1115 	;	E_xr1116 	;	E_xr1117 	;	E_xr1118 	;	E_xr1019 	;	E_xr1020 	;	E_xr1021 	;	E_xr1022 	;	E_xr1023 	;	E_xr1024 	;	E_xr12_1 	;	E_xr12_2 	;	E_xr123 	;	E_xr124 	;	E_xr125 	;	E_xr126 	;	E_xr127 	;	E_xr128 	;	E_xr129 	;	E_xr1210 	;	E_xr1211 	;	E_xr1213 	;	E_xr1214 	;	E_xr1215 	;	E_xr1216 	;	E_xr1217 	;	E_xr1218 	;	E_xr1219 	;	E_xr1220 	;	E_xr1221 	;	E_xr1222 	;	E_xr1223 	;	E_xr1224 	;	E_xr131 	;	E_xr132 	;	E_xr133 	;	E_xr134 	;	E_xr135 	;	E_xr136 	;	E_xr137 	;	E_xr138 	;	E_xr139 	;	E_xr1310 	;	E_xr1311 	;	E_xr1312 	;	E_xr1314 	;	E_xr1315 	;	E_xr1316 	;	E_xr1317 	;	E_xr1318 	;	E_xr1319 	;	E_xr1320 	;	E_xr1321 	;	E_xr1322 	;	E_xr1323 	;	E_xr1324 	;	E_xr141 	;	E_xr142 	;	E_xr143 	;	E_xr144 	;	E_xr145 	;	E_xr146 	;	E_xr147 	;	E_xr148 	;	E_xr149 	;	E_xr1410 	;	E_xr1411 	;	E_xr1412 	;	E_xr1413 	;	E_xr1415 	;	E_xr1416 	;	E_xr1417 	;	E_xr1418 	;	E_xr1419 	;	E_xr1420 	;	E_xr1421 	;	E_xr1422 	;	E_xr1423 	;	E_xr1424 	;	E_xr151 	;	E_xr152 	;	E_xr153 	;	E_xr154 	;	E_xr155 	;	E_xr156 	;	E_xr157 	;	E_xr158 	;	E_xr159 	;	E_xr1510 	;	E_xr1511 	;	E_xr1512 	;	E_xr1513 	;	E_xr1514 	;	E_xr1516 	;	E_xr1517 	;	E_xr1518 	;	E_xr1519 	;	E_xr1520 	;	E_xr1521 	;	E_xr1522 	;	E_xr1523 	;	E_xr1524 	;	E_xr161 	;	E_xr162 	;	E_xr163 	;	E_xr164 	;	E_xr165 	;	E_xr166 	;	E_xr167 	;	E_xr168 	;	E_xr169 	;	E_xr1610 	;	E_xr1611 	;	E_xr1612 	;	E_xr1613 	;	E_xr1614 	;	E_xr1615 	;	E_xr1617 	;	E_xr1618 	;	E_xr1619 	;	E_xr1620 	;	E_xr1621 	;	E_xr1622 	;	E_xr1623 	;	E_xr1624 	;	E_xr171 	;	E_xr172 	;	E_xr173 	;	E_xr174 	;	E_xr175 	;	E_xr176 	;	E_xr177 	;	E_xr178 	;	E_xr179 	;	E_xr1710 	;	E_xr1711 	;	E_xr1712 	;	E_xr1713 	;	E_xr1714 	;	E_xr1715 	;	E_xr1716 	;	E_xr1718 	;	E_xr1719 	;	E_xr1720 	;	E_xr1721 	;	E_xr1722 	;	E_xr1723 	;	E_xr1724 	;	E_xr181 	;	E_xr182 	;	E_xr183 	;	E_xr184 	;	E_xr185 	;	E_xr186 	;	E_xr187 	;	E_xr188 	;	E_xr189 	;	E_xr1810 	;	E_xr1811 	;	E_xr1812 	;	E_xr1813 	;	E_xr1814 	;	E_xr1815 	;	E_xr1816 	;	E_xr1817 	;	E_xr1819 	;	E_xr1820 	;	E_xr1821 	;	E_xr1822 	;	E_xr1823 	;	E_xr1824 	;	E_xr191 	;	E_xr192 	;	E_xr193 	;	E_xr194 	;	E_xr195 	;	E_xr196 	;	E_xr197 	;	E_xr198 	;	E_xr199 	;	E_xr1910 	;	E_xr1911 	;	E_xr1912 	;	E_xr1913 	;	E_xr1914 	;	E_xr1915 	;	E_xr1916 	;	E_xr1917 	;	E_xr1918 	;	E_xr1920 	;	E_xr1921 	;	E_xr1922 	;	E_xr1923 	;	E_xr1924 	;	E_xr201 	;	E_xr202 	;	E_xr203 	;	E_xr204 	;	E_xr205 	;	E_xr206 	;	E_xr207 	;	E_xr208 	;	E_xr209 	;	E_xr2010 	;	E_xr2011 	;	E_xr2012 	;	E_xr2013 	;	E_xr2014 	;	E_xr2015 	;	E_xr2016 	;	E_xr2017 	;	E_xr2018 	;	E_xr2019 	;	E_xr2021 	;	E_xr2022 	;	E_xr2023 	;	E_xr2024 	;	E_xr21_1 	;	E_xr21_2 	;	E_xr21_3 	;	E_xr214 	;	E_xr215 	;	E_xr216 	;	E_xr217 	;	E_xr218 	;	E_xr219 	;	E_xr2110 	;	E_xr2111 	;	E_xr2112 	;	E_xr2113 	;	E_xr2114 	;	E_xr2115 	;	E_xr2116 	;	E_xr2117 	;	E_xr2118 	;	E_xr2119 	;	E_xr2120 	;	E_xr2122 	;	E_xr2123 	;	E_xr2124 	;	E_xr22_1 	;	E_xr22_2 	;	E_xr22_3 	;	E_xr22_4 	;	E_xr225 	;	E_xr226 	;	E_xr227 	;	E_xr228 	;	E_xr229 	;	E_xr2210 	;	E_xr2211 	;	E_xr2212 	;	E_xr2213 	;	E_xr2214 	;	E_xr2215 	;	E_xr2216 	;	E_xr2217 	;	E_xr2218 	;	E_xr2219 	;	E_xr2220 	;	E_xr2221 	;	E_xr2223 	;	E_xr2224 	;	E_xr231 	;	E_xr232 	;	E_xr233 	;	E_xr234 	;	E_xr235 	;	E_xr236 	;	E_xr237 	;	E_xr238 	;	E_xr239 	;	E_xr2310 	;	E_xr2311 	;	E_xr2312 	;	E_xr2313 	;	E_xr2314 	;	E_xr2315 	;	E_xr2316 	;	E_xr2317 	;	E_xr2318 	;	E_xr2319 	;	E_xr2320 	;	E_xr2321 	;	E_xr2322 	;	E_xr2324 	;	E_xr241 	;	E_xr242 	;	E_xr243 	;	E_xr244 	;	E_xr245 	;	E_xr246 	;	E_xr247 	;	E_xr248 	;	E_xr249 	;	E_xr2410 	;	E_xr2411 	;	E_xr2412 	;	E_xr2413 	;	E_xr2414 	;	E_xr2415 	;	E_xr2416 	;	E_xr2417 	;	E_xr2418 	;	E_xr2419 	;	E_xr2420 	;	E_xr2421 	;	E_xr2422 	;	E_xr2423  ];
    
E_xi = [ E_xi12 ;	E_xi13 	;	E_xi14 	;	E_xi15 	;	E_xi16 	;	E_xi17 	;	E_xi18 	;	E_xi19 	;	E_xi110 	;	E_xi1_11 	;	E_xi1_12 	;	E_xi1_13 	;	E_xi1_14 	;	E_xi1_15 	;	E_xi1_16 	;	E_xi1_17 	;	E_xi1_18 	;	E_xi1_19 	;	E_xi1_20 	;	E_xi1_21 	;	E_xi1_22 	;	E_xi1_23 	;	E_xi1_24 	;	E_xi21 	;	E_xi23 	;	E_xi24 	;	E_xi25 	;	E_xi26 	;	E_xi27 	;	E_xi28 	;	E_xi29 	;	E_xi210 	;	E_xi2_11 	;	E_xi2_12 	;	E_xi2_13 	;	E_xi2_14 	;	E_xi2_15 	;	E_xi2_16 	;	E_xi2_17 	;	E_xi2_18 	;	E_xi2_19 	;	E_xi2_20 	;	E_xi2_21 	;	E_xi2_22 	;	E_xi2_23 	;	E_xi2_24 	;	E_xi31 	;	E_xi32 	;	E_xi34 	;	E_xi35 	;	E_xi36 	;	E_xi37 	;	E_xi38 	;	E_xi39 	;	E_xi310 	;	E_xi311 	;	E_xi312 	;	E_xi313 	;	E_xi314 	;	E_xi315 	;	E_xi316 	;	E_xi317 	;	E_xi318 	;	E_xi319 	;	E_xi320 	;	E_xi321 	;	E_xi322 	;	E_xi323 	;	E_xi324 	;	E_xi41 	;	E_xi42 	;	E_xi43 	;	E_xi45 	;	E_xi46 	;	E_xi47 	;	E_xi48 	;	E_xi49 	;	E_xi410 	;	E_xi411 	;	E_xi412 	;	E_xi413 	;	E_xi414 	;	E_xi415 	;	E_xi416 	;	E_xi417 	;	E_xi418 	;	E_xi419 	;	E_xi420 	;	E_xi421 	;	E_xi422 	;	E_xi423 	;	E_xi424 	;	E_xi51 	;	E_xi52 	;	E_xi53 	;	E_xi54 	;	E_xi56 	;	E_xi57 	;	E_xi58 	;	E_xi59 	;	E_xi510 	;	E_xi511 	;	E_xi512 	;	E_xi513 	;	E_xi514 	;	E_xi515 	;	E_xi516 	;	E_xi517 	;	E_xi518 	;	E_xi519 	;	E_xi520 	;	E_xi521 	;	E_xi522 	;	E_xi523 	;	E_xi524 	;	E_xi61 	;	E_xi62 	;	E_xi63 	;	E_xi64 	;	E_xi65 	;	E_xi67 	;	E_xi68 	;	E_xi69 	;	E_xi610 	;	E_xi611 	;	E_xi612 	;	E_xi613 	;	E_xi614 	;	E_xi615 	;	E_xi616 	;	E_xi617 	;	E_xi618 	;	E_xi619 	;	E_xi620 	;	E_xi621 	;	E_xi622 	;	E_xi623 	;	E_xi624 	;	E_xi71 	;	E_xi72 	;	E_xi73 	;	E_xi74 	;	E_xi75 	;	E_xi76 	;	E_xi78 	;	E_xi79 	;	E_xi710 	;	E_xi711 	;	E_xi712 	;	E_xi713 	;	E_xi714 	;	E_xi715 	;	E_xi716 	;	E_xi717 	;	E_xi718 	;	E_xi719 	;	E_xi720 	;	E_xi721 	;	E_xi722 	;	E_xi723 	;	E_xi724 	;	E_xi81 	;	E_xi82 	;	E_xi83 	;	E_xi84 	;	E_xi85 	;	E_xi86 	;	E_xi87 	;	E_xi89 	;	E_xi810 	;	E_xi811 	;	E_xi812 	;	E_xi813 	;	E_xi814 	;	E_xi815 	;	E_xi816 	;	E_xi817 	;	E_xi818 	;	E_xi819 	;	E_xi820 	;	E_xi821 	;	E_xi822 	;	E_xi823 	;	E_xi824 	;	E_xi91 	;	E_xi92 	;	E_xi93 	;	E_xi94 	;	E_xi95 	;	E_xi96 	;	E_xi97 	;	E_xi98 	;	E_xi910 	;	E_xi911 	;	E_xi912 	;	E_xi913 	;	E_xi914 	;	E_xi915 	;	E_xi916 	;	E_xi917 	;	E_xi918 	;	E_xi919 	;	E_xi920 	;	E_xi921 	;	E_xi922 	;	E_xi923 	;	E_xi924 	;	E_xi101 	;	E_xi102 	;	E_xi103 	;	E_xi104 	;	E_xi105 	;	E_xi106 	;	E_xi107 	;	E_xi108 	;	E_xi109 	;	E_xi1011 	;	E_xi1012 	;	E_xi1013 	;	E_xi1014 	;	E_xi1015 	;	E_xi1016 	;	E_xi1017 	;	E_xi1018 	;	E_xi10_19 	;	E_xi10_20 	;	E_xi10_21 	;	E_xi10_22 	;	E_xi10_23 	;	E_xi10_24 	;	E_xi11_1 	;	E_xi11_2 	;	E_xi113 	;	E_xi114 	;	E_xi115 	;	E_xi116 	;	E_xi117 	;	E_xi118 	;	E_xi119 	;	E_xi1110 	;	E_xi1112 	;	E_xi1113 	;	E_xi1114 	;	E_xi1115 	;	E_xi1116 	;	E_xi1117 	;	E_xi1118 	;	E_xi1019 	;	E_xi1020 	;	E_xi1021 	;	E_xi1022 	;	E_xi1023 	;	E_xi1024 	;	E_xi12_1 	;	E_xi12_2 	;	E_xi123 	;	E_xi124 	;	E_xi125 	;	E_xi126 	;	E_xi127 	;	E_xi128 	;	E_xi129 	;	E_xi1210 	;	E_xi1211 	;	E_xi1213 	;	E_xi1214 	;	E_xi1215 	;	E_xi1216 	;	E_xi1217 	;	E_xi1218 	;	E_xi1219 	;	E_xi1220 	;	E_xi1221 	;	E_xi1222 	;	E_xi1223 	;	E_xi1224 	;	E_xi131 	;	E_xi132 	;	E_xi133 	;	E_xi134 	;	E_xi135 	;	E_xi136 	;	E_xi137 	;	E_xi138 	;	E_xi139 	;	E_xi1310 	;	E_xi1311 	;	E_xi1312 	;	E_xi1314 	;	E_xi1315 	;	E_xi1316 	;	E_xi1317 	;	E_xi1318 	;	E_xi1319 	;	E_xi1320 	;	E_xi1321 	;	E_xi1322 	;	E_xi1323 	;	E_xi1324 	;	E_xi141 	;	E_xi142 	;	E_xi143 	;	E_xi144 	;	E_xi145 	;	E_xi146 	;	E_xi147 	;	E_xi148 	;	E_xi149 	;	E_xi1410 	;	E_xi1411 	;	E_xi1412 	;	E_xi1413 	;	E_xi1415 	;	E_xi1416 	;	E_xi1417 	;	E_xi1418 	;	E_xi1419 	;	E_xi1420 	;	E_xi1421 	;	E_xi1422 	;	E_xi1423 	;	E_xi1424 	;	E_xi151 	;	E_xi152 	;	E_xi153 	;	E_xi154 	;	E_xi155 	;	E_xi156 	;	E_xi157 	;	E_xi158 	;	E_xi159 	;	E_xi1510 	;	E_xi1511 	;	E_xi1512 	;	E_xi1513 	;	E_xi1514 	;	E_xi1516 	;	E_xi1517 	;	E_xi1518 	;	E_xi1519 	;	E_xi1520 	;	E_xi1521 	;	E_xi1522 	;	E_xi1523 	;	E_xi1524 	;	E_xi161 	;	E_xi162 	;	E_xi163 	;	E_xi164 	;	E_xi165 	;	E_xi166 	;	E_xi167 	;	E_xi168 	;	E_xi169 	;	E_xi1610 	;	E_xi1611 	;	E_xi1612 	;	E_xi1613 	;	E_xi1614 	;	E_xi1615 	;	E_xi1617 	;	E_xi1618 	;	E_xi1619 	;	E_xi1620 	;	E_xi1621 	;	E_xi1622 	;	E_xi1623 	;	E_xi1624 	;	E_xi171 	;	E_xi172 	;	E_xi173 	;	E_xi174 	;	E_xi175 	;	E_xi176 	;	E_xi177 	;	E_xi178 	;	E_xi179 	;	E_xi1710 	;	E_xi1711 	;	E_xi1712 	;	E_xi1713 	;	E_xi1714 	;	E_xi1715 	;	E_xi1716 	;	E_xi1718 	;	E_xi1719 	;	E_xi1720 	;	E_xi1721 	;	E_xi1722 	;	E_xi1723 	;	E_xi1724 	;	E_xi181 	;	E_xi182 	;	E_xi183 	;	E_xi184 	;	E_xi185 	;	E_xi186 	;	E_xi187 	;	E_xi188 	;	E_xi189 	;	E_xi1810 	;	E_xi1811 	;	E_xi1812 	;	E_xi1813 	;	E_xi1814 	;	E_xi1815 	;	E_xi1816 	;	E_xi1817 	;	E_xi1819 	;	E_xi1820 	;	E_xi1821 	;	E_xi1822 	;	E_xi1823 	;	E_xi1824 	;	E_xi191 	;	E_xi192 	;	E_xi193 	;	E_xi194 	;	E_xi195 	;	E_xi196 	;	E_xi197 	;	E_xi198 	;	E_xi199 	;	E_xi1910 	;	E_xi1911 	;	E_xi1912 	;	E_xi1913 	;	E_xi1914 	;	E_xi1915 	;	E_xi1916 	;	E_xi1917 	;	E_xi1918 	;	E_xi1920 	;	E_xi1921 	;	E_xi1922 	;	E_xi1923 	;	E_xi1924 	;	E_xi201 	;	E_xi202 	;	E_xi203 	;	E_xi204 	;	E_xi205 	;	E_xi206 	;	E_xi207 	;	E_xi208 	;	E_xi209 	;	E_xi2010 	;	E_xi2011 	;	E_xi2012 	;	E_xi2013 	;	E_xi2014 	;	E_xi2015 	;	E_xi2016 	;	E_xi2017 	;	E_xi2018 	;	E_xi2019 	;	E_xi2021 	;	E_xi2022 	;	E_xi2023 	;	E_xi2024 	;	E_xi21_1 	;	E_xi21_2 	;	E_xi21_3 	;	E_xi214 	;	E_xi215 	;	E_xi216 	;	E_xi217 	;	E_xi218 	;	E_xi219 	;	E_xi2110 	;	E_xi2111 	;	E_xi2112 	;	E_xi2113 	;	E_xi2114 	;	E_xi2115 	;	E_xi2116 	;	E_xi2117 	;	E_xi2118 	;	E_xi2119 	;	E_xi2120 	;	E_xi2122 	;	E_xi2123 	;	E_xi2124 	;	E_xi22_1 	;	E_xi22_2 	;	E_xi22_3 	;	E_xi22_4 	;	E_xi225 	;	E_xi226 	;	E_xi227 	;	E_xi228 	;	E_xi229 	;	E_xi2210 	;	E_xi2211 	;	E_xi2212 	;	E_xi2213 	;	E_xi2214 	;	E_xi2215 	;	E_xi2216 	;	E_xi2217 	;	E_xi2218 	;	E_xi2219 	;	E_xi2220 	;	E_xi2221 	;	E_xi2223 	;	E_xi2224 	;	E_xi231 	;	E_xi232 	;	E_xi233 	;	E_xi234 	;	E_xi235 	;	E_xi236 	;	E_xi237 	;	E_xi238 	;	E_xi239 	;	E_xi2310 	;	E_xi2311 	;	E_xi2312 	;	E_xi2313 	;	E_xi2314 	;	E_xi2315 	;	E_xi2316 	;	E_xi2317 	;	E_xi2318 	;	E_xi2319 	;	E_xi2320 	;	E_xi2321 	;	E_xi2322 	;	E_xi2324 	;	E_xi241 	;	E_xi242 	;	E_xi243 	;	E_xi244 	;	E_xi245 	;	E_xi246 	;	E_xi247 	;	E_xi248 	;	E_xi249 	;	E_xi2410 	;	E_xi2411 	;	E_xi2412 	;	E_xi2413 	;	E_xi2414 	;	E_xi2415 	;	E_xi2416 	;	E_xi2417 	;	E_xi2418 	;	E_xi2419 	;	E_xi2420 	;	E_xi2421 	;	E_xi2422 	;	E_xi2423  ];

E_yr = [ E_yr12 ;	E_yr13 	;	E_yr14 	;	E_yr15 	;	E_yr16 	;	E_yr17 	;	E_yr18 	;	E_yr19 	;	E_yr110 	;	E_yr1_11 	;	E_yr1_12 	;	E_yr1_13 	;	E_yr1_14 	;	E_yr1_15 	;	E_yr1_16 	;	E_yr1_17 	;	E_yr1_18 	;	E_yr1_19 	;	E_yr1_20 	;	E_yr1_21 	;	E_yr1_22 	;	E_yr1_23 	;	E_yr1_24 	;	E_yr21 	;	E_yr23 	;	E_yr24 	;	E_yr25 	;	E_yr26 	;	E_yr27 	;	E_yr28 	;	E_yr29 	;	E_yr210 	;	E_yr2_11 	;	E_yr2_12 	;	E_yr2_13 	;	E_yr2_14 	;	E_yr2_15 	;	E_yr2_16 	;	E_yr2_17 	;	E_yr2_18 	;	E_yr2_19 	;	E_yr2_20 	;	E_yr2_21 	;	E_yr2_22 	;	E_yr2_23 	;	E_yr2_24 	;	E_yr31 	;	E_yr32 	;	E_yr34 	;	E_yr35 	;	E_yr36 	;	E_yr37 	;	E_yr38 	;	E_yr39 	;	E_yr310 	;	E_yr311 	;	E_yr312 	;	E_yr313 	;	E_yr314 	;	E_yr315 	;	E_yr316 	;	E_yr317 	;	E_yr318 	;	E_yr319 	;	E_yr320 	;	E_yr321 	;	E_yr322 	;	E_yr323 	;	E_yr324 	;	E_yr41 	;	E_yr42 	;	E_yr43 	;	E_yr45 	;	E_yr46 	;	E_yr47 	;	E_yr48 	;	E_yr49 	;	E_yr410 	;	E_yr411 	;	E_yr412 	;	E_yr413 	;	E_yr414 	;	E_yr415 	;	E_yr416 	;	E_yr417 	;	E_yr418 	;	E_yr419 	;	E_yr420 	;	E_yr421 	;	E_yr422 	;	E_yr423 	;	E_yr424 	;	E_yr51 	;	E_yr52 	;	E_yr53 	;	E_yr54 	;	E_yr56 	;	E_yr57 	;	E_yr58 	;	E_yr59 	;	E_yr510 	;	E_yr511 	;	E_yr512 	;	E_yr513 	;	E_yr514 	;	E_yr515 	;	E_yr516 	;	E_yr517 	;	E_yr518 	;	E_yr519 	;	E_yr520 	;	E_yr521 	;	E_yr522 	;	E_yr523 	;	E_yr524 	;	E_yr61 	;	E_yr62 	;	E_yr63 	;	E_yr64 	;	E_yr65 	;	E_yr67 	;	E_yr68 	;	E_yr69 	;	E_yr610 	;	E_yr611 	;	E_yr612 	;	E_yr613 	;	E_yr614 	;	E_yr615 	;	E_yr616 	;	E_yr617 	;	E_yr618 	;	E_yr619 	;	E_yr620 	;	E_yr621 	;	E_yr622 	;	E_yr623 	;	E_yr624 	;	E_yr71 	;	E_yr72 	;	E_yr73 	;	E_yr74 	;	E_yr75 	;	E_yr76 	;	E_yr78 	;	E_yr79 	;	E_yr710 	;	E_yr711 	;	E_yr712 	;	E_yr713 	;	E_yr714 	;	E_yr715 	;	E_yr716 	;	E_yr717 	;	E_yr718 	;	E_yr719 	;	E_yr720 	;	E_yr721 	;	E_yr722 	;	E_yr723 	;	E_yr724 	;	E_yr81 	;	E_yr82 	;	E_yr83 	;	E_yr84 	;	E_yr85 	;	E_yr86 	;	E_yr87 	;	E_yr89 	;	E_yr810 	;	E_yr811 	;	E_yr812 	;	E_yr813 	;	E_yr814 	;	E_yr815 	;	E_yr816 	;	E_yr817 	;	E_yr818 	;	E_yr819 	;	E_yr820 	;	E_yr821 	;	E_yr822 	;	E_yr823 	;	E_yr824 	;	E_yr91 	;	E_yr92 	;	E_yr93 	;	E_yr94 	;	E_yr95 	;	E_yr96 	;	E_yr97 	;	E_yr98 	;	E_yr910 	;	E_yr911 	;	E_yr912 	;	E_yr913 	;	E_yr914 	;	E_yr915 	;	E_yr916 	;	E_yr917 	;	E_yr918 	;	E_yr919 	;	E_yr920 	;	E_yr921 	;	E_yr922 	;	E_yr923 	;	E_yr924 	;	E_yr101 	;	E_yr102 	;	E_yr103 	;	E_yr104 	;	E_yr105 	;	E_yr106 	;	E_yr107 	;	E_yr108 	;	E_yr109 	;	E_yr1011 	;	E_yr1012 	;	E_yr1013 	;	E_yr1014 	;	E_yr1015 	;	E_yr1016 	;	E_yr1017 	;	E_yr1018 	;	E_yr10_19 	;	E_yr10_20 	;	E_yr10_21 	;	E_yr10_22 	;	E_yr10_23 	;	E_yr10_24 	;	E_yr11_1 	;	E_yr11_2 	;	E_yr113 	;	E_yr114 	;	E_yr115 	;	E_yr116 	;	E_yr117 	;	E_yr118 	;	E_yr119 	;	E_yr1110 	;	E_yr1112 	;	E_yr1113 	;	E_yr1114 	;	E_yr1115 	;	E_yr1116 	;	E_yr1117 	;	E_yr1118 	;	E_yr1019 	;	E_yr1020 	;	E_yr1021 	;	E_yr1022 	;	E_yr1023 	;	E_yr1024 	;	E_yr12_1 	;	E_yr12_2 	;	E_yr123 	;	E_yr124 	;	E_yr125 	;	E_yr126 	;	E_yr127 	;	E_yr128 	;	E_yr129 	;	E_yr1210 	;	E_yr1211 	;	E_yr1213 	;	E_yr1214 	;	E_yr1215 	;	E_yr1216 	;	E_yr1217 	;	E_yr1218 	;	E_yr1219 	;	E_yr1220 	;	E_yr1221 	;	E_yr1222 	;	E_yr1223 	;	E_yr1224 	;	E_yr131 	;	E_yr132 	;	E_yr133 	;	E_yr134 	;	E_yr135 	;	E_yr136 	;	E_yr137 	;	E_yr138 	;	E_yr139 	;	E_yr1310 	;	E_yr1311 	;	E_yr1312 	;	E_yr1314 	;	E_yr1315 	;	E_yr1316 	;	E_yr1317 	;	E_yr1318 	;	E_yr1319 	;	E_yr1320 	;	E_yr1321 	;	E_yr1322 	;	E_yr1323 	;	E_yr1324 	;	E_yr141 	;	E_yr142 	;	E_yr143 	;	E_yr144 	;	E_yr145 	;	E_yr146 	;	E_yr147 	;	E_yr148 	;	E_yr149 	;	E_yr1410 	;	E_yr1411 	;	E_yr1412 	;	E_yr1413 	;	E_yr1415 	;	E_yr1416 	;	E_yr1417 	;	E_yr1418 	;	E_yr1419 	;	E_yr1420 	;	E_yr1421 	;	E_yr1422 	;	E_yr1423 	;	E_yr1424 	;	E_yr151 	;	E_yr152 	;	E_yr153 	;	E_yr154 	;	E_yr155 	;	E_yr156 	;	E_yr157 	;	E_yr158 	;	E_yr159 	;	E_yr1510 	;	E_yr1511 	;	E_yr1512 	;	E_yr1513 	;	E_yr1514 	;	E_yr1516 	;	E_yr1517 	;	E_yr1518 	;	E_yr1519 	;	E_yr1520 	;	E_yr1521 	;	E_yr1522 	;	E_yr1523 	;	E_yr1524 	;	E_yr161 	;	E_yr162 	;	E_yr163 	;	E_yr164 	;	E_yr165 	;	E_yr166 	;	E_yr167 	;	E_yr168 	;	E_yr169 	;	E_yr1610 	;	E_yr1611 	;	E_yr1612 	;	E_yr1613 	;	E_yr1614 	;	E_yr1615 	;	E_yr1617 	;	E_yr1618 	;	E_yr1619 	;	E_yr1620 	;	E_yr1621 	;	E_yr1622 	;	E_yr1623 	;	E_yr1624 	;	E_yr171 	;	E_yr172 	;	E_yr173 	;	E_yr174 	;	E_yr175 	;	E_yr176 	;	E_yr177 	;	E_yr178 	;	E_yr179 	;	E_yr1710 	;	E_yr1711 	;	E_yr1712 	;	E_yr1713 	;	E_yr1714 	;	E_yr1715 	;	E_yr1716 	;	E_yr1718 	;	E_yr1719 	;	E_yr1720 	;	E_yr1721 	;	E_yr1722 	;	E_yr1723 	;	E_yr1724 	;	E_yr181 	;	E_yr182 	;	E_yr183 	;	E_yr184 	;	E_yr185 	;	E_yr186 	;	E_yr187 	;	E_yr188 	;	E_yr189 	;	E_yr1810 	;	E_yr1811 	;	E_yr1812 	;	E_yr1813 	;	E_yr1814 	;	E_yr1815 	;	E_yr1816 	;	E_yr1817 	;	E_yr1819 	;	E_yr1820 	;	E_yr1821 	;	E_yr1822 	;	E_yr1823 	;	E_yr1824 	;	E_yr191 	;	E_yr192 	;	E_yr193 	;	E_yr194 	;	E_yr195 	;	E_yr196 	;	E_yr197 	;	E_yr198 	;	E_yr199 	;	E_yr1910 	;	E_yr1911 	;	E_yr1912 	;	E_yr1913 	;	E_yr1914 	;	E_yr1915 	;	E_yr1916 	;	E_yr1917 	;	E_yr1918 	;	E_yr1920 	;	E_yr1921 	;	E_yr1922 	;	E_yr1923 	;	E_yr1924 	;	E_yr201 	;	E_yr202 	;	E_yr203 	;	E_yr204 	;	E_yr205 	;	E_yr206 	;	E_yr207 	;	E_yr208 	;	E_yr209 	;	E_yr2010 	;	E_yr2011 	;	E_yr2012 	;	E_yr2013 	;	E_yr2014 	;	E_yr2015 	;	E_yr2016 	;	E_yr2017 	;	E_yr2018 	;	E_yr2019 	;	E_yr2021 	;	E_yr2022 	;	E_yr2023 	;	E_yr2024 	;	E_yr21_1 	;	E_yr21_2 	;	E_yr21_3 	;	E_yr214 	;	E_yr215 	;	E_yr216 	;	E_yr217 	;	E_yr218 	;	E_yr219 	;	E_yr2110 	;	E_yr2111 	;	E_yr2112 	;	E_yr2113 	;	E_yr2114 	;	E_yr2115 	;	E_yr2116 	;	E_yr2117 	;	E_yr2118 	;	E_yr2119 	;	E_yr2120 	;	E_yr2122 	;	E_yr2123 	;	E_yr2124 	;	E_yr22_1 	;	E_yr22_2 	;	E_yr22_3 	;	E_yr22_4 	;	E_yr225 	;	E_yr226 	;	E_yr227 	;	E_yr228 	;	E_yr229 	;	E_yr2210 	;	E_yr2211 	;	E_yr2212 	;	E_yr2213 	;	E_yr2214 	;	E_yr2215 	;	E_yr2216 	;	E_yr2217 	;	E_yr2218 	;	E_yr2219 	;	E_yr2220 	;	E_yr2221 	;	E_yr2223 	;	E_yr2224 	;	E_yr231 	;	E_yr232 	;	E_yr233 	;	E_yr234 	;	E_yr235 	;	E_yr236 	;	E_yr237 	;	E_yr238 	;	E_yr239 	;	E_yr2310 	;	E_yr2311 	;	E_yr2312 	;	E_yr2313 	;	E_yr2314 	;	E_yr2315 	;	E_yr2316 	;	E_yr2317 	;	E_yr2318 	;	E_yr2319 	;	E_yr2320 	;	E_yr2321 	;	E_yr2322 	;	E_yr2324 	;	E_yr241 	;	E_yr242 	;	E_yr243 	;	E_yr244 	;	E_yr245 	;	E_yr246 	;	E_yr247 	;	E_yr248 	;	E_yr249 	;	E_yr2410 	;	E_yr2411 	;	E_yr2412 	;	E_yr2413 	;	E_yr2414 	;	E_yr2415 	;	E_yr2416 	;	E_yr2417 	;	E_yr2418 	;	E_yr2419 	;	E_yr2420 	;	E_yr2421 	;	E_yr2422 	;	E_yr2423  ];

E_yi = [ E_yi12 ;	E_yi13 	;	E_yi14 	;	E_yi15 	;	E_yi16 	;	E_yi17 	;	E_yi18 	;	E_yi19 	;	E_yi110 	;	E_yi1_11 	;	E_yi1_12 	;	E_yi1_13 	;	E_yi1_14 	;	E_yi1_15 	;	E_yi1_16 	;	E_yi1_17 	;	E_yi1_18 	;	E_yi1_19 	;	E_yi1_20 	;	E_yi1_21 	;	E_yi1_22 	;	E_yi1_23 	;	E_yi1_24 	;	E_yi21 	;	E_yi23 	;	E_yi24 	;	E_yi25 	;	E_yi26 	;	E_yi27 	;	E_yi28 	;	E_yi29 	;	E_yi210 	;	E_yi2_11 	;	E_yi2_12 	;	E_yi2_13 	;	E_yi2_14 	;	E_yi2_15 	;	E_yi2_16 	;	E_yi2_17 	;	E_yi2_18 	;	E_yi2_19 	;	E_yi2_20 	;	E_yi2_21 	;	E_yi2_22 	;	E_yi2_23 	;	E_yi2_24 	;	E_yi31 	;	E_yi32 	;	E_yi34 	;	E_yi35 	;	E_yi36 	;	E_yi37 	;	E_yi38 	;	E_yi39 	;	E_yi310 	;	E_yi311 	;	E_yi312 	;	E_yi313 	;	E_yi314 	;	E_yi315 	;	E_yi316 	;	E_yi317 	;	E_yi318 	;	E_yi319 	;	E_yi320 	;	E_yi321 	;	E_yi322 	;	E_yi323 	;	E_yi324 	;	E_yi41 	;	E_yi42 	;	E_yi43 	;	E_yi45 	;	E_yi46 	;	E_yi47 	;	E_yi48 	;	E_yi49 	;	E_yi410 	;	E_yi411 	;	E_yi412 	;	E_yi413 	;	E_yi414 	;	E_yi415 	;	E_yi416 	;	E_yi417 	;	E_yi418 	;	E_yi419 	;	E_yi420 	;	E_yi421 	;	E_yi422 	;	E_yi423 	;	E_yi424 	;	E_yi51 	;	E_yi52 	;	E_yi53 	;	E_yi54 	;	E_yi56 	;	E_yi57 	;	E_yi58 	;	E_yi59 	;	E_yi510 	;	E_yi511 	;	E_yi512 	;	E_yi513 	;	E_yi514 	;	E_yi515 	;	E_yi516 	;	E_yi517 	;	E_yi518 	;	E_yi519 	;	E_yi520 	;	E_yi521 	;	E_yi522 	;	E_yi523 	;	E_yi524 	;	E_yi61 	;	E_yi62 	;	E_yi63 	;	E_yi64 	;	E_yi65 	;	E_yi67 	;	E_yi68 	;	E_yi69 	;	E_yi610 	;	E_yi611 	;	E_yi612 	;	E_yi613 	;	E_yi614 	;	E_yi615 	;	E_yi616 	;	E_yi617 	;	E_yi618 	;	E_yi619 	;	E_yi620 	;	E_yi621 	;	E_yi622 	;	E_yi623 	;	E_yi624 	;	E_yi71 	;	E_yi72 	;	E_yi73 	;	E_yi74 	;	E_yi75 	;	E_yi76 	;	E_yi78 	;	E_yi79 	;	E_yi710 	;	E_yi711 	;	E_yi712 	;	E_yi713 	;	E_yi714 	;	E_yi715 	;	E_yi716 	;	E_yi717 	;	E_yi718 	;	E_yi719 	;	E_yi720 	;	E_yi721 	;	E_yi722 	;	E_yi723 	;	E_yi724 	;	E_yi81 	;	E_yi82 	;	E_yi83 	;	E_yi84 	;	E_yi85 	;	E_yi86 	;	E_yi87 	;	E_yi89 	;	E_yi810 	;	E_yi811 	;	E_yi812 	;	E_yi813 	;	E_yi814 	;	E_yi815 	;	E_yi816 	;	E_yi817 	;	E_yi818 	;	E_yi819 	;	E_yi820 	;	E_yi821 	;	E_yi822 	;	E_yi823 	;	E_yi824 	;	E_yi91 	;	E_yi92 	;	E_yi93 	;	E_yi94 	;	E_yi95 	;	E_yi96 	;	E_yi97 	;	E_yi98 	;	E_yi910 	;	E_yi911 	;	E_yi912 	;	E_yi913 	;	E_yi914 	;	E_yi915 	;	E_yi916 	;	E_yi917 	;	E_yi918 	;	E_yi919 	;	E_yi920 	;	E_yi921 	;	E_yi922 	;	E_yi923 	;	E_yi924 	;	E_yi101 	;	E_yi102 	;	E_yi103 	;	E_yi104 	;	E_yi105 	;	E_yi106 	;	E_yi107 	;	E_yi108 	;	E_yi109 	;	E_yi1011 	;	E_yi1012 	;	E_yi1013 	;	E_yi1014 	;	E_yi1015 	;	E_yi1016 	;	E_yi1017 	;	E_yi1018 	;	E_yi10_19 	;	E_yi10_20 	;	E_yi10_21 	;	E_yi10_22 	;	E_yi10_23 	;	E_yi10_24 	;	E_yi11_1 	;	E_yi11_2 	;	E_yi113 	;	E_yi114 	;	E_yi115 	;	E_yi116 	;	E_yi117 	;	E_yi118 	;	E_yi119 	;	E_yi1110 	;	E_yi1112 	;	E_yi1113 	;	E_yi1114 	;	E_yi1115 	;	E_yi1116 	;	E_yi1117 	;	E_yi1118 	;	E_yi1019 	;	E_yi1020 	;	E_yi1021 	;	E_yi1022 	;	E_yi1023 	;	E_yi1024 	;	E_yi12_1 	;	E_yi12_2 	;	E_yi123 	;	E_yi124 	;	E_yi125 	;	E_yi126 	;	E_yi127 	;	E_yi128 	;	E_yi129 	;	E_yi1210 	;	E_yi1211 	;	E_yi1213 	;	E_yi1214 	;	E_yi1215 	;	E_yi1216 	;	E_yi1217 	;	E_yi1218 	;	E_yi1219 	;	E_yi1220 	;	E_yi1221 	;	E_yi1222 	;	E_yi1223 	;	E_yi1224 	;	E_yi131 	;	E_yi132 	;	E_yi133 	;	E_yi134 	;	E_yi135 	;	E_yi136 	;	E_yi137 	;	E_yi138 	;	E_yi139 	;	E_yi1310 	;	E_yi1311 	;	E_yi1312 	;	E_yi1314 	;	E_yi1315 	;	E_yi1316 	;	E_yi1317 	;	E_yi1318 	;	E_yi1319 	;	E_yi1320 	;	E_yi1321 	;	E_yi1322 	;	E_yi1323 	;	E_yi1324 	;	E_yi141 	;	E_yi142 	;	E_yi143 	;	E_yi144 	;	E_yi145 	;	E_yi146 	;	E_yi147 	;	E_yi148 	;	E_yi149 	;	E_yi1410 	;	E_yi1411 	;	E_yi1412 	;	E_yi1413 	;	E_yi1415 	;	E_yi1416 	;	E_yi1417 	;	E_yi1418 	;	E_yi1419 	;	E_yi1420 	;	E_yi1421 	;	E_yi1422 	;	E_yi1423 	;	E_yi1424 	;	E_yi151 	;	E_yi152 	;	E_yi153 	;	E_yi154 	;	E_yi155 	;	E_yi156 	;	E_yi157 	;	E_yi158 	;	E_yi159 	;	E_yi1510 	;	E_yi1511 	;	E_yi1512 	;	E_yi1513 	;	E_yi1514 	;	E_yi1516 	;	E_yi1517 	;	E_yi1518 	;	E_yi1519 	;	E_yi1520 	;	E_yi1521 	;	E_yi1522 	;	E_yi1523 	;	E_yi1524 	;	E_yi161 	;	E_yi162 	;	E_yi163 	;	E_yi164 	;	E_yi165 	;	E_yi166 	;	E_yi167 	;	E_yi168 	;	E_yi169 	;	E_yi1610 	;	E_yi1611 	;	E_yi1612 	;	E_yi1613 	;	E_yi1614 	;	E_yi1615 	;	E_yi1617 	;	E_yi1618 	;	E_yi1619 	;	E_yi1620 	;	E_yi1621 	;	E_yi1622 	;	E_yi1623 	;	E_yi1624 	;	E_yi171 	;	E_yi172 	;	E_yi173 	;	E_yi174 	;	E_yi175 	;	E_yi176 	;	E_yi177 	;	E_yi178 	;	E_yi179 	;	E_yi1710 	;	E_yi1711 	;	E_yi1712 	;	E_yi1713 	;	E_yi1714 	;	E_yi1715 	;	E_yi1716 	;	E_yi1718 	;	E_yi1719 	;	E_yi1720 	;	E_yi1721 	;	E_yi1722 	;	E_yi1723 	;	E_yi1724 	;	E_yi181 	;	E_yi182 	;	E_yi183 	;	E_yi184 	;	E_yi185 	;	E_yi186 	;	E_yi187 	;	E_yi188 	;	E_yi189 	;	E_yi1810 	;	E_yi1811 	;	E_yi1812 	;	E_yi1813 	;	E_yi1814 	;	E_yi1815 	;	E_yi1816 	;	E_yi1817 	;	E_yi1819 	;	E_yi1820 	;	E_yi1821 	;	E_yi1822 	;	E_yi1823 	;	E_yi1824 	;	E_yi191 	;	E_yi192 	;	E_yi193 	;	E_yi194 	;	E_yi195 	;	E_yi196 	;	E_yi197 	;	E_yi198 	;	E_yi199 	;	E_yi1910 	;	E_yi1911 	;	E_yi1912 	;	E_yi1913 	;	E_yi1914 	;	E_yi1915 	;	E_yi1916 	;	E_yi1917 	;	E_yi1918 	;	E_yi1920 	;	E_yi1921 	;	E_yi1922 	;	E_yi1923 	;	E_yi1924 	;	E_yi201 	;	E_yi202 	;	E_yi203 	;	E_yi204 	;	E_yi205 	;	E_yi206 	;	E_yi207 	;	E_yi208 	;	E_yi209 	;	E_yi2010 	;	E_yi2011 	;	E_yi2012 	;	E_yi2013 	;	E_yi2014 	;	E_yi2015 	;	E_yi2016 	;	E_yi2017 	;	E_yi2018 	;	E_yi2019 	;	E_yi2021 	;	E_yi2022 	;	E_yi2023 	;	E_yi2024 	;	E_yi21_1 	;	E_yi21_2 	;	E_yi21_3 	;	E_yi214 	;	E_yi215 	;	E_yi216 	;	E_yi217 	;	E_yi218 	;	E_yi219 	;	E_yi2110 	;	E_yi2111 	;	E_yi2112 	;	E_yi2113 	;	E_yi2114 	;	E_yi2115 	;	E_yi2116 	;	E_yi2117 	;	E_yi2118 	;	E_yi2119 	;	E_yi2120 	;	E_yi2122 	;	E_yi2123 	;	E_yi2124 	;	E_yi22_1 	;	E_yi22_2 	;	E_yi22_3 	;	E_yi22_4 	;	E_yi225 	;	E_yi226 	;	E_yi227 	;	E_yi228 	;	E_yi229 	;	E_yi2210 	;	E_yi2211 	;	E_yi2212 	;	E_yi2213 	;	E_yi2214 	;	E_yi2215 	;	E_yi2216 	;	E_yi2217 	;	E_yi2218 	;	E_yi2219 	;	E_yi2220 	;	E_yi2221 	;	E_yi2223 	;	E_yi2224 	;	E_yi231 	;	E_yi232 	;	E_yi233 	;	E_yi234 	;	E_yi235 	;	E_yi236 	;	E_yi237 	;	E_yi238 	;	E_yi239 	;	E_yi2310 	;	E_yi2311 	;	E_yi2312 	;	E_yi2313 	;	E_yi2314 	;	E_yi2315 	;	E_yi2316 	;	E_yi2317 	;	E_yi2318 	;	E_yi2319 	;	E_yi2320 	;	E_yi2321 	;	E_yi2322 	;	E_yi2324 	;	E_yi241 	;	E_yi242 	;	E_yi243 	;	E_yi244 	;	E_yi245 	;	E_yi246 	;	E_yi247 	;	E_yi248 	;	E_yi249 	;	E_yi2410 	;	E_yi2411 	;	E_yi2412 	;	E_yi2413 	;	E_yi2414 	;	E_yi2415 	;	E_yi2416 	;	E_yi2417 	;	E_yi2418 	;	E_yi2419 	;	E_yi2420 	;	E_yi2421 	;	E_yi2422 	;	E_yi2423  ];

%% aqui somamos as componentes (elementos das colunas) em cada dire��o e elevamos ao quadrado

Exr = sum(E_xr).^2;

Exi = sum(E_xi).^2;

Eyr = sum(E_yr).^2;

Eyi = sum(E_yi).^2;

%% aqui retira-se o m�dulo dos quadrados:

% V/m -> kV/cm = 10^-3/10^2
% V/m -> V/cm = *10^-2
Erms = ((Exr + Exi + Eyr + Eyi).^(1/2))*(10^-2);
thetad = rad2deg(theta); %converte o �ngulo de radianos para graus

%vetores para plot da curva refer�ncia condutor 1
xp5 = [ 0.382086461	,	1.28294072	,	2.63266891	,	3.982397101	,	5.332125291	,	6.681853482	,	8.480455604	,	9.830183795	,	11.62878592	,	13.42738804	,	15.22599016	,	17.02459229	,	18.82319441	,	20.62179653	,	22.86927259	,	24.66787471	,	26.46570023	,	28.26352576	,	29.61170075	,	31.85840021	,	34.10432306	,	36.35024592	,	38.59539218	,	41.28941237	,	43.53378203	,	46.22702563	,	48.47139529	,	51.16308568	,	53.40667875	,	55.64949521	,	57.89231167	,	60.13435154	,	62.37639141	,	65.0673052	,	67.30934507	,	69.55138493	,	71.34377427	,	73.1369402	,	74.92932954	,	76.27362154	,	78.06601087	,	79.85840021	,	81.65078954	,	82.99430494	,	84.33782035	,	86.13020968	,	87.92182242	,	89.26533782	,	90.60807662	,	92.39968936	,	93.74242816	,	95.5340409	,	96.8767797	,	98.21951851	,	100.0103546	,	101.8011908	,	103.143153	,	104.4851152	,	105.8270774	,	107.1698162	,	108.5117784	,	109.4056433	,	110.7483821	,	112.0903443	,	112.9842092	,	114.3261714	,	115.6681336	,	117.0100958	,	118.352058	,	119.6940202	,	121.0359824	,	122.3787212	,	123.7206834	,	125.0626456	,	126.4046078	,	127.74657	,	129.0885322	,	130.4304944	,	131.7724566	,	132.6663215	,	134.0090603	,	134.9029252	,	136.245664	,	137.5876262	,	138.4807145	,	139.8226767	,	141.1654155	,	142.5081543	,	143.8508931	,	144.7455346	,	146.0874968	,	147.429459	,	148.7721978	,	150.11416	,	151.4568988	,	152.7996376	,	154.1423764	,	155.4851152	,	156.827854	,	158.1705928	,	159.512555	,	160.8552938	,	162.1980326	,	163.5407714	,	164.8835102	,	166.226249	,	168.0170852	,	169.3590474	,	171.1506601	,	172.9414962	,	174.7323324	,	176.0750712	,	177.8659073	,	179.6567435	,	181.4475796	,	183.2391923	,	185.0308051	,	186.8231944	,	188.6155837	,	190.4071965	,	192.1995858	,	194.4400725	,	196.6805592	,	198.4729485	,	200.7142118	,	202.9562516	,	205.1975149	,	207.4387782	,	209.680818	,	212.3717318	,	214.6137717	,	216.8565881	,	219.1001812	,	221.3437743	,	223.5873673	,	226.2806109	,	228.9738545	,	231.6678747	,	233.9137976	,	236.1597204	,	238.8545172	,	241.1012167	,	243.3486927	,	245.1457417	,	247.3924411	,	249.1918198	,	250.9896454	,	253.2371214	,	255.0365001	,	256.3870049	,	258.1848304	,	259.9842092	,	261.7835879	,	263.5837432	,	265.8327725	,	267.6321512	,	269.4307533	,	270.7812581	,	272.1317629	,	273.9311416	,	275.7305203	,	277.5306756	,	278.8811804	,	280.2316852	,	282.0318405	,	283.3823453	,	284.7328501	,	286.5330054	,	287.8835102	,	289.6836655	,	291.4838209	,	293.2839762	,	294.6352576	,	295.9857624	,	297.3370437	,	298.6875485	,	300.0388299	,	301.3901113	,	302.7413927	,	304.541548	,	306.3417033	,	307.6929847	,	309.0442661	,	310.3955475	,	311.7460523	,	313.0973337	,	313.9989645	,	315.3494693	,	316.7007507	,	318.0520321	,	319.4033135	,	320.3049443	,	321.2065752	,	322.55708	,	323.9083614	,	325.2596428	,	326.6109242	,	327.9622055	,	329.3134869	,	330.6647683	,	332.0160497	,	333.3681077	,	335.168263	,	336.5195444	,	337.8700492	,	339.220554	,	340.5718354	,	341.9231167	,	343.2736215	,	344.6249029	,	345.9754077	,	347.3259125	,	348.6764173	,	350.0269221	,	351.3774269	,	352.7279317	,	354.0784364	,	355.4289412	,	356.7786694	,	358.1291742	,	359.479679 ];
yp5 = [ 0.030934352	,	0.030844611	,	0.030754983	,	0.030665355	,	0.030575727	,	0.030486099	,	0.030396583	,	0.030306955	,	0.030217439	,	0.030127923	,	0.030038407	,	0.029948891	,	0.029859375	,	0.029769859	,	0.029680456	,	0.02959094	,	0.029523915	,	0.02945689	,	0.029412244	,	0.029345332	,	0.02930091	,	0.029256489	,	0.029234559	,	0.029212741	,	0.029213301	,	0.029213974	,	0.029214535	,	0.029260191	,	0.029283243	,	0.029328786	,	0.029374329	,	0.029442363	,	0.029510398	,	0.029578544	,	0.029646579	,	0.029714613	,	0.029805026	,	0.029872948	,	0.029963362	,	0.030031172	,	0.030121585	,	0.030211998	,	0.030302412	,	0.030392713	,	0.030483014	,	0.030573427	,	0.030686332	,	0.030776633	,	0.030889425	,	0.03100233	,	0.031115122	,	0.031228027	,	0.031340819	,	0.031453611	,	0.031589007	,	0.031724402	,	0.031859686	,	0.031994969	,	0.032130253	,	0.032243045	,	0.032378329	,	0.032491009	,	0.032603801	,	0.032739084	,	0.032851765	,	0.032987048	,	0.033122332	,	0.033257615	,	0.033392898	,	0.033528182	,	0.033663465	,	0.033776258	,	0.033911541	,	0.034046825	,	0.034182108	,	0.034317391	,	0.034452675	,	0.034587958	,	0.034723242	,	0.034835922	,	0.034948714	,	0.035061394	,	0.035174187	,	0.03530947	,	0.035444641	,	0.035579925	,	0.035692717	,	0.03580551	,	0.035918302	,	0.036008491	,	0.036143774	,	0.036279058	,	0.03639185	,	0.036527133	,	0.036639926	,	0.036752718	,	0.03686551	,	0.036978303	,	0.037091095	,	0.037203887	,	0.037339171	,	0.037451963	,	0.037564755	,	0.037677548	,	0.03779034	,	0.037903132	,	0.038038528	,	0.038173811	,	0.038286716	,	0.038422111	,	0.038557507	,	0.038670299	,	0.038805695	,	0.038941091	,	0.039076486	,	0.039189391	,	0.039302295	,	0.039392709	,	0.039483122	,	0.039596026	,	0.03968644	,	0.039799456	,	0.039912473	,	0.040002886	,	0.040093412	,	0.040161446	,	0.040251972	,	0.040342497	,	0.040410532	,	0.040478678	,	0.040546712	,	0.040592256	,	0.040615308	,	0.04063836	,	0.040661412	,	0.040662085	,	0.040662758	,	0.04064094	,	0.040596518	,	0.040552097	,	0.040507788	,	0.040440875	,	0.040351471	,	0.040306938	,	0.040240025	,	0.040128018	,	0.040060993	,	0.039971589	,	0.039859582	,	0.039747463	,	0.039680438	,	0.039568431	,	0.039456424	,	0.039321926	,	0.03918754	,	0.039075533	,	0.038986017	,	0.038873898	,	0.038761778	,	0.038649771	,	0.038537764	,	0.038403266	,	0.038291147	,	0.038179028	,	0.038044529	,	0.03793241	,	0.037820291	,	0.037685793	,	0.037573673	,	0.037439175	,	0.037304677	,	0.037170179	,	0.037035568	,	0.036923449	,	0.036788839	,	0.036676719	,	0.036542109	,	0.036407498	,	0.036272888	,	0.03613839	,	0.036003892	,	0.035869281	,	0.035734671	,	0.03560006	,	0.035487941	,	0.035353331	,	0.035241099	,	0.03512898	,	0.03499437	,	0.034859759	,	0.034725149	,	0.034612917	,	0.034500686	,	0.034388567	,	0.034253956	,	0.034119346	,	0.033984736	,	0.033850125	,	0.033715515	,	0.033580904	,	0.033446294	,	0.033289192	,	0.033154694	,	0.033020084	,	0.032907964	,	0.032795845	,	0.032661235	,	0.032526624	,	0.032414505	,	0.032279895	,	0.032167775	,	0.032055656	,	0.031943537	,	0.031831418	,	0.031719298	,	0.031607179	,	0.03149506	,	0.031382941	,	0.031293313	,	0.031181193	,	0.031069074 ];

%vetores para plot da curva refer�ncia condutor 2

figure();
plot(thetad,Erms,'LineWidth',4.5);
ylabel('Campo el�trico (V/cm)')
xlabel('�ngulo (�)')
%title('Campo el�trico superficial condutor cinco') %t�tulo (trocar ao mudar condutor avaliado)
grid
hold on
% yy5 = smooth(xp5,yp5,'rloess'); %comando suavizar a curva refer�ncia condutor 5
% plot(xp5,yy5,'--','LineWidth',2.5); %plot curva refer�ncia (trocar vari�veis ao mudar condutor avaliado)
plot(xp5,yp5,'--','LineWidth',2.5); %plot curva refer�ncia (trocar vari�veis ao mudar condutor avaliado)
xlim([0 360])
legend('Simulado', '(SARMA, 1969)')

h = figure();
plot(x(1:4), y(1:4),'o', 'color', 'blue')
hold on
plot(x(5:8), y(5:8),'*', 'color', 'blue')
hold on
plot(x(9:12), y(9:12),'x', 'color', 'blue')
xlabel('Faixa de passagem (m)')
ylabel('Altura do cabo (m)')
%title('Configura��o geom�trica dos condutores')
legend('Cabos condutores fase A','Cabos condutores fase B','Cabos condutores fase C');
%title('Configura��o geom�trica dos condutores')
grid
xlim([-20 20])
ylim([17 25])

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'caso5configsup','-dpdf','-r0')

Emax1 = max(Erms);
Emin1 = min(Erms);

Emax2 = max(yp5);
Emin2 = min(yp5);

err = abs((yp5 - Erms)/yp5)*100;
