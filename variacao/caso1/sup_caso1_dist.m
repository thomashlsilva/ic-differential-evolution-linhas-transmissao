% APROXIMA 1 E 2 m em x da configura��o original
clear all; close all; clc

e_0 = 8.854*(10^(-12));
r = 9.155*(10^-3); % raio do condutor (est� igual ao seu trabalho)
n = 3; % n�mero de condutores do sistema
nc = n*2; % n�mero de condutores total (com imagens)
ci = (2*n-1); % n�mero de cargas imagens

%% matriz P

% posi��o dos condutores reais do sistema
xr = [-2 0 2];
yr = [14.01 14.01 14.01];
yi = -yr;

x = [xr xr]; % posi��o x dos condutores (com imagens)
y = [yr -yr]; % posi��o y dos condutores (com imagens)

theta = linspace(0,2*pi,360); %gera a superf�cie do condutor em 360 pontos

xcj = x(3); % posi��o x do centro condutor (para trocar o condutor avaliado deve-se mudar este valor)
ycj = y(3); % posi��o y do centro condutor

xf = r.*cos(theta) + xcj; % eixo x do ponto de avalia��o
yf = r.*sin(theta) + ycj; % eixo y do ponto de avalia��o

%% 

% posi��o dos condutores reais do sistema
xr1 = [-1 0 1];
yr1 = [14.01 14.01 14.01];
yi1 = -yr1;

x1 = [xr1 xr1]; % posi��o x dos condutores (com imagens)
y1 = [yr1 -yr1]; % posi��o y dos condutores (com imagens)

xcj1 = x1(3); % posi��o x do centro condutor (para trocar o condutor avaliado deve-se mudar este valor)
ycj1 = y1(3); % posi��o y do centro condutor

xf1 = r.*cos(theta) + xcj1; % eixo x do ponto de avalia��o
yf1 = r.*sin(theta) + ycj1; % eixo y do ponto de avalia��o

%%

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

V = 145*10^3;

%tens�o condutor 1 fase a
V_ra = V/sqrt(3);
V_ia = 0;

%tens�o condutor 2 fase b
V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

%tens�o condutor 3 fase c
V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

Vf = [ V_ra+V_ia ; V_rb+V_ib ; V_rc+V_ic ];

%% C�lculo densidade de carga

rho = P\Vf;

% busco cada uma das posicoes de rho real e imaginario por fase

% rho condutor 1 fase a
rho_ra = real(rho(1));
rho_ia = imag(rho(1));

% rho condutor 2 fase b
rho_rb = real(rho(2));
rho_ib = imag(rho(2));

% rho condutor 3 fase c
rho_rc = real(rho(3));
rho_ic = imag(rho(3));

%% Dist�ncia entre condutores

%sendo: i -> o n�mero do primeiro condutor e j -> o n�mero do segundo condutor
%exemplo: 1,2 � a dist�ncia entre condutor 1 e o condutor 2, os pr�ximos
%c�lculos seguem esta mesma l�gica
%sendo 4 a imagem do condutor 1, 5 a imagem do condutor 2 e 6 a imagem do
%condutor 3

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

%% E_xr componente x real campo el�trico condutor 2 fase b assim segue:

%busca o valor da carga imagem no condutor 1 real, valor do eixo x dos
%pontos de avalia��o e as posi��es geradas pelas matrizes posx e posy das
%cargas imagens

E_xr12 = (-rho_rb/(2*pi*e_0)).*((xf - posx(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_xr13 = (-rho_rc/(2*pi*e_0)).*((xf - posx(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_xr14 = (rho_ra/(2*pi*e_0)).*((xf - posx(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_xr15 = (rho_rb/(2*pi*e_0)).*((xf - posx(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_xr16 = (rho_rc/(2*pi*e_0)).*((xf - posx(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));

E_xr21 = (-rho_ra/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xr23 = (-rho_rc/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xr24 = (rho_ra/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xr25 = (rho_rb/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xr26 = (rho_rc/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));

E_xr31 = (-rho_ra/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xr32 = (-rho_rb/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xr34 = (rho_ra/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xr35 = (rho_rb/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xr36 = (rho_rc/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));

E_xr41 = (-rho_ra/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xr42 = (-rho_rb/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xr43 = (-rho_rc/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xr45 = (rho_rb/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xr46 = (rho_rc/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));

E_xr51 = (-rho_ra/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xr52 = (-rho_rb/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xr53 = (-rho_rc/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xr54 = (rho_ra/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xr56 = (rho_rc/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));

E_xr61 = (-rho_ra/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xr62 = (-rho_rb/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xr63 = (-rho_rc/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xr64 = (rho_ra/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xr65 = (rho_rb/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));


%% E_xi componente x imaginario campo el�trico condutor 2 fase b assim segue:

E_xi12 = (-rho_ib/(2*pi*e_0)).*((xf - posx(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_xi13 = (-rho_ic/(2*pi*e_0)).*((xf - posx(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_xi14 = (rho_ia/(2*pi*e_0)).*((xf - posx(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_xi15 = (rho_ib/(2*pi*e_0)).*((xf - posx(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_xi16 = (rho_ic/(2*pi*e_0)).*((xf - posx(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));

E_xi21 = (-rho_ia/(2*pi*e_0)).*((xf - posx(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_xi23 = (-rho_ic/(2*pi*e_0)).*((xf - posx(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_xi24 = (rho_ia/(2*pi*e_0)).*((xf - posx(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_xi25 = (rho_ib/(2*pi*e_0)).*((xf - posx(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_xi26 = (rho_ic/(2*pi*e_0)).*((xf - posx(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));

E_xi31 = (-rho_ia/(2*pi*e_0)).*((xf - posx(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_xi32 = (-rho_ib/(2*pi*e_0)).*((xf - posx(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_xi34 = (rho_ia/(2*pi*e_0)).*((xf - posx(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_xi35 = (rho_ib/(2*pi*e_0)).*((xf - posx(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_xi36 = (rho_ic/(2*pi*e_0)).*((xf - posx(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));

E_xi41 = (-rho_ia/(2*pi*e_0)).*((xf - posx(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_xi42 = (-rho_ib/(2*pi*e_0)).*((xf - posx(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_xi43 = (-rho_ic/(2*pi*e_0)).*((xf - posx(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_xi45 = (rho_ib/(2*pi*e_0)).*((xf - posx(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_xi46 = (rho_ic/(2*pi*e_0)).*((xf - posx(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));

E_xi51 = (-rho_ia/(2*pi*e_0)).*((xf - posx(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_xi52 = (-rho_ib/(2*pi*e_0)).*((xf - posx(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_xi53 = (-rho_ic/(2*pi*e_0)).*((xf - posx(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_xi54 = (rho_ia/(2*pi*e_0)).*((xf - posx(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_xi56 = (rho_ic/(2*pi*e_0)).*((xf - posx(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));

E_xi61 = (-rho_ia/(2*pi*e_0)).*((xf - posx(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_xi62 = (-rho_ib/(2*pi*e_0)).*((xf - posx(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_xi63 = (-rho_ic/(2*pi*e_0)).*((xf - posx(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_xi64 = (rho_ia/(2*pi*e_0)).*((xf - posx(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_xi65 = (rho_ib/(2*pi*e_0)).*((xf - posx(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));

%% E_yr componente y real campo el�trico condutor 2 fase b

E_yr12 = (-rho_rb/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yr13 = (-rho_rc/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yr14 = (rho_ra/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yr15 = (rho_rb/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yr16 = (rho_rc/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));

E_yr21 = (-rho_ra/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yr23 = (-rho_rc/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yr24 = (rho_ra/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yr25 = (rho_rb/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yr26 = (rho_rc/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));

E_yr31 = (-rho_ra/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yr32 = (-rho_rb/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yr34 = (rho_ra/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yr35 = (rho_rb/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yr36 = (rho_rc/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));

E_yr41 = (-rho_ra/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yr42 = (-rho_rb/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yr43 = (-rho_rc/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yr45 = (rho_rb/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yr46 = (rho_rc/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));

E_yr51 = (-rho_ra/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yr52 = (-rho_rb/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yr53 = (-rho_rc/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yr54 = (rho_ra/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yr56 = (rho_rc/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));

E_yr61 = (-rho_ra/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yr62 = (-rho_rb/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yr63 = (-rho_rc/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yr64 = (rho_ra/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yr65 = (rho_rb/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));


%% E_yi21 componente y imaginario campo el�trico condutor 2 fase b

E_yi12 = (-rho_ib/(2*pi*e_0)).*((yf - posy(1,2))./((xf - posx(1,2)).^2 + (yf - posy(1,2)).^2));
E_yi13 = (-rho_ic/(2*pi*e_0)).*((yf - posy(1,3))./((xf - posx(1,3)).^2 + (yf - posy(1,3)).^2));
E_yi14 = (rho_ia/(2*pi*e_0)).*((yf - posy(1,4))./((xf - posx(1,4)).^2 + (yf - posy(1,4)).^2));
E_yi15 = (rho_ib/(2*pi*e_0)).*((yf - posy(1,5))./((xf - posx(1,5)).^2 + (yf - posy(1,5)).^2));
E_yi16 = (rho_ic/(2*pi*e_0)).*((yf - posy(1,6))./((xf - posx(1,6)).^2 + (yf - posy(1,6)).^2));

E_yi21 = (-rho_ia/(2*pi*e_0)).*((yf - posy(2,1))./((xf - posx(2,1)).^2 + (yf - posy(2,1)).^2));
E_yi23 = (-rho_ic/(2*pi*e_0)).*((yf - posy(2,3))./((xf - posx(2,3)).^2 + (yf - posy(2,3)).^2));
E_yi24 = (rho_ia/(2*pi*e_0)).*((yf - posy(2,4))./((xf - posx(2,4)).^2 + (yf - posy(2,4)).^2));
E_yi25 = (rho_ib/(2*pi*e_0)).*((yf - posy(2,5))./((xf - posx(2,5)).^2 + (yf - posy(2,5)).^2));
E_yi26 = (rho_ic/(2*pi*e_0)).*((yf - posy(2,6))./((xf - posx(2,6)).^2 + (yf - posy(2,6)).^2));

E_yi31 = (-rho_ia/(2*pi*e_0)).*((yf - posy(3,1))./((xf - posx(3,1)).^2 + (yf - posy(3,1)).^2));
E_yi32 = (-rho_ib/(2*pi*e_0)).*((yf - posy(3,2))./((xf - posx(3,2)).^2 + (yf - posy(3,2)).^2));
E_yi34 = (rho_ia/(2*pi*e_0)).*((yf - posy(3,4))./((xf - posx(3,4)).^2 + (yf - posy(3,4)).^2));
E_yi35 = (rho_ib/(2*pi*e_0)).*((yf - posy(3,5))./((xf - posx(3,5)).^2 + (yf - posy(3,5)).^2));
E_yi36 = (rho_ic/(2*pi*e_0)).*((yf - posy(3,6))./((xf - posx(3,6)).^2 + (yf - posy(3,6)).^2));

E_yi41 = (-rho_ia/(2*pi*e_0)).*((yf - posy(4,1))./((xf - posx(4,1)).^2 + (yf - posy(4,1)).^2));
E_yi42 = (-rho_ib/(2*pi*e_0)).*((yf - posy(4,2))./((xf - posx(4,2)).^2 + (yf - posy(4,2)).^2));
E_yi43 = (-rho_ic/(2*pi*e_0)).*((yf - posy(4,3))./((xf - posx(4,3)).^2 + (yf - posy(4,3)).^2));
E_yi45 = (rho_ib/(2*pi*e_0)).*((yf - posy(4,5))./((xf - posx(4,5)).^2 + (yf - posy(4,5)).^2));
E_yi46 = (rho_ic/(2*pi*e_0)).*((yf - posy(4,6))./((xf - posx(4,6)).^2 + (yf - posy(4,6)).^2));

E_yi51 = (-rho_ia/(2*pi*e_0)).*((yf - posy(5,1))./((xf - posx(5,1)).^2 + (yf - posy(5,1)).^2));
E_yi52 = (-rho_ib/(2*pi*e_0)).*((yf - posy(5,2))./((xf - posx(5,2)).^2 + (yf - posy(5,2)).^2));
E_yi53 = (-rho_ic/(2*pi*e_0)).*((yf - posy(5,3))./((xf - posx(5,3)).^2 + (yf - posy(5,3)).^2));
E_yi54 = (rho_ia/(2*pi*e_0)).*((yf - posy(5,4))./((xf - posx(5,4)).^2 + (yf - posy(5,4)).^2));
E_yi56 = (rho_ic/(2*pi*e_0)).*((yf - posy(5,6))./((xf - posx(5,6)).^2 + (yf - posy(5,6)).^2));

E_yi61 = (-rho_ia/(2*pi*e_0)).*((yf - posy(6,1))./((xf - posx(6,1)).^2 + (yf - posy(6,1)).^2));
E_yi62 = (-rho_ib/(2*pi*e_0)).*((yf - posy(6,2))./((xf - posx(6,2)).^2 + (yf - posy(6,2)).^2));
E_yi63 = (-rho_ic/(2*pi*e_0)).*((yf - posy(6,3))./((xf - posx(6,3)).^2 + (yf - posy(6,3)).^2));
E_yi64 = (rho_ia/(2*pi*e_0)).*((yf - posy(6,4))./((xf - posx(6,4)).^2 + (yf - posy(6,4)).^2));
E_yi65 = (rho_ib/(2*pi*e_0)).*((yf - posy(6,5))./((xf - posx(6,5)).^2 + (yf - posy(6,5)).^2));

%% aqui em cada matriz dessa tem as componentes geradas na dire��o x e y, pela parte real e imagin�ria de cada carga el�trica:

E_xr = [ E_xr12 ; E_xr13 ; E_xr14 ; E_xr15 ; E_xr16 ; E_xr21 ; E_xr23 ; E_xr24 ; E_xr25 ; E_xr26 ; E_xr31 ; E_xr32 ; E_xr34 ; E_xr35 ; E_xr36 ; E_xr41 ; E_xr42 ; E_xr43 ; E_xr45 ; E_xr46 ; E_xr51 ; E_xr52 ; E_xr53 ; E_xr54 ; E_xr56 ; E_xr61 ; E_xr62 ; E_xr63 ; E_xr64 ; E_xr65 ];

E_ximag = [ E_xi12 ; E_xi13 ; E_xi14 ; E_xi15 ; E_xi16 ; E_xi21 ; E_xi23 ; E_xi24 ; E_xi25 ; E_xi26 ; E_xi31 ; E_xi32 ; E_xi34 ; E_xi35 ; E_xi36 ; E_xi41 ; E_xi42 ; E_xi43 ; E_xi45 ; E_xi46 ; E_xi51 ; E_xi52 ; E_xi53 ; E_xi54 ; E_xi56 ; E_xi61 ; E_xi62 ; E_xi63 ; E_xi64 ; E_xi65 ];

E_yr = [ E_yr12 ; E_yr13 ; E_yr14 ; E_yr15 ; E_yr16 ; E_yr21 ; E_yr23 ; E_yr24 ; E_yr25 ; E_yr26 ; E_yr31 ; E_yr32 ; E_yr34 ; E_yr35 ; E_yr36 ; E_yr41 ; E_yr42 ; E_yr43 ; E_yr45 ; E_yr46 ; E_yr51 ; E_yr52 ; E_yr53 ; E_yr54 ; E_yr56 ; E_yr61 ; E_yr62 ; E_yr63 ; E_yr64 ; E_yr65 ];

E_yimag = [ E_yi12 ; E_yi13 ; E_yi14 ; E_yi15 ; E_yi16 ; E_yi21 ; E_yi23 ; E_yi24 ; E_yi25 ; E_yi26 ; E_yi31 ; E_yi32 ; E_yi34 ; E_yi35 ; E_yi36 ; E_yi41 ; E_yi42 ; E_yi43 ; E_yi45 ; E_yi46 ; E_yi51 ; E_yi52 ; E_yi53 ; E_yi54 ; E_yi56 ; E_yi61 ; E_yi62 ; E_yi63 ; E_yi64 ; E_yi65 ];


%%

Exr = sum(E_xr).^2;

Eximag = sum(E_ximag).^2;

Eyr = sum(E_yr).^2;

Eyimag = sum(E_yimag).^2;

%% aqui retira-se o m�dulo dos quadrados:

% V/m -> kV/cm = 10^-3/10^2
Erms = ((Exr+Eximag+Eyr+Eyimag).^(1/2))*(10^-5);
thetad = rad2deg(theta); %converte o �ngulo de radianos para graus

Emed = mean(Erms);
Emax = max(Erms);
Kirreg = Emax/Emed;

%% matriz P

P1 = zeros(n,n);
for i = 1:n
     for j = 1:n
        if(i==j)
           P1(i,j) = (1/(2*pi*e_0))*log((4*yr1(i))/(2*r));
        else
           P1(i,j) = (1/(2*pi*e_0))*log(sqrt((xr1(j)-xr1(i))^2+(yi1(j)-yr1(i))^2)/(sqrt((xr1(j)-xr1(i))^2+(yr1(j)-yr1(i))^2)));
        end
     end
end

%% C�lculo tens�o por fase

V1 = 145*10^3;

%tens�o condutor 1 fase a
V1_ra = V1/sqrt(3);
V1_ia = 0;

%tens�o condutor 2 fase b
V1_rb = V1*(cos(2*pi/3))/sqrt(3);
V1_ib = 1i*V1*(sin(2*pi/3))/sqrt(3);

%tens�o condutor 3 fase c
V1_rc = V1*(cos(-2*pi/3))/sqrt(3);
V1_ic = 1i*V1*(sin(-2*pi/3))/sqrt(3);

Vf1 = [ V1_ra+V1_ia ; V1_rb+V1_ib ; V1_rc+V1_ic ];

%% C�lculo densidade de carga

rho1 = P1\Vf1;

% busco cada uma das posicoes de rho real e imaginario por fase

% rho condutor 1 fase a
rho1_ra = real(rho1(1));
rho1_ia = imag(rho1(1));

% rho condutor 2 fase b
rho1_rb = real(rho1(2));
rho1_ib = imag(rho1(2));

% rho condutor 3 fase c
rho1_rc = real(rho1(3));
rho1_ic = imag(rho1(3));

%% Dist�ncia entre condutores

%sendo: i -> o n�mero do primeiro condutor e j -> o n�mero do segundo condutor
%exemplo: 1,2 � a dist�ncia entre condutor 1 e o condutor 2, os pr�ximos
%c�lculos seguem esta mesma l�gica
%sendo 4 a imagem do condutor 1, 5 a imagem do condutor 2 e 6 a imagem do
%condutor 3

distance1 = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j
            distance1(i,j) = 0; %como n�o existe dist�ncia no centro do pr�prio condutor ent�o � zero
        else
            distance1(i,j) = sqrt(((x1(j)-x1(i))^2)+((y1(j)-y1(i))^2)); %c�lculo dist�ncia entre dois pontos (entre centro de dois condutores)
        end
    end
end

%% C�lculo de delta

delta1 = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j
            delta1(i,j) = 0; %delta para o centro do pr�prio condutor � zero
        else
            delta1(i,j) = (r^2)./(distance1(i,j)); %c�lculo de delta utilizando a f�rmula e a dist�ncia calculada acima
        end
    end
end

%% C�lculo das posi��es em x das cargas imagens
 
phix1 = zeros(nc,nc);
posx1 = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j || x1(j)==x1(i) %c�lculo de posi��o para cargas posicionadas no centro dos condutores
            posx1(i,j) = x1(i);
        elseif y1(j)==y1(i) %c�lculo para cargas de mesma altura
            if x1(j) > x1(i)
                posx1(i,j) = x1(i) + delta1(i,j);
            elseif x1(i) > x1(j)
                posx1(i,j) = x1(i) - delta1(i,j);
            end
        elseif x1(i)~=x1(j) && y1(j) > y1(i)
            if x1(j) > x1(i)
                phix1(i,j) = asin((y1(j)-y1(i))/(distance1(i,j)));
                posx1(i,j) = x1(i) + delta1(i,j)*cos(phix1(i,j));
            elseif x1(i) > x1(j)
                phix1(i,j) = asin((y1(j)-y1(i))/(distance1(i,j)));
                posx1(i,j) = x1(i) - delta1(i,j)*cos(phix1(i,j));
            end
        elseif x1(i)~=x1(j) &&  y1(i) > y1(j)
            if x1(j) > x1(i)
                phix1(i,j) = acos((y1(i)-y1(j))/(distance1(i,j)));
                posx1(i,j) = x1(i) + delta1(i,j)*sin(phix1(i,j));
            elseif x1(i) > x1(j)
                phix1(i,j) = acos((y1(i)-y1(j))/(distance1(i,j)));
                posx1(i,j) = x1(i) - delta1(i,j)*sin(phix1(i,j));
            end
        end
    end
end

%% C�lculo das posi��es em y das cargas imagens

phiy1 = zeros(nc,nc);
posy1 = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        if i==j || y1(j)==y1(i)
            posy1(i,j) = y1(i);
        elseif x1(i)==x1(j)
            if y1(j) > y1(i) || y1(i) < 0 && y1(j) < 0
                posy1(i,j) = y1(i) + delta1(i,j);
            elseif y1(i) > y1(j)
                posy1(i,j) = y1(i) - delta1(i,j);
            end
        elseif x1(i)~=x1(j) && y1(j) > y1(i)
            if  y1(i) > 0 || y1(i) < 0 && y1(j) < 0
                phiy1(i,j) = asin((y1(j)-y1(i))/(distance1(i,j)));
                posy1(i,j) = y1(i) + delta1(i,j)*sin(phiy1(i,j));
            elseif y1(i) < 0 && y1(j) > 0
                phiy1(i,j) = asin((y1(j)-y1(i))/(distance1(i,j)));
                posy1(i,j) = y1(i) - delta1(i,j)*sin(phiy1(i,j));
            end
        elseif x1(i)~=x1(j) && y1(i) > y1(j)
            phiy1(i,j) = acos((y1(i)-y1(j))/(distance1(i,j)));
            posy1(i,j) = y1(i) - delta1(i,j)*cos(phiy1(i,j));
        end
    end
end


%% E1_xr componente x real campo el�trico condutor 2 fase b assim segue:

%busca o valor da carga imagem no condutor 1 real, valor do eixo x dos
%pontos de avalia��o e as posi��es geradas pelas matrizes posx e posy das
%cargas imagens

E1_xr12 = (-rho1_rb/(2*pi*e_0)).*((xf1 - posx1(1,2))./((xf1 - posx1(1,2)).^2 + (yf1 - posy1(1,2)).^2));
E1_xr13 = (-rho1_rc/(2*pi*e_0)).*((xf1 - posx1(1,3))./((xf1 - posx1(1,3)).^2 + (yf1 - posy1(1,3)).^2));
E1_xr14 = (rho1_ra/(2*pi*e_0)).*((xf1 - posx1(1,4))./((xf1 - posx1(1,4)).^2 + (yf1 - posy1(1,4)).^2));
E1_xr15 = (rho1_rb/(2*pi*e_0)).*((xf1 - posx1(1,5))./((xf1 - posx1(1,5)).^2 + (yf1 - posy1(1,5)).^2));
E1_xr16 = (rho1_rc/(2*pi*e_0)).*((xf1 - posx1(1,6))./((xf1 - posx1(1,6)).^2 + (yf1 - posy1(1,6)).^2));

E1_xr21 = (-rho1_ra/(2*pi*e_0)).*((xf1 - posx1(2,1))./((xf1 - posx1(2,1)).^2 + (yf1 - posy1(2,1)).^2));
E1_xr23 = (-rho1_rc/(2*pi*e_0)).*((xf1 - posx1(2,3))./((xf1 - posx1(2,3)).^2 + (yf1 - posy1(2,3)).^2));
E1_xr24 = (rho1_ra/(2*pi*e_0)).*((xf1 - posx1(2,4))./((xf1 - posx1(2,4)).^2 + (yf1 - posy1(2,4)).^2));
E1_xr25 = (rho1_rb/(2*pi*e_0)).*((xf1 - posx1(2,5))./((xf1 - posx1(2,5)).^2 + (yf1 - posy1(2,5)).^2));
E1_xr26 = (rho1_rc/(2*pi*e_0)).*((xf1 - posx1(2,6))./((xf1 - posx1(2,6)).^2 + (yf1 - posy1(2,6)).^2));

E1_xr31 = (-rho1_ra/(2*pi*e_0)).*((xf1 - posx1(3,1))./((xf1 - posx1(3,1)).^2 + (yf1 - posy1(3,1)).^2));
E1_xr32 = (-rho1_rb/(2*pi*e_0)).*((xf1 - posx1(3,2))./((xf1 - posx1(3,2)).^2 + (yf1 - posy1(3,2)).^2));
E1_xr34 = (rho1_ra/(2*pi*e_0)).*((xf1 - posx1(3,4))./((xf1 - posx1(3,4)).^2 + (yf1 - posy1(3,4)).^2));
E1_xr35 = (rho1_rb/(2*pi*e_0)).*((xf1 - posx1(3,5))./((xf1 - posx1(3,5)).^2 + (yf1 - posy1(3,5)).^2));
E1_xr36 = (rho1_rc/(2*pi*e_0)).*((xf1 - posx1(3,6))./((xf1 - posx1(3,6)).^2 + (yf1 - posy1(3,6)).^2));

E1_xr41 = (-rho1_ra/(2*pi*e_0)).*((xf1 - posx1(4,1))./((xf1 - posx1(4,1)).^2 + (yf1 - posy1(4,1)).^2));
E1_xr42 = (-rho1_rb/(2*pi*e_0)).*((xf1 - posx1(4,2))./((xf1 - posx1(4,2)).^2 + (yf1 - posy1(4,2)).^2));
E1_xr43 = (-rho1_rc/(2*pi*e_0)).*((xf1 - posx1(4,3))./((xf1 - posx1(4,3)).^2 + (yf1 - posy1(4,3)).^2));
E1_xr45 = (rho1_rb/(2*pi*e_0)).*((xf1 - posx1(4,5))./((xf1 - posx1(4,5)).^2 + (yf1 - posy1(4,5)).^2));
E1_xr46 = (rho1_rc/(2*pi*e_0)).*((xf1 - posx1(4,6))./((xf1 - posx1(4,6)).^2 + (yf1 - posy1(4,6)).^2));

E1_xr51 = (-rho1_ra/(2*pi*e_0)).*((xf1 - posx1(5,1))./((xf1 - posx1(5,1)).^2 + (yf1 - posy1(5,1)).^2));
E1_xr52 = (-rho1_rb/(2*pi*e_0)).*((xf1 - posx1(5,2))./((xf1 - posx1(5,2)).^2 + (yf1 - posy1(5,2)).^2));
E1_xr53 = (-rho1_rc/(2*pi*e_0)).*((xf1 - posx1(5,3))./((xf1 - posx1(5,3)).^2 + (yf1 - posy1(5,3)).^2));
E1_xr54 = (rho1_ra/(2*pi*e_0)).*((xf1 - posx1(5,4))./((xf1 - posx1(5,4)).^2 + (yf1 - posy1(5,4)).^2));
E1_xr56 = (rho1_rc/(2*pi*e_0)).*((xf1 - posx1(5,6))./((xf1 - posx1(5,6)).^2 + (yf1 - posy1(5,6)).^2));

E1_xr61 = (-rho1_ra/(2*pi*e_0)).*((xf1 - posx1(6,1))./((xf1 - posx1(6,1)).^2 + (yf1 - posy1(6,1)).^2));
E1_xr62 = (-rho1_rb/(2*pi*e_0)).*((xf1 - posx1(6,2))./((xf1 - posx1(6,2)).^2 + (yf1 - posy1(6,2)).^2));
E1_xr63 = (-rho1_rc/(2*pi*e_0)).*((xf1 - posx1(6,3))./((xf1 - posx1(6,3)).^2 + (yf1 - posy1(6,3)).^2));
E1_xr64 = (rho1_ra/(2*pi*e_0)).*((xf1 - posx1(6,4))./((xf1 - posx1(6,4)).^2 + (yf1 - posy1(6,4)).^2));
E1_xr65 = (rho1_rb/(2*pi*e_0)).*((xf1 - posx1(6,5))./((xf1 - posx1(6,5)).^2 + (yf1 - posy1(6,5)).^2));


%% E1_xi componente x imaginario campo el�trico condutor 2 fase b assim segue:

E1_xi12 = (-rho1_ib/(2*pi*e_0)).*((xf1 - posx1(1,2))./((xf1 - posx1(1,2)).^2 + (yf1 - posy1(1,2)).^2));
E1_xi13 = (-rho1_ic/(2*pi*e_0)).*((xf1 - posx1(1,3))./((xf1 - posx1(1,3)).^2 + (yf1 - posy1(1,3)).^2));
E1_xi14 = (rho1_ia/(2*pi*e_0)).*((xf1 - posx1(1,4))./((xf1 - posx1(1,4)).^2 + (yf1 - posy1(1,4)).^2));
E1_xi15 = (rho1_ib/(2*pi*e_0)).*((xf1 - posx1(1,5))./((xf1 - posx1(1,5)).^2 + (yf1 - posy1(1,5)).^2));
E1_xi16 = (rho1_ic/(2*pi*e_0)).*((xf1 - posx1(1,6))./((xf1 - posx1(1,6)).^2 + (yf1 - posy1(1,6)).^2));

E1_xi21 = (-rho1_ia/(2*pi*e_0)).*((xf1 - posx1(2,1))./((xf1 - posx1(2,1)).^2 + (yf1 - posy1(2,1)).^2));
E1_xi23 = (-rho1_ic/(2*pi*e_0)).*((xf1 - posx1(2,3))./((xf1 - posx1(2,3)).^2 + (yf1 - posy1(2,3)).^2));
E1_xi24 = (rho1_ia/(2*pi*e_0)).*((xf1 - posx1(2,4))./((xf1 - posx1(2,4)).^2 + (yf1 - posy1(2,4)).^2));
E1_xi25 = (rho1_ib/(2*pi*e_0)).*((xf1 - posx1(2,5))./((xf1 - posx1(2,5)).^2 + (yf1 - posy1(2,5)).^2));
E1_xi26 = (rho1_ic/(2*pi*e_0)).*((xf1 - posx1(2,6))./((xf1 - posx1(2,6)).^2 + (yf1 - posy1(2,6)).^2));

E1_xi31 = (-rho1_ia/(2*pi*e_0)).*((xf1 - posx1(3,1))./((xf1 - posx1(3,1)).^2 + (yf1 - posy1(3,1)).^2));
E1_xi32 = (-rho1_ib/(2*pi*e_0)).*((xf1 - posx1(3,2))./((xf1 - posx1(3,2)).^2 + (yf1 - posy1(3,2)).^2));
E1_xi34 = (rho1_ia/(2*pi*e_0)).*((xf1 - posx1(3,4))./((xf1 - posx1(3,4)).^2 + (yf1 - posy1(3,4)).^2));
E1_xi35 = (rho1_ib/(2*pi*e_0)).*((xf1 - posx1(3,5))./((xf1 - posx1(3,5)).^2 + (yf1 - posy1(3,5)).^2));
E1_xi36 = (rho1_ic/(2*pi*e_0)).*((xf1 - posx1(3,6))./((xf1 - posx1(3,6)).^2 + (yf1 - posy1(3,6)).^2));

E1_xi41 = (-rho1_ia/(2*pi*e_0)).*((xf1 - posx1(4,1))./((xf1 - posx1(4,1)).^2 + (yf1 - posy1(4,1)).^2));
E1_xi42 = (-rho1_ib/(2*pi*e_0)).*((xf1 - posx1(4,2))./((xf1 - posx1(4,2)).^2 + (yf1 - posy1(4,2)).^2));
E1_xi43 = (-rho1_ic/(2*pi*e_0)).*((xf1 - posx1(4,3))./((xf1 - posx1(4,3)).^2 + (yf1 - posy1(4,3)).^2));
E1_xi45 = (rho1_ib/(2*pi*e_0)).*((xf1 - posx1(4,5))./((xf1 - posx1(4,5)).^2 + (yf1 - posy1(4,5)).^2));
E1_xi46 = (rho1_ic/(2*pi*e_0)).*((xf1 - posx1(4,6))./((xf1 - posx1(4,6)).^2 + (yf1 - posy1(4,6)).^2));

E1_xi51 = (-rho1_ia/(2*pi*e_0)).*((xf1 - posx1(5,1))./((xf1 - posx1(5,1)).^2 + (yf1 - posy1(5,1)).^2));
E1_xi52 = (-rho1_ib/(2*pi*e_0)).*((xf1 - posx1(5,2))./((xf1 - posx1(5,2)).^2 + (yf1 - posy1(5,2)).^2));
E1_xi53 = (-rho1_ic/(2*pi*e_0)).*((xf1 - posx1(5,3))./((xf1 - posx1(5,3)).^2 + (yf1 - posy1(5,3)).^2));
E1_xi54 = (rho1_ia/(2*pi*e_0)).*((xf1 - posx1(5,4))./((xf1 - posx1(5,4)).^2 + (yf1 - posy1(5,4)).^2));
E1_xi56 = (rho1_ic/(2*pi*e_0)).*((xf1 - posx1(5,6))./((xf1 - posx1(5,6)).^2 + (yf1 - posy1(5,6)).^2));

E1_xi61 = (-rho1_ia/(2*pi*e_0)).*((xf1 - posx1(6,1))./((xf1 - posx1(6,1)).^2 + (yf1 - posy1(6,1)).^2));
E1_xi62 = (-rho1_ib/(2*pi*e_0)).*((xf1 - posx1(6,2))./((xf1 - posx1(6,2)).^2 + (yf1 - posy1(6,2)).^2));
E1_xi63 = (-rho1_ic/(2*pi*e_0)).*((xf1 - posx1(6,3))./((xf1 - posx1(6,3)).^2 + (yf1 - posy1(6,3)).^2));
E1_xi64 = (rho1_ia/(2*pi*e_0)).*((xf1 - posx1(6,4))./((xf1 - posx1(6,4)).^2 + (yf1 - posy1(6,4)).^2));
E1_xi65 = (rho1_ib/(2*pi*e_0)).*((xf1 - posx1(6,5))./((xf1 - posx1(6,5)).^2 + (yf1 - posy1(6,5)).^2));

%% E1_yr componente y real campo el�trico condutor 2 fase b

E1_yr12 = (-rho1_rb/(2*pi*e_0)).*((yf1 - posy1(1,2))./((xf1 - posx1(1,2)).^2 + (yf1 - posy1(1,2)).^2));
E1_yr13 = (-rho1_rc/(2*pi*e_0)).*((yf1 - posy1(1,3))./((xf1 - posx1(1,3)).^2 + (yf1 - posy1(1,3)).^2));
E1_yr14 = (rho1_ra/(2*pi*e_0)).*((yf1 - posy1(1,4))./((xf1 - posx1(1,4)).^2 + (yf1 - posy1(1,4)).^2));
E1_yr15 = (rho1_rb/(2*pi*e_0)).*((yf1 - posy1(1,5))./((xf1 - posx1(1,5)).^2 + (yf1 - posy1(1,5)).^2));
E1_yr16 = (rho1_rc/(2*pi*e_0)).*((yf1 - posy1(1,6))./((xf1 - posx1(1,6)).^2 + (yf1 - posy1(1,6)).^2));

E1_yr21 = (-rho1_ra/(2*pi*e_0)).*((yf1 - posy1(2,1))./((xf1 - posx1(2,1)).^2 + (yf1 - posy1(2,1)).^2));
E1_yr23 = (-rho1_rc/(2*pi*e_0)).*((yf1 - posy1(2,3))./((xf1 - posx1(2,3)).^2 + (yf1 - posy1(2,3)).^2));
E1_yr24 = (rho1_ra/(2*pi*e_0)).*((yf1 - posy1(2,4))./((xf1 - posx1(2,4)).^2 + (yf1 - posy1(2,4)).^2));
E1_yr25 = (rho1_rb/(2*pi*e_0)).*((yf1 - posy1(2,5))./((xf1 - posx1(2,5)).^2 + (yf1 - posy1(2,5)).^2));
E1_yr26 = (rho1_rc/(2*pi*e_0)).*((yf1 - posy1(2,6))./((xf1 - posx1(2,6)).^2 + (yf1 - posy1(2,6)).^2));

E1_yr31 = (-rho1_ra/(2*pi*e_0)).*((yf1 - posy1(3,1))./((xf1 - posx1(3,1)).^2 + (yf1 - posy1(3,1)).^2));
E1_yr32 = (-rho1_rb/(2*pi*e_0)).*((yf1 - posy1(3,2))./((xf1 - posx1(3,2)).^2 + (yf1 - posy1(3,2)).^2));
E1_yr34 = (rho1_ra/(2*pi*e_0)).*((yf1 - posy1(3,4))./((xf1 - posx1(3,4)).^2 + (yf1 - posy1(3,4)).^2));
E1_yr35 = (rho1_rb/(2*pi*e_0)).*((yf1 - posy1(3,5))./((xf1 - posx1(3,5)).^2 + (yf1 - posy1(3,5)).^2));
E1_yr36 = (rho1_rc/(2*pi*e_0)).*((yf1 - posy1(3,6))./((xf1 - posx1(3,6)).^2 + (yf1 - posy1(3,6)).^2));

E1_yr41 = (-rho1_ra/(2*pi*e_0)).*((yf1 - posy1(4,1))./((xf1 - posx1(4,1)).^2 + (yf1 - posy1(4,1)).^2));
E1_yr42 = (-rho1_rb/(2*pi*e_0)).*((yf1 - posy1(4,2))./((xf1 - posx1(4,2)).^2 + (yf1 - posy1(4,2)).^2));
E1_yr43 = (-rho1_rc/(2*pi*e_0)).*((yf1 - posy1(4,3))./((xf1 - posx1(4,3)).^2 + (yf1 - posy1(4,3)).^2));
E1_yr45 = (rho1_rb/(2*pi*e_0)).*((yf1 - posy1(4,5))./((xf1 - posx1(4,5)).^2 + (yf1 - posy1(4,5)).^2));
E1_yr46 = (rho1_rc/(2*pi*e_0)).*((yf1 - posy1(4,6))./((xf1 - posx1(4,6)).^2 + (yf1 - posy1(4,6)).^2));

E1_yr51 = (-rho1_ra/(2*pi*e_0)).*((yf1 - posy1(5,1))./((xf1 - posx1(5,1)).^2 + (yf1 - posy1(5,1)).^2));
E1_yr52 = (-rho1_rb/(2*pi*e_0)).*((yf1 - posy1(5,2))./((xf1 - posx1(5,2)).^2 + (yf1 - posy1(5,2)).^2));
E1_yr53 = (-rho1_rc/(2*pi*e_0)).*((yf1 - posy1(5,3))./((xf1 - posx1(5,3)).^2 + (yf1 - posy1(5,3)).^2));
E1_yr54 = (rho1_ra/(2*pi*e_0)).*((yf1 - posy1(5,4))./((xf1 - posx1(5,4)).^2 + (yf1 - posy1(5,4)).^2));
E1_yr56 = (rho1_rc/(2*pi*e_0)).*((yf1 - posy1(5,6))./((xf1 - posx1(5,6)).^2 + (yf1 - posy1(5,6)).^2));

E1_yr61 = (-rho1_ra/(2*pi*e_0)).*((yf1 - posy1(6,1))./((xf1 - posx1(6,1)).^2 + (yf1 - posy1(6,1)).^2));
E1_yr62 = (-rho1_rb/(2*pi*e_0)).*((yf1 - posy1(6,2))./((xf1 - posx1(6,2)).^2 + (yf1 - posy1(6,2)).^2));
E1_yr63 = (-rho1_rc/(2*pi*e_0)).*((yf1 - posy1(6,3))./((xf1 - posx1(6,3)).^2 + (yf1 - posy1(6,3)).^2));
E1_yr64 = (rho1_ra/(2*pi*e_0)).*((yf1 - posy1(6,4))./((xf1 - posx1(6,4)).^2 + (yf1 - posy1(6,4)).^2));
E1_yr65 = (rho1_rb/(2*pi*e_0)).*((yf1 - posy1(6,5))./((xf1 - posx1(6,5)).^2 + (yf1 - posy1(6,5)).^2));


%% E1_yi21 componente y imaginario campo el�trico condutor 2 fase b

E1_yi12 = (-rho1_ib/(2*pi*e_0)).*((yf1 - posy1(1,2))./((xf1 - posx1(1,2)).^2 + (yf1 - posy1(1,2)).^2));
E1_yi13 = (-rho1_ic/(2*pi*e_0)).*((yf1 - posy1(1,3))./((xf1 - posx1(1,3)).^2 + (yf1 - posy1(1,3)).^2));
E1_yi14 = (rho1_ia/(2*pi*e_0)).*((yf1 - posy1(1,4))./((xf1 - posx1(1,4)).^2 + (yf1 - posy1(1,4)).^2));
E1_yi15 = (rho1_ib/(2*pi*e_0)).*((yf1 - posy1(1,5))./((xf1 - posx1(1,5)).^2 + (yf1 - posy1(1,5)).^2));
E1_yi16 = (rho1_ic/(2*pi*e_0)).*((yf1 - posy1(1,6))./((xf1 - posx1(1,6)).^2 + (yf1 - posy1(1,6)).^2));

E1_yi21 = (-rho1_ia/(2*pi*e_0)).*((yf1 - posy1(2,1))./((xf1 - posx1(2,1)).^2 + (yf1 - posy1(2,1)).^2));
E1_yi23 = (-rho1_ic/(2*pi*e_0)).*((yf1 - posy1(2,3))./((xf1 - posx1(2,3)).^2 + (yf1 - posy1(2,3)).^2));
E1_yi24 = (rho1_ia/(2*pi*e_0)).*((yf1 - posy1(2,4))./((xf1 - posx1(2,4)).^2 + (yf1 - posy1(2,4)).^2));
E1_yi25 = (rho1_ib/(2*pi*e_0)).*((yf1 - posy1(2,5))./((xf1 - posx1(2,5)).^2 + (yf1 - posy1(2,5)).^2));
E1_yi26 = (rho1_ic/(2*pi*e_0)).*((yf1 - posy1(2,6))./((xf1 - posx1(2,6)).^2 + (yf1 - posy1(2,6)).^2));

E1_yi31 = (-rho1_ia/(2*pi*e_0)).*((yf1 - posy1(3,1))./((xf1 - posx1(3,1)).^2 + (yf1 - posy1(3,1)).^2));
E1_yi32 = (-rho1_ib/(2*pi*e_0)).*((yf1 - posy1(3,2))./((xf1 - posx1(3,2)).^2 + (yf1 - posy1(3,2)).^2));
E1_yi34 = (rho1_ia/(2*pi*e_0)).*((yf1 - posy1(3,4))./((xf1 - posx1(3,4)).^2 + (yf1 - posy1(3,4)).^2));
E1_yi35 = (rho1_ib/(2*pi*e_0)).*((yf1 - posy1(3,5))./((xf1 - posx1(3,5)).^2 + (yf1 - posy1(3,5)).^2));
E1_yi36 = (rho1_ic/(2*pi*e_0)).*((yf1 - posy1(3,6))./((xf1 - posx1(3,6)).^2 + (yf1 - posy1(3,6)).^2));

E1_yi41 = (-rho1_ia/(2*pi*e_0)).*((yf1 - posy1(4,1))./((xf1 - posx1(4,1)).^2 + (yf1 - posy1(4,1)).^2));
E1_yi42 = (-rho1_ib/(2*pi*e_0)).*((yf1 - posy1(4,2))./((xf1 - posx1(4,2)).^2 + (yf1 - posy1(4,2)).^2));
E1_yi43 = (-rho1_ic/(2*pi*e_0)).*((yf1 - posy1(4,3))./((xf1 - posx1(4,3)).^2 + (yf1 - posy1(4,3)).^2));
E1_yi45 = (rho1_ib/(2*pi*e_0)).*((yf1 - posy1(4,5))./((xf1 - posx1(4,5)).^2 + (yf1 - posy1(4,5)).^2));
E1_yi46 = (rho1_ic/(2*pi*e_0)).*((yf1 - posy1(4,6))./((xf1 - posx1(4,6)).^2 + (yf1 - posy1(4,6)).^2));

E1_yi51 = (-rho1_ia/(2*pi*e_0)).*((yf1 - posy1(5,1))./((xf1 - posx1(5,1)).^2 + (yf1 - posy1(5,1)).^2));
E1_yi52 = (-rho1_ib/(2*pi*e_0)).*((yf1 - posy1(5,2))./((xf1 - posx1(5,2)).^2 + (yf1 - posy1(5,2)).^2));
E1_yi53 = (-rho1_ic/(2*pi*e_0)).*((yf1 - posy1(5,3))./((xf1 - posx1(5,3)).^2 + (yf1 - posy1(5,3)).^2));
E1_yi54 = (rho1_ia/(2*pi*e_0)).*((yf1 - posy1(5,4))./((xf1 - posx1(5,4)).^2 + (yf1 - posy1(5,4)).^2));
E1_yi56 = (rho1_ic/(2*pi*e_0)).*((yf1 - posy1(5,6))./((xf1 - posx1(5,6)).^2 + (yf1 - posy1(5,6)).^2));

E1_yi61 = (-rho1_ia/(2*pi*e_0)).*((yf1 - posy1(6,1))./((xf1 - posx1(6,1)).^2 + (yf1 - posy1(6,1)).^2));
E1_yi62 = (-rho1_ib/(2*pi*e_0)).*((yf1 - posy1(6,2))./((xf1 - posx1(6,2)).^2 + (yf1 - posy1(6,2)).^2));
E1_yi63 = (-rho1_ic/(2*pi*e_0)).*((yf1 - posy1(6,3))./((xf1 - posx1(6,3)).^2 + (yf1 - posy1(6,3)).^2));
E1_yi64 = (rho1_ia/(2*pi*e_0)).*((yf1 - posy1(6,4))./((xf1 - posx1(6,4)).^2 + (yf1 - posy1(6,4)).^2));
E1_yi65 = (rho1_ib/(2*pi*e_0)).*((yf1 - posy1(6,5))./((xf1 - posx1(6,5)).^2 + (yf1 - posy1(6,5)).^2));

%% aqui em cada matriz dessa tem as componentes geradas na dire��o x e y, pela parte real e imagin�ria de cada carga el�trica:

E1_xr = [ E1_xr12 ; E1_xr13 ; E1_xr14 ; E1_xr15 ; E1_xr16 ; E1_xr21 ; E1_xr23 ; E1_xr24 ; E1_xr25 ; E1_xr26 ; E1_xr31 ; E1_xr32 ; E1_xr34 ; E1_xr35 ; E1_xr36 ; E1_xr41 ; E1_xr42 ; E1_xr43 ; E1_xr45 ; E1_xr46 ; E1_xr51 ; E1_xr52 ; E1_xr53 ; E1_xr54 ; E1_xr56 ; E1_xr61 ; E1_xr62 ; E1_xr63 ; E1_xr64 ; E1_xr65 ];

E1_ximag = [ E1_xi12 ; E1_xi13 ; E1_xi14 ; E1_xi15 ; E1_xi16 ; E1_xi21 ; E1_xi23 ; E1_xi24 ; E1_xi25 ; E1_xi26 ; E1_xi31 ; E1_xi32 ; E1_xi34 ; E1_xi35 ; E1_xi36 ; E1_xi41 ; E1_xi42 ; E1_xi43 ; E1_xi45 ; E1_xi46 ; E1_xi51 ; E1_xi52 ; E1_xi53 ; E1_xi54 ; E1_xi56 ; E1_xi61 ; E1_xi62 ; E1_xi63 ; E1_xi64 ; E1_xi65 ];

E1_yr = [ E1_yr12 ; E1_yr13 ; E1_yr14 ; E1_yr15 ; E1_yr16 ; E1_yr21 ; E1_yr23 ; E1_yr24 ; E1_yr25 ; E1_yr26 ; E1_yr31 ; E1_yr32 ; E1_yr34 ; E1_yr35 ; E1_yr36 ; E1_yr41 ; E1_yr42 ; E1_yr43 ; E1_yr45 ; E1_yr46 ; E1_yr51 ; E1_yr52 ; E1_yr53 ; E1_yr54 ; E1_yr56 ; E1_yr61 ; E1_yr62 ; E1_yr63 ; E1_yr64 ; E1_yr65 ];

E1_yimag = [ E1_yi12 ; E1_yi13 ; E1_yi14 ; E1_yi15 ; E1_yi16 ; E1_yi21 ; E1_yi23 ; E1_yi24 ; E1_yi25 ; E1_yi26 ; E1_yi31 ; E1_yi32 ; E1_yi34 ; E1_yi35 ; E1_yi36 ; E1_yi41 ; E1_yi42 ; E1_yi43 ; E1_yi45 ; E1_yi46 ; E1_yi51 ; E1_yi52 ; E1_yi53 ; E1_yi54 ; E1_yi56 ; E1_yi61 ; E1_yi62 ; E1_yi63 ; E1_yi64 ; E1_yi65 ];


%%

E1xr = sum(E1_xr).^2;

E1ximag = sum(E1_ximag).^2;

E1yr = sum(E1_yr).^2;

E1yimag = sum(E1_yimag).^2;

%% aqui retira-se o m�dulo dos quadrados:

% V/m -> kV/cm = 10^-3/10^2
E1rms = ((E1xr+E1ximag+E1yr+E1yimag).^(1/2))*(10^-5);

fs = 0.82; %fator de superf�cie (constante)
pa = 760; %press�o atmosf�rica [mmHg]
Ecrit = fcn_supcrit(r*10^(2));
%Ecrit = 18.11*fs*pa*(1+(0.54187/sqrt(r*pa)))*(10^-3);
E1_crit = linspace(Ecrit,Ecrit,360);

Emed1 = mean(E1rms);
Emax1 = max(E1rms);
Kirreg1 = Emax/Emed;

%%

E1rms1 = [15.0083309467484 15.0083101486078 15.0082673830908 15.0082026636604 15.0081160107551 15.0080074517881 15.0078770211326 15.0077247601120 15.0075507169857 15.0073549469312 15.0071375120273 15.0068984812280 15.0066379303438 15.0063559420099 15.0060526056625 15.0057280175035 15.0053822804683 15.0050155041901 15.0046278049642 15.0042193056991 15.0037901358863 15.0033404315413 15.0028703351681 15.0023799957016 15.0018695684622 15.0013392150952 15.0007891035176 15.0002194078629 14.9996303084157 14.9990219915514 14.9983946496743 14.9977484811464 14.9970836902192 14.9964004869691 14.9956990872184 14.9949797124634 14.9942425898016 14.9934879518449 14.9927160366519 14.9919270876355 14.9911213534864 14.9902990880842 14.9894605504149 14.9886060044752 14.9877357191980 14.9868499683411 14.9859490304074 14.9850331885525 14.9841027304760 14.9831579483406 14.9821991386583 14.9812266022080 14.9802406439140 14.9792415727641 14.9782297016914 14.9772053474800 14.9761688306520 14.9751204753631 14.9740606093031 14.9729895635719 14.9719076725843 14.9708152739514 14.9697127083788 14.9686003195372 14.9674784539736 14.9663474609771 14.9652076924779 14.9640595029233 14.9629032491756 14.9617392903833 14.9605679878709 14.9593897050200 14.9582048071632 14.9570136614491 14.9558166367339 14.9546141034689 14.9534064335743 14.9521940003258 14.9509771782283 14.9497563429108 14.9485318709943 14.9473041399839 14.9460735281447 14.9448404143840 14.9436051781320 14.9423681992266 14.9411298577952 14.9398905341350 14.9386506085933 14.9374104614584 14.9361704728322 14.9349310225223 14.9336924899175 14.9324552538834 14.9312196926348 14.9299861836264 14.9287551034473 14.9275268276899 14.9263017308460 14.9250801862002 14.9238625657066 14.9226492398857 14.9214405777108 14.9202369464990 14.9190387118012 14.9178462372977 14.9166598846854 14.9154800135781 14.9143069813964 14.9131411432623 14.9119828519037 14.9108324575359 14.9096903077806 14.9085567475458 14.9074321189440 14.9063167611755 14.9052110104476 14.9041151998663 14.9030296593496 14.9019547155294 14.9008906916577 14.8998379075152 14.8987966793230 14.8977673196543 14.8967501373420 14.8957454373975 14.8947535209201 14.8937746850163 14.8928092227124 14.8918574228801 14.8909195701538 14.8899959448456 14.8890868228794 14.8881924757019 14.8873131702210 14.8864491687188 14.8856007287942 14.8847681032817 14.8839515401852 14.8831512826168 14.8823675687192 14.8816006316101 14.8808506993204 14.8801179947235 14.8794027354819 14.8787051339879 14.8780253973102 14.8773637271272 14.8767203196861 14.8760953657405 14.8754890505074 14.8749015536132 14.8743330490454 14.8737837051101 14.8732536843865 14.8727431436785 14.8722522339837 14.8717811004438 14.8713298823117 14.8708987129131 14.8704877196123 14.8700970237789 14.8697267407536 14.8693769798237 14.8690478441880 14.8687394309328 14.8684518310102 14.8681851292064 14.8679394041261 14.8677147281716 14.8675111675219 14.8673287821143 14.8671676256342 14.8670277454962 14.8669091828338 14.8668119724901 14.8667361430062 14.8666817166162 14.8666487092408 14.8666371304816 14.8666469836234 14.8666782656288 14.8667309671396 14.8668050724834 14.8669005596734 14.8670174004163 14.8671555601201 14.8673149979017 14.8674956665994 14.8676975127853 14.8679204767764 14.8681644926558 14.8684294882852 14.8687153853244 14.8690220992547 14.8693495394013 14.8696976089508 14.8700662049846 14.8704552185031 14.8708645344507 14.8712940317530 14.8717435833430 14.8722130561977 14.8727023113719 14.8732112040357 14.8737395835136 14.8742872933189 14.8748541712063 14.8754400492032 14.8760447536605 14.8766681052999 14.8773099192576 14.8779700051343 14.8786481670510 14.8793442036938 14.8800579083770 14.8807890690867 14.8815374685480 14.8823028842809 14.8830850886600 14.8838838489702 14.8846989274786 14.8855300814977 14.8863770634422 14.8872396209085 14.8881174967365 14.8890104290780 14.8899181514782 14.8908403929354 14.8917768779887 14.8927273267805 14.8936914551477 14.8946689746876 14.8956595928450 14.8966630129953 14.8976789345206 14.8987070528976 14.8997470597800 14.9007986430893 14.9018614870972 14.9029352725205 14.9040196766069 14.9051143732276 14.9062190329678 14.9073333232217 14.9084569082915 14.9095894494746 14.9107306051729 14.9118800309743 14.9130373797698 14.9142023018363 14.9153744449569 14.9165534545030 14.9177389735527 14.9189306429873 14.9201281015971 14.9213309861921 14.9225389317013 14.9237515712870 14.9249685364488 14.9261894571344 14.9274139618484 14.9286416777665 14.9298722308386 14.9311052459128 14.9323403468411 14.9335771565849 14.9348152973447 14.9360543906666 14.9372940575515 14.9385339185824 14.9397735940289 14.9410127039717 14.9422508684105 14.9434877073910 14.9447228411137 14.9459558900541 14.9471864750815 14.9484142175770 14.9496387395491 14.9508596637513 14.9520766138012 14.9532892143029 14.9544970909552 14.9556998706829 14.9568971817418 14.9580886538438 14.9592739182707 14.9604526079988 14.9616243578118 14.9627888044079 14.9639455865344 14.9650943450929 14.9662347232579 14.9673663665857 14.9684889231430 14.9696020436078 14.9707053813927 14.9717985927449 14.9728813368798 14.9739532760690 14.9750140757681 14.9760634047180 14.9771009350627 14.9781263424431 14.9791393061200 14.9801395090729 14.9811266381030 14.9821003839449 14.9830604413604 14.9840065092548 14.9849382907582 14.9858554933454 14.9867578289176 14.9876450139141 14.9885167693926 14.9893728211366 14.9902128997463 14.9910367407169 14.9918440845487 14.9926346768207 14.9934082682838 14.9941646149432 14.9949034781448 14.9956246246510 14.9963278267312 14.9970128622280 14.9976795146427 14.9983275732069 14.9989568329520 14.9995670947875 15.0001581655635 15.0007298581377 15.0012819914449 15.0018143905570 15.0023268867420 15.0028193175280 15.0032915267561 15.0037433646326 15.0041746877873 15.0045853593181 15.0049752488461 15.0053442325497 15.0056921932221 15.0060190202993 15.0063246099080 15.0066088648991 15.0068716948813 15.0071130162526 15.0073327522329 15.0075308328863 15.0077071951514 15.0078617828578 15.0079945467505 15.0081054445054 15.0081944407449 15.0082615070503 15.0083066219704 15.0083297710339 15.0083309467484];
E1rms2 = [16.1618417473238 16.1618496985500 16.1618568432931 16.1618631804730 16.1618687099941 16.1618734327436 16.1618773505877 16.1618804663672 16.1618827838915 16.1618843079308 16.1618850442104 16.1618849993979 16.1618841810965 16.1618825978293 16.1618802590296 16.1618771750267 16.1618733570295 16.1618688171111 16.1618635681934 16.1618576240253 16.1618509991687 16.1618437089724 16.1618357695574 16.1618271977886 16.1618180112589 16.1618082282596 16.1617978677565 16.1617869493689 16.1617754933388 16.1617635205052 16.1617510522804 16.1617381106153 16.1617247179719 16.1617108972997 16.1616966719973 16.1616820658881 16.1616671031864 16.1616518084632 16.1616362066227 16.1616203228586 16.1616041826304 16.1615878116257 16.1615712357286 16.1615544809829 16.1615375735697 16.1615205397558 16.1615034058720 16.1614861982821 16.1614689433345 16.1614516673438 16.1614343965451 16.1614171570723 16.1613999749093 16.1613828758719 16.1613658855641 16.1613490293525 16.1613323323270 16.1613158192719 16.1612995146409 16.1612834425102 16.1612676265628 16.1612520900501 16.1612368557699 16.1612219460207 16.1612073825977 16.1611931867432 16.1611793791323 16.1611659798372 16.1611530083157 16.1611404833707 16.1611284231351 16.1611168450442 16.1611057658255 16.1610952014585 16.1610851671647 16.1610756773928 16.1610667457929 16.1610583852025 16.1610506076242 16.1610434242221 16.1610368452949 16.1610308802734 16.1610255377026 16.1610208252318 16.1610167496035 16.1610133166483 16.1610105312770 16.1610083974730 16.1610069182840 16.1610060958277 16.1610059312774 16.1610064248701 16.1610075758971 16.1610093827169 16.1610118427447 16.1610149524607 16.1610187074246 16.1610231022625 16.1610281306828 16.1610337854917 16.1610400585888 16.1610469409869 16.1610544228192 16.1610624933527 16.1610711410010 16.1610803533426 16.1610901171302 16.1611004183150 16.1611112420602 16.1611225727581 16.1611343940582 16.1611466888719 16.1611594394162 16.1611726272145 16.1611862331397 16.1612002374198 16.1612146196796 16.1612293589537 16.1612444337238 16.1612598219405 16.1612755010495 16.1612914480233 16.1613076393905 16.1613240512677 16.1613406593841 16.1613574391211 16.1613743655349 16.1613914133948 16.1614085572088 16.1614257712647 16.1614430296589 16.1614603063260 16.1614775750813 16.1614948096422 16.1615119836759 16.1615290708178 16.1615460447207 16.1615628790791 16.1615795476633 16.1615960243621 16.1616122832034 16.1616282983953 16.1616440443649 16.1616594957793 16.1616746275841 16.1616894150382 16.1617038337477 16.1617178596865 16.1617314692427 16.1617446392356 16.1617573469575 16.1617695701971 16.1617812872695 16.1617924770467 16.1618031189849 16.1618131931475 16.1618226802419 16.1618315616343 16.1618398193801 16.1618474362481 16.1618543957420 16.1618606821246 16.1618662804369 16.1618711765225 16.1618753570426 16.1618788094981 16.1618815222461 16.1618834845152 16.1618846864222 16.1618851189873 16.1618847741458 16.1618836447613 16.1618817246355 16.1618790085201 16.1618754921241 16.1618711721217 16.1618660461583 16.1618601128575 16.1618533718226 16.1618458236413 16.1618374698866 16.1618283131172 16.1618183568773 16.1618076056934 16.1617960650730 16.1617837414976 16.1617706424199 16.1617567762552 16.1617421523737 16.1617267810914 16.1617106736608 16.1616938422578 16.1616762999704 16.1616580607834 16.1616391395660 16.1616195520541 16.1615993148332 16.1615784453205 16.1615569617477 16.1615348831370 16.1615122292847 16.1614890207355 16.1614652787612 16.1614410253363 16.1614162831138 16.1613910754004 16.1613654261243 16.1613393598187 16.1613129015845 16.1612860770649 16.1612589124175 16.1612314342816 16.1612036697477 16.1611756463325 16.1611473919352 16.1611189348200 16.1610903035681 16.1610615270553 16.1610326344177 16.1610036550152 16.1609746183917 16.1609455542531 16.1609164924271 16.1608874628210 16.1608584954002 16.1608296201437 16.1608008670096 16.1607722659093 16.1607438466572 16.1607156389514 16.1606876723249 16.1606599761232 16.1606325794605 16.1606055111897 16.1605787998719 16.1605524737352 16.1605265606476 16.1605010880777 16.1604760830718 16.1604515722128 16.1604275815955 16.1604041367933 16.1603812628281 16.1603589841383 16.1603373245517 16.1603163072626 16.1602959547947 16.1602762889862 16.1602573309486 16.1602391010598 16.1602216189208 16.1602049033543 16.1601889723600 16.1601738431102 16.1601595319205 16.1601460542326 16.1601334246002 16.1601216566646 16.1601107631448 16.1601007558188 16.1600916455104 16.1600834420761 16.1600761543970 16.1600697903588 16.1600643568562 16.1600598597774 16.1600563039876 16.1600536933422 16.1600520306723 16.1600513177761 16.1600515554295 16.1600527433729 16.1600548803213 16.1600579639555 16.1600619909367 16.1600669569014 16.1600728564690 16.1600796832504 16.1600874298553 16.1600960878983 16.1601056480092 16.1601160998441 16.1601274321033 16.1601396325338 16.1601526879569 16.1601665842717 16.1601813064786 16.1601968386928 16.1602131641713 16.1602302653274 16.1602481237419 16.1602667202036 16.1602860347198 16.1603060465446 16.1603267341947 16.1603480754907 16.1603700475685 16.1603926269168 16.1604157893909 16.1604395102640 16.1604637642311 16.1604885254565 16.1605137675954 16.1605394638342 16.1605655869060 16.1605921091395 16.1606190024839 16.1606462385393 16.1606737885965 16.1607016236641 16.1607297145126 16.1607580316938 16.1607865455905 16.1608152264396 16.1608440443783 16.1608729694638 16.1609019717229 16.1609310211841 16.1609600879003 16.1609891420016 16.1610181537190 16.1610470934223 16.1610759316528 16.1611046391624 16.1611331869384 16.1611615462486 16.1611896886659 16.1612175861042 16.1612452108528 16.1612725356029 16.1612995334869 16.1613261781031 16.1613524435453 16.1613783044382 16.1614037359628 16.1614287138833 16.1614532145796 16.1614772150697 16.1615006930345 16.1615236268497 16.1615459956023 16.1615677791211 16.1615889579925 16.1616095135881 16.1616294280802 16.1616486844654 16.1616672665807 16.1616851591215 16.1617023476588 16.1617188186553 16.1617345594760 16.1617495584064 16.1617638046600 16.1617772883930 16.1617900007114 16.1618019336805 16.1618130803321 16.1618234346711 16.1618329916803 16.1618417473238];
E1rms3 = [14.8666393781072 14.8666402409309 14.8666625348031 14.8667062532818 14.8667713836030 14.8668579066905 14.8669657971553 14.8670950233043 14.8672455471481 14.8674173244104 14.8676103045414 14.8678244307261 14.8680596399052 14.8683158627851 14.8685930238614 14.8688910414347 14.8692098276328 14.8695492884336 14.8699093236928 14.8702898271610 14.8706906865244 14.8711117834205 14.8715529934798 14.8720141863520 14.8724952257455 14.8729959694568 14.8735162694109 14.8740559717033 14.8746149166359 14.8751929387608 14.8757898669267 14.8764055243204 14.8770397285126 14.8776922915138 14.8783630198163 14.8790517144491 14.8797581710343 14.8804821798327 14.8812235258119 14.8819819886937 14.8827573430199 14.8835493582089 14.8843577986220 14.8851824236204 14.8860229876446 14.8868792402618 14.8877509262469 14.8886377856539 14.8895395538742 14.8904559617247 14.8913867355069 14.8923315970985 14.8932902640108 14.8942624494859 14.8952478625618 14.8962462081631 14.8972571871753 14.8982804965312 14.8993158293027 14.9003628747718 14.9014213185322 14.9024908425694 14.9035711253590 14.9046618419415 14.9057626640373 14.9068732601197 14.9079932955217 14.9091224325230 14.9102603304592 14.9114066458070 14.9125610322887 14.9137231409700 14.9148926203715 14.9160691165537 14.9172522732289 14.9184417318703 14.9196371318091 14.9208381103449 14.9220443028448 14.9232553428642 14.9244708622417 14.9256904912183 14.9269138585432 14.9281405915837 14.9293703164365 14.9306026580427 14.9318372403000 14.9330736861747 14.9343116178131 14.9355506566645 14.9367904235854 14.9380305389642 14.9392706228283 14.9405102949721 14.9417491750614 14.9429868827561 14.9442230378360 14.9454572603015 14.9466891705001 14.9479183892514 14.9491445379523 14.9503672387051 14.9515861144312 14.9528007889913 14.9540108873022 14.9552160354592 14.9564158608485 14.9576099922724 14.9587980600627 14.9599796961977 14.9611545344278 14.9623222103771 14.9634823616811 14.9646346280836 14.9657786515713 14.9669140764705 14.9680405495799 14.9691577202729 14.9702652406202 14.9713627654997 14.9724499527077 14.9735264630720 14.9745919605654 14.9756461124166 14.9766885892164 14.9777190650316 14.9787372175077 14.9797427279798 14.9807352815734 14.9817145673168 14.9826802782426 14.9836321114830 14.9845697683846 14.9854929545939 14.9864013801704 14.9872947596704 14.9881728122572 14.9890352617859 14.9898818368988 14.9907122711246 14.9915263029567 14.9923236759509 14.9931041388148 14.9938674454837 14.9946133552114 14.9953416326522 14.9960520479445 14.9967443767791 14.9974184004884 14.9980739061128 14.9987106864807 14.9993285402761 14.9999272721064 15.0005066925750 15.0010666183446 15.0016068721956 15.0021272830987 15.0026276862632 15.0031079231998 15.0035678417749 15.0040072962627 15.0044261473974 15.0048242624182 15.0052015151219 15.0055577858997 15.0058929617826 15.0062069364837 15.0064996104275 15.0067708907912 15.0070206915364 15.0072489334370 15.0074555441059 15.0076404580253 15.0078036165641 15.0079449680013 15.0080644675441 15.0081620773419 15.0082377665007 15.0082915110939 15.0083232941683 15.0083331057546 15.0083209428665 15.0082868095015 15.0082307166434 15.0081526822541 15.0080527312692 15.0079308955888 15.0077872140650 15.0076217324898 15.0074345035784 15.0072255869480 15.0069950491019 15.0067429634012 15.0064694100393 15.0061744760163 15.0058582551092 15.0055208478315 15.0051623614065 15.0047829097261 15.0043826133077 15.0039615992585 15.0035200012248 15.0030579593490 15.0025756202188 15.0020731368167 15.0015506684673 15.0010083807767 15.0004464455851 14.9998650408961 14.9992643508204 14.9986445655131 14.9980058811034 14.9973484996284 14.9966726289667 14.9959784827600 14.9952662803483 14.9945362466823 14.9937886122564 14.9930236130269 14.9922414903309 14.9914424907973 14.9906268662724 14.9897948737318 14.9889467751839 14.9880828375931 14.9872033327821 14.9863085373378 14.9853987325271 14.9844742041894 14.9835352426524 14.9825821426217 14.9816152030954 14.9806347272514 14.9796410223525 14.9786343996464 14.9776151742547 14.9765836650725 14.9755401946594 14.9744850891371 14.9734186780767 14.9723412943950 14.9712532742411 14.9701549568874 14.9690466846158 14.9679288026072 14.9668016588334 14.9656656039346 14.9645209911178 14.9633681760257 14.9622075166396 14.9610393731466 14.9598641078425 14.9586820849955 14.9574936707447 14.9562992329754 14.9550991412027 14.9538937664584 14.9526834811666 14.9514686590318 14.9502496749167 14.9490269047259 14.9478007252866 14.9465715142343 14.9453396498862 14.9441055111354 14.9428694773262 14.9416319281288 14.9403932434377 14.9391538032466 14.9379139875259 14.9366741761183 14.9354347486107 14.9341960842276 14.9329585617060 14.9317225591916 14.9304884541150 14.9292566230813 14.9280274417572 14.9268012847589 14.9255785255368 14.9243595362646 14.9231446877279 14.9219343492194 14.9207288884191 14.9195286712984 14.9183340619995 14.9171454227360 14.9159631136818 14.9147874928736 14.9136189161012 14.9124577367953 14.9113043059416 14.9101589719692 14.9090220806519 14.9078939750028 14.9067749951889 14.9056654784207 14.9045657588648 14.9034761675363 14.9023970322249 14.9013286773788 14.9002714240296 14.8992255896927 14.8981914882863 14.8971694300288 14.8961597213672 14.8951626648836 14.8941785592087 14.8932076989452 14.8922503745771 14.8913068724031 14.8903774744378 14.8894624583538 14.8885620973884 14.8876766602838 14.8868064111968 14.8859516096405 14.8851125104102 14.8842893635021 14.8834824140642 14.8826919023155 14.8819180634873 14.8811611277562 14.8804213201853 14.8796988606571 14.8789939638251 14.8783068390436 14.8776376903207 14.8769867162603 14.8763541100062 14.8757400591975 14.8751447459118 14.8745683466176 14.8740110321311 14.8734729675680 14.8729543122984 14.8724552199087 14.8719758381582 14.8715163089369 14.8710767682346 14.8706573460989 14.8702581666081 14.8698793478271 14.8695210017903 14.8691832344589 14.8688661457016 14.8685698292660 14.8682943727531 14.8680398575933 14.8678063590292 14.8675939460899 14.8674026815788 14.8672326220518 14.8670838178060 14.8669563128652 14.8668501449681 14.8667653455591 14.8667019397772 14.8666599464551 14.8666393781072];

Emed2 = mean(E1rms2);
Emax2 = max(E1rms2);
Kirreg2 = Emax2/Emed2;

h = figure();
plot(thetad,E1rms3,'-','LineWidth',3.5);
ylabel('Campo el�trico (kV/cm)')
xlabel('�ngulo (�)')
grid
hold on
plot(thetad,Erms,':','LineWidth',3.5);
hold on
plot(thetad,E1rms,'-.','LineWidth',3.5);
hold on
plot(thetad,E1_crit,'--','LineWidth',3.5);
xlim([0 360])
% ylim([13.5 15.5])
legend('Original', 'Aprox. 1 m', 'Aprox. 2 m', 'Campo el�trico cr�tico')
% title('Campo el�trico superficial condutor um')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'solo_caso1_dist','-dpdf','-r0')