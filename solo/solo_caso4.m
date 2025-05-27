clear all; close all; clc

n = 12; % número de condutores
e_0 = 8.854*(10^(-12));
r = 15.980*(10^-3); % raio do condutor

x = [-7.975 -7.025 -7.025 -7.975 -0.475 0.475 0.475 -0.475 7.025 7.975 7.975 7.025];
y = [18.45 18.45 17.5 17.5 25.95 25.95 25 25 18.45 18.45 17.5 17.5];
yi = -y;

P = zeros(n,n);
for i = 1:n
     for j = 1:n
        if(i==j)
           P(i,j) = (1/(2*pi*e_0))*log((4*y(i))/(2*r));
        else
           P(i,j) = (1/(2*pi*e_0))*log(sqrt((x(j)-x(i))^2+(yi(j)-y(i))^2)/(sqrt((x(j)-x(i))^2+(y(j)-y(i))^2)));
        end
     end
end

%%

V = 500*10^3;

%todos os cabos da fase A ( Tensão da Fase A com angulo de 0 grau),
%Todos os cabos da fase B ( Tensão da Fase B com angulo de +120 graus),
%Todos os cabos da Fase C ( Tensão da Fase C com angulo de -120 graus),

%v(k,1) = V/sqrt(3);
%v(k+b,1) = V*(cos(2*pi/3)+1i*sin(2*pi/3))/sqrt(3);
%v(k+2*b,1) = V*(cos(-2*pi/3)+1i*sin(-2*pi/3))/sqrt(3);

V_ra = V/sqrt(3);
V_ia = 0;

V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

V = [ V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ];

%% AP : Funciona destas duas maneiras e a ultima é mais rápida a barra invertida indica inversão de P:
%rho = inv(P)*V;

rho = P\V;

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

% x é o vetor de pontos onde se calcula o somatório dos campos
% (contribuição dos cabos de cada uma das fases e das imagens dos cabos

% os limites deste vetor devem ser maiores do que os limites dados pelas
% posições extremas dos cabos:

xf = linspace(-25,25,179);

% especifica a altura da linha de interesse
yf = 1;

%% AP: aqui tinha faltado o . , ele indica que a operação deve ser vetorial ou seja pega o vetor x e subtrai x1,x(2),x(3) e etc:

%E_xr1 componente x real campo elétrico fase a cabo 1
%E_xr2 componente x real campo elétrico fase a cabo 2
%E_xr3 componente x real campo elétrico fase a cabo 3
%E_xr4 componente x real campo elétrico fase a cabo 4

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));
E_xr4 = (rho_r4/(2*pi*e_0)).*((xf - x(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (xf - x(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));

%E_xr5 componente x real campo elétrico fase b cabo 5
%E_xr6 componente x real campo elétrico fase b cabo 6
%E_xr7 componente x real campo elétrico fase b cabo 7
%E_xr8 componente x real campo elétrico fase b cabo 8

E_xr5 = (rho_r5/(2*pi*e_0)).*((xf - x(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (xf - x(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_xr6 = (rho_r6/(2*pi*e_0)).*((xf - x(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (xf - x(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));
E_xr7 = (rho_r7/(2*pi*e_0)).*((xf - x(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (xf - x(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_xr8 = (rho_r8/(2*pi*e_0)).*((xf - x(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (xf - x(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));

%E_xr9 componente x real campo elétrico fase c cabo 9
%E_xr10 componente x real campo elétrico fase c cabo 10
%E_xr11 componente x real campo elétrico fase c cabo 11
%E_xr12 componente x real campo elétrico fase c cabo 12

E_xr9 = (rho_r9/(2*pi*e_0)).*((xf - x(9))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (xf - x(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));
E_xr10 = (rho_r10/(2*pi*e_0)).*((xf - x(10))./(((xf - x(10)).^2) + (yf - y(10)).^2) - (xf - x(10))./(((xf - x(10)).^2) + (yf + y(10)).^2));
E_xr11 = (rho_r11/(2*pi*e_0)).*((xf - x(11))./(((xf - x(11)).^2) + (yf - y(11)).^2) - (xf - x(11))./(((xf - x(11)).^2) + (yf + y(11)).^2));
E_xr12 = (rho_r12/(2*pi*e_0)).*((xf - x(12))./(((xf - x(12)).^2) + (yf - y(12)).^2) - (xf - x(12))./(((xf - x(12)).^2) + (yf + y(12)).^2));
%%

%E_xi1 componente x imaginario campo elétrico fase a cabo 1
%E_xi2 componente x imaginario campo elétrico fase a cabo 2
%E_xi3 componente x imaginario campo elétrico fase a cabo 3
%E_xi4 componente x imaginario campo elétrico fase a cabo 4

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));
E_xi4 = (rho_i4/(2*pi*e_0)).*((xf - x(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (xf - x(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));

%E_xi5 componente x imaginario campo elétrico fase b cabo 5
%E_xi6 componente x imaginario campo elétrico fase b cabo 6
%E_xi7 componente x imaginario campo elétrico fase b cabo 7
%E_xi8 componente x imaginario campo elétrico fase b cabo 8

E_xi5 = (rho_i5/(2*pi*e_0)).*((xf - x(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (xf - x(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_xi6 = (rho_i6/(2*pi*e_0)).*((xf - x(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (xf - x(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));
E_xi7 = (rho_i7/(2*pi*e_0)).*((xf - x(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (xf - x(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_xi8 = (rho_i8/(2*pi*e_0)).*((xf - x(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (xf - x(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));

%E_xi9 componente x imaginario campo elétrico fase c cabo 9
%E_xi10 componente x imaginario campo elétrico fase c cabo 10
%E_xi11 componente x imaginario campo elétrico fase c cabo 11
%E_xi12 componente x imaginario campo elétrico fase c cabo 12

E_xi9 = (rho_i9/(2*pi*e_0)).*((xf - x(7))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (xf - x(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));
E_xi10 = (rho_i10/(2*pi*e_0)).*((xf - x(10))./(((xf - x(10)).^2) + (yf - y(10)).^2) - (xf - x(10))./(((xf - x(10)).^2) + (yf + y(10)).^2));
E_xi11 = (rho_i11/(2*pi*e_0)).*((xf - x(11))./(((xf - x(11)).^2) + (yf - y(11)).^2) - (xf - x(11))./(((xf - x(11)).^2) + (yf + y(11)).^2));
E_xi12 = (rho_i12/(2*pi*e_0)).*((xf - x(12))./(((xf - x(12)).^2) + (yf - y(12)).^2) - (xf - x(12))./(((xf - x(12)).^2) + (yf + y(12)).^2));

%%

%E_yr1 componente y real campo elétrico fase a cabo 1
%E_yr2 componente y real campo elétrico fase a cabo 2
%E_yr3 componente y real campo elétrico fase a cabo 3
%E_yr4 componente y real campo elétrico fase a cabo 4

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (yf + y(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (yf + y(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));
E_yr4 = (rho_r4/(2*pi*e_0)).*((yf - y(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (yf + y(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));

%E_yr5 componente y real campo elétrico fase b cabo 5
%E_yr6 componente y real campo elétrico fase b cabo 6
%E_yr7 componente y real campo elétrico fase b cabo 7
%E_yr8 componente y real campo elétrico fase b cabo 8

E_yr5 = (rho_r5/(2*pi*e_0)).*((yf - y(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (yf + y(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_yr6 = (rho_r6/(2*pi*e_0)).*((yf - y(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (yf + y(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));
E_yr7 = (rho_r7/(2*pi*e_0)).*((yf - y(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (yf + y(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_yr8 = (rho_r8/(2*pi*e_0)).*((yf - y(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (yf + y(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));

%E_yr9 componente y real campo elétrico fase c cabo 9
%E_yr10 componente y real campo elétrico fase c cabo 10
%E_yr11 componente y real campo elétrico fase c cabo 11
%E_yr12 componente y real campo elétrico fase c cabo 12

E_yr9 = (rho_r9/(2*pi*e_0)).*((yf - y(9))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (yf + y(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));
E_yr10 = (rho_r10/(2*pi*e_0)).*((yf - y(10))./(((xf - x(10)).^2) + (yf - y(10)).^2) - (yf + y(10))./(((xf - x(10)).^2) + (yf + y(10)).^2));
E_yr11 = (rho_r11/(2*pi*e_0)).*((yf - y(11))./(((xf - x(11)).^2) + (yf - y(11)).^2) - (yf + y(11))./(((xf - x(11)).^2) + (yf + y(11)).^2));
E_yr12 = (rho_r12/(2*pi*e_0)).*((yf - y(12))./(((xf - x(12)).^2) + (yf - y(12)).^2) - (yf + y(12))./(((xf - x(12)).^2) + (yf + y(12)).^2));

%%

%E_yi1 componente y imaginario campo elétrico fase a cabo 1
%E_yi2 componente y imaginario campo elétrico fase a cabo 2
%E_yi3 componente y imaginario campo elétrico fase a cabo 3
%E_yi4 componente y imaginario campo elétrico fase a cabo 4

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (yf + y(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (yf + y(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));
E_yi4 = (rho_i4/(2*pi*e_0)).*((yf - y(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (yf + y(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));

%E_yi5 componente y imaginario campo elétrico fase b cabo 5
%E_yi6 componente y imaginario campo elétrico fase b cabo 6
%E_yi7 componente y imaginario campo elétrico fase b cabo 7
%E_yi8 componente y imaginario campo elétrico fase b cabo 8

E_yi5 = (rho_i5/(2*pi*e_0)).*((yf - y(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (yf + y(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_yi6 = (rho_i6/(2*pi*e_0)).*((yf - y(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (yf + y(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));
E_yi7 = (rho_i7/(2*pi*e_0)).*((yf - y(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (yf + y(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_yi8 = (rho_i8/(2*pi*e_0)).*((yf - y(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (yf + y(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));

%E_yi9 componente y imaginario campo elétrico fase c cabo 9
%E_yi10 componente y imaginario campo elétrico fase c cabo 10
%E_yi11 componente y imaginario campo elétrico fase c cabo 11
%E_yi12 componente y imaginario campo elétrico fase c cabo 12

E_yi9 = (rho_i9/(2*pi*e_0)).*((yf - y(9))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (yf + y(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));
E_yi10 = (rho_i10/(2*pi*e_0)).*((yf - y(10))./(((xf - x(10)).^2) + (yf - y(10)).^2) - (yf + y(10))./(((xf - x(10)).^2) + (yf + y(10)).^2));
E_yi11 = (rho_i11/(2*pi*e_0)).*((yf - y(11))./(((xf - x(11)).^2) + (yf - y(11)).^2) - (yf + y(11))./(((xf - x(11)).^2) + (yf + y(11)).^2));
E_yi12 = (rho_i12/(2*pi*e_0)).*((yf - y(12))./(((xf - x(12)).^2) + (yf - y(12)).^2) - (yf + y(12))./(((xf - x(12)).^2) + (yf + y(12)).^2));

%% AP: aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr1;  E_xr2 ;  E_xr3 ; E_xr4 ; E_xr5 ; E_xr6 ; E_xr7 ; E_xr8 ; E_xr9 ; E_xr10 ; E_xr11 ; E_xr12 ];

E_ximag = [ E_xi1;  E_xi2 ;  E_xi3 ; E_xi4 ; E_xi5 ; E_xi6 ; E_xi7 ; E_xi8 ; E_xi9 ; E_xi10 ; E_xi11 ; E_xi12 ];

E_yr = [ E_yr1;  E_yr2 ;  E_yr3 ; E_yr4 ; E_yr5 ; E_yr6 ; E_yr7 ; E_yr8 ; E_yr9 ; E_yr10 ; E_yr11 ; E_yr12 ];

E_yimag = [ E_yi1;  E_yi2 ;  E_yi3 ; E_yi4 ; E_yi5 ; E_yi6 ; E_yi7 ; E_yi8 ; E_yi9 ; E_yi10 ; E_yi11 ; E_yi12 ];


%% AP: aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

Exr = sum(E_xr).^2;

Eximag = sum(E_ximag).^2;

Eyr = sum(E_yr).^2;

Eyimag = sum(E_yimag).^2;

%% AP: aqui retira-se o módulo dos quadrados:

Erms = (Exr+Eximag+Eyr+Eyimag).^(1/2);

%% Tem outro erro aqui tb - é tudo soma:

x_p = [-25.00	,	-24.80	,	-24.53	,	-24.30	,	-24.06	,	-23.75	,	-23.48	,	-23.22	,	-22.99	,	-22.77	,	-22.53	,	-22.25	,	-21.97	,	-21.71	,	-21.55	,	-21.27	,	-20.96	,	-20.72	,	-20.48	,	-20.22	,	-20.03	,	-19.72	,	-19.44	,	-19.18	,	-18.97	,	-18.73	,	-18.50	,	-18.22	,	-17.97	,	-17.72	,	-17.48	,	-17.18	,	-16.91	,	-16.63	,	-16.42	,	-16.14	,	-15.95	,	-15.67	,	-15.39	,	-15.14	,	-14.95	,	-14.69	,	-14.37	,	-14.14	,	-13.91	,	-13.67	,	-13.47	,	-13.14	,	-12.88	,	-12.59	,	-12.22	,	-11.85	,	-11.49	,	-11.14	,	-10.77	,	-10.38	,	-10.03	,	-9.68	,	-9.33	,	-8.87	,	-8.59	,	-8.36	,	-8.08	,	-7.81	,	-7.57	,	-7.33	,	-7.12	,	-6.91	,	-6.65	,	-6.32	,	-6.10	,	-5.86	,	-5.62	,	-5.38	,	-5.11	,	-4.83	,	-4.59	,	-4.35	,	-4.10	,	-3.78	,	-3.54	,	-3.33	,	-3.11	,	-2.85	,	-2.60	,	-2.29	,	-2.05	,	-1.82	,	-1.54	,	-1.25	,	-1.02	,	-0.68	,	-0.38	,	-0.01	,	0.38	,	0.78	,	1.02	,	1.29	,	1.64	,	1.90	,	2.18	,	2.41	,	2.60	,	2.82	,	3.02	,	3.21	,	3.35	,	3.50	,	3.68	,	3.85	,	4.04	,	4.29	,	4.50	,	4.74	,	4.98	,	5.30	,	5.54	,	5.77	,	6.04	,	6.26	,	6.50	,	6.82	,	7.10	,	7.30	,	7.56	,	7.75	,	8.05	,	8.35	,	8.63	,	8.86	,	9.34	,	9.77	,	10.13	,	10.51	,	10.95	,	11.35	,	11.77	,	12.21	,	12.64	,	13.15	,	13.51	,	13.84	,	14.28	,	14.64	,	14.96	,	15.29	,	15.51	,	15.83	,	16.02	,	16.25	,	16.51	,	16.78	,	17.14	,	17.43	,	17.76	,	18.00	,	18.46	,	18.68	,	18.92	,	19.15	,	19.39	,	19.71	,	19.98	,	20.22	,	20.46	,	20.67	,	20.90	,	21.19	,	21.48	,	21.90	,	22.30	,	22.58	,	22.95	,	23.23	,	23.59	,	23.98	,	24.29	,	24.54	,	25.00];
y_p = [2216.42	,	2259.20	,	2298.01	,	2335.82	,	2375.62	,	2417.41	,	2457.21	,	2501.99	,	2545.77	,	2586.57	,	2631.34	,	2677.11	,	2719.90	,	2763.68	,	2805.47	,	2850.25	,	2894.03	,	2938.81	,	2980.60	,	3030.35	,	3078.11	,	3124.88	,	3175.62	,	3220.40	,	3262.19	,	3310.95	,	3359.70	,	3404.48	,	3447.26	,	3492.04	,	3533.83	,	3577.61	,	3620.40	,	3669.15	,	3709.95	,	3750.75	,	3789.55	,	3830.35	,	3867.16	,	3907.96	,	3944.78	,	3974.63	,	4009.45	,	4040.30	,	4068.16	,	4097.01	,	4121.89	,	4145.77	,	4164.68	,	4192.54	,	4218.41	,	4238.31	,	4250.25	,	4258.21	,	4249.25	,	4233.33	,	4213.43	,	4190.55	,	4166.67	,	4135.82	,	4100.00	,	4075.12	,	4030.35	,	3996.52	,	3949.75	,	3912.94	,	3866.17	,	3814.43	,	3763.68	,	3710.95	,	3655.22	,	3596.52	,	3530.85	,	3466.17	,	3398.51	,	3332.84	,	3265.17	,	3200.50	,	3131.84	,	3069.15	,	3000.50	,	2937.81	,	2876.12	,	2816.42	,	2752.74	,	2701.99	,	2648.26	,	2602.49	,	2554.73	,	2522.89	,	2486.07	,	2456.22	,	2441.29	,	2425.37	,	2443.28	,	2478.11	,	2510.95	,	2544.78	,	2593.53	,	2634.33	,	2686.07	,	2730.85	,	2772.64	,	2816.42	,	2864.18	,	2908.96	,	2957.71	,	3001.49	,	3054.23	,	3099.00	,	3145.77	,	3201.49	,	3267.16	,	3329.85	,	3389.55	,	3464.18	,	3524.88	,	3596.52	,	3653.23	,	3708.96	,	3762.69	,	3815.42	,	3871.14	,	3920.90	,	3963.68	,	4002.49	,	4041.29	,	4078.11	,	4116.92	,	4145.77	,	4181.59	,	4211.44	,	4234.33	,	4247.26	,	4258.21	,	4244.28	,	4225.37	,	4206.47	,	4176.62	,	4138.81	,	4105.97	,	4068.16	,	4012.44	,	3969.65	,	3932.84	,	3879.10	,	3845.27	,	3804.48	,	3772.64	,	3729.85	,	3687.06	,	3628.36	,	3571.64	,	3522.89	,	3476.12	,	3425.37	,	3359.70	,	3309.95	,	3261.19	,	3212.44	,	3175.62	,	3121.89	,	3074.13	,	3028.36	,	2978.61	,	2929.85	,	2889.05	,	2847.26	,	2801.49	,	2721.89	,	2662.19	,	2614.43	,	2542.79	,	2491.04	,	2431.34	,	2379.60	,	2324.88	,	2288.06	,	2225.37];

figure();
plot(xf,Erms,'LineWidth',4.5)
xlabel('Faixa de passagem (m)')
ylabel('Campo elétrico (V/m)')
%title('Campo elétrico ao nível do solo quatro condutores por fase')
grid
% hold on
% yy2 = smooth(x_p,y_p,'rloess');
% plot(x_p,yy2,'--','LineWidth',2.5);
% xlim([-25 25])
% ylim([1500 5000])
% legend('Simulado', '(PAGANOTTI, 2012)')

% h = figure;
% plot(1:10);
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'filename','-dpdf','-r0')


h = figure();
plot(x(1:4), y(1:4),'o', 'color', 'blue')
hold on
plot(x(5:8), y(5:8),'*', 'color', 'blue')
hold on
plot(x(9:12), y(9:12),'x', 'color', 'blue')
xlabel('Faixa de passagem (m)')
ylabel('Altura do cabo (m)')
%title('Configuração geométrica dos condutores')
legend('Cabos condutores fase A','Cabos condutores fase B','Cabos condutores fase C');
grid
xlim([-8.5 8.5])
ylim([17 27])

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'caso4config','-dpdf','-r0')

err = abs((y_p - Erms)/y_p)*100;