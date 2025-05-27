clear all; close all; clc

% 10.478 m entre fase ao centro e extremos
% espaçamento 0.228 m
% nível de tensão adotado: 500 kV
% altura 16.530 m e 16.758 acima do solo
% raio do cabo 0.001437 m
% campo elétrico ao nível do solo:
% x=-25m a x = 25m e y = 1m em 100 pontos diferentes

%A seguir será tratado o caso para 3 cabos por fase

%e_0 é a permissividade do vácuo [F/m]
%P é a matriz de coeficientes potenciais
%C é a matriz de capacitâncias, sendo a inversa de P
%V_r são componentes reais das tensões de fase para cada condutor
%V_i são componentes imaginárias das tensões de fase para cada condutor
%rho são as densidades de carga para cada condutor

n = 9; % número de condutores
e_0 = 8.854*(10^(-12));
r = 1.437*(10^(-2)); % raio do condutor

x = [-10.478 -10.25 -10.021 -0.228 0 0.228 10.021 10.25 10.478];
y = [16.530 16.758 16.530 16.530 16.758 16.530 16.530 16.758 16.530];
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

V = [ V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ];

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

% x é o vetor de pontos onde se calcula o somatório dos campos
% (contribuição dos cabos de cada uma das fases e das imagens dos cabos

% os limites deste vetor devem ser maiores do que os limites dados pelas
% posições extremas dos cabos:

xf = linspace(-25,25,161);

% especifica a altura da linha de interesse
yf = 1;

%% AP: aqui tinha faltado o . , ele indica que a operação deve ser vetorial ou seja pega o vetor x e subtrai x1,x(2),x(3) e etc:

%E_xr1 componente x real campo elétrico fase a cabo 1
%E_xr2 componente x real campo elétrico fase a cabo 2
%E_xr3 componente x real campo elétrico fase a cabo 3

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));

%E_xr4 componente x real campo elétrico fase b cabo 4
%E_xr5 componente x real campo elétrico fase b cabo 5
%E_xr6 componente x real campo elétrico fase b cabo 6

E_xr4 = (rho_r4/(2*pi*e_0)).*((xf - x(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (xf - x(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));
E_xr5 = (rho_r5/(2*pi*e_0)).*((xf - x(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (xf - x(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_xr6 = (rho_r6/(2*pi*e_0)).*((xf - x(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (xf - x(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));

%E_xr7 componente x real campo elétrico fase c cabo 7
%E_xr8 componente x real campo elétrico fase c cabo 8
%E_xr9 componente x real campo elétrico fase c cabo 9

E_xr7 = (rho_r7/(2*pi*e_0)).*((xf - x(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (xf - x(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_xr8 = (rho_r8/(2*pi*e_0)).*((xf - x(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (xf - x(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));
E_xr9 = (rho_r9/(2*pi*e_0)).*((xf - x(9))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (xf - x(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));

%%

%E_xi1 componente x imaginario campo elétrico fase a cabo 1
%E_xi2 componente x imaginario campo elétrico fase a cabo 2
%E_xi3 componente x imaginario campo elétrico fase a cabo 3

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));

%E_xi4 componente x imaginario campo elétrico fase b cabo 4
%E_xi5 componente x imaginario campo elétrico fase b cabo 5
%E_xi6 componente x imaginario campo elétrico fase b cabo 6

E_xi4 = (rho_i4/(2*pi*e_0)).*((xf - x(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (xf - x(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));
E_xi5 = (rho_i5/(2*pi*e_0)).*((xf - x(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (xf - x(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_xi6 = (rho_i6/(2*pi*e_0)).*((xf - x(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (xf - x(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));

%E_xi4 componente x imaginario campo elétrico fase c cabo 7
%E_xi5 componente x imaginario campo elétrico fase c cabo 8
%E_xi6 componente x imaginario campo elétrico fase c cabo 9

E_xi7 = (rho_i7/(2*pi*e_0)).*((xf - x(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (xf - x(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_xi8 = (rho_i8/(2*pi*e_0)).*((xf - x(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (xf - x(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));
E_xi9 = (rho_i9/(2*pi*e_0)).*((xf - x(9))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (xf - x(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));

%%

%E_yr1 componente y real campo elétrico fase a cabo 1
%E_yr2 componente y real campo elétrico fase a cabo 2
%E_yr3 componente y real campo elétrico fase b cabo 3

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (yf + y(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (yf + y(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));

%E_yr4 componente y real campo elétrico fase b cabo 4
%E_yr5 componente y real campo elétrico fase b cabo 5
%E_yr6 componente y real campo elétrico fase b cabo 6

E_yr4 = (rho_r4/(2*pi*e_0)).*((yf - y(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (yf + y(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));
E_yr5 = (rho_r5/(2*pi*e_0)).*((yf - y(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (yf + y(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_yr6 = (rho_r6/(2*pi*e_0)).*((yf - y(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (yf + y(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));

%E_yr7 componente y real campo elétrico fase c cabo 7
%E_yr8 componente y real campo elétrico fase c cabo 8
%E_yr9 componente y real campo elétrico fase c cabo 9

E_yr7 = (rho_r7/(2*pi*e_0)).*((yf - y(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (yf + y(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_yr8 = (rho_r8/(2*pi*e_0)).*((yf - y(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (yf + y(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));
E_yr9 = (rho_r9/(2*pi*e_0)).*((yf - y(9))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (yf + y(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));

%%

%E_yi1 componente y imaginario campo elétrico fase a cabo 1
%E_yi2 componente y imaginario campo elétrico fase a cabo 2
%E_yi3 componente y imaginario campo elétrico fase a cabo 3

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y(2))./(((xf - x(2)).^2) + (yf - y(2)).^2) - (yf + y(2))./(((xf - x(2)).^2) + (yf + y(2)).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y(3))./(((xf - x(3)).^2) + (yf - y(3)).^2) - (yf + y(3))./(((xf - x(3)).^2) + (yf + y(3)).^2));

%E_yi4 componente y imaginario campo elétrico fase b cabo 4
%E_yi5 componente y imaginario campo elétrico fase b cabo 5
%E_yi6 componente y imaginario campo elétrico fase b cabo 6

E_yi4 = (rho_i4/(2*pi*e_0)).*((yf - y(4))./(((xf - x(4)).^2) + (yf - y(4)).^2) - (yf + y(4))./(((xf - x(4)).^2) + (yf + y(4)).^2));
E_yi5 = (rho_i5/(2*pi*e_0)).*((yf - y(5))./(((xf - x(5)).^2) + (yf - y(5)).^2) - (yf + y(5))./(((xf - x(5)).^2) + (yf + y(5)).^2));
E_yi6 = (rho_i6/(2*pi*e_0)).*((yf - y(6))./(((xf - x(6)).^2) + (yf - y(6)).^2) - (yf + y(6))./(((xf - x(6)).^2) + (yf + y(6)).^2));

%E_yi7 componente y imaginario campo elétrico fase c cabo 7
%E_yi8 componente y imaginario campo elétrico fase c cabo 8
%E_yi9 componente y imaginario campo elétrico fase c cabo 9

E_yi7 = (rho_i7/(2*pi*e_0)).*((yf - y(7))./(((xf - x(7)).^2) + (yf - y(7)).^2) - (yf + y(7))./(((xf - x(7)).^2) + (yf + y(7)).^2));
E_yi8 = (rho_i8/(2*pi*e_0)).*((yf - y(8))./(((xf - x(8)).^2) + (yf - y(8)).^2) - (yf + y(8))./(((xf - x(8)).^2) + (yf + y(8)).^2));
E_yi9 = (rho_i9/(2*pi*e_0)).*((yf - y(9))./(((xf - x(9)).^2) + (yf - y(9)).^2) - (yf + y(9))./(((xf - x(9)).^2) + (yf + y(9)).^2));

%% AP: aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr1; E_xr2 ; E_xr3 ; E_xr4 ; E_xr5 ; E_xr6 ; E_xr7 ; E_xr8 ; E_xr9 ];

E_ximag = [ E_xi1; E_xi2 ; E_xi3 ; E_xi4 ; E_xi5 ; E_xi6 ; E_xi7 ; E_xi8 ; E_xi9  ];

E_yr = [ E_yr1; E_yr2 ; E_yr3 ; E_yr4 ; E_yr5 ; E_yr6 ; E_yr7 ; E_yr8 ; E_yr9  ];

E_yimag = [ E_yi1; E_yi2 ; E_yi3 ; E_yi4 ; E_yi5 ; E_yi6 ; E_yi7 ; E_yi8 ; E_yi9 ];


%% AP: aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

Exr = sum(E_xr).^2;

Eximag = sum(E_ximag).^2;

Eyr = sum(E_yr).^2;

Eyimag = sum(E_yimag).^2;

%% AP: aqui retira-se o módulo dos quadrados:

Erms = (Exr+Eximag+Eyr+Eyimag).^(1/2);

%% Tem outro erro aqui tb - é tudo soma:

x_p = [-25.00	,	-24.77	,	-24.50	,	-24.27	,	-24.02	,	-23.61	,	-23.37	,	-23.13	,	-22.91	,	-22.66	,	-22.46	,	-22.25	,	-21.99	,	-21.74	,	-21.51	,	-21.28	,	-21.01	,	-20.78	,	-20.49	,	-20.25	,	-19.99	,	-19.75	,	-19.43	,	-19.18	,	-18.88	,	-18.58	,	-18.28	,	-17.93	,	-17.65	,	-17.35	,	-16.97	,	-16.62	,	-16.28	,	-15.89	,	-15.47	,	-15.11	,	-14.69	,	-14.26	,	-13.84	,	-13.40	,	-13.02	,	-12.61	,	-12.19	,	-11.79	,	-11.35	,	-10.92	,	-10.65	,	-10.36	,	-10.09	,	-9.81	,	-9.51	,	-9.25	,	-9.00	,	-8.75	,	-8.55	,	-8.40	,	-8.21	,	-7.99	,	-7.75	,	-7.56	,	-7.35	,	-7.14	,	-6.97	,	-6.77	,	-6.50	,	-6.18	,	-5.95	,	-5.64	,	-5.44	,	-5.14	,	-4.83	,	-4.56	,	-4.31	,	-3.96	,	-3.67	,	-3.38	,	-3.07	,	-2.82	,	-2.36	,	-1.89	,	-1.45	,	-1.02	,	-0.57	,	-0.27	,	0.03	,	0.28	,	0.71	,	1.19	,	1.67	,	1.94	,	2.27	,	2.71	,	3.02	,	3.33	,	3.64	,	3.99	,	4.33	,	4.59	,	4.83	,	5.07	,	5.33	,	5.52	,	5.80	,	5.94	,	6.23	,	6.46	,	6.84	,	7.08	,	7.31	,	7.53	,	7.77	,	8.03	,	8.31	,	8.59	,	8.79	,	9.06	,	9.37	,	9.61	,	9.85	,	10.13	,	10.49	,	10.82	,	11.26	,	11.70	,	12.07	,	12.42	,	12.84	,	13.35	,	13.88	,	14.49	,	14.95	,	15.42	,	15.93	,	16.39	,	16.79	,	17.27	,	17.55	,	17.92	,	18.37	,	18.63	,	18.93	,	19.18	,	19.42	,	19.81	,	20.13	,	20.40	,	20.67	,	20.90	,	21.27	,	21.64	,	21.95	,	22.27	,	22.62	,	22.92	,	23.22	,	23.53	,	23.87	,	24.16	,	24.43	,	24.68	,	25.00];
y_p = [2417.91	,	2460.20	,	2500.00	,	2549.75	,	2595.36	,	2655.06	,	2687.40	,	2727.20	,	2764.51	,	2808.46	,	2849.09	,	2887.23	,	2926.20	,	2964.34	,	3007.46	,	3052.24	,	3097.84	,	3140.96	,	3185.74	,	3223.05	,	3276.12	,	3322.55	,	3367.33	,	3410.45	,	3452.74	,	3511.61	,	3548.09	,	3596.19	,	3637.65	,	3679.93	,	3729.68	,	3767.83	,	3805.97	,	3850.75	,	3883.08	,	3909.62	,	3935.32	,	3959.37	,	3971.81	,	3981.76	,	3970.98	,	3949.42	,	3922.89	,	3893.03	,	3854.06	,	3811.77	,	3775.29	,	3733.00	,	3688.23	,	3641.79	,	3583.75	,	3545.61	,	3495.85	,	3447.76	,	3409.62	,	3361.53	,	3314.26	,	3264.51	,	3220.56	,	3172.47	,	3123.55	,	3070.48	,	3023.22	,	2974.30	,	2903.81	,	2832.50	,	2771.14	,	2710.61	,	2649.25	,	2572.14	,	2507.46	,	2438.64	,	2378.94	,	2317.58	,	2272.80	,	2212.27	,	2157.55	,	2102.82	,	2040.63	,	1982.59	,	1942.79	,	1916.25	,	1893.86	,	1883.91	,	1875.62	,	1883.08	,	1901.33	,	1928.69	,	1966.00	,	2000.00	,	2042.29	,	2096.19	,	2139.30	,	2201.49	,	2277.78	,	2334.16	,	2398.84	,	2458.54	,	2519.07	,	2568.82	,	2626.04	,	2683.25	,	2747.10	,	2805.97	,	2871.48	,	2931.18	,	2995.85	,	3061.36	,	3120.23	,	3191.54	,	3247.10	,	3308.46	,	3359.87	,	3412.11	,	3461.86	,	3516.58	,	3564.68	,	3608.62	,	3655.06	,	3710.61	,	3756.22	,	3805.14	,	3852.40	,	3885.57	,	3919.57	,	3941.96	,	3961.86	,	3975.12	,	3961.03	,	3936.98	,	3908.79	,	3878.94	,	3840.80	,	3792.70	,	3736.32	,	3683.25	,	3643.45	,	3591.21	,	3529.85	,	3485.07	,	3436.98	,	3393.86	,	3352.40	,	3293.53	,	3251.24	,	3196.52	,	3145.11	,	3099.50	,	3044.78	,	2989.22	,	2934.49	,	2881.43	,	2824.21	,	2766.17	,	2708.13	,	2658.37	,	2601.99	,	2557.21	,	2509.95	,	2464.34	,	2417.91];

figure();
plot(xf,Erms,'LineWidth',4.5)
xlabel('Faixa de passagem (m)')
ylabel('Campo elétrico (V/m)')
grid
% hold on
% yy2 = smooth(x_p,y_p,'rloess');
% plot(x_p,yy2,'--','LineWidth',2.5);
% xlim([-25 25])
% ylim([1500 4500])
% legend('Simulado', '(PAGANOTTI, 2012)')

h = figure();
plot(x(1:3), y(1:3),'o', 'color', 'blue')
hold on
plot(x(4:6), y(4:6),'*', 'color', 'blue')
hold on
plot(x(7:9), y(7:9),'x', 'color', 'blue')
xlabel('Faixa de passagem (m)')
ylabel('Altura do cabo (m)')
%title('Configuração geométrica dos condutores')
legend('Cabos condutores fase A','Cabos condutores fase B','Cabos condutores fase C');
grid
xlim([-15 15])
ylim([16.5 16.85])

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'caso3config','-dpdf','-r0')

% máximos para comparar

Emax1 = max(Erms);
Emin1 = min(Erms);

Emax2 = max(y_p);
Emin2 = min(y_p);

%% AP: Aqui não entendi pois pra cima se tem apenas um cabo por fase e aqui parece que vc quer trabalhar com 3 cabos por fase:


%=================================
%código para gerar a figura representativa dos cabos
% 
% cabo1a = -10.48;
% cabo2a = -10.25;
% cabo3a = -10.02;
% cabo1b = -0.23;
% cabo2b = 0;
% cabo3b = 0.23;
% cabo1c = 10.02;
% cabo2c = 10.25;
% cabo3c = 10.48;
% ya = 16.53;
% yb = 16.78;
% 
% figure();
% plot([cabo1a cabo2a cabo3a], [ya yb ya],'o');
% hold on
% plot([cabo1b cabo2b cabo3b], [ya yb ya],'o');
% hold on
% plot([cabo1c cabo2c cabo3c], [ya yb ya],'o');
% axis([-15 15 0 20]);
% grid
% legend('Cabos da fase A','Cabos da fase B','Cabos da fase C');
% xlabel('Faixa de Passagem (m)')
% ylabel('Altura do cabo (m)')
% 
% %=====================================
% 
% figure();
% subplot(2,1,1);
% plot(x,Erms)
% ylabel('Campo Elétrico (V/m)')
% grid
% 
% subplot(2,1,2); 
% plot([x(1) x(2) x(3)],[y(1) y(2) y(3)],'*')
% legend('Cabos condutores equivalentes por fase');
% xlabel('Faixa de Passagem (m)')
% ylabel('Altura do cabo (m)')
% grid

% H_k_1 = 16.530;
% H_k_2 = 16.758;
% S_kl_1 = 0.228; 
% S_kl_2 = 0.456;
% S_kl_3 = 9.793;
% S_kl_4 = 10.021; % = 6
% S_kl_5 = 10.25; % = 7, 9
% S_kl_8 = 10.478; % = 10
% S_kl_11 = 10.706;
% S_kl_12 = 20.042; 
% S_kl_13 = 20.271; % = 15
% S_kl_14 = 20.499; % = 16, 18
% S_kl_17 = 20.728; % = 19
% S_kl_20 = 20.956;

% P = zeros(n,n);
% for i = 1:n
%      for j = 1:n
%         % i == j -> 1
%         if(i==1 && j==1  || i==3 && j==3 || i==4 && j==4 || i==6 && j==6 || i==7 && j==7 || i==9 && j==9)
%            P(i,j) = (1/(2*pi*e_0))*log((4*H_k_1)/(2*r));
%         % i == j -> 2
%         elseif(i==2 && j==2  || i==5 && j==5 || i==8 && j==8)
%            P(i,j) = (1/(2*pi*e_0))*log((4*H_k_2)/(2*r));
%         %1 *
%         elseif(i==1 && j==2 || i==2 && j==1 || i==2 && j==3  || i==3 && j==2  || i==4 && j==5  || i==5 && j==4 || i==5 && j==6  || i==6 && j==5 || i==7 && j==8  || i==8 && j==7 || i==8 && j==9  || i==9 && j==8)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_1^2+((H_k_2*2)^2)))/S_kl_1);
%         %2
%         elseif(i==1 && j==3  || i==3 && j==1 || i==4 && j==6 || i==6 && j==4 || i==7 && j==9 || i==9 && j==7)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_2^2+((H_k_1*2)^2)))/S_kl_2);
%         %3
%         elseif(i==3 && j==4 || i==4 && j==3 || i==6 && j==7 || i==7 && j==6)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_3^2+((H_k_1*2)^2)))/S_kl_3);
%         %4 *
%         elseif(i==2 && j==4 || i==4 && j==2 || i==5 && j==7 || i==7 && j==5)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_4^2+((H_k_2*2)^2)))/S_kl_4);
%         %5
%         elseif(i==1 && j==4 || i==4 && j==1 || i==4 && j==7 || i==7 && j==4)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_5^2+((H_k_1*2)^2)))/S_kl_5);
%         %6 = 4 *
%         elseif(i==3 && j==5 || i==5 && j==3 || i==6 && j==8 || i==8 && j==6)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_4^2+((H_k_2*2)^2)))/S_kl_4);
%         %7 = 5 *
%         elseif(i==2 && j==5 || i==5 && j==2 || i==5 && j==8 || i==8 && j==5)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_5^2+((H_k_2*2)^2)))/S_kl_5);
%         %8 *
%         elseif(i==1 && j==5 || i==5 && j==1 || i==4 && j==8 || i==8 && j==4)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_8^2+((H_k_2*2)^2)))/S_kl_8);
%         %9 = 5
%         elseif(i==3 && j==6 || i==6 && j==3 || i==6 && j==9 || i==9 && j==6)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_5^2+((H_k_1*2)^2)))/S_kl_5);
%         %10 = 8 *
%         elseif(i==2 && j==6 || i==6 && j==2 || i==5 && j==9 || i==9 && j==5)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_8^2+((H_k_2*2)^2)))/S_kl_8);
%         %11
%         elseif(i==1 && j==6 || i==6 && j==1 || i==4 && j==9 || i==9 && j==4)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_11^2+((H_k_1*2)^2)))/S_kl_11);
%         %12
%         elseif(i==3 && j==7 || i==7 && j==3)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_12^2+((H_k_1*2)^2)))/S_kl_12);
%         %13 *
%         elseif(i==2 && j==7 || i==7 && j==2)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_13^2+((H_k_2*2)^2)))/S_kl_13);
%         %14
%         elseif(i==1 && j==7 || i==7 && j==1)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_14^2+((H_k_1*2)^2)))/S_kl_14);
%         %15 = 13 *
%         elseif(i==3 && j==8 || i==8 && j==3)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_13^2+((H_k_2*2)^2)))/S_kl_13);
%         %16 = 14 *
%         elseif(i==2 && j==8 || i==8 && j==2)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_14^2+((H_k_2*2)^2)))/S_kl_14);
%         %17 *
%         elseif(i==1 && j==8 || i==8 && j==1)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_17^2+((H_k_2*2)^2)))/S_kl_17);
%         %18 =14
%         elseif(i==3 && j==9 || i==9 && j==3)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_14^2+((H_k_1*2)^2)))/S_kl_14);
%         %19 = 17 *
%         elseif(i==2 && j==9 || i==9 && j==2)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_17^2+((H_k_2*2)^2)))/S_kl_17);
%         %20
%         elseif(i==1 && j==9 || i==9 && j==1)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_20^2+((H_k_1*2)^2)))/S_kl_20);
%         end
%     end
% end

err = abs((y_p - Erms)/y_p)*100;