clear all; close all; clc

% 10.478 m entre fase mais distante e centro
% espaçamento 0.228 m
% nível de tensão adotado: 345 kV
% altura 14.290 m acima do solo
% raio do cabo 0.001437 m
% campo elétrico ao nível do solo:
% x=-25m a x = 25m e y = 1m em 100 pontos diferentes

%A seguir será tratado o caso para 2 cabos por fase

n = 6; % número de condutores
e_0 = 8.854*(10^(-12));
r = 14.37*(10^(-3)); % raio do condutor

x = [-10.478 -10.25 -0.228 0 10.022 10.25];
y = [14.290 14.290 14.290 14.290 14.290 14.290];
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

V = 345*10^3;

%todos os cabos da fase A ( Tensão da Fase A com angulo de 0 grau),
%Todos os cabos da fase B ( Tensão da Fase B com angulo de +120 graus),
%Todos os cabos da Fase C ( Tensão da Fase C com angulo de -120 graus),

V_ra = V/sqrt(3);
V_ia = 0;

V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ];

%% AP : Funciona destas duas maneiras e a ultima é mais rápida a barra invertida indica inversão de P:
%rho = inv(P)*V;

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

% x é o vetor de pontos onde se calcula o somatório dos campos
% (contribuição dos cabos de cada uma das fases e das imagens dos cabos

% os limites deste vetor devem ser maiores do que os limites dados pelas
% posições extremas dos cabos:

xf = linspace(-25,25,124);

% especifica a altura da linha de interesse
yf = 1;


%% AP: aqui tinha faltado o . , ele indica que a operação deve ser vetorial ou seja pega o vetor x e subtrai x1,x2,x3 e etc:

%E_xr1 componente x real campo elétrico fase a cabo 1
%E_xr2 componente x real campo elétrico fase a cabo 2
%E_xr3 componente x real campo elétrico fase b cabo 3

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y(1)).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y(1)).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y(1)).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y(1)).^2));

%E_xr4 componente x real campo elétrico fase b cabo 4
%E_xr5 componente x real campo elétrico fase c cabo 5
%E_xr6 componente x real campo elétrico fase c cabo 6

E_xr4 = (rho_r4/(2*pi*e_0)).*((xf - x(4))./(((xf - x(4)).^2) + (yf - y(1)).^2) - (xf - x(4))./(((xf - x(4)).^2) + (yf + y(1)).^2));
E_xr5 = (rho_r5/(2*pi*e_0)).*((xf - x(5))./(((xf - x(5)).^2) + (yf - y(1)).^2) - (xf - x(5))./(((xf - x(5)).^2) + (yf + y(1)).^2));
E_xr6 = (rho_r6/(2*pi*e_0)).*((xf - x(6))./(((xf - x(6)).^2) + (yf - y(1)).^2) - (xf - x(6))./(((xf - x(6)).^2) + (yf + y(1)).^2));

%%

%E_xi1 componente x imaginario campo elétrico fase a cabo 1
%E_xi2 componente x imaginario campo elétrico fase a cabo 2
%E_xi3 componente x imaginario campo elétrico fase b cabo 3

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y(1)).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y(1)).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y(1)).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y(1)).^2));

%E_xi4 componente x imaginario campo elétrico fase b cabo 4
%E_xi5 componente x imaginario campo elétrico fase c cabo 5
%E_xi6 componente x imaginario campo elétrico fase c cabo 6

E_xi4 = (rho_i4/(2*pi*e_0)).*((xf - x(4))./(((xf - x(4)).^2) + (yf - y(1)).^2) - (xf - x(4))./(((xf - x(4)).^2) + (yf + y(1)).^2));
E_xi5 = (rho_i5/(2*pi*e_0)).*((xf - x(5))./(((xf - x(5)).^2) + (yf - y(1)).^2) - (xf - x(5))./(((xf - x(5)).^2) + (yf + y(1)).^2));
E_xi6 = (rho_i6/(2*pi*e_0)).*((xf - x(6))./(((xf - x(6)).^2) + (yf - y(1)).^2) - (xf - x(6))./(((xf - x(6)).^2) + (yf + y(1)).^2));

%%

%E_yr1 componente y real campo elétrico fase a cabo 1
%E_yr2 componente y real campo elétrico fase a cabo 2
%E_yr3 componente y real campo elétrico fase b cabo 3

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y(1))./(((xf - x(2)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(2)).^2) + (yf + y(1)).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y(1))./(((xf - x(3)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(3)).^2) + (yf + y(1)).^2));

%E_yr4 componente y real campo elétrico fase b cabo 4
%E_yr5 componente y real campo elétrico fase c cabo 5
%E_yr6 componente y real campo elétrico fase c cabo 6

E_yr4 = (rho_r4/(2*pi*e_0)).*((yf - y(1))./(((xf - x(4)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(4)).^2) + (yf + y(1)).^2));
E_yr5 = (rho_r5/(2*pi*e_0)).*((yf - y(1))./(((xf - x(5)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(5)).^2) + (yf + y(1)).^2));
E_yr6 = (rho_r6/(2*pi*e_0)).*((yf - y(1))./(((xf - x(6)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(6)).^2) + (yf + y(1)).^2));

%%

%E_yi1 componente y imaginario campo elétrico fase a cabo 1
%E_yi2 componente y imaginario campo elétrico fase a cabo 2
%E_yi3 componente y imaginario campo elétrico fase b cabo 3

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y(1))./(((xf - x(1)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(1)).^2) + (yf + y(1)).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y(1))./(((xf - x(2)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(2)).^2) + (yf + y(1)).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y(1))./(((xf - x(3)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(3)).^2) + (yf + y(1)).^2));

%E_yi4 componente y imaginario campo elétrico fase b cabo 4
%E_yi5 componente y imaginario campo elétrico fase c cabo 5
%E_yi6 componente y imaginario campo elétrico fase c cabo 6

E_yi4 = (rho_i4/(2*pi*e_0)).*((yf - y(1))./(((xf - x(4)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(4)).^2) + (yf + y(1)).^2));
E_yi5 = (rho_i5/(2*pi*e_0)).*((yf - y(1))./(((xf - x(5)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(5)).^2) + (yf + y(1)).^2));
E_yi6 = (rho_i6/(2*pi*e_0)).*((yf - y(1))./(((xf - x(6)).^2) + (yf - y(1)).^2) - (yf + y(1))./(((xf - x(6)).^2) + (yf + y(1)).^2));


%% AP: aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr1;  E_xr2 ;  E_xr3 ; E_xr4 ; E_xr5 ; E_xr6 ];

E_ximag = [ E_xi1;  E_xi2 ;  E_xi3 ; E_xi4 ; E_xi5 ; E_xi6 ];

E_yr = [ E_yr1;  E_yr2 ;  E_yr3 ; E_yr4 ; E_yr5 ; E_yr6 ];

E_yimag = [ E_yi1;  E_yi2 ;  E_yi3 ; E_yi4 ; E_yi5 ; E_yi6 ];


%% AP: aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

Exr = sum(E_xr).^2;

Eximag = sum(E_ximag).^2;

Eyr = sum(E_yr).^2;

Eyimag = sum(E_yimag).^2;

%% AP: aqui retira-se o módulo dos quadrados:

Erms = (Exr+Eximag+Eyr+Eyimag).^(1/2);

%% Tem outro erro aqui tb - é tudo soma:

x_p = [-25.00	,	-24.77	,	-24.52	,	-24.10	,	-23.76	,	-23.43	,	-23.08	,	-22.74	,	-22.45	,	-22.07	,	-21.76	,	-21.46	,	-21.15	,	-20.86	,	-20.53	,	-20.27	,	-19.90	,	-19.52	,	-19.20	,	-18.93	,	-18.57	,	-18.14	,	-17.81	,	-17.49	,	-17.11	,	-16.76	,	-16.34	,	-16.00	,	-15.14	,	-14.70	,	-14.21	,	-13.71	,	-13.25	,	-12.77	,	-12.32	,	-11.56	,	-11.03	,	-10.63	,	-10.21	,	-9.82	,	-9.48	,	-9.19	,	-8.79	,	-8.38	,	-8.10	,	-7.79	,	-7.48	,	-7.22	,	-6.89	,	-6.61	,	-6.31	,	-5.96	,	-5.73	,	-5.43	,	-5.07	,	-4.78	,	-4.41	,	-4.03	,	-3.61	,	-3.16	,	-2.66	,	-2.10	,	-1.52	,	-0.92	,	-0.40	,	0.28	,	0.93	,	1.62	,	2.23	,	2.78	,	3.21	,	3.59	,	4.01	,	4.49	,	4.79	,	5.20	,	5.61	,	5.85	,	6.23	,	6.55	,	6.89	,	7.15	,	7.45	,	7.79	,	8.22	,	8.54	,	8.97	,	9.49	,	9.98	,	10.40	,	11.01	,	11.49	,	12.09	,	12.58	,	13.06	,	13.62	,	14.18	,	14.73	,	15.28	,	15.77	,	16.16	,	16.56	,	16.93	,	17.42	,	17.75	,	18.06	,	18.44	,	18.86	,	19.11	,	19.44	,	19.89	,	20.29	,	20.59	,	20.95	,	21.41	,	21.71	,	22.04	,	22.49	,	22.86	,	23.17	,	23.55	,	24.04	,	24.52	,	25.00];

y_p = [1574.63	,	1611.44	,	1653.23	,	1705.97	,	1743.78	,	1787.56	,	1842.29	,	1890.05	,	1930.85	,	1978.61	,	2021.39	,	2067.16	,	2109.95	,	2154.73	,	2213.43	,	2258.21	,	2313.93	,	2366.67	,	2412.44	,	2461.19	,	2524.88	,	2577.61	,	2619.40	,	2676.12	,	2721.89	,	2765.67	,	2817.41	,	2866.17	,	2946.77	,	2982.59	,	3015.42	,	3039.30	,	3054.23	,	3059.20	,	3051.24	,	3022.39	,	2989.55	,	2949.75	,	2906.97	,	2859.20	,	2810.45	,	2761.69	,	2708.96	,	2634.33	,	2589.55	,	2539.80	,	2476.12	,	2419.40	,	2354.73	,	2307.96	,	2258.21	,	2195.52	,	2148.76	,	2097.01	,	2047.26	,	2002.49	,	1948.76	,	1907.96	,	1867.16	,	1824.38	,	1789.55	,	1761.69	,	1741.79	,	1731.84	,	1724.88	,	1726.87	,	1739.80	,	1757.71	,	1786.57	,	1825.37	,	1868.16	,	1901.00	,	1942.79	,	1995.52	,	2046.27	,	2097.01	,	2155.72	,	2214.43	,	2290.05	,	2342.79	,	2403.48	,	2443.28	,	2511.94	,	2577.61	,	2643.28	,	2688.06	,	2751.74	,	2825.37	,	2891.04	,	2945.77	,	2993.53	,	3018.41	,	3038.31	,	3045.27	,	3037.31	,	3015.42	,	2982.59	,	2939.80	,	2891.04	,	2844.28	,	2799.50	,	2751.74	,	2692.04	,	2630.35	,	2585.57	,	2540.80	,	2481.09	,	2418.41	,	2370.65	,	2325.87	,	2263.18	,	2202.49	,	2154.73	,	2096.02	,	2035.32	,	1986.57	,	1939.80	,	1879.10	,	1822.39	,	1777.61	,	1725.87	,	1665.17	,	1603.48	,	1547.76];

figure();
plot(xf,Erms,'LineWidth',4.5)
xlabel('Faixa de passagem (m)')
ylabel('Campo elétrico (V/m)')
%title('Campo elétrico ao nível do solo dois condutores por fase')
grid
hold on
yy2 = smooth(x_p,y_p,'rloess');
plot(x_p,yy2,'--','LineWidth',2.5);
% plot(x_p,y_p,'--','LineWidth',2.5);
xlim([-25 25])
ylim([500 3500])
legend('Simulado', '(PAGANOTTI, 2012)')

h = figure();
plot(x(1:2), y(1:2),'o', 'color', 'blue')
hold on
plot(x(3:4), y(3:4),'*', 'color', 'blue')
hold on
plot(x(5:6), y(5:6),'x', 'color', 'blue')
xlabel('Faixa de passagem (m)')
ylabel('Altura do cabo (m)')
%title('Configuração geométrica dos condutores')
legend('Cabos condutores fase A','Cabos condutores fase B','Cabos condutores fase C');
grid
xlim([-15 15])
ylim([13 15.5])

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'caso2config','-dpdf','-r0')

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
% plot([x1 x2 x3],[y1 y2 y3],'*')
% legend('Cabos condutores equivalentes por fase');
% xlabel('Faixa de Passagem (m)')
% ylabel('Altura do cabo (m)')
% grid

% S_kl_1 = 0.228; % 11,21,34,43,56,65 -> espaçamento
% S_kl_4 = 10.478; % 14,41,36,63 -> cabo mais distante do cabo em x = 0
% S_kl_2 = S_kl_4 - 2*S_kl_1; % 23,32,45,54
% S_kl_3 = S_kl_1 + S_kl_2; % 13,31,24,42,35,53,46,64
% S_kl_5 = S_kl_3 + S_kl_2; % 25,52
% S_kl_6 = S_kl_4 + S_kl_2; % 15,51,26,62
% S_kl_7 = S_kl_4 + S_kl_3; % 16,61

% P = zeros(n,n);
% 
% for i = 1:n
%      for j = 1:n
%         if(i==j)
%            P(i,j) = (1/(2*pi*e_0))*log((4*y)/(2*r));   
%         elseif(i==1 && j==2 || i==2 && j==1 || i==3 && j==4  || i==4 && j==3  || i==5 && j==6  || i==6 && j==5)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_1^2+((y*2)^2)))/S_kl_1);
%         elseif(i==2 && j==3  || i==3 && j==2 || i==4 && j==5 || i==5 && j==4)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_2^2+((y*2)^2)))/S_kl_2);
%         elseif(i==1 && j==3 || i==3 && j==1 || i==2 && j==4 || i==4 && j==2 || i==3 && j==5 || i==5 && j==3 || i==4 && j==6 || i==6 && j==4)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_3^2+((y*2)^2)))/S_kl_3);
%         elseif(i==1 && j==4 || i==4 && j==1 || i==3 && j==6 || i==6 && j==3)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_4^2+((y*2)^2)))/S_kl_4);
%         elseif(i==2 && j==5 || i==5 && j==2)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_5^2+((y*2)^2)))/S_kl_5);
%         elseif(i==1 && j==5 || i==5 && j==1 || i==2 && j==6 || i==6 && j==2)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_6^2+((y*2)^2)))/S_kl_6);
%         elseif(i==1 && j==6 || i==6 && j==1)
%            P(i,j) = (1/(2*pi*e_0))*log((sqrt(S_kl_7^2+((y*2)^2)))/S_kl_7);
%         end
%     end
% end

err = abs((yy2.' - Erms)/yy2.')*100;