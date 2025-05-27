clear all; close all; clc

% Example of lateral profile of electric field at
% ground. 525-kV line with flat configuration, 10-m phase
% spacing, 10.6-m height above ground, regular
% 3-conductor bundles, 3.3-cm diameter, 45-cm spacing.

n = 3; % número de subcondutores = número de fases
e_0 = 8.854*(10^(-12));
s = 0.45; % espaçamento
d = 0.033; % diâmetro de um condutor
d_b = s/(sin(pi/n)); % diâmetro do feixe
d_eq = d_b*(((n*d)/d_b))^(1/n); % diâmetro equivalente

x = [-10 0 10]; % posição em x dos cabos equivalentes
y = 10.6; % posição em y dos cabos equivalentes
yi = -y; % posição da imagem de y

P = zeros(n,n); % matriz P
for i = 1:n
     for j = 1:n
        if(i==j)
           P(i,j) = (1/(2*pi*e_0))*log((4*y)/(d_eq));
        else
           P(i,j) = (1/(2*pi*e_0))*log(sqrt((x(j)-x(i))^2+(yi-y)^2)/(sqrt((x(j)-x(i))^2)));
        end
     end
end


%%

V = 500*10^3; % tensão adotada, se colocar 505kv fica certo

V_ra = V/sqrt(3);
V_ia = 0;

V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

Vf = [ V_ra+V_ia ; V_rb+V_ib ; V_rc+V_ic ];

%%

rho = P\Vf;

%cabo eq 1 fase a
rho_r1 = real(rho(1)); 
rho_i1 = imag(rho(1));

%cabo eq 2 fase b
rho_r2 = real(rho(2));
rho_i2 = imag(rho(2));

%cabo eq 3 fase c
rho_r3 = real(rho(3));
rho_i3 = imag(rho(3));

% posição da faixa de passagem:

xf = linspace(0,60,163); 

yf = 0;

%% AP: aqui tinha faltado o . , ele indica que a operação deve ser vetorial ou seja pega o vetor x e subtrai x1,x(2),x(3) e etc:

%E_xr1 componente x real campo elétrico fase a cabo 1
%E_xr2 componente x real campo elétrico fase b cabo 2
%E_xr3 componente x real campo elétrico fase c cabo 3

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y).^2));

%%

%E_xi1 componente x imaginario campo elétrico fase a cabo 1
%E_xi2 componente x imaginario campo elétrico fase b cabo 2
%E_xi3 componente x imaginario campo elétrico fase c cabo 3

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y).^2));

%%

%E_yr1 componente y real campo elétrico fase a cabo 1
%E_yr2 componente y real campo elétrico fase b cabo 2
%E_yr3 componente y real campo elétrico fase c cabo 3

E_yr1 = (rho_r1/(2*pi*e_0))*((yf - y)./(((xf - x(1)).^2) + (yf - y)^2) - (yf + y)./(((xf - x(1)).^2) + (yf + y)^2));
E_yr2 = (rho_r2/(2*pi*e_0))*((yf - y)./(((xf - x(2)).^2) + (yf - y)^2) - (yf + y)./(((xf - x(2)).^2) + (yf + y)^2));
E_yr3 = (rho_r3/(2*pi*e_0))*((yf - y)./(((xf - x(3)).^2) + (yf - y)^2) - (yf + y)./(((xf - x(3)).^2) + (yf + y)^2));

%%

%E_yi1 componente y imaginario campo elétrico fase a cabo 1
%E_yi2 componente y imaginario campo elétrico fase b cabo 2
%E_yi3 componente y imaginario campo elétrico fase c cabo 3

E_yi1 = (rho_i1/(2*pi*e_0))*((yf - y)./(((xf - x(1)).^2) + (yf - y)^2) - (yf + y)./(((xf - x(1)).^2) + (yf + y)^2));
E_yi2 = (rho_i2/(2*pi*e_0))*((yf - y)./(((xf - x(2)).^2) + (yf - y)^2) - (yf + y)./(((xf - x(2)).^2) + (yf + y)^2));
E_yi3 = (rho_i3/(2*pi*e_0))*((yf - y)./(((xf - x(3)).^2) + (yf - y)^2) - (yf + y)./(((xf - x(3)).^2) + (yf + y)^2));


%% AP: aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr1; E_xr2;  E_xr3];

E_ximag = [ E_xi1; E_xi2;  E_xi3];

E_yr = [ E_yr1; E_yr2;  E_yr3];

E_yimag = [ E_yi1; E_yi2;  E_yi3];


%% AP: aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

Exr = (E_xr1 + E_xr2 + E_xr3).^2;

Eximag = (E_xi1 + E_xi2 + E_xi3).^2;

Eyr = (E_yr1 + E_yr2 + E_yr3).^2;

Eyimag = (E_yi1 + E_yi2 + E_yi3).^2;

%% AP: aqui retira-se o módulo dos quadrados:

Erms = (Exr+Eximag+Eyr+Eyimag).^(1/2);

%%

%vetores para a curva do livro

x_p = [0.06	,	0.30	,	0.56	,	0.85	,	1.12	,	1.41	,	1.73	,	2.06	,	2.39	,	2.73	,	3.12	,	3.48	,	3.89	,	4.25	,	4.57	,	4.93	,	5.29	,	5.60	,	5.90	,	6.18	,	6.46	,	6.79	,	7.06	,	7.35	,	7.65	,	8.04	,	8.35	,	8.65	,	8.94	,	9.26	,	9.62	,	10.03	,	10.45	,	10.91	,	11.31	,	11.74	,	12.18	,	12.56	,	12.99	,	13.32	,	13.68	,	14.04	,	14.44	,	14.86	,	15.15	,	15.54	,	15.89	,	16.21	,	16.57	,	16.87	,	17.09	,	17.38	,	17.68	,	18.03	,	18.34	,	18.68	,	19.06	,	19.29	,	19.59	,	19.92	,	20.26	,	20.56	,	20.94	,	21.24	,	21.64	,	22.03	,	22.38	,	22.71	,	23.02	,	23.30	,	23.60	,	23.85	,	24.15	,	24.38	,	24.71	,	25.03	,	25.38	,	25.69	,	25.99	,	26.35	,	26.71	,	27.00	,	27.37	,	27.70	,	27.97	,	28.31	,	28.73	,	29.09	,	29.54	,	29.91	,	30.29	,	30.62	,	31.03	,	31.40	,	31.75	,	32.15	,	32.45	,	32.84	,	33.20	,	33.56	,	33.97	,	34.32	,	34.71	,	35.03	,	35.41	,	35.76	,	36.12	,	36.56	,	36.93	,	37.32	,	37.77	,	38.15	,	38.63	,	39.02	,	39.40	,	39.83	,	40.22	,	40.62	,	41.03	,	41.38	,	41.78	,	42.14	,	42.54	,	43.03	,	43.31	,	43.74	,	44.09	,	44.46	,	44.85	,	45.29	,	45.68	,	46.11	,	46.55	,	46.90	,	47.31	,	47.63	,	48.05	,	48.55	,	49.07	,	49.48	,	49.93	,	50.38	,	50.86	,	51.36	,	51.84	,	52.29	,	52.76	,	53.25	,	53.74	,	54.13	,	54.58	,	55.00	,	55.48	,	55.97	,	56.46	,	56.97	,	57.52	,	58.02	,	58.51	,	58.94	,	59.30	,	59.64	,	59.89];
y_p = [5874.58	,	5874.58	,	5868.63	,	5862.69	,	5844.89	,	5785.97	,	5716.05	,	5629.83	,	5561.79	,	5516.89	,	5483.45	,	5472.35	,	5477.90	,	5539.29	,	5612.74	,	5727.64	,	5821.25	,	5958.51	,	6092.82	,	6223.85	,	6409.45	,	6600.58	,	6776.78	,	6943.59	,	7150.65	,	7371.35	,	7568.12	,	7754.41	,	7913.16	,	8066.97	,	8215.45	,	8307.53	,	8392.13	,	8434.76	,	8451.87	,	8443.31	,	8409.16	,	8349.72	,	8190.52	,	8058.80	,	7913.16	,	7754.41	,	7614.27	,	7416.30	,	7260.16	,	7049.93	,	6866.63	,	6674.56	,	6481.29	,	6299.99	,	6161.11	,	5988.77	,	5815.36	,	5641.25	,	5461.27	,	5308.51	,	5118.36	,	4980.23	,	4806.71	,	4639.24	,	4500.34	,	4352.35	,	4196.45	,	4070.81	,	3897.26	,	3716.01	,	3590.17	,	3472.11	,	3357.94	,	3267.32	,	3172.71	,	3083.96	,	3000.74	,	2919.76	,	2820.88	,	2739.20	,	2649.12	,	2556.82	,	2477.76	,	2396.28	,	2317.48	,	2239.00	,	2160.99	,	2089.93	,	2017.11	,	1952.76	,	1882.81	,	1824.59	,	1755.68	,	1694.51	,	1630.50	,	1572.10	,	1512.72	,	1471.89	,	1429.27	,	1386.48	,	1339.53	,	1299.43	,	1259.24	,	1219.07	,	1187.37	,	1147.17	,	1111.69	,	1089.39	,	1053.57	,	1024.10	,	995.45	,	961.74	,	932.00	,	905.93	,	876.14	,	849.05	,	816.98	,	794.93	,	768.79	,	749.56	,	727.86	,	707.50	,	684.93	,	670.51	,	650.43	,	633.52	,	617.67	,	594.95	,	581.24	,	566.70	,	553.09	,	539.80	,	526.83	,	512.10	,	498.28	,	485.32	,	471.27	,	463.22	,	450.26	,	440.78	,	431.07	,	416.89	,	404.00	,	394.69	,	386.38	,	375.96	,	362.86	,	352.35	,	344.58	,	333.93	,	325.25	,	316.79	,	307.93	,	302.36	,	293.90	,	287.43	,	279.67	,	271.85	,	267.20	,	259.20	,	252.21	,	245.65	,	239.51	,	234.46	,	230.69	,	226.07	,	223.79];

figure();
plot(xf,Erms,'LineWidth',4.5)
xlabel('Faixa de passagem (m)')
ylabel('Campo elétrico (V/m)')
%title('Campo elétrico ao nível do solo um condutor equivalente por fase')
grid
hold on
% yy2 = smooth(x_p,y_p,'rloess');
% plot(x_p,yy2,'--','LineWidth',2.5);
plot(x_p,y_p,'--','LineWidth',2.5);
xlim([0 60])
ylim([100 10000])
legend('Simulado', '(EPRI, 2005)')

h = figure();
plot(x(1), y,'o', 'color', 'blue')
hold on
plot(x(2), y,'*', 'color', 'blue')
hold on
plot(x(3), y,'x', 'color', 'blue')
xlabel('Faixa de passagem (m)')
ylabel('Altura do cabo (m)')
%title('Configuração geométrica dos condutores')
legend('Cabo condutor fase A','Cabo condutor fase B','Cabo condutor fase C');
%title('Configuração geométrica dos condutores')
grid
xlim([-10.5 10.5])
ylim([9 12])

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'caso5config','-dpdf','-r0')

err = abs((y_p - Erms)/y_p)*100;