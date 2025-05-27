clear all; close all; clc

% 3 m entre fases
% nível de tensão adotado: 138 kV
% altura 14.01 m acima do solo
% raio do cabo 0.001437 m
% campo elétrico ao nível do solo:
% x=-25m a x = 25m e y = 1m em 100 pontos diferentes

%A seguir será tratado o caso para 1 cabo por fase

n = 3; % numero de condutores
e_0 = 8.854*(10^(-12));
r = 9.155*(10^-3); % raio do condutor

x = [-3 0 3];
y = 14.01;
yi = -y;

P = zeros(n,n);
for i = 1:n
     for j = 1:n
        if(i==j)
           P(i,j) = (1/(2*pi*e_0))*log((4*y)/(2*r));
        else
           P(i,j) = (1/(2*pi*e_0))*log(sqrt((x(j)-x(i))^2+(yi-y)^2)/(sqrt((x(j)-x(i))^2)));
        end
     end
end

%%

V = 138*10^3;

%todos os cabos da fase A ( Tensão da Fase A com angulo de 0 grau),
%Todos os cabos da fase B ( Tensão da Fase B com angulo de +120 graus),
%Todos os cabos da Fase C ( Tensão da Fase C com angulo de -120 graus),

V_ra = V/sqrt(3);
V_ia = 0;

V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

V = [ V_ra+V_ia ; V_rb+V_ib; V_rc+V_ic ];

%% AP : Funciona destas duas maneiras e a ultima é mais rápida a barra invertida indica inversão de P:
%rho = inv(P)*V;

rho = P\V;

%1 = fase a, 2 = fase b, 3 = fase c

%busco cada uma das posicoes de rho real e imaginario (ver indice)

rho_r1 = real(rho(1));
rho_i1 = imag(rho(1));

rho_r2 = real(rho(2));
rho_i2 = imag(rho(2));

rho_r3 = real(rho(3));
rho_i3 = imag(rho(3));

% xf é o vetor de pontos onde se calcula o somatório dos campos
% (contribuição dos cabos de cada uma das fases e das imagens dos cabos

% os limites deste vetor devem ser maiores do que os limites dados pelas
% posições extremas dos cabos:

xf = linspace(-25,25,174);

% especifica a altura da linha de interesse
yf = 1;


%% AP: aqui tinha faltado o . , ele indica que a operação deve ser vetorial ou seja pega o vetor x e subtrai x1,x2,x3 e etc:

%E_xr1 componente x campo elétrico fase a
%E_xr2 componente x campo elétrico fase b
%E_xr3 componente x campo elétrico fase c

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y).^2));

%E_xi1 componente x imaginario campo elétrico fase a
%E_xi2 componente x imaginario campo elétrico fase b
%E_xi3 componente x imaginario campo elétrico fase c

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1))./(((xf - x(1)).^2) + (yf - y).^2) - (xf - x(1))./(((xf - x(1)).^2) + (yf + y).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(2))./(((xf - x(2)).^2) + (yf - y).^2) - (xf - x(2))./(((xf - x(2)).^2) + (yf + y).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(3))./(((xf - x(3)).^2) + (yf - y).^2) - (xf - x(3))./(((xf - x(3)).^2) + (yf + y).^2));

%E_yr1 componente y real campo elétrico fase a
%E_yr2 componente y real campo elétrico fase b
%E_yr3 componente y real campo elétrico fase c

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y)./(((xf - x(1)).^2) + (yf - y).^2) - (yf + y)./(((xf - x(1)).^2) + (yf + y).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y)./(((xf - x(2)).^2) + (yf - y).^2) - (yf + y)./(((xf - x(2)).^2) + (yf + y).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y)./(((xf - x(3)).^2) + (yf - y).^2) - (yf + y)./(((xf - x(3)).^2) + (yf + y).^2));

%E_yi1 componente y imaginario campo elétrico fase a
%E_yi2 componente y imaginario campo elétrico fase b
%E_yi3 componente y imaginario campo elétrico fase c

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y)./(((xf - x(1)).^2) + (yf - y).^2) - (yf + y)./(((xf - x(1)).^2) + (yf + y).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y)./(((xf - x(2)).^2) + (yf - y).^2) - (yf + y)./(((xf - x(2)).^2) + (yf + y).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y)./(((xf - x(3)).^2) + (yf - y).^2) - (yf + y)./(((xf - x(3)).^2) + (yf + y).^2));


%% AP: aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E_xr = [ E_xr1;  E_xr2;  E_xr3];

E_ximag = [ E_xi1;  E_xi2;  E_xi3];

E_yr = [ E_yr1;  E_yr2;  E_yr3];

E_yimag = [ E_yi1;  E_yi2;  E_yi3];


%% AP: aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

Exr = sum(E_xr).^2;

Eximag = sum(E_ximag).^2;

Eyr = sum(E_yr).^2;

Eyimag = sum(E_yimag).^2;

%% AP: aqui retira-se o módulo dos quadrados:

Erms = (Exr+Eximag+Eyr+Eyimag).^(1/2);


%% Tem outro erro aqui tb - é tudo soma:

x_p = [-25.00	,	-24.68	,	-24.39	,	-24.08	,	-23.88	,	-23.62	,	-23.43	,	-23.14	,	-22.92	,	-22.72	,	-22.48	,	-22.20	,	-21.88	,	-21.56	,	-21.28	,	-20.96	,	-20.76	,	-20.45	,	-20.26	,	-19.99	,	-19.71	,	-19.43	,	-19.07	,	-18.93	,	-18.57	,	-18.27	,	-17.89	,	-17.64	,	-17.39	,	-17.09	,	-16.83	,	-16.46	,	-16.27	,	-16.02	,	-15.81	,	-15.58	,	-15.28	,	-15.02	,	-14.77	,	-14.45	,	-14.22	,	-13.90	,	-13.63	,	-13.42	,	-13.15	,	-12.98	,	-12.72	,	-12.46	,	-12.18	,	-11.93	,	-11.49	,	-11.14	,	-10.70	,	-10.35	,	-9.82	,	-9.33	,	-8.80	,	-8.35	,	-7.81	,	-7.36	,	-6.97	,	-6.75	,	-6.36	,	-5.85	,	-5.57	,	-5.29	,	-5.05	,	-4.74	,	-4.57	,	-4.37	,	-4.18	,	-3.98	,	-3.78	,	-3.59	,	-3.39	,	-3.24	,	-3.10	,	-2.97	,	-2.82	,	-2.55	,	-2.34	,	-2.15	,	-2.06	,	-1.87	,	-1.73	,	-1.53	,	-1.38	,	-1.22	,	-0.99	,	-0.84	,	-0.45	,	-0.19	,	0.20	,	0.46	,	0.73	,	0.98	,	1.17	,	1.26	,	1.38	,	1.54	,	1.66	,	1.85	,	2.01	,	2.13	,	2.26	,	2.34	,	2.52	,	2.64	,	2.78	,	2.92	,	3.06	,	3.15	,	3.29	,	3.43	,	3.57	,	3.67	,	3.82	,	4.03	,	4.23	,	4.38	,	4.59	,	4.80	,	5.02	,	5.39	,	5.67	,	5.96	,	6.24	,	6.63	,	6.92	,	7.40	,	7.94	,	8.36	,	8.87	,	9.53	,	10.08	,	10.54	,	11.02	,	11.39	,	11.81	,	12.13	,	12.46	,	12.83	,	13.19	,	13.46	,	13.75	,	14.09	,	14.36	,	14.59	,	14.95	,	15.27	,	15.53	,	15.88	,	16.14	,	16.44	,	16.82	,	17.20	,	17.57	,	17.97	,	18.41	,	18.81	,	19.11	,	19.50	,	19.99	,	20.50	,	20.96	,	21.51	,	21.98	,	22.54	,	23.06	,	23.46	,	23.96	,	24.35	,	24.72	,	25.00];
y_p = [140.61	,	144.06	,	147.51	,	152.55	,	154.41	,	157.06	,	159.19	,	163.17	,	166.22	,	170.07	,	173.12	,	176.83	,	180.55	,	186.12	,	190.10	,	194.21	,	198.06	,	202.84	,	206.55	,	211.99	,	216.77	,	222.47	,	227.91	,	231.49	,	237.86	,	243.03	,	250.60	,	255.24	,	261.08	,	267.45	,	272.89	,	279.65	,	284.96	,	290.40	,	295.31	,	301.94	,	307.38	,	315.07	,	320.78	,	327.01	,	331.39	,	338.69	,	344.66	,	350.90	,	355.67	,	361.51	,	365.62	,	372.65	,	378.23	,	385.52	,	391.63	,	396.80	,	403.43	,	410.46	,	415.77	,	420.81	,	424.93	,	426.12	,	422.27	,	417.36	,	412.45	,	408.08	,	401.44	,	389.10	,	379.95	,	370.66	,	359.25	,	348.51	,	340.15	,	331.79	,	320.91	,	312.69	,	300.75	,	289.87	,	281.91	,	271.69	,	260.02	,	253.38	,	240.91	,	225.66	,	218.23	,	206.29	,	193.68	,	185.59	,	172.99	,	160.38	,	151.49	,	139.02	,	124.03	,	109.83	,	98.96	,	90.33	,	90.20	,	98.96	,	109.70	,	121.24	,	129.34	,	139.55	,	151.49	,	159.72	,	168.34	,	175.90	,	186.12	,	195.80	,	205.09	,	215.84	,	222.74	,	231.23	,	238.92	,	250.46	,	259.49	,	266.25	,	274.08	,	284.03	,	290.80	,	298.36	,	305.79	,	313.48	,	320.51	,	328.08	,	338.96	,	346.65	,	358.06	,	371.19	,	384.86	,	393.88	,	401.18	,	407.55	,	413.12	,	418.69	,	422.80	,	425.85	,	422.80	,	418.56	,	412.72	,	406.88	,	399.05	,	393.35	,	386.32	,	379.55	,	372.26	,	362.84	,	355.27	,	348.24	,	342.94	,	336.97	,	329.00	,	323.43	,	317.06	,	308.71	,	302.60	,	295.97	,	290.40	,	282.70	,	274.08	,	266.78	,	259.49	,	250.33	,	241.58	,	234.41	,	229.24	,	220.61	,	213.18	,	203.90	,	195.94	,	187.84	,	180.15	,	172.45	,	165.29	,	158.52	,	153.08	,	148.18	,	144.06	,	140.75];


figure();
plot(xf,Erms,'LineWidth',4.5)
xlabel('Faixa de passagem (m)')
ylabel('Campo elétrico (V/m)')
%title('Campo elétrico ao nível do solo um condutor por fase')
grid
hold on
yy2 = smooth(x_p,y_p,'rloess');
plot(x_p,yy2,'--','LineWidth',2.5);
% plot(x_p,y_p,'--','LineWidth',2.5);
xlim([-25 25])
ylim([50 450])
legend('Simulado', '(PAGANOTTI, 2012)')
set(gca,'YTick', 25:50:475, 'YTickLabel', 25:50:475);

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

grid
xlim([-3.5 3.5])
ylim([13 14.5])

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'caso1config','-dpdf','-r0')


err = abs((yy2.' - Erms)/yy2.')*100;