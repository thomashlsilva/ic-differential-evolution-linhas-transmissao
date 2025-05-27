% faixa de passagem = 1.5 m
% norma campo elétrico ao nível do solo = 8.33 kV/m
% ACIMA e ABAIXO 2 m em y da configuração original

clear all; close all; clc

var = 2;
e_0 = 8.854*(10^(-12));
%faixa de passagem
xf = linspace(-25,25,100);
yf = 1.5;

%% ACIMA

n1 = 3; % numero de condutores
r1 = 9.155*(10^-3); % raio do condutor

x1 = [-3 0 3];
y1 = [14.01 14.01 14.01]+var;
yi1 = -y1;

%% ABAIXO

n2 = 3;
r2 = 9.155*(10^-3);

x2 = [-3 0 3];
y2 = [14.01 14.01 14.01]-var;
yi2 = -y2;


%%

P1 = zeros(n1,n1);
for i = 1:n1
     for j = 1:n1
        if(i==j)
           P1(i,j) = (1/(2*pi*e_0))*log((4*y1(i))/(2*r1));
        else
           P1(i,j) = (1/(2*pi*e_0))*log(sqrt((x1(j)-x1(i))^2+(yi1(j)-y1(i))^2)/(sqrt((x1(j)-x1(i))^2+(y1(j)-y1(i))^2)));
        end
     end
end


%%

V1 = 138*10^3;

V1_ra = V1/sqrt(3);
V1_ia = 0;

V1_rb = V1*(cos(2*pi/3))/sqrt(3);
V1_ib = 1i*V1*(sin(2*pi/3))/sqrt(3);

V1_rc = V1*(cos(-2*pi/3))/sqrt(3);
V1_ic = 1i*V1*(sin(-2*pi/3))/sqrt(3);

Vf1 = [ V1_ra+V1_ia ; V1_rb+V1_ib; V1_rc+V1_ic ];

%%

rho1 = P1\Vf1;

rho1_r1 = real(rho1(1));
rho1_i1 = imag(rho1(1));

rho1_r2 = real(rho1(2));
rho1_i2 = imag(rho1(2));

rho1_r3 = real(rho1(3));
rho1_i3 = imag(rho1(3));

%%

E1_xr1 = (rho1_r1/(2*pi*e_0)).*((xf - x1(1))./(((xf - x1(1)).^2) + (yf - y1(1)).^2) - (xf - x1(1))./(((xf - x1(1)).^2) + (yf + y1(1)).^2));
E1_xr2 = (rho1_r2/(2*pi*e_0)).*((xf - x1(2))./(((xf - x1(2)).^2) + (yf - y1(2)).^2) - (xf - x1(2))./(((xf - x1(2)).^2) + (yf + y1(2)).^2));
E1_xr3 = (rho1_r3/(2*pi*e_0)).*((xf - x1(3))./(((xf - x1(3)).^2) + (yf - y1(3)).^2) - (xf - x1(3))./(((xf - x1(3)).^2) + (yf + y1(3)).^2));

E1_xi1 = (rho1_i1/(2*pi*e_0)).*((xf - x1(1))./(((xf - x1(1)).^2) + (yf - y1(1)).^2) - (xf - x1(1))./(((xf - x1(1)).^2) + (yf + y1(1)).^2));
E1_xi2 = (rho1_i2/(2*pi*e_0)).*((xf - x1(2))./(((xf - x1(2)).^2) + (yf - y1(2)).^2) - (xf - x1(2))./(((xf - x1(2)).^2) + (yf + y1(2)).^2));
E1_xi3 = (rho1_i3/(2*pi*e_0)).*((xf - x1(3))./(((xf - x1(3)).^2) + (yf - y1(3)).^2) - (xf - x1(3))./(((xf - x1(3)).^2) + (yf + y1(3)).^2));

E1_yr1 = (rho1_r1/(2*pi*e_0)).*((yf - y1(1))./(((xf - x1(1)).^2) + (yf - y1(1)).^2) - (yf + y1(1))./(((xf - x1(1)).^2) + (yf + y1(1)).^2));
E1_yr2 = (rho1_r2/(2*pi*e_0)).*((yf - y1(2))./(((xf - x1(2)).^2) + (yf - y1(2)).^2) - (yf + y1(2))./(((xf - x1(2)).^2) + (yf + y1(2)).^2));
E1_yr3 = (rho1_r3/(2*pi*e_0)).*((yf - y1(3))./(((xf - x1(3)).^2) + (yf - y1(3)).^2) - (yf + y1(3))./(((xf - x1(3)).^2) + (yf + y1(3)).^2));

E1_yi1 = (rho1_i1/(2*pi*e_0)).*((yf - y1(1))./(((xf - x1(1)).^2) + (yf - y1(1)).^2) - (yf + y1(1))./(((xf - x1(1)).^2) + (yf + y1(1)).^2));
E1_yi2 = (rho1_i2/(2*pi*e_0)).*((yf - y1(2))./(((xf - x1(2)).^2) + (yf - y1(2)).^2) - (yf + y1(2))./(((xf - x1(2)).^2) + (yf + y1(2)).^2));
E1_yi3 = (rho1_i3/(2*pi*e_0)).*((yf - y1(3))./(((xf - x1(3)).^2) + (yf - y1(3)).^2) - (yf + y1(3))./(((xf - x1(3)).^2) + (yf + y1(3)).^2));


%%

E1_xr = [ E1_xr1;  E1_xr2;  E1_xr3];

E1_ximag = [ E1_xi1;  E1_xi2;  E1_xi3];

E1_yr = [ E1_yr1;  E1_yr2;  E1_yr3];

E1_yimag = [ E1_yi1;  E1_yi2;  E1_yi3];


%%

E1xr = sum(E1_xr).^2;

E1ximag = sum(E1_ximag).^2;

E1yr = sum(E1_yr).^2;

E1yimag = sum(E1_yimag).^2;

E1rms = (E1xr+E1ximag+E1yr+E1yimag).^(1/2);

%%

P2 = zeros(n2,n2);
for i = 1:n2
     for j = 1:n2
        if(i==j)
           P2(i,j) = (1/(2*pi*e_0))*log((4*y2(i))/(2*r2));
        else
           P2(i,j) = (1/(2*pi*e_0))*log(sqrt((x2(j)-x2(i))^2+(yi2(j)-y2(i))^2)/(sqrt((x2(j)-x2(i))^2+(y2(j)-y2(i))^2)));
        end
     end
end


%%

V2 = 138*10^3;

%todos os cabos da fase A ( Tensão da Fase A com angulo de 0 grau),
%Todos os cabos da fase B ( Tensão da Fase B com angulo de +120 graus),
%Todos os cabos da Fase C ( Tensão da Fase C com angulo de -120 graus),

V2_ra = V2/sqrt(3);
V2_ia = 0;

V2_rb = V2*(cos(2*pi/3))/sqrt(3);
V2_ib = 1i*V2*(sin(2*pi/3))/sqrt(3);

V2_rc = V2*(cos(-2*pi/3))/sqrt(3);
V2_ic = 1i*V2*(sin(-2*pi/3))/sqrt(3);

Vf2 = [ V2_ra+V2_ia ; V2_rb+V2_ib ; V2_rc+V2_ic ];

%% AP : Funciona destas duas maneiras e a ultima é mais rápida a barra invertida indica inversão de P:
%rho = inv(P)*V;

rho2 = P2\Vf2;

%busco cada uma das posicoes de rho real e imaginario (ver indice)

%cabo 1 fase a
rho2_r1 = real(rho2(1)); 
rho2_i1 = imag(rho2(1));

%cabo 2 fase a
rho2_r2 = real(rho2(2));
rho2_i2 = imag(rho2(2));

%cabo 3 fase b
rho2_r3 = real(rho2(3));
rho2_i3 = imag(rho2(3));


%% 

E2_xr1 = (rho2_r1/(2*pi*e_0)).*((xf - x2(1))./(((xf - x2(1)).^2) + (yf - y2(1)).^2) - (xf - x2(1))./(((xf - x2(1)).^2) + (yf + y2(1)).^2));
E2_xr2 = (rho2_r2/(2*pi*e_0)).*((xf - x2(2))./(((xf - x2(2)).^2) + (yf - y2(2)).^2) - (xf - x2(2))./(((xf - x2(2)).^2) + (yf + y2(2)).^2));
E2_xr3 = (rho2_r3/(2*pi*e_0)).*((xf - x2(3))./(((xf - x2(3)).^2) + (yf - y2(3)).^2) - (xf - x2(3))./(((xf - x2(3)).^2) + (yf + y2(3)).^2));

E2_xi1 = (rho2_i1/(2*pi*e_0)).*((xf - x2(1))./(((xf - x2(1)).^2) + (yf - y2(1)).^2) - (xf - x2(1))./(((xf - x2(1)).^2) + (yf + y2(1)).^2));
E2_xi2 = (rho2_i2/(2*pi*e_0)).*((xf - x2(2))./(((xf - x2(2)).^2) + (yf - y2(2)).^2) - (xf - x2(2))./(((xf - x2(2)).^2) + (yf + y2(2)).^2));
E2_xi3 = (rho2_i3/(2*pi*e_0)).*((xf - x2(3))./(((xf - x2(3)).^2) + (yf - y2(3)).^2) - (xf - x2(3))./(((xf - x2(3)).^2) + (yf + y2(3)).^2));

E2_yr1 = (rho2_r1/(2*pi*e_0)).*((yf - y2(1))./(((xf - x2(1)).^2) + (yf - y2(1)).^2) - (yf + y2(1))./(((xf - x2(1)).^2) + (yf + y2(1)).^2));
E2_yr2 = (rho2_r2/(2*pi*e_0)).*((yf - y2(2))./(((xf - x2(2)).^2) + (yf - y2(2)).^2) - (yf + y2(2))./(((xf - x2(2)).^2) + (yf + y2(2)).^2));
E2_yr3 = (rho2_r3/(2*pi*e_0)).*((yf - y2(3))./(((xf - x2(3)).^2) + (yf - y2(3)).^2) - (yf + y2(3))./(((xf - x2(3)).^2) + (yf + y2(3)).^2));

E2_yi1 = (rho2_i1/(2*pi*e_0)).*((yf - y2(1))./(((xf - x2(1)).^2) + (yf - y2(1)).^2) - (yf + y2(1))./(((xf - x2(1)).^2) + (yf + y2(1)).^2));
E2_yi2 = (rho2_i2/(2*pi*e_0)).*((yf - y2(2))./(((xf - x2(2)).^2) + (yf - y2(2)).^2) - (yf + y2(2))./(((xf - x2(2)).^2) + (yf + y2(2)).^2));
E2_yi3 = (rho2_i3/(2*pi*e_0)).*((yf - y2(3))./(((xf - x2(3)).^2) + (yf - y2(3)).^2) - (yf + y2(3))./(((xf - x2(3)).^2) + (yf + y2(3)).^2));


%% AP: aqui em cada matriz dessa tem as componentes geradas na direção x e y, pela parte real e imaginária de cada carga elétrica:

E2_xr = [ E2_xr1;  E2_xr2 ;  E2_xr3 ];

E2_ximag = [ E2_xi1;  E2_xi2 ;  E2_xi3 ];

E2_yr = [ E2_yr1;  E2_yr2 ;  E2_yr3 ];

E2_yimag = [ E2_yi1;  E2_yi2 ;  E2_yi3 ];


%% AP: aqui somamos as componentes (elementos das colunas) em cada direção e elevamos ao quadrado

E2xr = sum(E2_xr).^2;

E2ximag = sum(E2_ximag).^2;

E2yr = sum(E2_yr).^2;

E2yimag = sum(E2_yimag).^2;

E2rms = (E2xr+E2ximag+E2yr+E2yimag).^(1/2);


%%

Enorma = linspace(8330,8330,100);

E1rms_or = [140.031939621531 145.894621221566 152.035036081088 158.464201122930 165.192890369772 172.231460945521 179.589642645991 187.276285318538 195.299057659792 203.664090425933 212.375556489434 221.435179722773 230.841664415196 240.590036927823 250.670891687262 261.069534562247 271.765018350675 282.729067752343 293.924895078252 305.305913343671 316.814360619363 328.379858877165 339.917942322717 351.328604522582 362.494930497401 373.281899096317 383.535461747088 393.082024959432 401.728484021902 409.262971794682 415.456496344934 420.065640899434 422.836485573882 423.509879433140 421.828141893168 417.543205204892 410.426129990279 400.277845278470 386.940904739468 370.312048945524 350.355484615013 327.117158750224 300.741174845565 271.491463205632 239.786378569131 206.265139646051 171.934580243388 138.524074424894 109.348493894436 90.7916851582907 90.7916851582906 109.348493894435 138.524074424893 171.934580243387 206.265139646051 239.786378569130 271.491463205631 300.741174845565 327.117158750223 350.355484615012 370.312048945524 386.940904739467 400.277845278470 410.426129990278 417.543205204892 421.828141893168 423.509879433139 422.836485573882 420.065640899434 415.456496344934 409.262971794682 401.728484021902 393.082024959431 383.535461747087 373.281899096317 362.494930497401 351.328604522582 339.917942322716 328.379858877164 316.814360619363 305.305913343670 293.924895078252 282.729067752343 271.765018350674 261.069534562246 250.670891687262 240.590036927822 230.841664415196 221.435179722773 212.375556489434 203.664090425932 195.299057659792 187.276285318538 179.589642645990 172.231460945521 165.192890369772 158.464201122929 152.035036081087 145.894621221566 140.031939621531];


h = figure();
plot(xf,E1rms_or,'-','LineWidth',3.5)
xlabel('Faixa de passagem (m)')
ylabel('Campo elétrico (V/m)')
grid
hold on
plot(xf,E1rms,':','LineWidth',3.5)
hold on
plot(xf,E2rms,'-.','LineWidth',3.5)
hold on
plot(xf,Enorma,'--','LineWidth',3.5)
legend('Original','Elevação', 'Descensão', 'Limite da norma')
set(gca, 'YScale', 'log')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'solo_caso1_altura','-dpdf','-r0')

xr1 = [-3 0 3];
yr1 = [14.01 14.01 14.01];

% xr2 = [-10.478 -10.25 -0.228 0 10.022 10.25];
% yr2 = [14.290 14.290 14.290 14.290 14.290 14.290];
% 
% xr3 = [-10.478 -10.25 -10.021 -0.228 0 0.228 10.021 10.25 10.478];
% yr3 = [16.530 16.758 16.530 16.530 16.758 16.530 16.530 16.758 16.530];
% 
% xr4 = [-7.975 -7.025 -7.025 -7.975 -0.475 0.475 0.475 -0.475 7.025 7.975 7.975 7.025];
% yr4 = [18.45 18.45 17.5 17.5 25.95 25.95 25 25 18.45 18.45 17.5 17.5];

g = figure();
plot(xr1,yr1,'o','color','b')
xlabel('Faixa de passagem (m)')
ylabel('Altura dos cabos (m)')
grid
hold on
plot(x1,y1,'*','color','r')
hold on
plot(x2,y2,'x','color','r')
legend('Original','Elevação', 'Descensão')
xlim([-3.5 3.5])
ylim([11.5 16.5])

set(g,'Units','Inches');
pos = get(g,'Position');
set(g,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(g,'solo_caso1_altura_config','-dpdf','-r0')
