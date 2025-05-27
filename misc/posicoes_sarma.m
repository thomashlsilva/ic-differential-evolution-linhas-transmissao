clear all; close all; clc

% medidas em metros
h1 = 19.812; % altura do solo até o centro do feixe
c1 = 15.24; % distância entre a origem e o centro do feixe
e1 = 0.4572; % espaçamento entre condutores
e2 = e1/2; 
d = 0.035052; % diâmetro do condutor

x1 = -(c1+e2);
y1 = h1-e2;

x2 = -(c1+e2);
y2 = h1+e2;

x3 = -(c1-e2);
y3 = h1-e2;

x4 = -(c1-e2);
y4 = h1+e2;

x5 = -e2;
y5 = h1-e2;

x6 = -e2;
y6 = h1+e2;

x7 = e2;
y7 = h1-e2;

x8 = e2;
y8 = h1+e2;

x9 = c1-e2;
y9 = h1-e2;

x10 = c1-e2;
y10 = h1+e2;

x11 = c1+e2;
y11 = h1-e2;

x12 = c1+e2;
y12 = h1+e2;

x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12];
y = [y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12];

yi = -y;
