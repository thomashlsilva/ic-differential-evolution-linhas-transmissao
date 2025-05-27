function [Ecrit] = fcn_supcrit(r)

%este raio deve entrar em cm
f_s = 0.82 ; %constante
h = 650; %altitude média da LT [m]
t = 50; %temperatura final do condutor da LT [ºC]

densidadeAr = 0.386 * (750-0.086*h)/(273+t);
Ecrit = 18.11 * f_s * densidadeAr * ( 1 + 0.54187/sqrt(r*densidadeAr)); % kV/cm
end