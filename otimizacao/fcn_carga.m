function [rho] = fcn_carga(V,x,y,r)

%cálculo de densidade de carga para as quatro configurações, condicionante
%depende da entrada do nível de tensão

e_0 = 8.854*(10^(-12)); %permissividade do ar

%tensão das fases
V_ra = V/sqrt(3);
V_ia = 0;

V_rb = V*(cos(2*pi/3))/sqrt(3);
V_ib = 1i*V*(sin(2*pi/3))/sqrt(3);

V_rc = V*(cos(-2*pi/3))/sqrt(3);
V_ic = 1i*V*(sin(-2*pi/3))/sqrt(3);

if size(x,2) == 6 %caso 1  
    Vf = [ V_ra+V_ia ; V_rb+V_ib; V_rc+V_ic ];

elseif size(x,2) == 12 %caso 2
    Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ];

elseif size(x,2) == 18 %caso 3
    Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ];

elseif size(x,2) == 24 %caso 4
    Vf = [ V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_ra+V_ia ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rb+V_ib ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ; V_rc+V_ic ];

end

yi = -y; %imagem de y (altura)
n = size(x,2)/2;

%matriz P
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

rho = P\Vf;