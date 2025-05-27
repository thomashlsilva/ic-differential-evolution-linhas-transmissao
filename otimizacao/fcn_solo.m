function [Erms] = fcn_solo(x,V,r)

%nota: este x entra com valores de x,y respectivamente na forma 1,:

e_0 = 8.854*(10^(-12)); %permissividade do ar

%linha de interesse
xf = linspace(-25,25,100);
%especifica a altura da linha de interesse
yf = 1.5;

y = x(1,(size(x,2)/2)+1:size(x,2));

%% Condicionante para qual configuração será realizado o cálculo do campo elétrico ao nível do solo 

if size(x,2) == 6 %caso 1

rho = fcn_carga(V,x,y,r); %entram valores para cálculo da densidade de carga
%1 = fase a, 2 = fase b, 3 = fase c
rho_r1 = real(rho(1));
rho_i1 = imag(rho(1));

rho_r2 = real(rho(2));
rho_i2 = imag(rho(2));

rho_r3 = real(rho(3));
rho_i3 = imag(rho(3));

%E_xr1 componente x campo elétrico fase a
%E_xr2 componente x campo elétrico fase b
%E_xr3 componente x campo elétrico fase c

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));

%E_xi1 componente x imaginario campo elétrico fase a
%E_xi2 componente x imaginario campo elétrico fase b
%E_xi3 componente x imaginario campo elétrico fase c

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));

%E_yr1 componente y real campo elétrico fase a
%E_yr2 componente y real campo elétrico fase b
%E_yr3 componente y real campo elétrico fase c

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));

%E_yi1 componente y imaginario campo elétrico fase a
%E_yi2 componente y imaginario campo elétrico fase b
%E_yi3 componente y imaginario campo elétrico fase c

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));

E_xr = [ E_xr1;  E_xr2;  E_xr3 ];
E_ximag = [ E_xi1;  E_xi2;  E_xi3 ];
E_yr = [ E_yr1;  E_yr2;  E_yr3 ];
E_yimag = [ E_yi1;  E_yi2;  E_yi3 ];

elseif size(x,2) == 12 %caso 2

rho = fcn_carga(V,x,y,r); %entram valores para cálculo da densidade de carga
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

%E_xr1 componente x real campo elétrico fase a cabo 1
%E_xr2 componente x real campo elétrico fase a cabo 2
%E_xr3 componente x real campo elétrico fase b cabo 3
%E_xr4 componente x real campo elétrico fase b cabo 4
%E_xr5 componente x real campo elétrico fase c cabo 5
%E_xr6 componente x real campo elétrico fase c cabo 6

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_xr4 = (rho_r4/(2*pi*e_0)).*((xf - x(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (xf - x(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_xr5 = (rho_r5/(2*pi*e_0)).*((xf - x(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (xf - x(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_xr6 = (rho_r6/(2*pi*e_0)).*((xf - x(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (xf - x(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));

%E_xi1 componente x imaginario campo elétrico fase a cabo 1
%E_xi2 componente x imaginario campo elétrico fase a cabo 2
%E_xi3 componente x imaginario campo elétrico fase b cabo 3
%E_xi4 componente x imaginario campo elétrico fase b cabo 4
%E_xi5 componente x imaginario campo elétrico fase c cabo 5
%E_xi6 componente x imaginario campo elétrico fase c cabo 6

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_xi4 = (rho_i4/(2*pi*e_0)).*((xf - x(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (xf - x(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_xi5 = (rho_i5/(2*pi*e_0)).*((xf - x(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (xf - x(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_xi6 = (rho_i6/(2*pi*e_0)).*((xf - x(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (xf - x(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));

%E_yr1 componente y real campo elétrico fase a cabo 1
%E_yr2 componente y real campo elétrico fase a cabo 2
%E_yr3 componente y real campo elétrico fase b cabo 3
%E_yr4 componente y real campo elétrico fase b cabo 4
%E_yr5 componente y real campo elétrico fase c cabo 5
%E_yr6 componente y real campo elétrico fase c cabo 6

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_yr4 = (rho_r4/(2*pi*e_0)).*((yf - y(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (yf + y(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_yr5 = (rho_r5/(2*pi*e_0)).*((yf - y(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (yf + y(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_yr6 = (rho_r6/(2*pi*e_0)).*((yf - y(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (yf + y(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));

%E_yi1 componente y imaginario campo elétrico fase a cabo 1
%E_yi2 componente y imaginario campo elétrico fase a cabo 2
%E_yi3 componente y imaginario campo elétrico fase b cabo 3
%E_yi4 componente y imaginario campo elétrico fase b cabo 4
%E_yi5 componente y imaginario campo elétrico fase c cabo 5
%E_yi6 componente y imaginario campo elétrico fase c cabo 6

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_yi4 = (rho_i4/(2*pi*e_0)).*((yf - y(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (yf + y(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_yi5 = (rho_i5/(2*pi*e_0)).*((yf - y(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (yf + y(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_yi6 = (rho_i6/(2*pi*e_0)).*((yf - y(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (yf + y(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));

E_xr = [ E_xr1;  E_xr2 ;  E_xr3 ; E_xr4 ; E_xr5 ; E_xr6 ];
E_ximag = [ E_xi1;  E_xi2 ;  E_xi3 ; E_xi4 ; E_xi5 ; E_xi6 ];
E_yr = [ E_yr1;  E_yr2 ;  E_yr3 ; E_yr4 ; E_yr5 ; E_yr6 ];
E_yimag = [ E_yi1;  E_yi2 ;  E_yi3 ; E_yi4 ; E_yi5 ; E_yi6 ];

elseif size(x,2) == 18

rho = fcn_carga(V,x,y,r); %entram valores para cálculo da densidade de carga
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

%E_xr1 componente x real campo elétrico fase a cabo 1
%E_xr2 componente x real campo elétrico fase a cabo 2
%E_xr3 componente x real campo elétrico fase a cabo 3
%E_xr4 componente x real campo elétrico fase b cabo 4
%E_xr5 componente x real campo elétrico fase b cabo 5
%E_xr6 componente x real campo elétrico fase b cabo 6
%E_xr7 componente x real campo elétrico fase c cabo 7
%E_xr8 componente x real campo elétrico fase c cabo 8
%E_xr9 componente x real campo elétrico fase c cabo 9

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_xr4 = (rho_r4/(2*pi*e_0)).*((xf - x(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (xf - x(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_xr5 = (rho_r5/(2*pi*e_0)).*((xf - x(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (xf - x(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_xr6 = (rho_r6/(2*pi*e_0)).*((xf - x(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (xf - x(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_xr7 = (rho_r7/(2*pi*e_0)).*((xf - x(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (xf - x(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_xr8 = (rho_r8/(2*pi*e_0)).*((xf - x(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (xf - x(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_xr9 = (rho_r9/(2*pi*e_0)).*((xf - x(1,9))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (xf - x(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));

%E_xi1 componente x imaginario campo elétrico fase a cabo 1
%E_xi2 componente x imaginario campo elétrico fase a cabo 2
%E_xi3 componente x imaginario campo elétrico fase a cabo 3
%E_xi4 componente x imaginario campo elétrico fase b cabo 4
%E_xi5 componente x imaginario campo elétrico fase b cabo 5
%E_xi6 componente x imaginario campo elétrico fase b cabo 6
%E_xi4 componente x imaginario campo elétrico fase c cabo 7
%E_xi5 componente x imaginario campo elétrico fase c cabo 8
%E_xi6 componente x imaginario campo elétrico fase c cabo 9

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_xi4 = (rho_i4/(2*pi*e_0)).*((xf - x(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (xf - x(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_xi5 = (rho_i5/(2*pi*e_0)).*((xf - x(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (xf - x(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_xi6 = (rho_i6/(2*pi*e_0)).*((xf - x(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (xf - x(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_xi7 = (rho_i7/(2*pi*e_0)).*((xf - x(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (xf - x(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_xi8 = (rho_i8/(2*pi*e_0)).*((xf - x(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (xf - x(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_xi9 = (rho_i9/(2*pi*e_0)).*((xf - x(1,9))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (xf - x(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));

%E_yr1 componente y real campo elétrico fase a cabo 1
%E_yr2 componente y real campo elétrico fase a cabo 2
%E_yr3 componente y real campo elétrico fase b cabo 3
%E_yr4 componente y real campo elétrico fase b cabo 4
%E_yr5 componente y real campo elétrico fase b cabo 5
%E_yr6 componente y real campo elétrico fase b cabo 6
%E_yr7 componente y real campo elétrico fase c cabo 7
%E_yr8 componente y real campo elétrico fase c cabo 8
%E_yr9 componente y real campo elétrico fase c cabo 9

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_yr4 = (rho_r4/(2*pi*e_0)).*((yf - y(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (yf + y(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_yr5 = (rho_r5/(2*pi*e_0)).*((yf - y(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (yf + y(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_yr6 = (rho_r6/(2*pi*e_0)).*((yf - y(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (yf + y(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_yr7 = (rho_r7/(2*pi*e_0)).*((yf - y(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (yf + y(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_yr8 = (rho_r8/(2*pi*e_0)).*((yf - y(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (yf + y(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_yr9 = (rho_r9/(2*pi*e_0)).*((yf - y(1,9))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (yf + y(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));


%E_yi1 componente y imaginario campo elétrico fase a cabo 1
%E_yi2 componente y imaginario campo elétrico fase a cabo 2
%E_yi3 componente y imaginario campo elétrico fase a cabo 3
%E_yi4 componente y imaginario campo elétrico fase b cabo 4
%E_yi5 componente y imaginario campo elétrico fase b cabo 5
%E_yi6 componente y imaginario campo elétrico fase b cabo 6
%E_yi7 componente y imaginario campo elétrico fase c cabo 7
%E_yi8 componente y imaginario campo elétrico fase c cabo 8
%E_yi9 componente y imaginario campo elétrico fase c cabo 9

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_yi4 = (rho_i4/(2*pi*e_0)).*((yf - y(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (yf + y(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_yi5 = (rho_i5/(2*pi*e_0)).*((yf - y(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (yf + y(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_yi6 = (rho_i6/(2*pi*e_0)).*((yf - y(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (yf + y(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_yi7 = (rho_i7/(2*pi*e_0)).*((yf - y(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (yf + y(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_yi8 = (rho_i8/(2*pi*e_0)).*((yf - y(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (yf + y(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_yi9 = (rho_i9/(2*pi*e_0)).*((yf - y(1,9))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (yf + y(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));

E_xr = [ E_xr1; E_xr2 ; E_xr3 ; E_xr4 ; E_xr5 ; E_xr6 ; E_xr7 ; E_xr8 ; E_xr9 ];
E_ximag = [ E_xi1; E_xi2 ; E_xi3 ; E_xi4 ; E_xi5 ; E_xi6 ; E_xi7 ; E_xi8 ; E_xi9  ];
E_yr = [ E_yr1; E_yr2 ; E_yr3 ; E_yr4 ; E_yr5 ; E_yr6 ; E_yr7 ; E_yr8 ; E_yr9  ];
E_yimag = [ E_yi1; E_yi2 ; E_yi3 ; E_yi4 ; E_yi5 ; E_yi6 ; E_yi7 ; E_yi8 ; E_yi9 ];

elseif size(x,2) == 24 %caso 4

rho = fcn_carga(V,x,y,r); %entram valores para cálculo da densidade de carga
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

%E_xr1 componente x real campo elétrico fase a cabo 1
%E_xr2 componente x real campo elétrico fase a cabo 2
%E_xr3 componente x real campo elétrico fase a cabo 3
%E_xr4 componente x real campo elétrico fase a cabo 4
%E_xr5 componente x real campo elétrico fase b cabo 5
%E_xr6 componente x real campo elétrico fase b cabo 6
%E_xr7 componente x real campo elétrico fase b cabo 7
%E_xr8 componente x real campo elétrico fase b cabo 8
%E_xr9 componente x real campo elétrico fase c cabo 9
%E_xr10 componente x real campo elétrico fase c cabo 10
%E_xr11 componente x real campo elétrico fase c cabo 11
%E_xr12 componente x real campo elétrico fase c cabo 12

E_xr1 = (rho_r1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xr2 = (rho_r2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xr3 = (rho_r3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_xr4 = (rho_r4/(2*pi*e_0)).*((xf - x(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (xf - x(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_xr5 = (rho_r5/(2*pi*e_0)).*((xf - x(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (xf - x(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_xr6 = (rho_r6/(2*pi*e_0)).*((xf - x(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (xf - x(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_xr7 = (rho_r7/(2*pi*e_0)).*((xf - x(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (xf - x(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_xr8 = (rho_r8/(2*pi*e_0)).*((xf - x(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (xf - x(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_xr9 = (rho_r9/(2*pi*e_0)).*((xf - x(1,9))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (xf - x(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));
E_xr10 = (rho_r10/(2*pi*e_0)).*((xf - x(1,10))./(((xf - x(1,10)).^2) + (yf - y(1,10)).^2) - (xf - x(1,10))./(((xf - x(1,10)).^2) + (yf + y(1,10)).^2));
E_xr11 = (rho_r11/(2*pi*e_0)).*((xf - x(1,11))./(((xf - x(1,11)).^2) + (yf - y(1,11)).^2) - (xf - x(1,11))./(((xf - x(1,11)).^2) + (yf + y(1,11)).^2));
E_xr12 = (rho_r12/(2*pi*e_0)).*((xf - x(1,12))./(((xf - x(1,12)).^2) + (yf - y(1,12)).^2) - (xf - x(1,12))./(((xf - x(1,12)).^2) + (yf + y(1,12)).^2));

%E_xi1 componente x imaginario campo elétrico fase a cabo 1
%E_xi2 componente x imaginario campo elétrico fase a cabo 2
%E_xi3 componente x imaginario campo elétrico fase a cabo 3
%E_xi4 componente x imaginario campo elétrico fase a cabo 4
%E_xi5 componente x imaginario campo elétrico fase b cabo 5
%E_xi6 componente x imaginario campo elétrico fase b cabo 6
%E_xi7 componente x imaginario campo elétrico fase b cabo 7
%E_xi8 componente x imaginario campo elétrico fase b cabo 8
%E_xi9 componente x imaginario campo elétrico fase c cabo 9
%E_xi10 componente x imaginario campo elétrico fase c cabo 10
%E_xi11 componente x imaginario campo elétrico fase c cabo 11
%E_xi12 componente x imaginario campo elétrico fase c cabo 12

E_xi1 = (rho_i1/(2*pi*e_0)).*((xf - x(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (xf - x(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_xi2 = (rho_i2/(2*pi*e_0)).*((xf - x(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (xf - x(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_xi3 = (rho_i3/(2*pi*e_0)).*((xf - x(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (xf - x(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_xi4 = (rho_i4/(2*pi*e_0)).*((xf - x(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (xf - x(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_xi5 = (rho_i5/(2*pi*e_0)).*((xf - x(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (xf - x(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_xi6 = (rho_i6/(2*pi*e_0)).*((xf - x(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (xf - x(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_xi7 = (rho_i7/(2*pi*e_0)).*((xf - x(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (xf - x(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_xi8 = (rho_i8/(2*pi*e_0)).*((xf - x(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (xf - x(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_xi9 = (rho_i9/(2*pi*e_0)).*((xf - x(1,7))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (xf - x(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));
E_xi10 = (rho_i10/(2*pi*e_0)).*((xf - x(1,10))./(((xf - x(1,10)).^2) + (yf - y(1,10)).^2) - (xf - x(1,10))./(((xf - x(1,10)).^2) + (yf + y(1,10)).^2));
E_xi11 = (rho_i11/(2*pi*e_0)).*((xf - x(1,11))./(((xf - x(1,11)).^2) + (yf - y(1,11)).^2) - (xf - x(1,11))./(((xf - x(1,11)).^2) + (yf + y(1,11)).^2));
E_xi12 = (rho_i12/(2*pi*e_0)).*((xf - x(1,12))./(((xf - x(1,12)).^2) + (yf - y(1,12)).^2) - (xf - x(1,12))./(((xf - x(1,12)).^2) + (yf + y(1,12)).^2));


%E_yr1 componente y real campo elétrico fase a cabo 1
%E_yr2 componente y real campo elétrico fase a cabo 2
%E_yr3 componente y real campo elétrico fase a cabo 3
%E_yr4 componente y real campo elétrico fase a cabo 4
%E_yr5 componente y real campo elétrico fase b cabo 5
%E_yr6 componente y real campo elétrico fase b cabo 6
%E_yr7 componente y real campo elétrico fase b cabo 7
%E_yr8 componente y real campo elétrico fase b cabo 8
%E_yr9 componente y real campo elétrico fase c cabo 9
%E_yr10 componente y real campo elétrico fase c cabo 10
%E_yr11 componente y real campo elétrico fase c cabo 11
%E_yr12 componente y real campo elétrico fase c cabo 12

E_yr1 = (rho_r1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yr2 = (rho_r2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yr3 = (rho_r3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_yr4 = (rho_r4/(2*pi*e_0)).*((yf - y(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (yf + y(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_yr5 = (rho_r5/(2*pi*e_0)).*((yf - y(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (yf + y(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_yr6 = (rho_r6/(2*pi*e_0)).*((yf - y(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (yf + y(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_yr7 = (rho_r7/(2*pi*e_0)).*((yf - y(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (yf + y(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_yr8 = (rho_r8/(2*pi*e_0)).*((yf - y(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (yf + y(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_yr9 = (rho_r9/(2*pi*e_0)).*((yf - y(1,9))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (yf + y(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));
E_yr10 = (rho_r10/(2*pi*e_0)).*((yf - y(1,10))./(((xf - x(1,10)).^2) + (yf - y(1,10)).^2) - (yf + y(1,10))./(((xf - x(1,10)).^2) + (yf + y(1,10)).^2));
E_yr11 = (rho_r11/(2*pi*e_0)).*((yf - y(1,11))./(((xf - x(1,11)).^2) + (yf - y(1,11)).^2) - (yf + y(1,11))./(((xf - x(1,11)).^2) + (yf + y(1,11)).^2));
E_yr12 = (rho_r12/(2*pi*e_0)).*((yf - y(1,12))./(((xf - x(1,12)).^2) + (yf - y(1,12)).^2) - (yf + y(1,12))./(((xf - x(1,12)).^2) + (yf + y(1,12)).^2));

%E_yi1 componente y imaginario campo elétrico fase a cabo 1
%E_yi2 componente y imaginario campo elétrico fase a cabo 2
%E_yi3 componente y imaginario campo elétrico fase a cabo 3
%E_yi4 componente y imaginario campo elétrico fase a cabo 4
%E_yi5 componente y imaginario campo elétrico fase b cabo 5
%E_yi6 componente y imaginario campo elétrico fase b cabo 6
%E_yi7 componente y imaginario campo elétrico fase b cabo 7
%E_yi8 componente y imaginario campo elétrico fase b cabo 8
%E_yi9 componente y imaginario campo elétrico fase c cabo 9
%E_yi10 componente y imaginario campo elétrico fase c cabo 10
%E_yi11 componente y imaginario campo elétrico fase c cabo 11
%E_yi12 componente y imaginario campo elétrico fase c cabo 12

E_yi1 = (rho_i1/(2*pi*e_0)).*((yf - y(1,1))./(((xf - x(1,1)).^2) + (yf - y(1,1)).^2) - (yf + y(1,1))./(((xf - x(1,1)).^2) + (yf + y(1,1)).^2));
E_yi2 = (rho_i2/(2*pi*e_0)).*((yf - y(1,2))./(((xf - x(1,2)).^2) + (yf - y(1,2)).^2) - (yf + y(1,2))./(((xf - x(1,2)).^2) + (yf + y(1,2)).^2));
E_yi3 = (rho_i3/(2*pi*e_0)).*((yf - y(1,3))./(((xf - x(1,3)).^2) + (yf - y(1,3)).^2) - (yf + y(1,3))./(((xf - x(1,3)).^2) + (yf + y(1,3)).^2));
E_yi4 = (rho_i4/(2*pi*e_0)).*((yf - y(1,4))./(((xf - x(1,4)).^2) + (yf - y(1,4)).^2) - (yf + y(1,4))./(((xf - x(1,4)).^2) + (yf + y(1,4)).^2));
E_yi5 = (rho_i5/(2*pi*e_0)).*((yf - y(1,5))./(((xf - x(1,5)).^2) + (yf - y(1,5)).^2) - (yf + y(1,5))./(((xf - x(1,5)).^2) + (yf + y(1,5)).^2));
E_yi6 = (rho_i6/(2*pi*e_0)).*((yf - y(1,6))./(((xf - x(1,6)).^2) + (yf - y(1,6)).^2) - (yf + y(1,6))./(((xf - x(1,6)).^2) + (yf + y(1,6)).^2));
E_yi7 = (rho_i7/(2*pi*e_0)).*((yf - y(1,7))./(((xf - x(1,7)).^2) + (yf - y(1,7)).^2) - (yf + y(1,7))./(((xf - x(1,7)).^2) + (yf + y(1,7)).^2));
E_yi8 = (rho_i8/(2*pi*e_0)).*((yf - y(1,8))./(((xf - x(1,8)).^2) + (yf - y(1,8)).^2) - (yf + y(1,8))./(((xf - x(1,8)).^2) + (yf + y(1,8)).^2));
E_yi9 = (rho_i9/(2*pi*e_0)).*((yf - y(1,9))./(((xf - x(1,9)).^2) + (yf - y(1,9)).^2) - (yf + y(1,9))./(((xf - x(1,9)).^2) + (yf + y(1,9)).^2));
E_yi10 = (rho_i10/(2*pi*e_0)).*((yf - y(1,10))./(((xf - x(1,10)).^2) + (yf - y(1,10)).^2) - (yf + y(1,10))./(((xf - x(1,10)).^2) + (yf + y(1,10)).^2));
E_yi11 = (rho_i11/(2*pi*e_0)).*((yf - y(1,11))./(((xf - x(1,11)).^2) + (yf - y(1,11)).^2) - (yf + y(1,11))./(((xf - x(1,11)).^2) + (yf + y(1,11)).^2));
E_yi12 = (rho_i12/(2*pi*e_0)).*((yf - y(1,12))./(((xf - x(1,12)).^2) + (yf - y(1,12)).^2) - (yf + y(1,12))./(((xf - x(1,12)).^2) + (yf + y(1,12)).^2));

E_xr = [ E_xr1 ;  E_xr2 ;  E_xr3 ; E_xr4 ; E_xr5 ; E_xr6 ; E_xr7 ; E_xr8 ; E_xr9 ; E_xr10 ; E_xr11 ; E_xr12 ];
E_ximag = [ E_xi1 ;  E_xi2 ;  E_xi3 ; E_xi4 ; E_xi5 ; E_xi6 ; E_xi7 ; E_xi8 ; E_xi9 ; E_xi10 ; E_xi11 ; E_xi12 ];
E_yr = [ E_yr1 ;  E_yr2 ;  E_yr3 ; E_yr4 ; E_yr5 ; E_yr6 ; E_yr7 ; E_yr8 ; E_yr9 ; E_yr10 ; E_yr11 ; E_yr12 ];
E_yimag = [ E_yi1 ;  E_yi2 ;  E_yi3 ; E_yi4 ; E_yi5 ; E_yi6 ; E_yi7 ; E_yi8 ; E_yi9 ; E_yi10 ; E_yi11 ; E_yi12 ];

end

Exr = sum(E_xr).^2;
Eximag = sum(E_ximag).^2;
Eyr = sum(E_yr).^2;
Eyimag = sum(E_yimag).^2;

Erms = (Exr+Eximag+Eyr+Eyimag).^(1/2);