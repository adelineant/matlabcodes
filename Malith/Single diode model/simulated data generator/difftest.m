clc,clear all,close all
syms Rs N Rsh integer
syms f(x)

%f(x) = f(x)*log(x);
%gradient(f, x)
%lambertw(0, 2)/(x*(lambertw(0, 2) + 1))
%double(subs(f,2))

syms y(I) dvdI(I) dvI2(I);

dvdI(I) = gradient(y,I)
dvdI2(I) = gradient(dvdI,I);

rad = ((((dvdI2))./((1 + dvdI.^2).^(3/2))));

grrad = gradient(rad,I)


syms dy(z) y(z) 

y(z) = z^2;
dy(z) = diff(y(z),z);
 gradient(dy(z),z)
 
 syms z(x) 
 z(x) = x^2;
 y(x) = 5*z(x);
 
 v = z(x) + y(z);
 gradient(v,x);
 
 
 
 syms V(Ireg) Ireg kb q Isc Rsh Rs N

 dvdI = ((N*kb/q)/((Isc + Ireg - V(Ireg)/Rsh + N*kb/(q*Rsh) + ((Ireg + Isc)*Rs)/Rsh))) +Rs;
 
 gradient(dvdI,Ireg)
 
 %-(N*kb*(Rs/Rsh - dVdi/Rsh + 1))/(q*(Ireg + Isc - V(Ireg)/Rsh + (Rs*(Ireg + Isc))/Rsh + (N*kb)/(Rsh*q))^2)
 
 %third der
 syms V(Ireg) Ireg kb q Isc Rsh Rs N
 dvdi2 = -(N*kb*(Rs/Rsh - gradient(V,Ireg)/Rsh + 1))/(q*(Ireg + Isc - V(Ireg)/Rsh + (Rs*(Ireg + Isc))/Rsh + (N*kb)/(Rsh*q))^2)
 dvdi3 = gradient(dvdi2,Ireg)
 
