clc,clear all, close all

V = 0:0.01:0.7;
Rs = 0.0364;
n = 1.484;
Rsh = 0.0538*1000;
Iph =0.7608;
I0 =0.3223*10^-6;
T = 306.15;

%y=SD_equation(V,I,Rs,Rsh,n)
Ical = zeros(size(V));

Ical(1)= fzero(@(I)SD_equation(V(1),I,Rs,Rsh,n,Iph,I0,T),-0.3);

for i=1:length(V)-1
    Ical(i+1) = fzero(@(I)SD_equation(V(i+1),I,Rs,Rsh,n,Iph,I0,T),1);
end
plot (V,-Ical)



function y=SD_equation(V,I,Rs,Rsh,n,Iph,I0,T)

q = 1.6012*10^(-19);
kb = 1.38*10^(-23);


y = I0.*(exp(q.*(V-Rs.*I)./(n.*kb.*T))-1)+(V-Rs.*I)./Rsh-Iph-I;
end