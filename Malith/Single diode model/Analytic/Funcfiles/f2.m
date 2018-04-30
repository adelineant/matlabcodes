function y=f2(I,V,I0, Rs, n, Rsh, Iph)
q = 1.6012E-19;
kb = 1.38E-23;
T = 306.15;
y = I0.*(exp(q.*(V-Rs.*-I)./(n.*kb.*T))-1)+(V-Rs.*-I)./Rsh-Iph+I;
end