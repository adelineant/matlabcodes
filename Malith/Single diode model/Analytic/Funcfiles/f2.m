function y=f2(I,V,Isc,Voc, Rs, n, Rsh, Iph)

  q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    T = 300;
    Vt = kb*T/q;
y = (Isc-Voc/Rsh)*exp(-Voc/(n*Vt)).*(exp(q.*(V-Rs.*-I)./(n.*kb.*T))-1)+(V-Rs.*-I)./Rsh-Iph+I;
end