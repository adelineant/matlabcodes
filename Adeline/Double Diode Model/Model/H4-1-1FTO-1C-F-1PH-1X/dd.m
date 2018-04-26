function y=dd(I,V, n1, n2, Rs, Rsh, I02, I01, Iph, Ns, T)
q=1.602E-19;
k=1.38065E-23;
a1=n1*Ns*k*T/q;
a2=n2*Ns*k*T/q;
Id2 = I02.*(exp((V+I.*Rs)/a2)-1);
Id1 = I01.*(exp((V+I.*Rs)/a1)-1);
y=Iph-Id1-Id2-(V+I.*Rs)/Rsh-I;
end