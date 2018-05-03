%major progress

Vm = V(bb);
Ireg = -Ireg(bb);
Isc = 0.0037;
N = 300;
Rsh = 100000;
Rs = 30;
kb = 1.38*10^(-23);
q = 1.6012*10^(-19);


first = 2*((Isc + Ireg - Vm/Rsh + N*kb/(q*Rsh) + ((Ireg + Isc)*Rs)/Rsh))/...
    (N*kb/q)+ 1/(Rsh +(Rs +Vm/Im)) 

second = 3*(((1 + (-Vm/Ireg)^2))^0.5)*(-Vm/Ireg)/((1-(Vm/Ireg))^1.5)