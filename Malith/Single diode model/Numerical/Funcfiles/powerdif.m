clc,clear all
N = 6.774560546875000e+02;
Rs = 32.660633313050130;
Rsh = 3.005510172434505e+03;
q = 1.6012e-19;
kb = 1.38e-23;
Isc = 0.0037;
Voc = 0.9694;
Vm = 0.723909000000000;
Im = 0.003083727759953;

x = fzero(@minN,N,[],Rs,Rsh,q,kb,Isc,Voc,Vm,Im)

function y = minN(N,Rs,Rsh,q,kb,Isc,Voc,Vm,Im)

Io = (Isc + ((Rs*Isc - Voc)/Rsh))*exp(-q*Voc/(N*kb))/...
    (1-exp(q*(Rs*Isc - Voc)/(N*kb)));
Iph = Io*(exp((q*Voc)/(N*kb)))+ Voc/Rsh;

y = -Iph + (Vm - Rs*Im)/Rsh + Vm/Rsh*(1- Rs*(Im/Vm)) + ...
    Io*(exp(q*(Vm-Rs*Im)/(N*kb)) -1 )...
    + (Vm*Io/N*kb)*(q - Rs*(Im/Vm))*exp(q*(Vm-Rs*Im)/(N*kb));
end

