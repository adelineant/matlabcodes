clc;clear all
b = fzero(@minthisshit,1000);
Rshcal = minthisshit(b);


function min = minthisshit(N)
    
    kb = 1.38e-23;
    q = 1.6012e-19;
    Voc = 0.969384;
    Isc = 0.003701;
    Rs0 = 32.66;
    Rsh0 = 3.005510172434505e+03;

    a = N*kb/q;
    c = (1 - exp(-q*Voc/(N*kb)));

    Rshcal = (((a/Rs0) - (Isc/c))/(-Voc/c + a))^(-1);

    %dvdi = (a/(Isc*(1/c -1) - 1/Rshcal*(Voc/c - 1/a)))

    Rshcal2 = (((a/Rsh0) - Isc*(1/c -1))/(a - Voc/c))^(-1);
    
    min = Rshcal2 - Rshcal;
end



