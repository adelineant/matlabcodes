clc;clear all    


N = 4.73*300;

   
    
    kb = 1.38e-23;
    q = 1.6012e-19;
    Voc = 1;
    Isc = 0.00200;
    Rs0 = 6.121173711272144;
    Rsh0 = 1.542555292209305e+03;

    a = N*kb/q;
    c = (1 - exp(-q*Voc/(N*kb)));

    Rshcal = (((a/Rs0) - (Isc/c))/(-Voc/c + a))^(-1);

    %dvdi = (a/(Isc*(1/c -1) - 1/Rshcal*(Voc/c - 1/a)))

    Rshcal2 = (((a/Rsh0) - Isc*(1/c -1))/(Voc+ a - Voc/c))^(-1);
