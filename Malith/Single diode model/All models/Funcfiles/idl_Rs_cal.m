function [n,Rs] = idl_Rs_cal(GradRs0,GradRsh0,Vt,Isc,Voc,Im,Vm)
    %n and Rs calculated by using the fact that dv/di = -V/I at
    %maxpowerpoint.
    %A method i came up with to improve it.
    a = (GradRs0 + Vm/Im)/Vt;
    b = (1/(Isc - Voc/GradRsh0) - 1/(Isc + Im - Vm/GradRsh0));
    n = a/b;
    Rs = (-Vm/Im) - (n*Vt)/(Isc + Im - Vm/GradRsh0);


end