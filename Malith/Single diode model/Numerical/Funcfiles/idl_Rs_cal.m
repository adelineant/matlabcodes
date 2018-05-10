function [n,Rs] = idl_Rs_cal(GradRs0,GradRsh0,Vt,Isc,Voc,Im,Vm)

    a = (GradRs0 + Vm/Im)/Vt;
    b = (1/(Isc - Voc/GradRsh0) - 1/(Isc + Im - Vm/GradRsh0));
    n = a/b;
    Rs = (-Vm/Im) - (n*Vt)/(Isc + Im - Vm/GradRsh0);


end