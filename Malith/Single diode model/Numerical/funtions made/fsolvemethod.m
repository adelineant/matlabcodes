clc;clear all;close all
fun = @parameters;
x0 = [30,500];
options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
x = fsolve(fun,x0,options)

function F = parameters(R)
Isc = 0.003701000000000;
q = 1.601200000000000e-19;
kb =  1.380000000000000e-23;
Im = 0.003083727759953;
Vm = 0.723909000000000;
Voc = 0.969384000000000;

Rs0 =  32.660633313050127;
Rsh0 = 3000;

alpha = ((Isc + (R(1)*Isc - Voc))/Rsh0)/...
    (1 - exp(q*((R(1)*Isc - Voc)/(R(2)*kb)))) + Voc/Rsh0;

F(1) = Rs0 - (R(1) + (R(2)*kb/q)/((alpha-(Voc-R(2)*kb/q)/Rsh0)));
F(2) = (Vm/Im) - (R(1) + (R(2)*kb/q)/((alpha-Im-(Vm - R(1)*(-Im) - R(2)*kb/q)/Rsh0)));
F(3) = Rsh0 - (R(1) + (R(2)*kb/q)/((alpha-Isc - (-R(1)*(-Isc) - R(2)*kb/q)/Rsh0)));
end