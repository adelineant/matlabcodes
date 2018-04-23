clear all; close all; clc

%Reading the data file and generating I and V vector
fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1)';
I = -A(:,2)';

%Known parameters
Rs=36;
k=1.38065E-23;
T=298;
q=1.602E-19;
Vt1=k*T/q;
Vt2=k*T/q;
a1=1;
a2=2
for Rs = 1:36
    Ns=1; %Is this one for only one solar cell?

    %Use point closest to Voc as Voc (Need to improve)
    Voc_index = find(abs(I)==min(abs(I-0)));
    Voc = V(Voc_index);

    %Use point closest to Isc as Isc (Need to improve)
    Isc_index = find(abs(V)==min(abs(V-0)));
    Isc = I(Isc_index);

    %Find Im and Vm from MPP tracking
    mpp = abs(I([Isc_index:Voc_index]).*V([Isc_index:Voc_index]));
    max_index = find(mpp==max(mpp));
    Im = I(max_index);
    Vm = V(max_index);

    %Plotting experimental data
    % plot (V,I)
    % hold on
    % plot (Voc, I(Voc_index), 'r*')
    % plot (V(Isc_index), Isc, 'r*')
    % plot (Vm, Im, 'b*')

    %Assumptions
    Xoc1 = exp(Voc/(a1*Ns*Vt1));
    Xoc2 = exp(Voc/(a2*Ns*Vt2));
    Xm1 = exp((Vm+Rs*Im)/(a1*Ns*Vt1));
    Xm2 = exp((Vm+Rs*Im)/(a2*Ns*Vt2));
    Xs1 = exp(Rs*Isc/(a1*Ns*Vt1));
    Xs2 = exp(Rs*Isc/(a2*Ns*Vt2));
    %Extracting parameters
    K=T^(2/5)/3.77;
    Is1 = (Voc*(Isc-Im)-Vm*Isc)/(Voc*(Xm1+K*Xm2)-Vm*(Xoc1+K*Xoc2));
    Is2 = K*Is1;
    Iph = (Voc*Im+Is1*(Voc*(Xm1+K*Xm2)-Vm*(Xoc1-K*Xoc2)))/(Voc-Vm);
    Rsh = (Vm + Im*Rs)/(Iph-Im-Is1*(Xm1-1)-Is2*(Xm2-1));

    I2=[];
    for i=1:length(V)
        I2=[I2 fzero(@(I)dd(I,V(i),a1, a2, Rs, Rsh, Is2, Is1, Iph, Ns, T),0)];
    end
    plot (V(1:1:Voc_index),I2(1:1:Voc_index), 'r')
    hold on
    plot (V(1:1:Voc_index),I(1:1:Voc_index), 'b')
end

