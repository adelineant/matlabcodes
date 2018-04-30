clear all; close all; clc

%Reading the data file and generating I and V vector
fileID = fopen('s75.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1)';
I = A(:,2)';

%Known parameters
k=1.38065E-23;
T=298.15;
q=1.602E-19;
Vt1=k*T/q;
Vt2=k*T/q;

%Input parameters
a1=1;%1.14
a2=2;%2.6
Ns=36;

%Initialising
Rs=0;
e=100;

while e>0.0001
   
    %Finding Voc (Need to include polynomial)
    Voc_index = find(abs(I)==min(abs(I-0)));
    Voc = V(Voc_index);

    %Finding Isc (Need to include polynomial)
    Isc_index = find(abs(V)==min(abs(V-0)));
    Isc = I(Isc_index);

    %Finding Im and Vm
    mpp = abs(I([Isc_index:Voc_index]).*V([Isc_index:Voc_index]));
    max_index = find(mpp==max(mpp));
    Im = I(max_index);
    Vm = V(max_index);

   % Assumptions
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
    
    %Determining when to stop loop
    Rscal = Vm/Im-1/(Is1/(a1*Ns*Vt1)*Xm1+Is2/(a2*Ns*Vt1)*Xm2+1/Rsh);
    e = Rs- Rscal;
    Rs = Rs+0.001;

end

% Plotting extracted parameters
    Ical=[];
    for i=1:length(V(1:1:Voc_index))
        Ical=[Ical fzero(@(I)dd(I,V(i),a1, a2, Rs, Rsh, Is2, Is1, Iph, Ns, T),0)];
    end
    plot (V(1:i),Ical(1:i), 'r')
    hold on
    plot (V(1:i),I(1:i), 'b')
    legend ('fitted', 'simulated')
    ylim([0 5])

%Calculating least squares mean error
error = mean((I(1:i)-Ical).^2);  

%Printing parameters
fprintf('DOUBLE DIODE PARAMETERS\n')
fprintf('a1: %.2d \n',a1)
fprintf('a2: %.2d \n',a2)
fprintf('Isc(A): %.2d \n',Isc)
fprintf('Voc(V): %.2d \n',Voc)
fprintf('Im,cal(A): %.2d \n',Im)
fprintf('Vm,cal(V): %.2d \n',Vm)
fprintf('Iph: %.2d \n',Iph)
fprintf('Is1: %.2d \n',Is1)
fprintf('Is2: %.2d \n',Is2)
fprintf('Rs(ohms): %.2d \n',Rs)
fprintf('Rsh(ohms): %.2d \n',Rsh)
fprintf('ERROR: %.2d \n',error)


    

