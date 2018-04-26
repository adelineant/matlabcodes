clear all; close all; clc

%Opening data file
fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1)';
I = -A(:,2)';

%Finding Isc & Voc
Isc_index = find(abs(V)==min(abs(V-0)));
Isc= I(Isc_index);
Voc_index = find(abs(I)==min(abs(I-0)));
Voc= V(Voc_index);

%Finding Im & Vm
mpp = abs(I([Isc_index:Voc_index]).*V([Isc_index:Voc_index]));
max_index = find(mpp==max(mpp));
Im = I(max_index);
Vm = V(max_index);

%Known Parameters
k=1.38E-23;
T= 306.15;
q = 1.6012E-19;
Vt = k*T/q;

%Extracting Rsh0 & Rs0 (Polynomial fitting)
Rsh0 = 3953.12;
Rs0 = 36.23;

%Extracting parameters
Rsh=Rsh0;
n=(Vm+Im.*Rs0-Voc)./(Vt.*(log(Isc-Vm./Rsh-Im)-log(Isc-Voc./Rsh)+Im./(Isc-Voc./Rsh)));
Is = (Isc-Voc/Rsh)*exp(-Voc/(n*Vt));
Rs = Rs0-n*Vt/Is*exp(-Voc/(n*Vt));
Iph = Isc*(1+Rs/Rsh)+Is*(exp(Isc*Rs/(n*Vt))-1);

%Plotting extracted parameters
 Ical=[];
for i=1:length(V(1:1:Voc_index))
    Ical=[Ical fzero(@(I)f2(I,V(i),Is, Rs, n, Rsh, Iph),0)];
end
plot (V(1:i),I(1:i), 'r');
hold on
plot (V(1:i),Ical(1:i), 'b');
legend ('experimental', 'fitted')
xlabel ('V')
ylabel ('I')

%Calculating least squares mean error
error = mean((I(1:i)-Ical).^2);  

%Printing parameters
fprintf('SINGLE DIODE ANALYTICAL PARAMETERS\n')
fprintf('Isc(A): %.2d \n',Isc)
fprintf('Voc(V): %.2d \n',Voc)
fprintf('Im,cal(A): %.2d \n',Im)
fprintf('Vm,cal(V): %.2d \n',Vm)
fprintf('Iph: %.2d \n',Iph)
fprintf('Is: %.2d \n',Is)
fprintf('Rs(ohms): %.2d \n',Rs)
fprintf('Rsh(ohms): %.2d \n',Rsh)
fprintf('ERROR: %.2d \n',error)




