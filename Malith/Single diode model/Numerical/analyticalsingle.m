clear all; close all; clc

fileID = fopen('simulateddata.txt','r');
%scan the data to A matrix. Should probably initialize this. There are only
%two columns but infite rows 
A = [fscanf(fileID,'%f',[2 Inf])]';
%close any opened file to prevent memory leaks
fclose(fileID);
%the first column SHOULD have voltage and the second column should have
%current
V = A(:,1);
I = A(:,2);

Isc_index = find(V==min(V(V>=0)));
Isc= (I(Isc_index));
Voc_index = find(I==min(I(I>=0)));
Voc= V(Voc_index);
Rs0 = -(V(Voc_index+1)-V(Voc_index))/(I(Voc_index+1)-I(Voc_index));
Rsh0 =-(V(Isc_index+1)-V(Isc_index))/(I(Isc_index+1)-I(Isc_index));

mpp = I.*V;
max_index = find(mpp==max(mpp));
Im = I(max_index);
Vm = V(max_index);

k=1.38*10^(-23);
T= 306.15;
q = 1.6012*10^(-19);
Vt = k*T/q;
Rsh=Rsh0;
n=(Vm+Im*Rs0-Voc)/(Vt*(log(Isc-Vm/Rsh0-Im)-log(Isc-Voc/Rsh)+Im/(Isc-Voc/Rsh0)));
Is = (Isc-Voc/Rsh)*exp(-Voc/(n*Vt));
Rs = Rs0-n*Vt/Is*exp(-Voc/(n*Vt));

Iph = Isc*(1+Rs/Rsh)+Is*(exp(Isc*Rs/(n*Vt))-1);

