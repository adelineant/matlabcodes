 clear all; close all; clc

%Opening data file
fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1);
I = -A(:,2);


%Defining the parameters in analytical single model. Use polynomial

Isc_index = find(abs(V)==min(abs(V-0)));
Isc= I(Isc_index);
Voc_index = find(abs(I)==min(abs(I-0)));
Voc= V(Voc_index);
%Find Im and Vm from MPP tracking
mpp = abs(I([1:1:Voc_index]).*V([1:1:Voc_index]));
max_index = find(mpp==max(mpp));
Im = I(max_index);
Vm = V(max_index);
%Applying Analytical Single Model
k=1.38E-23;
T= 306.15;
q = 1.6012E-19;

%[Rs0,Rsh0] = lineofbestfit(V,I);
%Rsh0 =-(V(Isc_index+1)-V(Isc_index))/(I(Isc_index+1)-I(Isc_index));
%Rs0 = -(V(Voc_index+1)-V(Voc_index))/(I(Voc_index+1)-I(Voc_index));
Rsh0 = 4.0242e+03;
Rs0 = 36.5178;
Vt = k*T/q;
Rsh=Rsh0;
n=(Vm+Im.*Rs0-Voc)./(Vt.*(log(Isc-Vm./Rsh-Im)-log(Isc-Voc./Rsh)+Im./(Isc-Voc./Rsh)));
Is = (Isc-Voc/Rsh)*exp(-Voc/(n*Vt));
Rs = Rs0-n*Vt/Is*exp(-Voc/(n*Vt));
Iph = Isc*(1+Rs/Rsh)+Is*(exp(Isc*Rs/(n*Vt))-1);

%Plotting using the extracted parameters
I2= fzero(@(I)f2(I,V(1),Is, Rs, n, Rsh, Iph),Isc);
for i=1:length(V)-1
    
   
        %I2=[I2 fzero(@(I)f2(I,V(i+1),Is, Rs, n, Rsh, Iph),I2(i)*2)]
    
        I2=[I2 fzero(@(I)f2(I,V(i+1),Is, Rs, n, Rsh, Iph),I2(i))];
        
        if(isnan(I2(i+1))== 1)
            break
        end
        
    
    
end
figure;
plot (V,I, 'ro');
hold on
plot (V(1:i),I2(1:i), 'b');

legend ('experimental', 'fitted')
xlabel ('V')
ylabel ('I')

