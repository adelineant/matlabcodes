 clear all; close all; clc


%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();


for fileiter = [1:1:length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    name = struArray{fileiter}.name;
    V = A(:,1);
    I = awgn(A(:,2),50,'measured');
    



    %Defining the parameters in analytical single model. Use polynomial


    %Applying Analytical Single Model
    q = 1.6012*10^(-19); 
    k = 1.38*10^(-23);
    T = 300;

    [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I_smooth,I_orignal_data] = lineofbestfit(V,I);

    Vt = k*T/q;
    Im = -Im;
    I = -I_smooth;

    
Rsh=Rsh0;
n=(Vm+Im.*Rs0-Voc)./(Vt.*(log(Isc-Vm./Rsh-Im)-log(Isc-Voc./Rsh)+Im./(Isc-Voc./Rsh)));
Is = (Isc-Voc/Rsh)*exp(-Voc/(n*Vt));
Rs = Rs0-n*Vt/Is*exp(-Voc/(n*Vt));
Iph = Isc*(1+Rs/Rsh)+Is*(exp(Isc*Rs/(n*Vt))-1);

%Plotting extracted parameters
 Ical=[];
 %options = optimset('TolX',1e-7,'TolFun',1e-7);
for i=1:length(V(1:1:Voc_index))
    Ical=[Ical fzero(@(I)f2(I,V(i), Isc,Voc,Rs, n, Rsh, Iph),0)];
end


   ffactor = abs(Vm*Im)/(Voc*Isc) * 100;
    
    fillfactor = sprintf('Fill Factor = %0.1f%%',ffactor);

     
     
     subplot(2,3,fileiter);
     name = sprintf('Case %d',fileiter);
     title(name);
     hold on
     plot (V(1:i),-I_orignal_data(1:i), '-r');
     hold on
     plot (V(1:i),Ical(1:i), 'b');
     plot(Vm,Im,'ko');
    
   text(0.10*Voc,Isc/2,fillfactor);
     %plot(Vm_fit,-Im_fit,'go')
     grid on 
     grid minor
     xlabel('Voltage (V)')
     ylabel('Current (A)')
     
     ylim([0 Isc*1.2])
     xlim([0 Voc*1.1])
     
     if fileiter == 1
         
        legend('Simulated data with noise','Model data','Simulated data max power point','location','southwest') 
         
     end
     
     RMSE = sqrt(mean((-I_orignal_data(1:i) - Ical(1:i)).^2));  % Root Mean Squared Error
     
     extractpara = [Rs,Rsh,n,RMSE];
     
     struArray{fileiter}.extract = extractpara;




   
end
hold off

for i = 1:1:fileiter
   
    sprintf('case %i',i) 
   Rs =  struArray{i}.extract(1)
   Rsh = struArray{i}.extract(2)
   n = struArray{i}.extract(3)
   RMSE = struArray{i}.extract(4)
    
end

cd(top);

function [struArray,top] = datagrab()

top = pwd;
cd Funcfiles
struArray = datareader(top);

cd(top);
cd Funcfiles

end
