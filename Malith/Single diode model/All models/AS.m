 clear all; close all; clc

%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();

for fileiter = [1:1:length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    %YOU NEED TO FIND A WAY A WAY TO MAKE THE PROGRAM KNOW WHAT SIGN
    I = A(:,2);
    



    %Defining the parameters in analytical single model. Use polynomial


    %Applying Analytical Single Model
    q = 1.6012*10^(-19); 
    k = 1.38*10^(-23);
    T = 300;

    [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I] = lineofbestfit(V,I);

    Vt = k*T/q;
    Im = -Im;
    I = -I;

    
Rsh=Rsh0;
n=(Vm+Im.*Rs0-Voc)./(Vt.*(log(Isc-Vm./Rsh-Im)-log(Isc-Voc./Rsh)+Im./(Isc-Voc./Rsh)));
Is = (Isc-Voc/Rsh)*exp(-Voc/(n*Vt));
Rs = Rs0-(n*Vt)/(Isc - Voc/Rsh);
Iph = Isc*(1+Rs/Rsh)+Is*(exp(Isc*Rs/(n*Vt))-1);


%Plotting extracted parameters
beta0 = [Rs,Rsh,n];
Ical = -fitter(beta0,V,Voc,Isc);


   ffactor = abs(Vm*Im)/(Voc*Isc) * 100;
    
     fillfactor = sprintf('Fill Factor = %0.1f%%',ffactor);
     subplot(4,2,fileiter);
     name = sprintf('Case %d',fileiter);
     title(name);
     hold on
     plot (V,I,'-r');
     hold on
     plot (V,Ical, 'b');
     plot(Vm,Im,'ko');
    
   text(0.50*Voc,Isc/3,fillfactor);
     %plot(Vm_fit,-Im_fit,'go')
     grid on 
     grid minor
     xlabel('Voltage (V)')
     ylabel('Current (A)')
     
     ylim([0 Isc*1.1])
     xlim([0 Voc*1.1])
     
     if fileiter == 1
         
        legend('Simulated Data with no noise','Model Data','Simulated data max power point','location','southwest') 
         
     end
     
     RMSE = sqrt(mean((I - Ical).^2));  % Root Mean Squared Error
     
     extractpara = [Rs,Rsh,n,RMSE];
     
     struArray{fileiter}.extract = extractpara;

     

   
end


cd(top);

function [struArray,top] = datagrab()

top = pwd;
cd Funcfiles
struArray = datareader(top);

cd(top);
cd Funcfiles

end

function Ireg = fitter(beta0,V,Voc,Isc)


Rs = beta0(1);
Rsh = beta0(2);
N = beta0(3)*300;

q = 1.6012E-19;

kb = 1.38*10^(-23);

Ireg = V/Rs - Rsh*(((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/(1-exp(q*(Rs*Isc - Voc)/(kb*N))) + ...
    (Rs*Voc)/Rsh + V)/(Rs*(Rs+Rsh))) + (kb*N)/(q*Rs)...
    *lambertw(((q*Rs)/(kb*N))*(Isc - Voc/(Rs+Rsh))*exp(-(q*Voc)/(kb*N))...
    /(1- exp(q*(Rs*Isc - Voc)/(kb*N)))*...
    exp((Rsh*q*((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))...
    /(1-exp(q*(Rs*Isc - Voc)/(kb*N)))+(Rs*Voc)/Rsh + V))/(N*kb*(Rs+Rsh))));

end
