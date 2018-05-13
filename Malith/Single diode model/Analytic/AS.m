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

     
     name = sprintf('Case %d',fileiter);
     title(name);
     subplot(1,4,fileiter);
     hold on
     plot (V(1:i),I(1:i), '-r');
     hold on
     plot (V(1:i),Ical(1:i), 'b');
     plot(Vm,Im,'ro');
    
     text(0.10*Voc,Isc/2,fillfactor);
     %plot(Vm_fit,-Im_fit,'go')
     

     
     legend('Simulated Data','Model Data','Simulated data max power point','location','best')
     
     grid on 
     grid minor
     xlabel('Voltage (V)')
     ylabel('Current (A)')
     hold off
     
    ylim([0 Isc*1.2])




   
end
cd(top);

function [struArray,top] = datagrab()

top = pwd;
cd Funcfiles
struArray = datareader(top);

cd(top);
cd Funcfiles

end
