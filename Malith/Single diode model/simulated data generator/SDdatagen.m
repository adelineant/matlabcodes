clc,clear all, close all
format long

Voc = 1;
Isc = (20/1000);
V = linspace(0,Voc*1.1,1000);

T = 300;
q = 1.6012*10^(-19);

kb = 1.38*10^(-23);


for Rs = [0.1]%{ 0.1*(Voc/Isc),0.4*(Voc/Isc),0.91*(Voc/Isc)]
    for Rsh = [200]
        for n = [2.5]%linspace(100,10000,10)
            
            I0 = (Isc - Voc/Rsh)*exp(-q*Voc/(n*T*kb))/(1 - exp(-q*Voc/(n*T*kb)));

            %addition = I0*(Voc/Rsh + exp(q*Voc/(N*kb)) - 1 );
            
            %thier = (Isc + (Rs*Isc - Voc)/Rsh)/(1- exp(q*(Rs*Isc - Voc)/(N*kb))) + Voc/Rsh;

            
           [Ireg] = SD_equation(V,Voc,Isc,Rs,Rsh,n,T);
           
         

           
            
            [Vm,Im] = mxpower(Ireg,V);
    
            
            plot(V,-Ireg,'-',Vm,-Im,'bo')
            hold on
            %plot(Vk,-Ik,'ro')
            
           
            xlabel('Voltage/ V');
            
            ylabel('Current/ A');
            
            
            grid on
            grid minor
           
            ylim([0 Inf])
     
            A = [V;Ireg];
            name = sprintf("Rs=%e Rsh=%e n=%f.txt",Rs,Rsh,n);
            fileID = fopen(name,'w');
            fprintf(fileID,'%12f %12f\r\n',A);
            fclose(fileID);
  
        end
    end
end


function [Ireg] = SD_equation(V,Voc,Isc,Rs,Rsh,n,T)

%I broke the terms in the long lambert equation to make it easy to read and
%follow
q = 1.6012*10^(-19);

kb = 1.38*10^(-23);

Ireg = V/Rs - Rsh*(((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/(1-exp(q*(Rs*Isc - Voc)/(kb*n*T))) + ...
    (Rs*Voc)/Rsh + V)/(Rs*(Rs+Rsh))) + (kb*n*T)/(q*Rs)...
    *lambertw(((q*Rs)/(kb*n*T))*(Isc - Voc/(Rs+Rsh))*exp(-(q*Voc)/(kb*n*T))...
    /(1- exp(q*(Rs*Isc - Voc)/(kb*n*T)))*...
    exp((Rsh*q*((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))...
    /(1-exp(q*(Rs*Isc - Voc)/(kb*n*T)))+(Rs*Voc)/Rsh + V))/(n*T*kb*(Rs+Rsh))));
end


 function [Vm,Im] = mxpower(I,V)
    %Find Im and Vm from MPP tracking
    %Calculating power for the reverse bias and forward bias measurements
    %may become an issue because indexes are swapped. e.g in forward I(1)
    %will contain values near Isc. But in reverse bias I(1) will be
    %negative or closer to zero
    %First if statement is for forward bias and the else is for negative
    Isc_index = find(abs(V)==min(abs(V-0)));
    
    Isc= I(Isc_index);
    
    if (Isc > 0)
        %if Isc is less than 0 multiple all the currents by zero
        disp('The current values supplied have been multiplied by -1');
        I = -I;
        
    end

    Voc_index = find(abs(I)==min(abs(I-0)));

    
    if (Voc_index > Isc_index )
        mapvect = Isc_index:1:Voc_index;
        mpp = abs(I(mapvect).*V(mapvect));
        max_index = find(mpp==max(mpp));
        Im = I(max_index + Isc_index);
        Vm = V(max_index + Isc_index);

    else 
        mapvect = Voc_index:1:Isc_index;
        mpp = abs(I(mapvect).*V(mapvect));
        max_index = find(mpp==max(mpp));
        Im = I(max_index + Voc_index);
        Vm = V(max_index +Voc_index);
 
    end
    
 
        
      


end 
