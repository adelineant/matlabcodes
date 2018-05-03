clc,clear all, close all

for Rs = [65];
    for Rsh = 100000%linspace(10000,20000,5)
        for N = [300]
            
            Voc = 0.9694;
            Isc = 0.0037;
            V = linspace(0,Voc*2.0,10000);
            figure ;
        
            %T = 306.15;

            %y=SD_equation(V,I,Rs,Rsh,n)
           [Ireg,Vm,Im,dvdI,dvdI2] = SD_equation(V,Voc,Isc,Rs,Rsh,N);
            
            %plot (V,-Ical)
           % hold on
   
            
            
            %plot(dVdI,Ireg'.');
            
            %plot(V,dVdI);
          % figure ;
           %rad = (abs(((dvdI2))./((1 + dvdI.^2).^(3/2))));
           %rad =  ((1 + dvdI.^2).^(3/2))./(abs(dvdI2));
           %not abs here. just to test the equation
           rad = ((((dvdI2))./((1 + dvdI.^2).^(3/2))));
           plot(V,rad);
           
            figure;
            plot(V,Ireg);
            hold on
            plot(Vm,Im,'o');
           
            xlabel('V');
            ylabel('I');
            format = sprintf("Rs=%d Rsh=%d N=%d",Rs,Rsh,N);
            legend(format,'powerpoint')
            ylim([0 Inf])
          
            bb = (find(min(rad) == rad));
           plot(V(bb),Ireg(bb),'*')
           
            Vm = V(bb);
            Ireg = -Ireg(bb);
            q = 1.6012*10^(-19);

            kb = 1.38*10^(-23);



first = 2*((Isc + Ireg - Vm/Rsh + N*kb/(q*Rsh) + ((Ireg + Isc)*Rs)/Rsh))/...
    (N*kb/q)+ 1/(Rsh +(Rs +Vm/Im)) 

second = 3*(((1 + (-Vm/Ireg)^2))^0.5)*(-Vm/Ireg)/((1-(Vm/Ireg))^1.5)
       
           
            %format = sprintf("Rs=%e_Rsh=%e_n=%f.txt",Rs,Rsh,N);
            

            %fileID = fopen(format,'w');
            %fprintf(fileID,'%6s %12s\r\n','V','Ical');
            %fprintf(fileID,'%6.2f %12.8f\r\n',A);
            %fclose(fileID);
        end
    end
end


function [Ireg,Vm,Im,dvdI,dvdI2] = SD_equation(V,Voc,Isc,Rs,Rsh,N)

%I broke the terms in the long lambert equation to make it easy to read and
%follow
q = 1.6012*10^(-19);

kb = 1.38*10^(-23);
t1 = V/Rs;

t2 = (Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/...
    (1-exp(q*(Rs*Isc - Voc)/(kb*N)));

t3 = (Rs*Voc)/Rsh + V;

t4 = (kb*N)/(q*Rs);

t5 = ((q*Rs)/(kb*N)) *...
    (Isc - Voc/(Rs+Rsh))*...
    exp(-(q*Voc)/(kb*N))/...
    (1- exp(q*(Rs*Isc - Voc)/(kb*N)));

Ireg = t1 - Rsh*((t2 + t3)/(Rs*(Rs+Rsh))) + t4*lambertw(t5 * exp((Rsh*q*(t2+t3))/(N*kb*(Rs+Rsh))));

a = (N*kb/q);
b = (Isc + Ireg - V/Rsh + N*kb/(q*Rsh) + ((Ireg + Isc)*Rs)/Rsh);
dvdI = (a./b) +Rs;
dvdI2 = -(a*(1+((Rs-dvdI)/Rsh))./b.^2);
%dvdI2 = -(N*kb*(Rs/Rsh - dvdI/Rsh + 1))./(q*(Ireg + Isc - V/Rsh + (Rs*(Ireg + Isc))/Rsh + (N*kb)/(Rsh*q)).^2);

    Isc_index = find(abs(V)==min(abs(V-0)));
    
    Isc= Ireg(Isc_index);
    
    if (Isc < 0)
        %if Isc is less than 0 multiple all the currents by zero
        disp('The current values supplied have been multiplied by -1');
        Ireg = -Ireg;
        Isc = -Isc;
        
    end
    Voc_index = find(abs(Ireg)==min(abs(Ireg-0)));

[Vm,Im] = mxpower(Voc_index,Isc_index,Ireg,V);





%[Vm_cal,Im_cal] = mxpower(Voc_index,Isc_index,Ireg,V);
%plot(Vm_cal,-Im_cal,'b*');
%

%syms dvdIequ(N,kb,q,Isc,Ireg,V,Rsh,Rs)

%dvdIequ = ((N*kb/q)./(Isc + Ireg - V/Rsh + N*kb/(q*Rsh) + (Ireg+ Isc)*Rs/Rsh)) +Rs;

%gradient(dvdIequ, Ireg);

end