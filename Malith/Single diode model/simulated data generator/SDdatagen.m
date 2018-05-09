%clc,clear all, close all
format long
for Rs = [245.982886107869];
    for Rsh = [2850.24830959519]
        for N = [60.6034974503790]%linspace(100,10000,10)
            
            
            
            Voc = 0.9694;
            Isc = 0.0037;
            V = linspace(0,Voc*2.0,1000);
            %figure ;
        
            %T = 306.15;

            %y=SD_equation(V,I,Rs,Rsh,n)
           [Ireg,Vm,Im,dvdI,dvdI2,dvdI3] = SD_equation(V,Voc,Isc,Rs,Rsh,N);
            
            %plot (V,-Ical)
           % hold on
           
           
   
            
            
            %plot(dVdI,Ireg'.');
            
            %plot(V,dVdI);
          % figure ;
           %rad = (abs(((dvdI2))./((1 + dvdI.^2).^(3/2))));
           %rad =  ((1 + dvdI.^2).^(3/2))./(abs(dvdI2));
           %not abs here. just to test the equation
           rad = ((((dvdI2))./((1 + dvdI.^2).^(3/2))));
           %plot(Ireg,rad);
           hold on
            %figure;
            plot(V,Ireg);
            hold on
            plot(Vm,Im,'o');
           
            xlabel('V');
            
            ylabel('I');
            
            fmt = sprintf("Rs=%d Rsh=%d N=%d",Rs,Rsh,N);
            
            legend(fmt,'powerpoint')
            
            ylim([0 Inf])
          
            bb = (min(rad) == rad);
            
            plot(V(bb),Ireg(bb),'*')
                            
            
            Isc_index = find(abs(V)==min(abs(V-0)));
            
            %fileID = fopen("lowshunt.txt",'w');
           %fprintf(fileID,'%6s %12s\r\n','V','Ical');
            
            A = [V;Ireg];
            figure
            plot(A(1,:),A(2,:))
            fileID = fopen("lowshunt.txt",'w');
            fprintf(fileID,'%12f %12f\r\n',A);
            fclose(fileID);
            
           % fun = @(x) (1 + 1/(Ireg(bb)^2) + 2/(-Ireg(bb)^3))/(x^2*((1+1/(Ireg(bb)^2))^2))  -  3*((1+(x/Ireg(bb))^2)^(-0.5))
            
           %fzero(fun,0.6)
       
           
            %format = sprintf("Rs=%e_Rsh=%e_n=%f.txt",Rs,Rsh,N);
            

            %fileID = fopen(format,'w');
            %fprintf(fileID,'%6s %12s\r\n','V','Ical');
            %fprintf(fileID,'%6.2f %12.8f\r\n',A);
            %fclose(fileID);
            Ireg(bb);
        end
    end
end


function [Ireg,Vm,Im,dvdI,dvdI2,dvdI3] = SD_equation(V,Voc,Isc,Rs,Rsh,N)

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
%dvdI2_comp = -(N*kb*(Rs/Rsh - dvdI/Rsh + 1))./(q*(Ireg + Isc - V/Rsh + (Rs*(Ireg + Isc))/Rsh + (N*kb)/(Rsh*q)).^2);
dvdI3 =  (2*N*kb*(Rs/Rsh - dvdI/Rsh + 1).^2)./(q*(Ireg + Isc - V/Rsh + (Rs*(Ireg + Isc))/Rsh + (N*kb)/(Rsh*q)).^3) + (N*kb*dvdI2)./(Rsh*q*(Ireg + Isc - V/Rsh + (Rs*(Ireg + Isc))/Rsh + (N*kb)/(Rsh*q)).^2);

%drad = dvdI3./(dvdI.^2 + 1).^(3/2) - (3*dvdI.*dvdI2.^2)./(dvdI.^2 + 1).^(5/2);
%plot(V,drad)


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