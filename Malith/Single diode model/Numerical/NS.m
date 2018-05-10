clc,clear all,close all
 
%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();


for fileiter = [1:1:length(struArray) ]%length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    %some I's will be negative or postive data it doesn't matter the
    %program will identify it
    A = struArray{fileiter}.data;
    V = A(:,1);
    I = A(:,2);
    %need to fix the signs of the current
    %phyical quantities
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    T = 300;
    Vt = kb*T/q;
    
    

    %The lineofbestfit must return postive gradients and Isc all the time
    %The I values will always be negative as it is returned as negative
    [GradRs0,GradRsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I] = lineofbestfit(V,I);

    [n,Rs] = idl_Rs_cal(GradRs0,GradRsh0,Vt,Isc,Voc,Im,Vm);
    
    %in some cases Rs may become negative and n may not be feasible.
    if Rs < 0
        
        Rs = 1*10^-5;
        
    end
    
    if n < 0 || n > 3
       
        n = 1;
    end
 
    beta0=[Rs,GradRsh0,n*300];

    
    %VERY IMPORTANT ALL THE I VALUES MUST BE NEGATIVE AT ISC FOR THIS TO
    %WORK
    
      b=lsqnonlin(@thisisfun,beta0,[],[],[],V,I,Voc,Isc,Vm,Im,kb,q,Voc_index,Isc_index);

      [min,Ical] = thisisfun(b,V,I,Voc,Isc,Vm,Im,kb,q,Voc_index,Isc_index);
    [Vm_cal,Im_cal] = mxpower(Voc_index,Isc_index,Ical,V);
    %plot both
    figure;
    plot(V,-Ical);
    hold on
    plot(V,-I);
    hold on
    %plot(Vm_cal,-Im_cal,'b*');
    plot(Vm,-Im,'r*');
    hold on
    plot(Vm_cal,-Im_cal,'o');




    %legend('Fitted IV curve','Actual IV data')
    ylim([0 Inf])
end
cd(top);

function [min,Ireg] = thisisfun(beta0,V,I,Voc,Isc,Vm,Im,kb,q,Voc_index,Isc_index)

Rs = beta0(1);
Rsh = beta0(2);
N = beta0(3);

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


min = Ireg - I;
end 

function [struArray,top] = datagrab()
%top is the current directory
top = pwd;
%go to Funcfiles
cd Funcfiles
%get the data
struArray = datareader(top);
%{return to back to Funcfiles since datareader will direct you to Datafiles
cd(top);
cd Funcfiles

end

