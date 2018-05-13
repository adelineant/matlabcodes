clc,clear all,close all
 
%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();



for fileiter = [1:1:length(struArray) ]%length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    %some I's will be negative or postive data it doesn't matter the
    %program will identify it
    A = struArray{fileiter}.data;
    name = struArray{fileiter}.name;
    V = A(:,1);
    I = -A(:,2);
    %need to fix the signs of the current
    %phyical quantities
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    T = 300;
    Vt = kb*T/q;
   
    %The lineofbestfit must return postive gradients and Isc all the time
    %The I values will always be negative as it is returned as negative
    [GradRs0,GradRsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I] = lineofbestfit(V,I);
    
    %incase you get more than one 

    [n,Rs] = idl_Rs_cal(GradRs0,GradRsh0,Vt,Isc,Voc,Im,Vm);
    

    %in some cases Rs may become negative and n may not be feasible.
    if Rs < 0
        
        Rs = 0.001;
        
    end
    
   if n < 1 
       
       % n = 1;
    end
 
    if GradRsh0 < 0 
       
        %GradRsh0 = Voc/Isc*1.5;
        
    end
  
 
    beta0=[Rs,GradRsh0,n];
    
    % options = optimset('Display','none',...
    %'TolX',1e-8,...
    %'TolFun',1e-8,...
    %'Algorithm','trust-region-reflective');

    b=lsqnonlin(@errorfuntion,beta0,[],[],[],V,I,Voc,Isc);
    I_fit = fitter(b,V,Voc,Isc);
    
    %VERY IMPORTANT ALL THE I VALUES MUST BE NEGATIVE AT ISC FOR THIS TO
    %WORK
    ffactor = abs(Vm*Im)/(Voc*Isc) * 100;
    
    fillfactor = sprintf('Fill Factor = %0.1f%%',ffactor);
     subplot(1,4,fileiter);
     
     name = sprintf('Case %d',fileiter);
     title(name);
     hold on
     plot(V,-I,'r-');
     plot(V,-I_fit,'b-');
     plot(Vm,-Im,'ro');
    
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



function [min] = errorfuntion(beta0,V,I,Voc,Isc)

Rs = beta0(1);
Rsh = beta0(2);
N = beta0(3)*300;

q = 1.6012*10^(-19);

kb = 1.38*10^(-23);

Ireg = V/Rs - Rsh*(((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/(1-exp(q*(Rs*Isc - Voc)/(kb*N))) + ...
    (Rs*Voc)/Rsh + V)/(Rs*(Rs+Rsh))) + (kb*N)/(q*Rs)...
    *lambertw(((q*Rs)/(kb*N))*(Isc - Voc/(Rs+Rsh))*exp(-(q*Voc)/(kb*N))...
    /(1- exp(q*(Rs*Isc - Voc)/(kb*N)))*...
    exp((Rsh*q*((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))...
    /(1-exp(q*(Rs*Isc - Voc)/(kb*N)))+(Rs*Voc)/Rsh + V))/(N*kb*(Rs+Rsh))));


min = Ireg - I;
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

