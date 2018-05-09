clc;clear all;close all
clc,clear all,close all
 

%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();


for fileiter = [1:1:length(struArray) ]%length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    I = A(:,2);
    %need to fix the signs of the current
    %phyical quantities
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    
    N = 1000 ;

   % b = (Isc - Voc/Rsh)/(1-exp((-q*Voc)/(N*kb)));  
   % alpha = b + Voc/Rsh;
   

       
    [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index] = lineofbestfit(V,-I);
    [a ,b] = Rshcalmeth(N,V,I,Voc,Isc,kb,q,Rs0,Rsh0)
    %logic = V < Vm;
    %b = lsqnonlin(@Rshcalmeth,N,[],[],[],V(logic),I(logic),Voc,Isc,kb,q,Rs0)

   
    %beta0=[1*10^-6,Rshcal,N];
    
    plot(V,-I);
    hold on
    
    b=lsqnonlin(@thisisfun,beta0,[],[],[],V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);
     [min,Ical] = thisisfun(b,V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);
        
        
          
          plot(V,-Ical);
   
% [Ireg,Vm1,Im1,dvdI1,dvdI21,dvdI31] = SD_equation(V,Voc,Isc,1*10^-6,Rshcal,N);
    % plot(V,Ireg);
    % ylim([0 Inf])
   
     
     %[Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index] = lineofbestfit(V,-I);
    %beta0=[1*10^-6,Rshcal,1000];
    % b=lsqnonlin(@thisisfun,beta0,[],[],[],V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);
    % [min,Ical] = thisisfun(b,V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);
       % plot(V,-I);


%legend('Fitted IV curve','Actual IV data')
ylim([0 Inf])
end
cd(top);

function [min,Ireg] = thisisfun(beta0,V,I,Voc,Isc,Vm,Im,kb,q,Rs,Rsh,Voc_index,Isc_index)

N = beta0(3);
Rs = beta0(1);
Rsh = beta0(2);


%I broke the terms in the long lambert equation to make it easy to read and
%follow

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

function [Rshcal , Rshcal2] = Rshcalmeth(N,V,I,Voc,Isc,kb,q,Rs0,Rsh0)

a = N*kb/q;
c = (1 - exp(-q*Voc/(N*kb)));

Rshcal = (((a/Rs0) - (Isc/c))/(-Voc/c + a))^(-1)
    
%dvdi = (a/(Isc*(1/c -1) - 1/Rshcal*(Voc/c - 1/a)))

Rshcal2 = (((a/Rsh0) - Isc*(1/c -1))/(1/a - Voc/a))^(-1)
    
    


end

