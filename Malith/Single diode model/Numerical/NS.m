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
    
    
  

    
    [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index] = lineofbestfit(V,-I);

 
    
    
    
        %Rs0 = Voc/Isc *0.90;
    beta0=[1*10^-6,3.005510172434505e+03,900];
     

     
     b=lsqnonlin(@thisisfun,beta0,[0,3.005510172434505e+03,100],[32,1.293616711954217e+18,1000],[],V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);
    [min,Ical] = thisisfun(b,V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);
     
     %[min,Ical] = thisisfun(beta0,V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);

%for now consider only points for -I < Im & -I > 0.

%b = fzero(@thisisfun,1000000,[],V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index)
 %[min,Ical] = thisisfun(b,V(),I(),Voc,Isc,Vm,Im,kb,q,Rs0,Rsh0,Voc_index,Isc_index);

[Vm_cal,Im_cal] = mxpower(Voc_index,Isc_index,-Ical,V);
%plot both
figure;
plot(V,-Ical);
hold on
plot(V,-I);
hold on
%plot(Vm_cal,-Im_cal,'b*');
plot(Vm,Im,'r*');
hold on
plot(Vm_cal,Im_cal,'o');




%legend('Fitted IV curve','Actual IV data')
ylim([0 Inf])
end
cd(top);

function [min,Ireg] = thisisfun(beta0,V,I,Voc,Isc,Vm,Im,kb,q,Rs,Rsh,Voc_index,Isc_index)


Rs = beta0(1);
Rsh = beta0(2);
N = beta0(3);

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

