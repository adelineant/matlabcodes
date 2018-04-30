clc,clear all,close all
 

%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();


for fileiter = [1:1:length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    I = A(:,2);
    %need to fix the signs of the current
    %phyical quantities
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    T = 306.15;
    
    [Rs0,Rsh0,Voc,Isc,Im,Vm] = lineofbestfit(V,-I);

    %beta 

     beta0=[1, Rs0,Rsh0];

%for now consider only points for -I < Im & -I > 0.

b=lsqnonlin(@thisisfun,beta0,[],[],[],V(-I < Im & -I > 0),I(-I < Im & -I > 0),Voc,Isc,T,kb,q);

%run with new parameters
[min,Ical] = thisisfun(b,V(),I(),Voc,Isc,T,kb,q);

%plot both
figure
plot(V,-Ical);
hold on
plot(V,-I);

legend('Fitted IV curve','Actual IV data')

end
cd(top);

function [min,Ireg] = thisisfun(beta0,V,I,Voc,Isc,T,kb,q)

n = beta0(1);
Rs = beta0(2);
Rsh = beta0(3);
n= 1;

%I broke the terms in the long lambert equation to make it easy to read and
%follow

t1 = V/Rs;

t2 = (Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/...
    (1-exp(q*(Rs*Isc - Voc)/(n*kb*T)));

t3 = (Rs*Voc)/Rsh + V;

t4 = (n*kb*T)/(q*Rs);

t5 = ((q*Rs)/(n*kb*T)) *...
    (Isc - Voc/(Rs+Rsh))*...
    exp(-(q*Voc)/(n*kb*T))/...
    (1- exp(q*(Rs*Isc - Voc)/(n*kb*T)));

Ireg = t1 - Rsh*((t2 + t3)/(Rs*(Rs+Rsh))) + t4*lambertw(t5 * exp((Rsh*q*(t2+t3))/(n*kb*T*(Rs+Rsh))));

min=Ireg-I;

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

