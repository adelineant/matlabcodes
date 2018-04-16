clc,clear all,close all
 
%open the file. You can test with the simulation file
fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
%scan the data to A matrix. Should probably initialize this. There are only
%two columns but infite rows 
A = [fscanf(fileID,'%f',[2 Inf])]';
%close any opened file to prevent memory leaks
fclose(fileID);
%the first column SHOULD have voltage and the second column should have
%current
V = A(:,1);
I = A(:,2);

%Isc = short circuit current
%Voc = open voltage
%q = charge density in si units
%kb boltzman constant in si units
%T in kelvin

%**very important the modulus of Isc should be input.
%Isc = 0.003701;
%Voc = 0.9687;

Isc_index = find(V==min(V(V>=0)));
Isc= abs(I(Isc_index));
Voc_index = find(I==min(I(I>=0)));
Voc= V(Voc_index);

q = 1.6012*10^(-19); 
kb = 1.38*10^(-23);
T = 306.15;

%%**Note if you want to test the simulated results 
%type in the simulateddate.text to fopen 
%set Isc to 0.708
%set Voc to 0.5728


%beta is a vector that has the changing paramters. So this model there are
%three parametern :- n (=idelity factor),Rs (= series resistance)
%Rsh (=shunts resistance) 
%therefroe beta0 = [n,Rs,Rsh]
beta0=[1, 12.4491,788];
%1.5 is a good starting point
%but find a method that can guess the n value.
%call the minimize least squares function
%the first input is the model which is the single diode lambert equation
%the second are the parameter array to be change
%the [] are bunch of other options. look at the documentation. I haven't
%change anything there. The other stuff are just inputs passed into the 
%funtion
b=lsqnonlin(@thisisfun,beta0,[],[],[],V,I,Voc,Isc,T,kb,q);
%the first output is not important of our purpose
%Ical is the regressed current we calculate the modified parameters
[min,Ical] = thisisfun(b,V,I,Voc,Isc,T,kb,q);
%plot both
plot(V,Ical);
hold on
plot(V,I);

legend('Fitted IV curve','Actual IV data')


function [min,Ireg] = thisisfun(beta0,V,I,Voc,Isc,T,kb,q)

n = beta0(1);
Rs = beta0(2);
Rsh = beta0(3);

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

%ingore this it is not useful
function y = myfun(b,V,I,Voc,Isc,T,kB,qchar)

Ical=-(-qchar*V+(-lambertw(qchar*b(2)*(Isc-(Voc-b(2)*Isc)/b(3))*exp(-qchar*Voc/(b(1) *kB*T))...
    *b(3)/(b(2)*b(1)*kB*T+b(3)*b(1)*kB*T)*exp(b(3)*qchar*(b(2)*...
    (Isc+b(2)*Isc/b(3))+V)/b(1)/kB/T/(b(3)+b(2))))+b(3)*qchar*(b(2)*...
    (Isc+b(2)*Isc/b(3))+V)/b(1)/kB/T/(b(3)+b(2)))*b(1)*kB*T)/qchar/b(2);
y=Ical-I;

end
