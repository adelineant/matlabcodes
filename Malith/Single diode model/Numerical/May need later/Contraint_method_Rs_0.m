
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
   
   plot(V,I);

       
    [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I] = lineofbestfit(V,I);
    [n_min] = fzero(@Rshcalmeth,300,[],V,Vm,Im,Voc,Isc,kb,q,Rs0,Rsh0);

ylim([0 Inf])
end
cd(top);

function [min,Ireg] = thisisfun(beta0,V,I,Voc,Isc,Vm,Im,kb,q,Rs,Rsh,Voc_index,Isc_index)

N = beta0(3);
Rs = beta0(1);
Rsh = beta0(2);


%I broke the terms in the long lambert equation to make it easy to read and
%follow
Ireg = V/Rs - Rsh*(((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/(1-exp(q*(Rs*Isc - Voc)/(kb*N))) + ...
    (Rs*Voc)/Rsh + V)/(Rs*(Rs+Rsh))) + (kb*N)/(q*Rs)...
    *lambertw(((q*Rs)/(kb*N))*(Isc - Voc/(Rs+Rsh))*exp(-(q*Voc)/(kb*N))...
    /(1- exp(q*(Rs*Isc - Voc)/(kb*N)))*...
    exp((Rsh*q*((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))...
    /(1-exp(q*(Rs*Isc - Voc)/(kb*N)))+(Rs*Voc)/Rsh + V))/(N*kb*(Rs+Rsh))));

min = Ireg - I;
end 

function [Ireg] = SD_equation(V,Voc,Isc,Rs,Rsh,N)

%I broke the terms in the long lambert equation to make it easy to read and
%follow
q = 1.6012*10^(-19);

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

function [minimmum] = Rshcalmeth(N,V,Vm,Im,Voc,Isc,kb,q,Rs0,Rsh0)

a = (N*kb)/q;
c = (1 - exp((-q*Voc)/(N*kb)));

D = Rs0;

Rshcal = (((a/D) - (Isc/c))/(-Voc/c + a))^(-1);

Pmax = Im*Vm;

[Ireg] = SD_equation(V,Voc,Isc,1e-8,Rshcal,N);


    Isc_index = find(abs(V)==min(abs(V-0)));
    
    Isc= Ireg(Isc_index);
    
    if (Isc > 0)
        %if Isc is less than 0 multiple all the currents by zero
        disp('The current values supplied have been multiplied by -1');
        Ireg = -Ireg;
        
    end

    Voc_index = find(abs(Ireg)==min(abs(Ireg-0)));
    
    Voc= V(Voc_index);
    %smooth the data
    %Ismooth = smoothdata(I,'sgolay');
        
   [Vm2,Im2] =  mxpower(Voc_index,Isc_index,Ireg,V);


    
minimmum = Pmax - Vm2*Im2;

%Rshcal2 = (((a/Rsh0) - Isc*(1/c -1))/(1/a - Voc/a))^(-1);


end





 

