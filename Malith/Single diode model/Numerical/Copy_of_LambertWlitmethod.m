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
    %I = awgn(A(:,2),50,'measured');
    I = A(:,2);
    
    %need to fix the signs of the current
    %phyical quantities
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    T = 300;
    Vt = kb*T/q;
   
    %The lineofbestfit must return postive gradients and Isc all the time
    %The I values will always be negative as it is returned as negative
    [GradRs0,GradRsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I_smooth,I_orignal_data] = lineofbestfit(V,I);
    
    %incase you get more than one 

    [n,Rs0] = idl_Rs_cal(GradRs0,GradRsh0,Vt,Isc,Voc,Im,Vm);
    

     
     %RMSE = sqrt(mean((I_smooth(span) - I_fit(span)).^2));  % Root Mean Squared Error
     
    % extractpara = [b(1),b(2),b(3),RMSE];
     
     %struArray{fileiter}.extract = extractpara;

     
    % if fileiter == 1
         
     %   legend('Simulated Data with no noise','Model Data','Simulated data max power point','location','southwest') 
         
     %end
      
end
%hold off
%legend('Simulated Data','Model Data','MMP Simulated Data','location','best')
%for i = 1:1:fileiter
 %  sprintf('case %i',i) 
  % Rs =  struArray{i}.extract(1)
   %Rsh = struArray{i}.extract(2)
   %n = struArray{i}.extract(3)
   %RMSE = struArray{i}.extract(4)

    
%end

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

