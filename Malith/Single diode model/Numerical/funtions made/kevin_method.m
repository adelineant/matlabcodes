%clc;clear all;close all


[struArray,top] = datagrab();
for fileiter = [1:1:28]%length(struArray)]
    
    
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    I = A(:,2);
    
    
   [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index] = lineofbestfit(V,-I);
   
    N = 300 ;
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);  
       
   % b = (Isc - Voc/Rsh)/(1-exp((-q*Voc)/(N*kb)));  
   % alpha = b + Voc/Rsh;
   
   c = N*kb/q;
   e = Isc/(1-exp((-q*Voc)/(N*kb)));
   f = (Voc)/(1-exp((-q*Voc)/(N*kb)))
   D = +Rs0;
    Rsh = ((c - e*D)/(D*c - D*f))^-1
    
   
  
end
cd(top);

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