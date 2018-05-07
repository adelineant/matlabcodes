clc;clear all;close all


[struArray,top] = datagrab();
for fileiter = [6:1:28]%length(struArray)]
    
    
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    I = A(:,2);
    I = smooth(V,I,0.05,'rloess');
    
   [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index] = lineofbestfit(V,-I);
   
   Vcal = (-I - Isc - (Isc/Voc).*V)./((-2/Voc)*(Isc/Voc));
   
  
   Ical = (-Isc/Voc).*Vcal + Isc;
   
  % logic = (Vcal > 0 & Vcal < Voc);
   
   %Vcal = Vcal(logic);
  % Ical = Ical(logic);
   
   d = ((V -Vcal).^2 - (I + Ical).^2).^0.5;

  %Ical = (Isc/Voc).*V + Isc*(1-((2/Voc).*V));
  
  %limit V and I to get the solution you want
  
  maxpoint = find (d == max(d((Vcal > 0 & Vcal < Voc))));
  
 
 
  
   plot(V,-I)
   hold on
  %plot(Vcal,Ical,'*')
   hold on
   plot(Vm,Im,'o')
   hold on 
   plot(Vcal,Ical,'*')
   plot(V(maxpoint), -I(maxpoint),'*')
   %plot(Vcal(maxpoint),Ical(maxpoint),'o')
   %figure;
   %plot(Vcal,real(d))
   %plot(Ical,Vcal)
   
   
   %f = @(Vm) Isc/Voc.*Vm + Isc*(1-((2/Voc).*Vm)) - Im;
   %plot(fzero(f,Voc),Im,'*')
   
   ylim([0 Inf])
    
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