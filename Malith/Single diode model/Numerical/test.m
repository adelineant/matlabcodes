clc;clear all;close all


[struArray,top] = datagrab();
for fileiter = [1:1:length(struArray)]
    
    
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    I = -A(:,2);
  %  I = smooth(V,I,0.05,'rloess');
    
   [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index] = lineofbestfit(V,-I);
   
   I = -I;
   
   %Vcal = (-I - Isc - (Isc/(Voc(2))).*V)./((-2/(Voc(2)))*(Isc/(Voc(2))));
   
  
   %Ical = ((-Isc/(Voc(2))).*Vcal + Isc);
   
  % logic = (Vcal > 0 & Vcal < Voc);
   
   %Vcal = Vcal(logic);
  % Ical = Ical(logic);
  %c = I - (Isc/Voc(1)).*V;
  
  %Ical = (1/2)*(Isc + c);
  
  %Vcal = (I - Isc)*(Voc(1)/-Isc);
  Voc = Voc(1)
  c = I - (Voc/Isc).*V;
  Vcal = (c - Isc)/(-(Isc/Voc +Voc/Isc));
  Ical = (-Isc/Voc)*Vcal + Isc;
  
  plot(Vcal,Ical)
   
   d = ((V -Vcal).^2 + (I - Ical).^2).^0.5;

  %Ical = (Isc/Voc).*V + Isc*(1-((2/Voc).*V));
  
  %limit V and I to get the solution you want
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
  
  maxpoint = find (d == max(d((Vcal > 0 & Vcal < Voc(1)))));
%  minRsh = -Vm/(Im - Isc);
  %N = (q/(kb))*((-Vm/Im) - (-Im/(Voc - Vm))^-1)/(1/(Isc + Im - Vm/minRsh) - 1/(Isc - Voc/minRsh))
 
 %N = (q/kb)*((2*Vm - Voc)/-Im)/(1/(Isc - Voc/minRsh) - 1/(Isc - Vm/minRsh - Im))
  %Im = -Im;

   plot(V,I,'b')
   hold on
  plot(Vcal,Ical,'r-')
   hold on
   %plot(Vm,Im,'bo')
   hold on 
   plot(V(maxpoint),I(maxpoint),'ro')
   %plot(Vcal,Ical,'*')
   
   plot(Vcal(maxpoint), Ical(maxpoint),'go')
   %ylim([0 Inf])
   %plot(Vcal(maxpoint),Ical(maxpoint),'o')
   %figure;
   %plot(Vcal,real(d))
   %plot(Ical,Vcal)
   
   
   %f = @(Vm) Isc/Voc.*Vm + Isc*(1-((2/Voc).*Vm)) - Im;
   %plot(fzero(f,Voc),Im,'*')
   
   ylim([0 Isc])
   %figure
    
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