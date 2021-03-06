clc;clear all;close all


[struArray,top] = datagrab();
for fileiter = [1:1:length(struArray)]
    
    
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    I = A(:,2);
  %  I = smooth(V,I,0.05,'rloess');
    
   [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I] = lineofbestfit(V,-I);
   

  Voc = Voc(1);
  Isc = Isc(1);
  c = -I + (Voc/Isc).*V;
  Vcal = (c + Isc)/((Isc/Voc +Voc/Isc));
  Ical = (Isc/Voc)*Vcal - Isc;
  
  plot(Vcal,Ical)
   
   d = ((V -Vcal).^2 + (I - Ical).^2).^0.5;

  %Ical = (Isc/Voc).*V + Isc*(1-((2/Voc).*V));
  
  %limit V and I to get the solution you want
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
  
  maxpoint = find (d == max(d((Vcal > 0 & Vcal < Voc(1)))));


   plot(V,I,'b')
   hold on
   plot(Vcal,Ical,'r-')
   hold on
   %plot(Vm,Im,'bo')
   hold on 
   plot(V(maxpoint),I(maxpoint),'ro')
   %plot(Vcal,Ical,'*')
   
   plot(Vcal(maxpoint), Ical(maxpoint),'go')

   %ylim([0 Isc])
   %figure
    
end

