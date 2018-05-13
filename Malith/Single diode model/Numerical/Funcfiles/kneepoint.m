function [Ik,Vk] = kneepoint(V,I,Voc,Isc)
  %not required as of now
  Voc = Voc(1);
  Isc = Isc(1);
  c = -I + (Voc/Isc).*V;
  Vcal = (c + Isc)/((Isc/Voc +Voc/Isc));
  Ical = (Isc/Voc)*Vcal - Isc;
  
   
  d = ((V -Vcal).^2 + (I - Ical).^2).^0.5;

   maxpoint = find (d == max(d((Vcal > 0 & Vcal < Voc(1)))));

   Vk = V(maxpoint);
   Ik= I(maxpoint);

end
