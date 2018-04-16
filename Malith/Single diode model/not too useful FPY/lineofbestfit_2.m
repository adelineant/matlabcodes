function [Rs0,Rsh0] = lineofbestfit(V,I)

%test if the user has inputted the is correct
if ( (isreal(V) || isreal(I)) == false) 
    
    error('the data contains imaginary numbers');
    
end
    Isc_index = find(abs(V)==min(abs(V-0)));
    Isc= I(Isc_index);
    Voc_index = find(abs(I)==min(abs(I-0)));
    Voc= V(Voc_index);
    negV = V(V < 0);

if isempty(negV) == 1
     
    disp('no negative voltage');
    negV = 0;


elseif (negV(1) < 0)
   
    sVvector = length(V);
    snegV = length(negV);
    nfrac = snegV/sVvector;
    

    
    if(nfrac > 0.5)
        
        disp('More than 50% of the voltage is  negative data')
        
    end
    
end


if (I(end) > I(1))

    disp('The current values supplied have been multiplied by -1');
    %I = -I;
end


%calculate Rs
Rs0 = Rs0fit(V,I,Voc_index,"Rs0");
Rsh0 = Rsh0fit(V,I,Isc_index,"Rsh0",Rs0,Voc_index);

%{
Rs0 = -(V(Voc_index+1)-V(Voc_index))/(I(Voc_index+1)-I(Voc_index));
Rsh0 =(-V(Isc_index+1)-V(Isc_index))/(I(Isc_index+1)-I(Isc_index));

%}


end


function grad = Rs0fit(V,I,z_index,type)

if (type == "Rs0")
    n = 2;
%%%continue from negV
    %take the Voc_index as the center and take 5% of the data and 10% of the
    %data and compute the difference of gradient
    toleV = 0.001;
    
    Vtol = find(V > (1-toleV)*V(z_index) & I < (1+toleV)*V(z_index));
    
    Vdatapoint = V(Vtol);
    Idatapoint = I(Vtol);

    [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
    Vfit = @(IData) Vpara(1)*IData.^2 + Vpara*IData;
    gradinteq =@(IData) 2*Vpara(1)*IData + Vpara(2);
   
    figure;

    plot(Vdatapoint,Idatapoint,'o')

    hold on

    plot(polyval(Vpara,Idatapoint),Idatapoint)

    grad = -gradinteq(0);

else
    %n = 1;
    %a = round(z_index*rn)
    %b = 2*z_index - (floor(z_index*rn))
    %take the Voc_index as the center and take 5% of the data and 10% of the
    %data and compute the difference of gradient
    
    %toleI = 0.01;
    
    %Itol = find(I > (1-toleI)*I(z_index) & I < (1+toleI)*I(z_index));

    
    %Vdatapoint = V(Itol);
   % Idatapoint = I(Itol);

    %[Vpara]= polyfit(Idatapoint,Vdatapoint,n);
    %Vfit = @(IData) Vpara(1)*IData + Vpara(2);
    %gradinteq =@(IData) (Vpara(1));
   
   
    %figure;

    %plot(Vdatapoint,Idatapoint,'o')

    %hold on

    %plot(polyval(Vpara,Idatapoint),Idatapoint)
    %Isc_Cal = fzero(Vfit,I(z_index));
    %grad = -gradinteq(Isc_Cal);
end

%Rsh0 =(-V(Isc_index+1)-V(Isc_index))/(I(Isc_index+1)-I(Isc_index));
%Rs0 = -(V(Voc_index+1)-V(Voc_index))/(I(Voc_index+1)-I(Voc_index));
%Rsh0 =(-V(Isc_index+1)-V(Isc_index))/(I(Isc_index+1)-I(Isc_index));


end

function grad = Rsh0fit(V,I,z_index,type,Rsh0,Voc_index)

    n = 2;
    mpp = abs(I([1:1:Voc_index]).*V([1:1:Voc_index]));
    max_index = find(mpp==max(mpp));
    toleI = 0.05;
    Itol = find(I > (1-toleI)*I(max_index) & I < (1+toleI)*I(max_index));
    Vdatapoint = V(Itol);
    Idatapoint = I(Itol);
    
    [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
    Vfit = @(IData) Vpara(1)*IData + Vpara(2);
    
       
    figure;

    plot(Vdatapoint,Idatapoint,'o')

    hold on

    plot(polyval(Vpara,Idatapoint),Idatapoint)
    %Isc_Cal = fzero(Vfit,I(z_index));
    %grad = -gradinteq(Isc_Cal);
    %gradinteq =@(IData) (Vpara(1));
    
    %n = 1;
    %a = round(z_index*rn)
    %b = 2*z_index - (floor(z_index*rn))
    %take the Voc_index as the center and take 5% of the data and 10% of the
    %data and compute the difference of gradient
    
    %toleI = 0.01;
    
    %Itol = find(I > (1-toleI)*I(z_index) & I < (1+toleI)*I(z_index));

    
    %Vdatapoint = V(Itol);
   % Idatapoint = I(Itol);

    %[Vpara]= polyfit(Idatapoint,Vdatapoint,n);
    %Vfit = @(IData) Vpara(1)*IData + Vpara(2);
    %gradinteq =@(IData) (Vpara(1));
   
   
    %figure;

    %plot(Vdatapoint,Idatapoint,'o')

    %hold on

    %plot(polyval(Vpara,Idatapoint),Idatapoint)
    %Isc_Cal = fzero(Vfit,I(z_index));
    %grad = -gradinteq(Isc_Cal);
end

%Rsh0 =(-V(Isc_index+1)-V(Isc_index))/(I(Isc_index+1)-I(Isc_index));
%Rs0 = -(V(Voc_index+1)-V(Voc_index))/(I(Voc_index+1)-I(Voc_index));
%Rsh0 =(-V(Isc_index+1)-V(Isc_index))/(I(Isc_index+1)-I(Isc_index));


