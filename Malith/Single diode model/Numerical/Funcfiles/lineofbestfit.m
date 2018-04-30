function [Rs0,Rsh0,Voc,Isc,Im,Vm] = lineofbestfit(V,I)

%test if the user has inputted the is correct
%check if the inputs are all real and contains no imaginary values
if ( (isreal(V) || isreal(I)) == false) 
    
    error('the data contains imaginary numbers');
    
end
%Find the Isc and Voc in the data that has been given

    Isc_index = find(abs(V)==min(abs(V-0)));
    
    Isc= I(Isc_index);
    if (Isc < 0)
        %if Isc is less than 0 multiple all the currents by zero
        disp('The current values supplied have been multiplied by -1');
        I = -I;
        
        
    end
    hold on
    plot(0,Isc,'*');
    Voc_index = find(abs(I)==min(abs(I-0)));
    
    Voc= V(Voc_index);
    %smooth the data
    Ismooth = smoothdata(I,'sgolay');
    


    %Find Im and Vm from MPP tracking
    %Calculating power for the reverse bias and forward bias measurements
    %may become an issue because indexes are swapped. e.g in forward I(1)
    %will contain values near Isc. But in reverse bias I(1) will be
    %negative or closer to zero
    %First if statement is for forward bias and the else is for negative
    if (Voc_index > Isc_index )
        mapvect = Isc_index:1:Voc_index;
        mpp = abs(Ismooth(mapvect).*V(mapvect));
        max_index = find(mpp==max(mpp));
        Im = Ismooth(max_index + Isc_index);
        Vm = V(max_index + Isc_index);
    
    else 
        mapvect = Voc_index:1:Isc_index;
        mpp = abs(Ismooth(mapvect).*V(mapvect));
        max_index = find(mpp==max(mpp));
        Im = Ismooth(max_index + Voc_index);
        Vm = V(max_index +Voc_index);
    end


%Here we call the the local function Rsfit that is a local function in the
%bestfit model

    Rs0 = Rsfit(V,Ismooth,"Rs0",Vm,Im,Isc,Voc);
    Rsh0 = Rsfit(V,Ismooth,"Rsh0",Vm,Im,Isc,Voc);


end


function grad = Rsfit(V,I,type,Vm,Im,Isc,Voc)


    if   strcmp(type,'Rsh0')
        %For Rsh0 data from Vm/2 is expected to be straight.
        %Not a bad assumption
        zlogic = (V < Vm/2);
        %mgrad is not used right now but will be implemented later
        mgrad = (Im - Isc)/(Vm - 0);
        %for Rsh0 straight line is ok. hence n = 1
        n = 1;
        
        Vdatapoint = V(zlogic);
        
        Idatapoint = I(zlogic);
        
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        
        %dV/dI = -Vpara(1)
        grad = -Vpara(1);
        
    else
        zlogic = (I < Im/2);
        
        mgrad = (0 - Im)/(Voc - Vm);

        Vdatapoint = V(zlogic);
        
        Idatapoint = I(zlogic);
        %ax^2 + bx + c when n = 2
        n =2;
       
       
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
    
        grad = -(2*Vpara(1)*0 + Vpara(2));

    end
    

end



