function [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index] = lineofbestfit(V,I)

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

    Voc_index = find(abs(I)==min(abs(I-0)));
    
    Voc= V(Voc_index);
    %smooth the data
    Ismooth = smoothdata(I,'sgolay');
        
   [Vm,Im] =  mxpower(Voc_index,Isc_index,Ismooth,V);


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



