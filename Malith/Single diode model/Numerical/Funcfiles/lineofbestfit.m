function [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I] = lineofbestfit(V,I)

%test if the user has inputted the is correct
%check if the inputs are all real and contains no imaginary values
if ( (isreal(V) || isreal(I)) == false) 
    
    error('the data contains imaginary numbers');
    
end
%Find the Isc and Voc in the data that has been given

    Isc_index = find(abs(V)==min(abs(V-0)));
    
    Isc= I(Isc_index);
    
    if (Isc > 0)
        %if Isc is less than 0 multiple all the currents by zero
        disp('The current values supplied have been multiplied by -1');
        I = -I;
        
    end

    Voc_index = find(abs(I)==min(abs(I-0)));
    
    Voc= V(Voc_index);
    %smooth the data
    %Ismooth = smoothdata(I,'sgolay');
        
   [Vm,Im] =  mxpower(Voc_index,Isc_index,I,V);


%Here we call the the local function Rsfit that is a local function in the
%bestfit model

    
    Rsh0 = Rshfit(V,I,Vm,Im,Isc);
    Rs0 = Rsfit(V,I,Im,Rsh0,Isc);
    
    %return absolute ISC
    
    Isc = abs(Isc);

end


function grad = Rsfit(V,I,Im,Rsh0,Isc)


        zlogic = (I > Im/2);
        


        Vdatapoint = V(zlogic);
        
        Idatapoint = I(zlogic);
        %ax^2 + bx + c when n = 2
        n =2;
       
       
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        
        gradient = @(x) (2*Vpara(1)*x + Vpara(2));
        
        %dVdI
    
        grad = gradient(0);
        

    

end

function grad = Rshfit(V,I,Vm,Im,Isc)



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
        grad = Vpara(1);
        
    

end 



