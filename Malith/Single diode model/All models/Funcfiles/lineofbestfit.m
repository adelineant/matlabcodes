function [Rs0,Rsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I_smooth,I_orignal_data,MPOC] = lineofbestfit(V,I)

%test if the user has inputted the is correct
%check if the inputs are all real and contains no imaginary values



if ( (isreal(V) || isreal(I)) == false) 
    
    error('the data contains imaginary numbers');
    
end
%Find the Isc and Voc in the data that has been given

    I_orignal_data = I;
    
   % plot(V,-I_orignal_data,'LineWidth',1.5)
    %hold on
       
    I_smooth = smoothdata(I_orignal_data,1,'rlowess',10);

    Isc_ind = find(abs(V)==min(abs(V-0)));
    
    Isc_index = Isc_ind(1);
    
    Isc= I_smooth(Isc_index);

    if (Isc > 0)
        %if Isc is less than 0 multiple all the currents by zero as the
        %equations do not work out properly otherwise
       
        I_orignal_data = -I;
        I_smooth = -I_smooth;
        Isc = -Isc;
    end
    

    Voc_ind = find(abs(I_smooth)==min(abs(I_smooth-0)));
    %incase two or more values exist
    Voc_index = Voc_ind(1);
    
    Voc= V(Voc_index);
    %smooth the data
    %Ismooth = smoothdata(I,'sgolay');
        
   [Vm,Im,MPOC] =  mxpower(Voc_index,Isc_index,I_smooth,V);


%Here we call the the local function Rsfit that is a local function in the
%bestfit model

    
    Rsh0 = Rshfit(V,I_smooth,Vm,Im,Isc,Voc);
    Rs0 = Rsfit(V,I_smooth,Im,Voc);
    
    %return absolute ISC
    
    Isc = abs(Isc);
    
   %{
    plot(Vm,-Im,'ko','MarkerFaceColor','black');
    ylim([0 Isc*1.2])
    xlim([0 Voc*1.1])
    legend('Simulated data case 6','2nd order polynomial 1','2nd order polynomial 2','MPP','location','southwest') 
    grid on
    grid minor
    title('Two 2nd order polynomials drawn on top of the case 6 simulated data');
      xlabel('Voltage (V)')
     ylabel('Current (A)')
     %}
end


function grad = Rsfit(V,I,Im,Voc)

        %you need to take points near Voc for this to make sense
        %
        zlogic = (I >= Im/2 & V <= Voc);
        


        Vdatapoint = V(zlogic);
        
        Idatapoint = I(zlogic);
        %ax^2 + bx + c when n = 2
        n =2;
       
       
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        
        gradient = @(x) (2*Vpara(1)*x + Vpara(2));
        
        %dVdI
    
        grad = gradient(0);
        %plot(polyval(Vpara,Idatapoint),-Idatapoint,'LineWidth',1.5)
    

end

function grad = Rshfit(V,I,Vm,Im,Isc,Voc)

zlogic = (V > 0 & V < Vm*0.8);
       
        %assume second order polynomials
        n = 2;
        
        Vdatapoint = V(zlogic);
        
        Idatapoint = I(zlogic);
        
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        
        
        %gradient = @(x) Vpara(1);
        gradient = @(x) (2*Vpara(1)*x + Vpara(2));
        
        grad = gradient(Isc);
       

        %line of code to test 
       if (grad < 0 || (abs(Im) > abs(Isc)))
           
           Im_temp = Im;
            if abs(Im) > abs(Isc)
                
                Im_temp = Isc*0.98;
                
            end
             
           
            grad = abs(Vm/(Isc - Im_temp))*1.2;
            
       end    
       % plot(polyval(Vpara,Idatapoint),-Idatapoint,'LineWidth',1.5)

       

end 



