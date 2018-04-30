function [Rs0,Rsh0,Voc,Isc,Im,Vm] = lineofbestfit(V,I)

%test if the user has inputted the is correct
%check if the inputs are all real and contains no imaginary values
if ( (isreal(V) || isreal(I)) == false) 
    
    error('the data contains imaginary numbers');
    
end
%Find the Isc and Voc in the data that has been given
%It may be redundant here as this was written in the main code
%But i think maybe it is better to pass this instead of being calculated
%here again
    Isc_index = find(abs(V)==min(abs(V-0)));
    Isc= I(Isc_index);
    Voc_index = find(abs(I)==min(abs(I-0)));
    Voc= V(Voc_index);
    Ismooth = smoothdata(I,'sgolay');

    %Find Im and Vm from MPP tracking
    mpp = abs(I([1:1:Voc_index]).*V([1:1:Voc_index]));
    max_index = find(mpp==max(mpp));
    Im = I(max_index);
    Vm = V(max_index);
   % plot(V,Ismooth,'b.',V,I,'r-')

% this block of code will check if there is any negative voltages
% this is common in experimental data but not so in the simulated data
% this code cause problems if we dont take notice of it

%negV is an logical array. It will be a vector containing all the negative
%number
%but is it useful to have the negative numbers? is one enough?
    negV = V(V < 0);
    
%if all the values are postive then the array will be empty and enter the
%loop. In his case the array only contains postive values
if isempty(negV) == 1
     
    disp('no negative voltage');
    %i defined it as 0. Maybe for future use
    negV = 0;

%just having else should be enough. There is no point of accounting for
%else if as the negV cannot be positive under any case. 

else
   
    sVvector = length(V);
    
    snegV = length(negV);
    
    nfrac = snegV/sVvector;
    
    %This is to check if the data is good or not. May be useful in the
    %future
    
    if(nfrac > 0.5)
        
        disp('More than 50% of the voltage is  negative data')
        
    end
    
end

%here we check of the current. The current at short circuit must be higher
%than the current at the end of the run. But one issue is that direction of
%measure will affect the sign of the data. One way to check this if the ISC
%current is less than 0
if (Isc < 0)

    disp('The current values supplied have been multiplied by -1');
    Ismooth = -Ismooth;
end


%Here we call the the local function Rsfit that is a local function in the
%bestfit model
Rs0 = Rsfit(V,Ismooth,Voc_index,"Rs0",Vm,Im,Isc,Voc);
Rsh0 = Rsfit(V,Ismooth,Isc_index,"Rsh0",Vm,Im,Isc,Voc);


end


function grad = Rsfit(V,I,z_index,type,Vm,Im,Isc,Voc)

%The if loop will determine which code to run.
%There may be some merit running it in a seperate code entierly
%rather than putting it inside the same function file

%For Rs0 it was found that a second order polynomial provided a good fit
%for the tested data
    n = 1;

    %we define a quantity called toleV which will be the permissible data
    %range allowed across V(z_index). V(z_index) is actually VoC. 
    
    tolerance = 0.01;
    
    %find the logical that meets this criteria
    %however, this line of code is subjected to failure if the input does
    %not contain data that meets the tolerance specified.
    %a while loop might need to be implemented
    %the while loop will need to consider the fit and ensure enough data
    %points
    
    %zlogic will be the logical that contains indices which meet the critea
    %which we are concerned with
    zlogic = 1;
    %define the R,P,grad and offset as follows otherwise it wont enter the loop
    R = 0;
    P = 1;
    grad = 0;
    offset = 0;
    %compare if the type is correct. If type is Rs0 then we allow the z to
    %be V which will find the Voc in the subsequent. If z = I we will be
    %find the Isc for Rsh0
    if   strcmp(type,'Rsh0')
        
        z = V(V < Vm/2);
        mgrad = (Im - Isc)/(Vm - 0);
       % plot(z,I(V < Vm/2));
        n = 1;
        %polyfitting the data
        Vdatapoint = z;
        Idatapoint = I(V < Vm/2);
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        %one issue is the gradient. It may be better to keep it as a line and
        %keep n = 2.
      % hold on
        %plot(polyval(Vpara,Idatapoint),Idatapoint,'k-')
        grad = -Vpara(1);
        
    else
        z = I(I < Im/2);
        mgrad = (0 - Im)/(Voc - Vm);
      %  figure;
       % plot(V(I < Im/2),z);
        Vdatapoint = V(I < Im/2);
        Idatapoint = z;
        n =2;
        %polyfitting the data
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
     
        %one issue is the gradient. It may be better to keep it as a line and
        %keep n = 2.
       % hold on
       % plot(polyval(Vpara,Idatapoint),Idatapoint,'k-')
        
        grad = -(2*Vpara(1)*0 + Vpara(2));

    end
    
    %take the data points to 1/2 of the data at Vm
    
    
    
    %{
    while (isscalar(zlogic)|| R <0.99 || P > 0.05 || grad <= 0)
        offset = offset + tolerance;
        zlogic = find(z > (1-offset)*z(z_index) & z < (1+offset)*z(z_index));

        %the data points that meet the criteria
        Vdatapoint = V(zlogic);
        Idatapoint = I(zlogic);

        %polyfitting the data
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        %one issue is the gradient. It may be better to keep it as a line and
        %keep n = 1.

        figure;

       % plot(Vdatapoint,Idatapoint,'o')

       % hold on

       % plot(polyval(Vpara,Idatapoint),Idatapoint,'k-')

        %i need to write code that will do a while loop on R while keeping vtol
        %as low as possible
        
        [R,P] = corrcoef(Vdatapoint,polyval(Vpara,Idatapoint));
        
        R = R(2,1);
        P = P(2,1);

        grad = -Vpara(1);
 
        %for now this will exit the Rsh0 loop
        if((length(zlogic)/length(V)) > 0.05)
            break 
        end
    end

    %}

end



