function [Rs0,Rsh0] = lineofbestfit(V,I)

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
    I = -I;
end


%Here we call the the local function Rsfit that is a local function in the
%bestfit model
Rs0 = Rsfit(V,I,Voc_index,"Rs0");
Rsh0 = Rsfit(V,I,Isc_index,"Rsh0");


end


function grad = Rsfit(V,I,z_index,type)

%The if loop will determine which code to run.
%There may be some merit running it in a seperate code entierly
%rather than putting it inside the same function file

%For Rs0 it was found that a second order polynomial provided a good fit
%for the tested data
if (type == "Rs0")
    n = 1;

    %we define a quantity called toleV which will be the permissible data
    %range allowed across V(z_index). V(z_index) is actually VoC. 
    
    toleV = 0.01;
    
    %find the logical that meets this criteria
    %however, this line of code is subjected to failure if the input does
    %not contain data that meets the tolerance specified.
    %a while loop might need to be implemented
    %the while loop will need to consider the fit and ensure enough data
    %points
    Vtol = 1;
    Wi = 1;
    R = 0;
    P = 1;
    grad = 0;
    while (isscalar(Vtol)&& R <0.9 && P > 0.5&& grad <= 0)
        toleV = toleV*Wi;
        Vtol = find(V > (1-toleV)*V(z_index) & V < (1+toleV)*V(z_index));

        %the data points that meet the criteria
        Vdatapoint = V(Vtol);
        Idatapoint = I(Vtol);

        %polyfitting the data
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        %one issue is the gradient. It may be better to keep it as a line and
        %keep n = 1.

        figure;

        plot(Vdatapoint,Idatapoint,'o')

        hold on

        plot(polyval(Vpara,Idatapoint),Idatapoint,'k-')

        %i need to write code that will do a while loop on R while keeping vtol
        %as low as possible
        
        [R,P] = corrcoef(Vdatapoint,polyval(Vpara,Idatapoint));
        
        R = R(2,1);
        P = P(2,1);

        grad = -Vpara(1);
    end

else
    n = 1;
    %will need to overhaul this code entirely
    %take the Voc_index as the center and take 5% of the data and 10% of the
    %data and compute the difference of gradient
    
    toleI = 0.01;
    
    Itol = 1;
    %{
    
    while (isscalar(Itol))
        toleV = toleV*Wi;
        Vtol = find(V > (1-toleV)*V(z_index) & V < (1+toleV)*V(z_index));

        %the data points that meet the criteria
        Vdatapoint = V(Vtol);
        Idatapoint = I(Vtol);

        %polyfitting the data
        [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
        %one issue is the gradient. It may be better to keep it as a line and
        %keep n = 1.

        figure;

        plot(Vdatapoint,Idatapoint,'o')

        hold on

        plot(polyval(Vpara,Idatapoint),Idatapoint,'k-')

        %i need to write code that will do a while loop on R while keeping vtol
        %as low as possible

        grad = -Vpara(1);
    end
    %}
    
    Itol = find(I > (1-toleI)*I(z_index) & I < (1+toleI)*I(z_index));

    
    Vdatapoint = V(Itol);
    Idatapoint = I(Itol);

    [Vpara]= polyfit(Idatapoint,Vdatapoint,n);
    Vfit = @(IData) Vpara(1)*IData + Vpara(2);
    gradinteq =@(IData) (Vpara(1));
   
   
    figure;

    plot(Vdatapoint,Idatapoint,'o')

    hold on

    plot(polyval(Vpara,Idatapoint),Idatapoint)
    Isc_Cal = fzero(Vfit,I(z_index));
    grad = -gradinteq(Isc_Cal);
end


end



