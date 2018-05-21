clc,clear all,close all
 
%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();



for fileiter = [1:1:length(struArray) ]%length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    %some I's will be negative or postive data it doesn't matter the
    %program will identify it
    A = struArray{fileiter}.data;
    name = struArray{fileiter}.name;
    V = A(:,1);
    %I = awgn(A(:,2),45,'measured');
    I = A(:,2);
    
    %need to fix the signs of the current
    %phyical quantities
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    T = 300;
    Vt = kb*T/q;
   
    %The lineofbestfit must return postive gradients and Isc all the time
    %The I values will always be negative as it is returned as negative
    [GradRs0,GradRsh0,Voc,Isc,Im,Vm,Voc_index,Isc_index,I_smooth,I_orignal_data] = lineofbestfit(V,I);
    
    %incase you get more than one 
    
     ffactor = abs(Vm*Im)/(Voc*Isc) * 100;
    %if ffactor < 26
       
     %   disp('Unable to chracterize')
      %  continue; 
   % end

    [n,Rs0] = idl_Rs_cal(GradRs0,GradRsh0,Vt,Isc,Voc,Im,Vm);
    
    
    
    Vt = kb*T/q;
   
    %I = -I; 
   % n=(Vm-Im.*Rs0-Voc)./(Vt.*(log(Isc-Vm./GradRsh0+Im)-log(Isc-Voc./GradRsh0)-Im./(Isc-Voc./GradRsh0)));
    %Is = (Isc-Voc/GradRsh0)*exp(-Voc/(n*Vt));
    %Rs = Rs0-n*Vt/Is*exp(-Voc/(n*Vt));
    

    %in some cases Rs may become negative and n may not be feasible.

    %beta0=[Rs0,GradRsh0,n];
    
    %Rsmax = (Vm/-Im) - (2*Vm - Voc)/(-Im/(Isc + Im) + log(Isc+ Im)/Isc)*(Isc + Im)
    
   
    
    beta0 = [Rs0,GradRsh0,n];
    
     options = optimset('Display','none',...
    'TolX',1e-10,...
    'TolFun',1e-10,... 
    'Algorithm','trust-region-reflective');
    
    if Voc_index > Isc_index
    
    span = Isc_index:1:Voc_index;
   
    else 
        
    span = Voc_index:1:Isc_index;
    
    end 
    %[0.0001,200,1],[0.7*(Voc/Isc),10^6,7]
    
    nmax = (GradRs0*Isc)/(Vt);
    
    if nmax < 7
        
    nmax = ceil(nmax) + 0.1;
     
    else 
        nmax  = 7;
    end
        
    Rsh_22 = (((Vt)/(-Vm/Im) - (Isc + Im))/(-Vm + Vt))^-1;
    
    lowb = [1e-4,Rsh_22,1];
    upb = [GradRs0*1.3,Inf,nmax];
    b=lsqnonlin(@errorfuntion,real(beta0),lowb,upb,options,V(span),I_smooth(span),Voc,Isc);
    I_fit = fitter(b,V,Voc,Isc);
    [GradRs0_fit,GradRsh0_fit,Voc_fit,Isc_fit,Im_fit,Vm_fit,Voc_index_fit,Isc_index_fit,I_smooth_fit,I_orignal_data_fit] = lineofbestfit(V,I_fit);
    
    
    %VERY IMPORTANT ALL THE I VALUES MUST BE NEGATIVE AT ISC FOR THIS TO
    %WORK
    poloo = subplot(2,2,fileiter);
     fillfactor = sprintf('Fill Factor = %0.1f%%',ffactor);
     %subplot(15,2,fileiter);
 
     name = sprintf('Case %d',fileiter);
     title(name);
     hold on
     plot(V,-I_orignal_data,'r-','LineWidth',1);
     %plot(V,-I_fit,'b-','LineWidth',1);
     %plot(Vm_fit,-Im_fit,'bo','MarkerFaceColor','blue')   
     plot(Vm,-Im,'ko','MarkerFaceColor','black');
     
    
  %  text(0.10*Voc,Isc/2,fillfactor);
     %plot(Vm_fit,-Im_fit,'go')
     
     grid on 
     grid minor
     xlabel('Voltage (V)')
     ylabel('Current (A)')
     
     ylim([0 Isc*1.2])
     xlim([0 Voc*1.1])
     
     RMSE = sqrt(mean((I_smooth(span) - I_fit(span)).^2));  % Root Mean Squared Error
     
     extractpara = [b(1),b(2),b(3),RMSE];
     
     struArray{fileiter}.extract = extractpara;

     
     if fileiter == 1
         legend('Simulated data with SR = 45','MPP','location','southwest') 
         
        %legend('Experimental Curve','Fitted Curve','Fitted MPP','Experimental MPP','location','northeast') 
         
     end
      
end

hold off
%legend('Simulated Data','Model Data','MMP Simulated Data','location','best')
for i = 1:1:fileiter
   sprintf('case %i',i) 
   Rs =  struArray{i}.extract(1)
   Rsh = struArray{i}.extract(2)
   n = struArray{i}.extract(3)
   RMSE = struArray{i}.extract(4)

    
end

cd(top);



function [min] = errorfuntion(beta0,V,I,Voc,Isc)

Rs = beta0(1);
Rsh = beta0(2);
N = beta0(3)*300;

q = 1.6012*10^(-19);

kb = 1.38*10^(-23);

Ireg = V/Rs - Rsh*(((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/(1-exp(q*(Rs*Isc - Voc)/(kb*N))) + ...
    (Rs*Voc)/Rsh + V)/(Rs*(Rs+Rsh))) + (kb*N)/(q*Rs)...
    *lambertw(((q*Rs)/(kb*N))*(Isc - Voc/(Rs+Rsh))*exp(-(q*Voc)/(kb*N))...
    /(1- exp(q*(Rs*Isc - Voc)/(kb*N)))*...
    exp((Rsh*q*((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))...
    /(1-exp(q*(Rs*Isc - Voc)/(kb*N)))+(Rs*Voc)/Rsh + V))/(N*kb*(Rs+Rsh))));


min = Ireg - I;
end 


function Ireg = fitter(beta0,V,Voc,Isc)


Rs = beta0(1);
Rsh = beta0(2);
N = beta0(3)*300;

q = 1.6012E-19;

kb = 1.38*10^(-23);

Ireg = V/Rs - Rsh*(((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))/(1-exp(q*(Rs*Isc - Voc)/(kb*N))) + ...
    (Rs*Voc)/Rsh + V)/(Rs*(Rs+Rsh))) + (kb*N)/(q*Rs)...
    *lambertw(((q*Rs)/(kb*N))*(Isc - Voc/(Rs+Rsh))*exp(-(q*Voc)/(kb*N))...
    /(1- exp(q*(Rs*Isc - Voc)/(kb*N)))*...
    exp((Rsh*q*((Rs*(Isc + ((Rs*Isc - Voc)/Rsh)))...
    /(1-exp(q*(Rs*Isc - Voc)/(kb*N)))+(Rs*Voc)/Rsh + V))/(N*kb*(Rs+Rsh))));

end

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

