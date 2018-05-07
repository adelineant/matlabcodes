 clear all; close all; clc


%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();


for fileiter = [1:1:length(struArray)]
    %Acess the structure and then store the data in A
    %save the column in I and V respectively
    A = struArray{fileiter}.data;
    V = A(:,1);
    %YOU NEED TO FIND A WAY A WAY TO MAKE THE PROGRAM KNOW WHAT SIGN
    I = -A(:,2);
    



    %Defining the parameters in analytical single model. Use polynomial


    %Applying Analytical Single Model
    k=1.38E-23;
    T= 306.15;
    q = 1.6012E-19;

    [Rs0,Rsh0,Voc,Isc,Im,Vm] = lineofbestfit(V,I);

    Vt = k*T/q;
    Rsh=Rsh0;
    n=(Vm+Im.*Rs0-Voc)./(Vt.*(log(Isc-Vm./Rsh-Im)-log(Isc-Voc./Rsh)+Im./(Isc-Voc./Rsh)));
    Is = (Isc-Voc/Rsh)*exp(-Voc/(n*Vt));
    Rs = Rs0-n*Vt/Is*exp(-Voc/(n*Vt));
    if Rs < 0
        
        Rs = 0;
        
    end
    Iph = Isc*(1+Rs/Rsh)+Is*(exp(Isc*Rs/(n*Vt))-1);

    %Plotting using the extracted parameters
    I2= fzero(@(I)f2(I,V(1),Is, Rs, n, Rsh, Iph),Isc);
    for i=1:length(V)-1


            %I2=[I2 fzero(@(I)f2(I,V(i+1),Is, Rs, n, Rsh, Iph),I2(i)*2)]

            I2=[I2 fzero(@(I)f2(I,V(i+1),Is, Rs, n, Rsh, Iph),I2(i))];

            if(isnan(I2(i+1))== 1)
                break
            end
    end
    figure;
    
    plot (V,I, 'r');
    hold on
    plot (V(1:i),I2(1:i), 'b');

    legend ('experimental', 'fitted')
    xlabel ('V')
    ylabel ('I')
    hold on
    plot(Vm,Im,'go');
    ylim([0 Inf])
    %grad = (Vm*Im - )/()
end
cd(top);

function [struArray,top] = datagrab()

top = pwd;
cd Funcfiles
struArray = datareader(top);

cd(top);
cd Funcfiles

end
