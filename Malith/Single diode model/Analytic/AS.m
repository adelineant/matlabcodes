 clear all; close all; clc


%open the files in the datareader. ma
[struArray,top] = datagrab();

%Opening data file
%fileID = fopen('H9-4-1FTO-1C-R-1PH-1X.txt','r');
%A = [fscanf(fileID,'%f',[2 Inf])]';

for fileiter = [1:1:length(struArray)]
    
    A = struArray{fileiter}.data;
    V = A(:,1);
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
    
end
cd(top);

function [struArray,top] = datagrab()

top = pwd;
cd Funcfiles
struArray = datareader(top);

cd(top);
cd Funcfiles

end
