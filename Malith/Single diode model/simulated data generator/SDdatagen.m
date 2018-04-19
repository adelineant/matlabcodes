clc,clear all, close all

for Rs = linspace(1e-8,1e-3,5)
    for Rsh = linspace(8,1,5)
        for n = linspace(1,3,6)
            
            V = 0:0.01:1;
            Iph =0.7608;
            I0 =0.3223*10^-6;
            T = 306.15;

            %y=SD_equation(V,I,Rs,Rsh,n)
            Ical = zeros(size(V));

            Ical(1)= fzero(@(I)SD_equation(V(1),I,Rs,Rsh,n,Iph,I0,T),-0.3);

            for i=1:length(V)-1
                Ical(i+1) = fzero(@(I)SD_equation(V(i+1),I,Rs,Rsh,n,Iph,I0,T),1);
            end
            
            %plot (V,-Ical)
           % hold on
            
            x = 0:.1:1;
            A = [V;-Ical];
            format = sprintf("Rs=%e_Rsh=%e_n=%f.txt",Rs,Rsh,n);
            

            fileID = fopen(format,'w');
            fprintf(fileID,'%6s %12s\r\n','V','Ical');
            fprintf(fileID,'%6.2f %12.8f\r\n',A);
            fclose(fileID);
        end
    end
end


function y=SD_equation(V,I,Rs,Rsh,n,Iph,I0,T)

q = 1.6012*10^(-19);
kb = 1.38*10^(-23);


y = I0.*(exp(q.*(V-Rs.*I)./(n.*kb.*T))-1)+(V-Rs.*I)./Rsh-Iph-I;
end