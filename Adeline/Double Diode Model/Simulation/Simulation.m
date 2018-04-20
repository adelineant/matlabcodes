clear all; close all; clc
% V=0:0.01:22;
% n1=1.26;
% n2=2.84;
% Rs=0.36;
% Rsh=233.46;
% I02=2.75E-8;
% I01=7.11E-8;
% Iph = 3.45;
% Ns=36;
% T=333;
% I=[];
% for i=1:length(V)
%     I=[I fzero(@(I)dd(I,V(i),n1, n2, Rs, Rsh, I02, I01, Iph, Ns, T),5)];
% end
% plot (V,I, 'b')
% hold on

fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1)';
I2 = -A(:,2)';

n1=1
n2=2
Rs=36
Rsh=108
I02=6.3E-19
I01=2.4E-19
Iph = 0.0024
Ns=1
T=303
I=[];
for i=1:length(V)
    I=[I fzero(@(I)dd(I,V(i),n1, n2, Rs, Rsh, I02, I01, Iph, Ns, T),0)];
end
plot (V,I, 'b')
hold on







