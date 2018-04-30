clear all; close all; clc

V=0:0.01:4;
%Constants
k=1.38E-23;
T= 298.15;
q = 1.6012E-19;

% %VARYING N
% Is= 7.72e-07;
% Rs= 1 ;
% n= 4;
% Rsh= 3.95e+03; 
% Iph=3.70e-03 ;
% for n=1:0.2:2
% %Plotting extracted parameters
%  I=[];
% for i=1:length(V)
%     I=[I fzero(@(I)f2(I,V(i),Is, Rs, n, Rsh, Iph),0)];
% end
% plot (V(1:i),I(1:i), 'b');
% hold on
% end
% ylim([0 max(I)])
% xlabel ('V')
% ylabel ('I')
% title ('Effect of varying ideality factor(n) on I-V curve')

%VARYING Rs
Is= 7.72e-07;
n=1;
Rsh= 3.95e+03; 
Iph=3.70e-03 ;
for Rs=0:5:20
%Plotting extracted parameters
 I=[];
for i=1:length(V)
    I=[I fzero(@(I)f2(I,V(i),Is, Rs, n, Rsh, Iph),0)];
end
plot (V(1:i),I(1:i));
hold on
end
ylim([0 max(I)])
xlabel ('Voltage (V)')
ylabel ('Current (A)')
title ('Effect of varying series resistant(Rs(ohms)) on I-V curve')
legend ('Rs=0', 'Rs=5', 'Rs=10', 'Rs=15', 'Rs=20')

% %VARYING Rsh
% Is= 7.72e-07;
% n=1;
% Rs = 1;
% Iph=3.70e-03 ;
% for Rsh = [100 500 1000 5000 10000]
% %Plotting extracted parameters
%  I=[];
% for i=1:length(V)
%     I=[I fzero(@(I)f2(I,V(i),Is, Rs, n, Rsh, Iph),0)];
% end
% plot (V(1:i),I(1:i));
% hold on
% end
% ylim([0 max(I)])
% xlabel ('Voltage(V)')
% ylabel ('Current (A)')
% title ('Effect of varying shunt resistant(Rsh(ohms)) on I-V curve')
% legend ('Rsh=100', 'Rs=500', 'Rs=1000', 'Rs=5000', 'Rs=10000')
