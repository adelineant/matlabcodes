clear all; close all; clc

V=0:0.01:4;
%Constants
k=1.38E-23;
T= 298.15;
q = 1.6012E-19;

%VARYING N
Is= 7.72e-07;
Rs= 1 ;
n= 4;
Rsh= 3.95e+03; 
Iph=3.70e-03 ;
for n=1:0.2:2
%Plotting extracted parameters
 I=[];
for i=1:length(V)
    I=[I fzero(@(I)f2(I,V(i),Is, Rs, n, Rsh, Iph),0)];
end
plot (V(1:i),I(1:i), 'b');
hold on
end
ylim([0 max(I)])
xlabel ('V')
ylabel ('I')
title ('Effect of varying ideality factor(n) on I-V curve')

