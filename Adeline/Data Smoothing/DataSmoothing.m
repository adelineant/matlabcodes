clear all; close all; clc

fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1)';
I = -A(:,2)';

%Use point closest to Voc as Voc (Need to improve)
Voc_index = find(abs(I)==min(abs(I-0)));
Voc = V(Voc_index);

%Use point closest to Isc as Isc (Need to improve)
Isc_index = find(abs(V)==min(abs(V-0)));
Isc = I(Isc_index);

B = smooth(I(1:1:Voc_index));
plot(V(1:1:Voc_index),I(1:1:Voc_index),'r*',V(1:1:Voc_index),B(1:1:Voc_index),'b')
legend('Original Data','Smoothed Data')