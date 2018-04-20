clear all; clc
fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1);
I = -A(:,2);

Voc_index = find(abs(I)==min(abs(I-0)));
b = Voc_index;

Rso2 = [];%Rs values for different ranges
for range = [10 20 30 40 50]; %CAN CHANGE
% Testing different no of data points taken left and right of centre
x=V(b); %Centre point of data to be fitted
y=I(b); %Centre point of data to be fitted
%Constructing data vector
for i=1:range
    x = [x V(b-i)];
    y = [y I(b-i)];
end
x=fliplr(x);
y=fliplr(y);
for i=1:range
    x= [x V(b+i)];
    y = [y I(b+i)];
end

%Fitting polynomial
n=1; %initialising
error = 100; %initialising
while (error>0.05)%Can change tolerance
    p=polyfit(x,y,n);
    y2=polyval(p,x);
    error = sum((y2-y).^2);
    n=n+1;
end
%plot (V,I)
%hold on
%plot (x,y2,'r')
%Extract Voc and Rso
Voc = roots(p);
%plot (Voc, 0, 'b*')
q=polyder(p);
Rso = -1/polyval(q,Voc);
Rso2 = [Rso2 Rso];
end
Rso2
meanRso = mean(Rso2)




