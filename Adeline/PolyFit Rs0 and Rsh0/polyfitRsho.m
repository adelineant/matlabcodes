clear all; close all; clc
fileID = fopen('H4-1-1FTO-1C-F-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1);
I = -A(:,2);

Isc_index = find(abs(V)==min(abs(V-0)));
a = Isc_index;

Rsho2 = [];
for range = [30 40 50 55]; %CAN CHANGE THIS ACCORDING TO DATA USED
    x=V(a);%Centre point of data to be fitted
    y=I(a);%centre point of data to be fitted
    for i=1:range
        x = [x V(a-i)];
        y = [y I(a-i)];
    end
    x=fliplr(x);
    y=fliplr(y);
    for i=1:range
        x= [x V(a+i)];
        y = [y I(a+i)];
    end

    %Fitting polynomial
    n=1; %initialising
    error = 100; %initialising
    while (error>0.1) % can change tolerance
        p=polyfit(x,y,n);
        y2=polyval(p,x);
        error = sum((y2-y).^2);
        n=n+1;
    end

    %Extract Isc & Rsho
    Isc = polyval(p,0);
    q=polyder(p);
    Rsho=-1/polyval(q,0);
    Rsho2 = [Rsho Rsho2];
end
Rsho2
meanRsh = mean(Rsho2)
plot (x,y, 'r*');
hold on
plot (x,y2, 'b');