clear all; close all; clc

%Reading the data file and generating I and V vector
fileID = fopen('sim2d.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1)';
I = A(:,2)';

%Centre point of data to be fitted around point closest to Isc
Isc_index = find(abs(V)==min(abs(V-0)));
a = Isc_index;

%Using a for loop for different ranges of data. Point closest to Isc is the centrepoint and
%the range is the no of points left and right to be included in the
%fitting.
Rsho2 = [];
Isc2 = [];
range1 = [2 3 4 5 6]; %Range of data to be fitted
for range = range1; 
    
    % Extracting the data to be fitted
        x=V(a); %x = V
        y=I(a); %y = I
%         %Extracting data left of centre point
%         for i=1:range
%             x = [x V(a-i)];
%             y = [y I(a-i)];
%         end
            % As data to the left of Isc is added to the end of the vector (i.e
            % the right), need to flip to corret this
            x=fliplr(x); 
            y=fliplr(y);
        % Extracting the data right of centre point
        for i=1:range
            x= [x V(a+i)];
            y = [y I(a+i)];
        end

    %Fitting polynomial to x  & y (or V & I) using inbuilt matlab function
    %polyfit(x, y, n) where x is indep variable, y is dependent variable
    %and n is the nth degree polynomial to be used. A while loop is used to
    %determine the n value which meets the error tolerance using least
    %squares error.
    n=1; %initialising. n is the nth degree polynomial
    error = 100; %initialising
    while (error>0.1) % can change tolerance
        p=polyfit(x,y,n);
        y2=polyval(p,x); % y2 returns the polyfit values for x
        error = sum((y2-y).^2); % using least squares error 
        n=n+1;
    end

    %Extract Isc & Rsho. Isc is when x=0 (or V=0). polyder evaluates the
    %derivate of a function given its coefficients p. 
    Isc = polyval(p,0);
    q=polyder(p); %coefficients of derivative function
    Rsho=-1/polyval(q,0); %evaluating derivative at Isc

    % Rsho array for different ranges
    Rsho2 = [Rsho2 Rsho];
    Isc2 = [Isc2 Isc];
end
%Tabulating range, Rsho, Isc
Rsho = Rsho2';
range = range1';
Isc = Isc2';
table (range, Rsho, Isc)
mean = mean(Rsho);
fprintf('The average Rsho2 value is %.2f',mean)
