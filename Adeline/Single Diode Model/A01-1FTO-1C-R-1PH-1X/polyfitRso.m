clear all; close all; clc

%Reading the data file and generating I and V vector
fileID = fopen('A01-1FTO-1C-R-1PH-1X.txt','r');
A = [fscanf(fileID,'%f',[2 Inf])]';
fclose(fileID);
V = A(:,1)';
I = -A(:,2)';
V = flip(V);
I = smoothdata(flip(I));

%Centre point of data to be fitted around point closest to Voc
Voc_index = find(abs(I)==min(abs(I-0)));
b = Voc_index;

%Using a for loop for different ranges of data. Point closest to Voc is the centrepoint and
%the range is the no of points left and right to be included in the
%fitting.
Rso2 = [];
Voc2 = [];
range1 = [1]; %Range of data to be fitted
for range = range1; 
    
    % Extracting the data to be fitted
        x=V(b); %x = V
        y=I(b); %y = I
        %Extracting data left of centre point
        for i=1:range
            x = [x V(b-i)];
            y = [y I(b-i)];
        end
            % As data to the left of Isc is added to the end of the vector (i.e
            % the right), need to flip to corret this
            x=fliplr(x); 
            y=fliplr(y);
        % Extracting the data right of centre point
        for i=1:range
            x= [x V(b+i)];
            y = [y I(b+i)];
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

    %Extract Voc & Rso. Voc is the root of the polynomial. polyder evaluates the
    %derivate of a function given its coefficients p. 
    Voc = roots(p);
    q=polyder(p);
    Rso = -1/polyval(q,Voc);
    
    %Rsho and Voc array for different ranges
    Rso2 = [Rso2 Rso];
    Voc2 = [Voc2 Voc];
end
%Tabulating range, Rsho, Isc
Rso = Rso2';
range = range1';
Voc = Voc2';
table (range, Rso, Voc)
mean = mean(Rso);
fprintf('The average Rso value is %.2f',mean)




