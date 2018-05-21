clc,clear all,close all
 
%open the files in the datareader.ma which returns a structure
[struArray,top] = datagrab();



for fileiter = [1:1:length(struArray) ]

    A = struArray{fileiter}.data;
    name = struArray{fileiter}.name;
    V = A(:,1)';
    I = awgn(A(:,2),45,'measured')';
    q = 1.6012*10^(-19); 
    kb = 1.38*10^(-23);
    T = 300;
    Vt = kb*T/q;
   
    name = sprintf('%s_SR_45.txt',name);
    A = [V;I];
       
    fileID = fopen(name,'w');
    fprintf(fileID,'%12f %12f\r\n',A);
    fclose(fileID);

      
end


function [struArray,top] = datagrab()
%top is the current directory
top = pwd;
%go to Funcfiles
cd Funcfiles
%get the data
struArray = datareader(top);
%{return to back to Funcfiles since datareader will direct you to Datafiles
cd(top);
cd Funcfiles

end

