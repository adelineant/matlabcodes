function [struArray] = datareader(topdir)

%Stop the directory you came from
    
    cd(topdir);
%Go to data files
    
    cd Datafiles;
%look for files with txt ending 
    fileList = dir('*.txt');
%make an array to contain the filenames    
    FileNames = string(ones(1,length(fileList)));
%define a cell     
    struArray = cell(1,length(fileList));
%into this cell pass a structure with data and the file name    
    for k=1:length(fileList)
        
    FileNames(k) =fileList(k).name;
    
    fileID = fopen(FileNames(k),'r');
    A = [fscanf(fileID,'%f',[2 Inf])]';
    fclose(fileID);
    s = struct('data',A,'name',FileNames(k));
    
    struArray{k} = s;
       
    end

end