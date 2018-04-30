function [struArray] = datareader(topdir)


    cd(topdir);
    cd Datafiles;
 
    fileList = dir('*.txt');
    FileNames = string(ones(1,length(fileList)));
    struArray = cell(1,length(fileList));
    for k=1:length(fileList)
        
    FileNames(k) =fileList(k).name;
    
    fileID = fopen(FileNames(k),'r');
    A = [fscanf(fileID,'%f',[2 Inf])]';
    fclose(fileID);
    s = struct('data',A,'name',FileNames(k));
    
    struArray{k} = s;
       
    end

end