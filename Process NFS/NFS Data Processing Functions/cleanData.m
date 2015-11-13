%%Script to read a csv file and remove lines containing non-data elements
%%Kevin Delgado

function [f_vec]=cleanData(file_in, file_out)

    f_vec=[];
    %if no output file is defined, the output file will be the input file
    if nargin < 2
       file_out = file_in;
    end
    %temporarily converts file to .txt
    [pathstr] = fileparts(file_in);
    if exist([pathstr,'\TEMP'],'dir')~=7
        mkdir([pathstr,'\TEMP'])
    end
    temp_file=[pathstr,'\TEMP\in.txt'];
    copyfile(file_in,temp_file); 

    %open the in file to be read and the out file to be written to 
    fid = fopen(temp_file, 'r+');
    fout = fopen(file_out, 'w');
    
    %check that the file has been opened
    if fid == -1
        error('Could not open input file');
    end
    
    %Until reaching the end of the file check if the current line is data
    %and if so, write it to the out file
%     curLine = fgetl(fid);
    while ~feof(fid)
        curLine = fgetl(fid);
        if isempty(strfind(curLine, '['))
            fprintf(fout, '%s\n', curLine);
        elseif ~isempty(findstr(curLine,'[Frequency]')) && ...
               ~isempty(findstr(curLine,'[Pol],1'))
            
            s=curLine(findstr(curLine,'[Frequency]')+12:end);
            f=sscanf(s,'%f');
            f_vec=cat(2,f_vec,f);
        end
    end
    
    
    fclose(fid);
    fclose(fout);
    
    return 
end