%%Script to read a csv file and remove lines containing non-data elements
%%Kevin Delgado & Jonah Gollub

function cleanData(file_in, file_out)

    data=[];
    f=[];
     %if no output file is defined, the output file will be the input file
    if nargin < 2
       file_out = [file_in(1:end-4),'.mat'];
    end
    
    %open the in file to be read and the out file to be written to
    fid = fopen(file_in, 'r+');    
    %check that the file has been opened
    if fid == -1
        error('Could not open input file\n');
    end
     fprintf('PROCESSING: ')
    %get parameters from preamble of file
for interation=1:10
    curLine = fgetl(fid);        
    if ~isempty(findstr(curLine,'[Beam count],'))
        s=curLine(findstr(curLine,'[Beam count],')+13:end);
        fNum=sscanf(s,'%f');
    elseif ~isempty(findstr(curLine,'[H Axis]'))
        
        index=findstr(curLine,'[Name],')+7;
        s=curLine(index:index+1);
        Haxis.Axis=sscanf(s,'%s');
        
        index=findstr(curLine,'[Unit],')+7;
        s=curLine(index:index+5);
        Haxis.Units=sscanf(s,'%s');
        
        index=findstr(curLine,'[Points],')+9;
        s=curLine(index:end);
        Haxis.Points=sscanf(s,'%u');
        
    elseif ~isempty(findstr(curLine,'[V Axis]'))
        
        index=findstr(curLine,'[Name],')+7;
        s=curLine(index:index+1);
        Vaxis.Axis=sscanf(s,'%s');
        
        index=findstr(curLine,'[Unit],')+7;
        s=curLine(index:index+5);
        Vaxis.Units=sscanf(s,'%s');
        
        index=findstr(curLine,'[Points],')+9;
        s=curLine(index:end);
        Vaxis.Points=sscanf(s,'%u');
        
    elseif ~isempty(findstr(curLine,'[Dual Pol],'))
        
        index=findstr(curLine,'[Dual Pol],')+11;
        s=curLine(index:end);
        logic=sscanf(s,'%s');
        if strcmp(logic,'Yes')
            polNum=2;
        else
            polNum=1;
        end
    end
end

if isempty(fNum) | isempty(Vaxis.Points) | isempty(Haxis.Points) | isempty(polNum)
    error('missing scan parameters in preamble');
end

%scan in data for each frequency
bufferSize=Vaxis.Points*Haxis.Points;

    %Until reaching the end of the file check if the current line is data
    %and if so, write it to the out file

while ~feof(fid)
    curLine=fgetl(fid);  %get first line of block      
    if ~isempty(findstr(curLine,'[Frequency]')) && ...
            ~isempty(findstr(curLine,'[Pol],1'))
        
        s=curLine(findstr(curLine,'[Frequency]')+12:end);
        fVal=sscanf(s,'%f');
        f=cat(1,f,fVal);
    end
   curLine=fgetl(fid); %get and ingnore second line of block 
   
   ScannedData = cell2mat(textscan(fid,'%f%f%f%f%f ',bufferSize,'Delimiter',','));
          data = cat(1,data,ScannedData(:,2:5));
  
          curLine=fgetl(fid); %get and ingnore second line of block
          
end
    
    fclose(fid);
    
    xNum=numel(unique(data(:,1)));
    yNum=numel(unique(data(:,2)));
    if   fNum~=numel(f) | Haxis.Points~=xNum | Vaxis.Points~=yNum
        error('number of scan positiosn or frequencies in scan setup and saved data are not equal')
    end
    
    fprintf('CSV DATA SORTED...')
%% convert to matrix form for matlab use

measurements=reshape(...
    10.^(data(:,3)/20).*exp(1.0j*(pi/180).*data(:,4)),...
    xNum,...
    yNum,...
    polNum,...
    fNum...
    );
X=reshape(...
    data(:,1),...
    xNum,...
    yNum,...
    polNum,...
    fNum...
    );
X=X(:,:,1,1);
X=X(:,:,1,1)*1000; 

Y=reshape(...
    data(:,2),...
    xNum,...
    yNum,...
    polNum,...
    fNum...
    );
Y=Y(:,:,1,1)*1000;

%permute to put into matlab format[xval yval pol freq] -> [yval xval freq pol]
X=permute(X,[2 1]);
Y=permute(Y,[2 1]);
Y=flip(Y,1);  
measurements=permute(measurements,[2 1 4 3]);
measurements=flip(measurements,1); %%!!!!! flipped because Y is flipped (delete to go back)   

fprintf('SAVING MATLAB DATA STRUCTURE...')
%save data
save(file_out,'measurements','f','fNum','X','Y','Vaxis','Haxis','polNum')
fprintf(['FINISHED: \n']);
    return 
end