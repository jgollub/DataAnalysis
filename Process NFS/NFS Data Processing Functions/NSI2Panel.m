%function takes in frequency processing points and any number of phase
%corrections to be made. Phase input will be removed (divided)

function [] = NSI2Panel(savePath,scanPath,f_NSI,f,polarizations,varargin)

%% Convert from NSI to Matlab
phase_correction=varargin;

% CSV format: [duration, horizontal(x), vertical(y), amplitude (dB), phase (degree)]

% Load data
csv=load(scanPath);

xNum=numel(unique(csv(:,2)));
yNum=numel(unique(csv(:,3)));
freq_num=size(csv,1)/(xNum*yNum*polarizations);
measurements=reshape(...
    10.^(csv(:,4)/20).*exp(1.0j*(pi/180).*csv(:,5)),...
    xNum,...
    yNum,...
    polarizations,...
    freq_num...
    );
X=reshape(...
    csv(:,2),...
    xNum,...
    yNum,...
    polarizations,...
    freq_num...
    );
X=X(:,:,1,1);
X=X(:,:,1,1)*1000; %!!!! fixed negative

Y=reshape(...
    csv(:,3),...
    xNum,...
    yNum,...
    polarizations,...
    freq_num...
    );
Y=Y(:,:,1,1)*1000;

%permute to put into matlab format[xval yval pol freq] -> [yval xval freq pol]
X=permute(X,[2 1]);
Y=permute(Y,[2 1]);
Y=flip(Y,1);  %!!!! flipped delete to go back to original
measurements=permute(measurements,[2 1 4 3]);
measurements=flip(measurements,1); %%!!!!! flipped because Y is flipped (delete to go back)

%check that frequency pts matchup


%% PHASE CORRECTIONS
[f_check, f_location]=ismember(f,f_NSI);
if ~all(f_check)
    error('Input frequencies do not match up with NSI data frequencies')
end


for el=1:numel(phase_correction)
    counter=0;
    for ii=f_location
        counter=counter+1;
        measurements_out(:,:,counter,:)=measurements(:,:,ii,:)/(phase_correction{el}(counter));
    end
end

measurements=measurements_out; %reset to known variable name

% PROBE CORRECTION
%     for ii=1:Freq_points
%     measurements(:,:,ii,:)=measurements(:,:,ii,:)/measured2(ii);
%     end
save(savePath,'X','Y','f','measurements');