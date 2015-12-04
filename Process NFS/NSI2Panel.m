%function takes in frequency processing points and any number of phase
%corrections to be made. Phase input will be removed (divided)

function [] = NSI2Panel(savePath,scanPath,f,varargin)
phase_correction=varargin;

%load file 
NSI=load(scanPath);

%Check that frequency pts matchup
[f_check, f_location]=ismember(f,NSI.f);
if ~all(f_check)
    error('Input frequencies do not match up with NSI data frequencies')
end

%Choose desired frequency pts
measurements=NSI.measurements(:,:,f_location,:);

% PHASE CORRECTIONS
%correct for each component (probe, connector, etc)
for el=1:numel(phase_correction)
    for ii=1:numel(f)
        measurements(:,:,ii,:)=measurements(:,:,ii,:)/(phase_correction{el}(ii));
    end
end

X=NSI.X;
Y=NSI.Y;

save(savePath,'X','Y','f','measurements');