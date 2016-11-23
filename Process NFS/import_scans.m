% Copyright 2015 Duke and Evolv Technologies, Inc.  Development of this software was
% supported in part by the U.S. Government under contract number HSHQDC-12-C-00049
%==================================================================================
% IMPORT NEAR FIELD SCAN FUNCTION
% data inputs: [data, polarization, measurement index,... ]
%
% Using the Near Field Scan (NFS) coordinate system (X,Y) assign the polarization to the data, 
% e.g. if your first measurement is for polarization along the Y-axis
% then,
%   panel=import(data,'Ey', 1) 
% or if you take two polarizations with the first being along x, then:
%   panel=import(data,'Ex', 1,'Ey',2) 
%
% Generally, this function removes the abiguity of the previous import panel function
% around how the magnetic dipoles are derived from the NFS.
% The Ex polarization results in magnetic dipoles oriented along the
% Z-axis (once translated into the imaging coordinate system (x=>Y,y=>Z)) and the
% Ey goes to magnetic dipoles oriented along the y-axis
%
%=================================================================================

function panel = import_scans(data,varargin)
%parse input
p=inputParser;
validationFcn=@(x) validateattributes(x,{'numeric'},{'integer','<=',2,'>',0});
addParameter(p,'Ex',NaN,validationFcn);
addParameter(p,'Ey',NaN,validationFcn);
addParameter(p,'Renormalize', 0);
parse(p,varargin{:})

% allow direct path input (but direct ok too)
if ischar(data)
    data = load(data);
end

%assign data to correct polarization    
if (p.Results.Ex==p.Results.Ey | (isnan(p.Results.Ey) & isnan(p.Results.Ex)));
    error('you have set Ey and Ez to the same value or neither have been assigned')
end

if ~isnan(p.Results.Ex)
    data_Ex=data;
    data_Ex.measurements=data_Ex.measurements(:,:,:,p.Results.Ex);
    panel.dipoles.z = [];
end

if ~isnan(p.Results.Ey)
    data_Ey=data;
    data_Ey.measurements=data_Ey.measurements(:,:,:,p.Results.Ey);
    panel.dipoles.y = [];
end

%units:
GHz=1e9;
mm=1e-3;
cm=1e-2;
in=2.54*cm;
u0 = 4*pi*1e-7;

%feed locations based on CAD file
 feed_positions      = [0, 0, 0];                  
 panel.feedLocs(:,1) = feed_positions(:,1); %set feedLocs by CAD positions
 panel.feedLocs(:,2) = feed_positions(:,2);
 panel.feedLocs(:,3) = feed_positions(:,3);
 panel.numfeeds = size(feed_positions,1); 
 
panel.y = data.X*mm; 
panel.z = data.Y*mm;
panel.f = data.f;

%default panel settings
panel.type='panel';

%default orientation
panel.u=[0; 1; 0]; 
panel.v=[0; 0; 1];

% data derivatives
panel.Size_y           = max(panel.y(:)) - min(panel.y(:));
panel.Size_z           = max(panel.z(:)) - min(panel.z(:));
panel.ElementSize_y    = abs(panel.y(2,2) - panel.y(1,1));
panel.ElementSize_z    = abs(panel.z(2,2) - panel.z(1,1));
A0 = panel.ElementSize_y * panel.ElementSize_z;
panel.fstart = panel.f(1);
panel.fstop  = panel.f(end);
panel.dims = size(panel.y);
panel.x = zeros(panel.dims);

%calculate magnetic dipole equivalents

    tmp_y = zeros([numel(panel.f), size(panel.y)]);
    tmp_z = zeros([numel(panel.f), size(panel.y)]);
    panel.dipoles.x=tmp_y;
    panel.dipoles.y=tmp_y;
    panel.dipoles.z=tmp_y;

    for fi = 1:length(panel.f)
        if ~isnan(p.Results.Ey)
            tmp_y(fi,:,:) = -2i*A0/(u0*2*pi*panel.f(fi))*data_Ey.measurements(:,:,fi);
        else
            tmp_y(fi,:,:) =0;
        end
        
        if ~isnan(p.Results.Ex)
            tmp_z(fi,:,:) = 2i*A0/(u0*2*pi*panel.f(fi))*data_Ex.measurements(:,:,fi); %FIXED NEG SIGN)
        else
            tmp_z(fi,:,:)=0;
        end
    end
    
    if ~isnan(p.Results.Ey)
        panel.dipoles.y = tmp_y; 
    end
    if  ~isnan(p.Results.Ex)
        panel.dipoles.z = tmp_z;
    end

%     panel.dipoles.x = zeros(size(panel.dipoles.y));

    
if p.Results.Renormalize
%     kappa = repmat(sqrt((power_calc(panel, 'self'))), [1 size(panel.x)]);
    kappa = sqrt(mean(power_calc(panel, 'self'))); %was 'nearest'
    panel.dipoles.x = panel.dipoles.x ./ kappa;
    panel.dipoles.y = panel.dipoles.y ./ kappa;
    panel.dipoles.z = panel.dipoles.z ./ kappa;
end

