% Copyright 2013 Evolv Technologies, Inc.  Development of this software was
% supported in part by the U.S. Government under contract number HSHQDC-12-C-00049
function panel = import_panel_nsi(data_m, data_x)
% if only one dataset is handed, just assume the off-axis component is zero


if nargin == 1
    data_x.measurements = 0*data_m.measurements;
end
renormalize = 1;

%units:
GHz=1e9;
mm=1e-3;
cm=1e-2;
in=2.54*cm;
u0 = 4*pi*1e-7;

% data import
if ischar(data_m)
    data_m = load(data_m);
end
if ischar(data_x)
    data_x = load(data_x);
end
panel.y = data_m.y; %neg sign duke to coordinate change
panel.z = data_m.z;
panel.f = data_m.f;

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

panel.feedLocs(:,1) = 0;
panel.numfeeds = length(data_m.measurements(1,1,1,:));
panel.dipoles.y = [];
panel.dipoles.z = [];
for i =1:panel.numfeeds
    for fi = 1:length(panel.f)
        tmp_y(fi,:,:) = -2i*A0/(u0*2*pi*panel.f(fi))*data_m.measurements(:,:,fi,i);  
        tmp_z(fi,:,:) = -2i*A0/(u0*2*pi*panel.f(fi))*data_x.measurements(:,:,fi,i); %(neg sign due to importing coordinates sys)
    end
    [~, mi] = max(sum(abs(tmp_y(:,:)).^2));
    panel.feedLocs(i,2) = panel.y(mi); %approximate the feed locations from the data maxima
    panel.feedLocs(i,3) = panel.z(mi);
    
    panel.dipoles.y = cat(1, panel.dipoles.y, tmp_y);
    panel.dipoles.z = cat(1, panel.dipoles.z, tmp_z);
end

panel.dipoles.x = zeros(size(panel.dipoles.y));

% if renormalize
% %     kappa = repmat(sqrt((power_calc(panel, 'self'))), [1 size(panel.x)]);
%     kappa = sqrt(mean(power_calc(panel, 'nearest')));
%     panel.dipoles.x = panel.dipoles.x ./ kappa;
%     panel.dipoles.y = panel.dipoles.y ./ kappa;
%     panel.dipoles.z = panel.dipoles.z ./ kappa;
% end