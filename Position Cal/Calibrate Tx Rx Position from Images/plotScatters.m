function [f_plot] = volume_plot3(f_raw, indices,imgDomain, domainGrid, varargin)
parser=inputParser;

defaultType='linear';
defaultThreshold=15;
addParameter(parser, 'Type', defaultType);
addParameter(parser, 'Threshold', defaultThreshold,@isnumeric);
% addParameter(parser, 'Indices',imgDomain)
parse(parser,varargin{:});
    
    f = zeros(size(imgDomain, 1), 1);
    f(indices) = abs(f_raw).^2/max(max(max(abs(f_raw).^2)));
    f_plot=reshape(f, domainGrid);
    
if strcmp(parser.Results.Type,'log')
    
    f_db       = 20*log10(abs(f_plot));
    thresh     = max(max(max(f_db))) -parser.Results.Threshold;
    %         f          = ones(size(f_plot))*thresh;
    %         f_plot = 20*log10(abs(f_plot));
    f_db(f_db  <thresh)= thresh;
    f_plot=f_db;
end

hf = vol3d('cdata',  f_plot,...
    'XData', [min(imgDomain(:,1)), max(imgDomain(:,1))],...
    'YData', [min(imgDomain(:,2)), max(imgDomain(:,2))],...
    'ZData', [min(imgDomain(:,3)), max(imgDomain(:,3))]...
    );

    xlim([min(imgDomain(indices,1)),max(imgDomain(indices,1))]);
    ylim([min(imgDomain(indices,2)),max(imgDomain(indices,2))]);
    zlim([min(imgDomain(indices,3)),max(imgDomain(indices,3))]);

view(3);
 axis equal;

end

