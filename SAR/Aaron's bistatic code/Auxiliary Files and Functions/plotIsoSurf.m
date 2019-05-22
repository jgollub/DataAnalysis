function [ fax ] = plotIsoSurf(plotting)
% fax = plotIsoSurf(recon,levels,alphas,colors)
%   recon should include x, y, and z vectors. It should also include
%   image_lin which contains the linear format of the 3D data set. The
%   number of levels/alphas/colors must be the same. Colors should be in
%   string format or 3vectors. The output is a useless figure handle.

[X,Y,Z] = meshgrid(plotting.x,plotting.y,plotting.z);
Image=20*log10(plotting.image_lin);

fax=figure;
for ii=1:length(plotting.levels)-1
    p = patch(isosurface(X,Y,Z,Image,plotting.levels(ii)));
    set(p,'FaceColor',plotting.colors{ii},'EdgeColor','none','FaceAlpha',plotting.alphas(ii));
end

p = patch(isosurface(X,Y,Z,Image,plotting.levels(end)));
set(p,'FaceColor',plotting.colors{end},'EdgeColor','none',...
    'FaceAlpha',plotting.alphas(end),'EdgeAlpha',0.1);

camlight(37,11)
lighting gouraud
box on, grid on
axis image

cellLeg=cell(numel(plotting.levels),1);
for ii=1:length(plotting.levels)-1; cellLeg{ii}=strcat(num2str(plotting.levels(ii)),' dB'); end
legend(cellLeg{1:3})

end

