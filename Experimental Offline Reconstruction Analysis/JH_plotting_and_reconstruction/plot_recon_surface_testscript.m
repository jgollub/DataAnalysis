%function plot_recon_surface(obj3D,Az,El,Z,threshold,upsample)
%plot a 3D scene in multipule range slices
%Inputs:
%obj3D - a 3D matrix where the first dimension is elevation, the second is azimuth, and the third is range
%Az,El - azimuth and elevation, 2D plaid matricies of the type produced by meshdrid()
%Z - 1D array of ranges
%threshold - [0,1], only display parts of the image greater than this fraction of the max value
%upsample - image upsampling factor
%range_indexes - 1D array of range indicies to plot

close all
obj3D = scene_data.obj_saved(4).reconstructed;
Az = scene_data.Az;
El = scene_data.El;
Z = scene_data.Z;
upsample = 3;
threshold = 0.0;

obj3D = flipdim(abs(obj3D),1);
%create an upsampled 3D image
objup = upsample_3d(obj3D,upsample);
objsum = sum(objup,3)./max(max(sum(objup,3)));

Az = upsample_image(Az,upsample);
El = upsample_image(El,upsample);
Z = linspace(min(Z),max(Z),round(upsample*length(Z)));

nAz = size(Az,2);
nEl = size(El,1);

%find max values along z of objup
[objmax ind] = max(objup,[],3);
objmax = objmax./max(max(objmax));
ind = repmat(ind,[1 1 length(Z)]);

depth = [];
for zn=1:length(Z)
depth = cat(3,depth,zn*ones(nEl,nAz));
end


Z = repmat(permute(Z,[1 3 2]),[nAz, nEl,1]);
Z(~(ind==depth)) = 0;
Z = sum(Z,3);

Z(objsum<=threshold) = NaN;

X = tan(Az).*Z;
Y = tan(El).*Z./cos(Az);


figure('Color',[0 0 0])
axes('Color',[0 0 0],'XColor',[1 1 1], 'YColor',[1 1 1], 'ZColor',[1 1 1])
%surface(X,Z,Y,objmax,'edgecolor','none')
surface(X,Z,Y,objmax,'AlphaData',objmax.^0.5,'FaceAlpha','interp','edgecolor','none')
%surface(X,Z,Y,objsum,'edgecolor','none')
%surface(X,Z,Y,objsum,'AlphaData',objsum,'FaceAlpha','interp','edgecolor','none')
axis equal
axis tight
colormap gray
