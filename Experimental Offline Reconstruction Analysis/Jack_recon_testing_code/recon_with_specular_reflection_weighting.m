clear
%close all
home

reg = 1E-3;
tol = 1E-4; %tolerance for stop criteria 1

specularity = -3;
upsample_rate = 3;
cmap = 'gray';

%filename = '2PaneTx_hornRx_Cross2';
%filename = 'SmileyFace_2';
%filename = 'Medium_SmileyFace';
%filename = 'Medium_SmileyFace_201pts';
filename = '\Scenes_with_kinect\Smiley_IFtest';

load(['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'.mat'])
H = scene_data.H;
g = scene_data.obj_saved(2).measurement;
F = scene_data.F;

Z = scene_data.Z;
Az = scene_data.Az;
El = scene_data.El;

%%secularity
Spec = specular_reflection_matrix(F,Az,El,Z,specularity);

H = H.*Spec;

%% TV minimization
S = svd(H);  
Hp = H./max(S);
gp = g./max(S);
lam = min(S)/max(S);
Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));

[obj obj_debias objfunc] = TwIST(gp,Hp,reg,'lambda',lam,'Initialization',0,'Phi',Phi,'Psi',Psi,'MaxiterA',2000,'StopCriterion',1,'ToleranceA',tol,'Verbose',0);

figure(3)
semilogy(abs(objfunc(2:end)-objfunc(1:(end-1))))
drawnow

obj = abs(obj);

% for zn=1:length(Z)
%     obj3D(:,:,zn) = reshape(obj((1:32^2)+(zn-1)*32^2),32,32);
% end

nz = length(Z);
nel = size(El,1);
naz = size(Az,2);
nae = nel*naz;
obj3D = zeros(nel,naz,nz);
for zn=1:nz
    objd = obj((1:nae)+nae*(zn-1));
    obj3D(:,:,zn) = reshape(objd,nel,naz);
end

%% plot
figure
%voxelimage_projections(robj(:,:,:,3),Az,El,Z)
plot_range_slices(obj3D,Az,El,Z,upsample_rate,1:length(Z),'y')
colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])

figure
[dum ind] = max(squeeze(sum(sum(obj3D,1),2)));
plot_range_slices(obj3D,Az,El,Z,upsample_rate,ind,'y')
axis equal;axis tight;colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])

figure
plot_recon_surface(obj3D,Az,El,Z,0,upsample_rate);


