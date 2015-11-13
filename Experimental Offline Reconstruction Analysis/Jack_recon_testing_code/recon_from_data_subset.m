clear
close all
home

reg = 6E-2;
tol = 1E-5; %tolerance for stop criteria 1

upsample_rate = 3;
cmap = 'gray';

filename = '2PaneTx_hornRx_Cross2';
%filename = 'Medium_SmileyFace';
load(['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'.mat'])
H = scene_data.H;
g = scene_data.obj_saved(1).measurement;

Z = scene_data.Z;
Az = scene_data.Az;
El = scene_data.El;

%% TV minimization
S = svd(H);  
Hp = H./max(S);
gp = g./max(S);
lam = min(S)/max(S);
Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));

[obj obj_debias objfunc] = TwIST(gp,Hp,reg,'lambda',lam,'Initialization',0,'Phi',Phi,'Psi',Psi,'MaxiterA',2000,'StopCriterion',1,'ToleranceA',tol,'Verbose',0);

figure(2)
semilogy(objfunc)
drawnow

obj = abs(obj);

for zn=1:length(Z)
    obj3D(:,:,zn) = reshape(obj((1:32^2)+(zn-1)*32^2),32,32);
end

%% plot
%this reconstructed image
h = figure(1);
%voxelimage_projections(robj(:,:,:,3),Az,El,Z)
plot_range_slices(obj3D,Az,El,Z,upsample_rate)
colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 15]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\all_ranges\','regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)

h = figure(1);
imagesc(Az(1,:)*180/pi,El(:,1)*180/pi,obj3D(:,:,5));
axis equal;axis tight;colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 12]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\','regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)

%% plot upsampled
%this reconstructed image
h = figure(1);
%voxelimage_projections(robj(:,:,:,3),Az,El,Z)
plot_range_slices(obj3D,Az,El,Z,upsample_rate)
colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 15]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\all_ranges\upsampled\','regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)

h = figure(1);
imagesc(Az(1,:)*180/pi,El(:,1)*180/pi,upsample_image(obj3D(:,:,5),upsample_rate));
axis equal;axis tight;colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 12]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\upsampled\','regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)

%% ***************************************************************************
%panel B
Ha = H((end/2+1):end,:);
ga = g((end/2+1):end);
%% TV minimization
S = svd(Ha);  
Hp = Ha./max(S);
gp = ga./max(S);
lam = min(S)/max(S);
Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));

[obj obj_debias objfunc] = TwIST(gp,Hp,reg,'lambda',lam,'Initialization',0,'Phi',Phi,'Psi',Psi,'MaxiterA',2000,'StopCriterion',1,'ToleranceA',tol,'Verbose',0);

figure(2)
semilogy(objfunc)
drawnow

obj = abs(obj);

for zn=1:length(Z)
    obj3D(:,:,zn) = reshape(obj((1:32^2)+(zn-1)*32^2),32,32);
end

%% plot
%this reconstructed image
h = figure(1);
%voxelimage_projections(robj(:,:,:,3),Az,El,Z)
plot_range_slices(obj3D,Az,El,Z,upsample_rate)
colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 15]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\all_ranges\','PanelB_regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)

h = figure(1);
imagesc(Az(1,:)*180/pi,El(:,1)*180/pi,obj3D(:,:,5));
axis equal;axis tight;colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 12]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\','PanelB_regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)

%% plot upsampled
%this reconstructed image
h = figure(1);
%voxelimage_projections(robj(:,:,:,3),Az,El,Z)
plot_range_slices(obj3D,Az,El,Z,upsample_rate)
colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 15]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\all_ranges\upsampled\','PanelB_regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)

h = figure(1);
imagesc(Az(1,:)*180/pi,El(:,1)*180/pi,upsample_image(obj3D(:,:,5),upsample_rate));
axis equal;axis tight;colormap(cmap)
title(['regularization: ' num2str(reg,'%1.2E')])
%set(h,'OuterPosition',[2270,500,365,500])
set(h, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 12]);
saveas(h,['C:\Users\John\Dropbox\Metaimager Project\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\',filename,'_figures\different_panels\upsampled\','PanelB_regularization_',num2str(reg,'%1.2E'),'.png'])
close(h)
