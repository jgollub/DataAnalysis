%% Load data
A2_data=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\SemiRigid 101pts 2014-04\scans at 1m\A2_F1_1m_NFS-20-Jun-2014.mat');
A2_data=A2_data.data;
A2_data.X=-1*A2_data.X*1e-3+.298;
A2_data.Y=A2_data.Y*1e-3;
A3_data=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\SemiRigid 101pts 2014-04\scans at 1m\A3_F1_1m_NFS-22-Jun-2014.mat');
A3_data=A3_data.data;
A3_data.X=-1*A3_data.X*1e-3+.298; %.298 offset because NFS was aligned with panel A2
A3_data.Y=A3_data.Y*1e-3;
%Near Field Scan (NFS) data location for panels
file_path='C:\Users\Jonah Gollub\Documents\MetaImagerData\SemiRigid_NFS_scans\';

%Load dummy data to get panel positions/rotations
expData=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\MetaImager Scenes\Semi-Rigid&Phase_stable\imager_16-Jun-2014 11-11-08.mat');
panelPosition       = expData.panelPosition;%+repmat([.05 0 0],12,1);
panelRotation       = expData.panelRotation;

%% plot data A2
freq=51;
try figure(A2_NFSat1m_plot)
catch
    A2_NFSat1m_plot=figure();
    scrsize=get(0,'ScreenSize')
    figsize=A2_NFSat1m_plot.Position(3:4);
    A2_NFSat1m_plot.Position(1:2)=[scrsize(3)-figsize(1),scrsize(4)-figsize(2)-80];
end
subplot(2,2,1)

imagesc(A2_data.X(1,:),A2_data.Y(:,1),abs(A2_data.measurements(:,:,freq)))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Measured Fields (Z-Pol)')

subplot(2,2,3)

imagesc(A2_data.X(1,:),A2_data.Y(:,1),angle(A2_data.measurements(:,:,freq)))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Measured Phase (Z-Pol)')

try figure(A3_NFSat1m_plot)
catch
    A3_NFSat1m_plot=figure();
    scrsize=get(0,'ScreenSize')
    figsize=A3_NFSat1m_plot.Position(3:4);
    A3_NFSat1m_plot.Position(1:2)=[scrsize(3)-figsize(1),scrsize(4)/2-figsize(2)-80];
end
subplot(2,2,1)

imagesc(A3_data.X(1,:),A3_data.Y(:,1),abs(A3_data.measurements(:,:,freq)))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Measured Fields (Z-Pol)')

subplot(2,2,3)
imagesc(A3_data.X(1,:),A3_data.Y(:,1),angle(A3_data.measurements(:,:,freq)))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Measured Phase (Z-Pol)')


%% import panel

panel=zeros(12,1).';
panel([2 3])=1;
panelData=find(panel);

try figure(NFS_plot)
catch
    NFS_plot=figure();
    scrsize=get(0,'ScreenSize')
    figsize=NFS_plot.Position(3:4);
    NFS_plot.Position(1:2)=[scrsize(3)/2-figsize(1),scrsize(4)-figsize(2)-80];
end
cla
for i=panelData
    %     try isstruct(panel{i});
    %     catch
    [pan_index,pan_letter]=ind2sub([4 3],i);
    filename=[char(64+pan_letter),num2str(pan_index)];
    
    FullPath=dir([file_path,filename,'.mat']);
    fprintf(['Loading file: ',FullPath.name,'\n']);
    
    %load NFS data
    load([file_path,FullPath.name]);
    
    importedPanel{i}= import_panel(data);
    importedPanel{i}.type='panel';
    
    fprintf('Convert loaded data to dipoles data: %.3f mins\n', toc/60);
    
    %offset panel appropriately
    importedPanel{i}=panel_offset(...
        importedPanel{i},...
        panelPosition(i,1),...
        panelPosition(i,2),...
        panelPosition(i,3)...
        );
    %rotatepanel
    importedPanel{i}=panel_rotate(...
        importedPanel{i},...
        panelRotation(i,1),...
        panelRotation(i,2),...
        panelRotation(i,3),...
        locate(importedPanel{i})...
        );
    
    %% plot effective dipoles of panels for frequency 1
    figure(NFS_plot)
    index_subplot=sub2ind([3 4],pan_letter,5-pan_index);
    subplot(4,3,index_subplot)
    imagesc(importedPanel{i}.y(1,:),...
        importedPanel{i}.z(:,1),...
        permute(abs(importedPanel{i}.dipoles.y(1,:,:)),[2 3 1]));
    xlabel('Y')
    ylabel('Z')
    axis xy
    set(gca,'xdir','reverse')
    title(['Panel ',...
        char(64+pan_letter),num2str(pan_index),...
        'Dipoles f(1)'])
    
end

%% generate fields
 [xmtx ymtx zmtx]=meshgrid(1,A2_data.X(1,:),A2_data.Y(:,1));
 xGrid=xmtx(:).';
 yGrid=ymtx(:).';
 zGrid=zmtx(:).';
 nx=numel(unique(xGrid));
  ny=numel(unique(yGrid));
   nz=numel(unique(zGrid));
 fprintf('cal Panels fields \n');
 fields_importedPanel=dipoles_to_fieldsEXP3({importedPanel{panelData}},cat(2,xGrid.',yGrid.',zGrid.'));
fields_importedPanel{1}.E=fields_importedPanel{1}.E(:,:,1:101);

A2_field_plot=reshape(fields_importedPanel{1}.E(3,:,freq),[ny nx nz]);
A2_field_plot=squeeze(permute(A2_field_plot,[3 1 2]));
figure(A2_NFSat1m_plot)
subplot(2,2,2)
imagesc(A2_data.X(1,:),A2_data.Y(:,1),abs(A2_field_plot))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Propagated Fields (Z-Pol)')

A2_field_plot=reshape(fields_importedPanel{1}.E(3,:,freq),[ny nx nz]);
A2_field_plot=squeeze(permute(A2_field_plot,[3 1 2]));
figure(A2_NFSat1m_plot)
subplot(2,2,4)
imagesc(A2_data.X(1,:),A2_data.Y(:,1),angle(A2_field_plot))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Propagated Phase (Z-Pol)')


A3_field_plot=reshape(fields_importedPanel{2}.E(3,:,freq),[ny nx nz]);
A3_field_plot=squeeze(permute(A3_field_plot,[3 1 2]));
figure(A3_NFSat1m_plot)
subplot(2,2,2)
imagesc(A3_data.X(1,:),A3_data.Y(:,1),abs(A3_field_plot))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Propagated Fields (Z-Pol)')

A3_field_plot=reshape(fields_importedPanel{2}.E(3,:,freq),[ny nx nz]);
A3_field_plot=squeeze(permute(A3_field_plot,[3 1 2]));
figure(A3_NFSat1m_plot)
subplot(2,2,4)
imagesc(A3_data.X(1,:),A3_data.Y(:,1),angle(A3_field_plot))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Propagated Phase (Z-Pol)')

%% compare A2

figure(101)

norm_prop_fields=abs(max(max(A2_field_plot)));
norm_measured_fields=abs(max(max(A2_data.measurements(:,:,freq))));

subplot(2,2,4)
imagesc(A2_data.X(1,:),A2_data.Y(:,1),abs((A2_field_plot/norm_prop_fields)...
    -(A2_data.measurements(:,:,freq)/norm_measured_fields)))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Propagated Phase (Z-Pol)')

subplot(2,2,4)
imagesc(A2_data.X(1,:),A2_data.Y(:,1),angle((A2_field_plot/norm_prop_fields)...
    -(A2_data.measurements(:,:,freq)/norm_measured_fields)))
axis equal
axis xy
axis tight
set(gca,'Xdir','Reverse')
xlabel('y(m)')
ylabel('z(m)')
title('Propagated Phase (Z-Pol)')

%% register
a2_register=figure()

fixed=abs(A2_field_plot);
moving=abs(A2_data.measurements(:,:,freq));
subplot(1,2,1)
imshowpair(fixed, moving,'falsecolor');
[optimizer, metric] = imregconfig('multimodal');
set(gca,'Ydir','normal');
title('Measured/Propagated (Overlaid)')

registration_A2=imregister(moving, fixed, 'affine',optimizer, metric);
subplot(1,2,2)
imshowpair(fixed, registration_A2,'scaling','independent');
set(gca,'Ydir','normal');
title('Measured/Propagated (Registered)')
suptitle('A2 Panel')
% a3
a3_register=figure()

fixed=abs(A3_field_plot);
moving=abs(A3_data.measurements(:,:,freq));
subplot(1,2,1)
imshowpair(fixed, moving,'falsecolor');
[optimizer, metric] = imregconfig('multimodal');

set(gca,'Ydir','normal');
title('Measured/Propagated (Overlaid)')

registration_A3=imregister(moving, fixed, 'affine',optimizer, metric);
subplot(1,2,2)
imshowpair(fixed, registration_A3,'scaling','independent');
set(gca,'Ydir','normal');
title('Measured/Propagated (Registered)')
suptitle('A3 Panel')
