close all;
clear all;
clc;

% Add path
addpath 'C:\Users\Jonah Gollub\Documents\code\Virtualizer\graphics functions';
addpath 'C:\Users\Jonah Gollub\Documents\code\Virtualizer\mexGpuFunctions';
freq=linspace(17.5e9,26.5e9,101);
numel_freq=numel(freq);
c=3e8;

% %loading cables response
% addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\NSI Calibration\';
% load cables_calibration2.csv;
% cables_response=10.^(cables_calibration2(1:4:404,4)/20).*exp(j*cables_calibration2(1:4:404,5)*pi/180);%10.^(cables_calibration(1:4:404,4)/20)

%direct measurement
addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\NSI Calibration\Direct VNA CAL of NSI cables'
cables_calibration=csvread('Phase_Cal_101pts.csv',8,0,[8,0,108,2]);
cables_response=[cables_calibration(:,2)+1i*cables_calibration(:,3)];

%load probe phase
load probephase.mat
addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\panel to mini panel\MiniMeNFS\';

%load panel response and correct its phase
load Large_panel.mat
data_panel=data;
for ff=1:101
     data.measurements(:,:,ff,:)=data.measurements(:,:,ff,:)*exp(-j*probe_phase_meas(ff));
end
panel=import_panel(data);
panel.type='panel';
clear data

%load NSI scan
 addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\printed dipole\';
numel_pol=2;
% load p1.csv;
% probe_nfs=p1;
load p1_noProbeCorrection.csv;
probe_nfs=p1_noProbeCorrection;

%reshape the data to be imported as a panel
data.f=linspace(17.5e9,26.5e9,101);

x_numel=numel(unique(probe_nfs(:,2)));
y_numel=numel(unique(probe_nfs(:,3)));
total_meas_pts=x_numel*y_numel;
data.y=reshape(probe_nfs(1:total_meas_pts,2),y_numel,x_numel);
data.x=reshape(probe_nfs(1:total_meas_pts,3),y_numel,x_numel);

data.measurements=reshape(10.^(probe_nfs(:,4)/20).*exp(j*(pi/180).*probe_nfs(:,5)),...
    x_numel,...
    y_numel,...
    numel_pol,...
    numel_freq...
    );

%permute to put into matlab format [yval xval freq pol]
data.measurements=permute(data.measurements, [1 2 4 3]);

for ii=1:101
data.measurements(:,:,ii,:)=data.measurements(:,:,ii,:)*exp(-j*probe_phase_meas(ii))/(cables_response(ii));
end
%only use polarization 1
data.measurements=data.measurements(:,:,:,1);

data_probe=data;
probe1 = import_panel_NSI_measured(data_probe);
probe1.type='panel';
probe2=probe1;
clear data

%place panel and probes at the correct location
panel=panel_offset(panel,.086,0,0);
probe2=panel_offset(probe2,0.075,-0.404,0);
probe1=panel_offset(probe1,0.075,0.395,0);

% Plot panels and "on elements"
figure
panel_plot(1,panel,probe1,probe2); axis image;
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');

%%loading data
 addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\MetaImager Scenes\Semi-Rigid\Printed Circuit Dipole Single Panel\';
load results_2probe_3cm_res_104cm.mat
g=savedata.g;

%region of interest
x_offset=1.05;
[locs,grid]=test_space2(x_offset,.04,0,.44,-.064,.2,.005);

%producing fields

srpr=dipoles_to_fieldsEXP3({probe1, probe2},locs.locs);
sr=dipoles_to_fieldsEXP3(panel,locs.locs);

%% Calculate H from the above fields

H_recon=makeH_faceted({srpr{:}},sr);
%% Recalculate target grids

%  f_est= (H_recon'*g);
 f_est= cgs((H_recon'*H_recon),(H_recon'*g), 1e-3, 5);
 
 f_est_reshaped=zeros([ grid(3) grid(1) grid(2)]);
for slice_index=1:grid(2)
    for slice_y=1:grid(3)
        f_est_tem(1:grid(1),slice_y)=abs(...
            f_est(((slice_index-1)*grid(1)+(slice_y-1)*grid(1)*grid(2)+1):...
            ((slice_index-1)*grid(1)+(slice_y-1)*grid(1)*grid(2)+grid(1) ) ) )...
            /max(max(abs(f_est)));
    end
    f_est_reshaped(:,:,slice_index)=permute(reshape(f_est_tem,grid(1), grid(3)),[2 1]);
    figure
    imagesc(locs.yrange,locs.zrange,f_est_reshaped(:,:,slice_index).^2,[0 1])
    set(gca,'Xdir','reverse')
    axis equal
    axis tight
    colorbar
end
f_est_plot=reshape(f_est,grid);
    f_est_plot=abs(f_est_plot./max(max(max((abs(f_est_plot))))));
figure
        vol3d('cdata',f_est_plot.^2,...
            'XData',locs.xmat,...
            'YData',locs.ymat,...
            'ZData',locs.zmat,...
            'Alpha',f_est_plot.^6 ...
              );
        
        axis equal;
        axis tight;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        view(-48,26);