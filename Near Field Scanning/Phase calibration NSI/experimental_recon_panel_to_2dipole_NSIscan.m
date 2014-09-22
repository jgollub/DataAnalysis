close all;
clear all;
clc;


% Add path
% addpath 'C:\Users\Jonah Gollub\Documents\code\Virtualizer\graphics functions';
% addpath 'C:\Users\Jonah Gollub\Documents\code\Virtualizer\mexGpuFunctions';
freq=linspace(17.5e9,26.5e9,101);
c=3e8;

%%loading cables response
addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\NSI Calibration\Direct VNA CAL of NSI cables'
load cables_calibration.csv;
cables_response=10.^(cables_calibration(1:4:404,4)/20).*exp(j*cables_calibration(1:4:404,5)*pi/180);%10.^(cables_calibration(1:4:404,4)/20)

%load probe phase
load probephase.mat
addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\panel to mini panel\MiniMeNFS\';

%load panel response and correct its phase
load Large_panel.mat
data_panel=data;
for ff=1:101
     data.measurements(:,:,ff,:)=data.measurements(:,:,ff,:)*exp(-j*probe_phase_meas(ff));
end
panel=import_panel(data);panel.type='panel';
clear data


%load NSI scan
 addpath 'C:\Users\sf158\Dropbox (Duke Electric & Comp)\MetaImager Data (1)\Near Field Scans\printed dipole\';

load p1.csv;
probe_nfs=p1;
% load p1_noProbeCorrection.csv;
% probe_nfs=p1_noProbeCorrection;

%reshape the data to be imported as a panel
data.f=linspace(17.5e9,26.5e9,101);
data.y=reshape(probe_nfs(1:61*61,2),61,61)';
data.z=reshape(probe_nfs(1:61*61,3),61,61)';
% pol 1
for ii=1:101
    data.measurements(:,:,ii)=10.^(reshape(probe_nfs(((2*ii-2)*61*61+1):61*61*(2*ii-1),4),61,61)'/20).*(exp(reshape(j*pi/180*probe_nfs(((2*ii-2)*61*61+1):61*61*(2*ii-1),5),61,61)'))...
        *exp(-j*probe_phase_meas(ii))/(cables_response(ii));
end
%pol 2
% for ii=1:101
%     data.measurements(:,:,ii)=10.^(reshape(probe_nfs(((2*ii-1)*61*61+1):61*61*(2*ii),4),61,61)/20).*(exp(reshape(j*pi/180*probe_nfs(((2*ii-1)*61*61+1):61*61*(2*ii),5),61,61)))...
%         *exp(-j*probe_phase_meas(ff))/(cables_response(ff));
% end

data_probe=data;
probe1 = import_panel_nsi(data_probe);
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
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)'); xlim([-.01 .01]);


%%loading data
addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\MetaImager Scenes\Semi-Rigid\Printed Circuit Dipole Single Panel\';
load results_2probe_3cm_res_104cm.mat
g=savedata.g;

%region of interest
x_offset=1.04;
[locs,grid]=test_space2(x_offset,.2,0,.44,-.064,.2,.01);

%producing fields
sr1=dipoles_to_fieldsEXP3(panel,locs.locs);
srpr1=dipoles_to_fieldsEXP3(probe1,locs.locs);
srpr2=dipoles_to_fieldsEXP3(probe2,locs.locs);

%plotting fields
% figure(3)
%  scatter3(srpr1.Xmat,srpr1.Ymat,srpr1.Zmat,50,20*log10(sqrt(abs(srpr1.E_x(11,:)).^2+abs(srpr1.E_z(11,:)).^2)),'fill'); axis image; 
%  xlim([x_offset-0.01 x_offset+0.01]);
%  xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
% figure(4)
%  scatter3(sr1.Xmat,sr1.Ymat,sr1.Zmat,50,20*log10(sqrt(abs(sr1.E_x(81,:)).^2+abs(sr1.E_z(81,:)).^2)),'fill'); axis image; 
%  xlim([x_offset-0.01 x_offset+0.01]);
%  xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');

%% Calculate H from the above fields


H_recon=makeH_faceted({srpr1{1},srpr2{1}},sr1);

% H_recon=makeH_faceted(sr1,{srpr2{1},srpr1{1}});

% H_recon=makeH_faceted({srpr1{1},srpr2{1}},sr1);

%% Recalculate target grids

 f_est= (H_recon'*g);
 %f_est= cgs((H_recon'*H_recon),(H_recon'*g), 1e-3, 5);
 
for slice_index=1:grid(2)
    for slice_z=1:grid(3)
        f_est_tem(1:grid(1),slice_z)=abs(f_est( ( (slice_index-1)*grid(1)+(slice_z-1)*grid(1)*grid(2)+1 ):  ((slice_index-1)*grid(1)+(slice_z-1)*grid(1)*grid(2)+grid(1) ) ) )       /max(max(abs(f_est)));
    end
    f_est_reshaped(:,:,slice_index)=reshape(f_est_tem,grid(1), grid(3));
    figure
    imagesc(f_est_reshaped(:,:,slice_index).^2)
    colorbar
end
f_est_plot=reshape(f_est,[grid(2) grid(1) grid(3)]);
    f_est_plot=abs(f_est_plot./max(max(max((abs(f_est_plot))))));
figure
        vol3d('cdata',f_est_plot.^2,...
            'XData',locs.xrange,...
            'YData',locs.yrange,...
            'ZData',locs.zrange,...
            'Alpha',f_est_plot.^6 ...
            );
        
        axis equal;
        axis tight;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        view(-48,26);