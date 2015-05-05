close all
clear all


c=3e8;
%direct measurement
addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\NSI Calibration\Direct VNA CAL of NSI cables'
cables_calibration=csvread('Phase_Cal_101pts.csv',8,0,[8,0,108,2]);
cables_response=[cables_calibration(:,2)+1i*cables_calibration(:,3)];

%load probe phase
load probephase.mat
addpath 'D:\Dropbox (Duke Electric & Comp)\MetaImager Data\panel to mini panel\MiniMeNFS\';
%load panel response and correct its phase

%load probe
addpath('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\printed dipole\Horizontal_Dipoles')
load P17.csv
 
numel_pol=1;

%reshape the data to be imported as a panel
data.f=linspace(17.5e9,26.5e9,101);
numel_freq=numel(data.f);

xGrid=unique(P17(:,2));
yGrid=unique(P17(:,3));
x_numel=numel(xGrid);
y_numel=numel(yGrid);

[data.x,data.y]=meshgrid(xGrid, yGrid);

data.measurements=reshape(10.^(P17(:,4)/20).*exp(j*(pi/180).*P17(:,5)),...
    x_numel,...
    y_numel,...
    numel_pol,...
    numel_freq...
    );

%permute to put into matlab format [yval xval freq pol]
data.measurements=permute(data.measurements, [2 1 4 3]);

for ii=1:101
data.measurements(:,:,ii,:)=data.measurements(:,:,ii,:)*exp(-j*probe_phase_meas(ii))/(cables_response(ii));
end
%only use polarization 1
data.measurements=data.measurements(:,:,:,1);

%plot printed dipole
figure;
for ii=1:numel_freq
imagesc(unique(data.x),unique(data.y),angle(data.measurements(:,:,ii)))
pause(.05)
end

