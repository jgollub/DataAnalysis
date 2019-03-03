
c=299792458;

probe      = load('C:\Users\Jonah Gollub\Documents\code\data Analysis\SAR\Calibration Data\rawHornCalMohammadreza.mat');
%  rawData    = (squeeze(squeeze(probe.data.measurements(1,1,:)./abs(probe.data.measurements(1,1,:))))); %?? why positive phase advance

data       = load('C:\Users\Jonah Gollub\Documents\code\data Analysis\SAR\Calibration Data\rawHorn2018-2-18.mat');
rawData    = data.raw./abs(data.raw);




cal_offset = 0% 0.197;
dist       = 0.431+cal_offset;
f          = probe.data.f.';

k=2*pi*f/c;

figure(1); 
subplot(2,3,1); cla; plot(f, real(rawData));
subplot(2,3,2); cla; plot(f, db(rawData));
subplot(2,3,3); cla; plot(f, unwrap(angle(rawData)));

%correct for free space propagation
calData   = rawData./(exp(-1j*k*dist));


subplot(2,3,4), cla; plot(f, real(calData));
subplot(2,3,5), cla; plot(f, db(calData));

phaseCalData   = unwrap(angle(calData))/2;
subplot(2,3,6), cla; plot(f, phaseCalData);


% fextra=linspace(f(1),f(end),10000);
% figure(2); cla;
% plot(fextra,(exp(-1j*(2*pi*fextra/c)*dist)));
% hold on;
% plot(f, (rawData));
% legend('derived','raw data')
% 
% 
% hold on;
% probe      = load('C:\Users\Jonah Gollub\Documents\code\data Analysis\SAR\Calibration Data\rawHornCalMohammadreza.mat');
% rawData    = (squeeze(squeeze(probe.data.measurements(1,1,:)./abs(probe.data.measurements(1,1,:))))); %?? why positive phase advance
% 
% figure(1); hold on;
% subplot(2,3,1);hold on;  plot(f, real(rawData),'-r');
% subplot(2,3,2); hold on; plot(f, db(rawData),'-r');
% subplot(2,3,3);hold on; plot(f, unwrap(angle(rawData)),'-r');
% 
% %correct for free space propagation
% calData   = rawData./(exp(-1j*k*dist));
% 
% subplot(2,3,4), hold on; plot(f, real(calData),'-r');
% subplot(2,3,5), hold on; plot(f, db(calData),'-r');
% 
% phaseCalData   = unwrap(angle(calData))/2;
% 
% subplot(2,3,6),hold on; plot(f, phaseCalData,'-r');