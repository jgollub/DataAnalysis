%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FMCW SAR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data and set parameters
close all
clear

%load data

%% test data: target at 50 cm
c=2.99792458*10^8;

%modulation
f0    = 24E9;           % Hz; 
BW    = 2E9;            % Hz; bandwidth
fc    = f0+BW/2;        % Hz;carrier frequency 
T     = 0.125E-3;%      % Sec, Sweep Time !!!!!
alpha =(BW/T);          % frequency ramp rate
Nt    = 1024;           % number of measurement points (mixed down)
% th   = 0:1/(2*(f0+BW)):T-1/(2*(f0+BW));  % fast time sampling (for simulator only)
th    = -T/2:1/(2*(BW+f0)):(T/2-1/(2*(BW+f0)));  % fast time sampling (for simulator only)
d_th  = mean(diff(th)); % time step (for simulator only)

%max range is less than d_ts*c?

% baseband sampling
d_ts = T/Nt;                    % baseband stepsize
ts   = linspace(-T/2,T/2-d_ts,Nt); % Baseband time sampling sampling
% ts   = linspace(-T/2,T/2-1/(2*(f0+BW)),Nt); % Baseband time sampling sampling
% d_ts = mean(diff(ts)); % Baseband time step 

% Simulation parameters

% Frequency ramp
ft =@(t) (fc.*rem(t,T/2) + (1/2)*alpha*rem(t,T/2).^2);   %chirp ramp
framp = @(t) (fc + alpha*rem(t,T/2));
st =exp(1.0j*2*pi*ft(th));       %chirp functon

% Received signal
% Delay time
tau = 2*5/c%2E-9 %ns
sr = exp(1.0j*2*pi*ft(th-tau)); %delayed signal
sm = sr.*conj(st);                  %mixed signal, note mixed in real domain

% figure
% test_sm = exp(1j*2*pi*(f0*tau+alpha*th*tau-(1/2)*alpha*tau^2));
% plot(th(sample_Nt),real(sm(sample_Nt))); hold on;
% plot(th(sample_Nt), test_sm(sample_Nt),'r');

%plot frequency ramp signal and and delayed return signal
figure(1); clf; 
subplot(3,1,1)
plot(th,framp(th),'-b'); hold on;
plot(th,framp(th-tau),'-r'); 
ylabel('Freq (Hz)');
xlabel('time (s)');
title('FMCW Signal Output & Return')

%plot mixed signal signal 
subplot(3,1,2) 
plot(th,real(st))
title('Generated Signal (Real Part)')
xlabel('Time (s)');
ylabel('Amplitude (s)');

subplot(3,1,3); cla; %plot out signal in frequency domain
th_freqs = (-1/(2*d_th):1/(T):1/(2*d_th)-1/T); 
st_FFT=ifftshift(fft(fftshift(st))); %convert to baseband & take fft
stb_FFT=ifftshift(fft(fftshift(exp(-1j*2*pi*fc*th).*st))); %convert to baseband & take fft

plot(th_freqs,abs(st_FFT),th_freqs,abs(stb_FFT)); hold on;
title('Baseband/Actual Signal, Frequency Domain')
xlabel('Freq (Hz)');
ylabel('Amplitude (a.u.)');

%plot mixed signal in time domain
% figure;plot(th(1:100),st(1:100),th(1:100),sr(1:100),th(1:100),sm(1:100))
figure(2); clf; 
subplot(2,1,1)
plot(th,real(sm));
title('Generated and Mixed Signal (Real Part)')
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');

%plot mixed signal in frequency domain
subplot(2,1,2)
plot(th_freqs,abs(ifftshift(d_th*fft(sm))));
title('Mixed Signal Fourier domain (Abs)')
xlabel('Freqs (Hz)');
ylabel('Amplitude (a.u.)');

sample_Nt=floor(1:length(th)/Nt:length(th)); %range 

figure(3); clf; %plot sampled signal
subplot(2,1,1)
plot(th(sample_Nt),sm(sample_Nt))
subplot(2,1,2)

ts_freqs=-1/(2*d_ts):1/T:1/(2*d_ts)-1/T;
sb=ifftshift(fft(sm(sample_Nt)))
plot(ts_freqs,abs(sb)); %frequency demand

%calc distance
[~,indx] = max(abs(sb(Nt/2+1:end)));
fb=ts_freqs(Nt/2+indx);

r=c*T*fb/(2*BW)
%beat frequency

% Filter
figure(32); clf;
% sm_filtered = lowpass(sm,.1,1/d_th); 
% plot(th, sm_filtered);
subplot(4,1,1);
 plot(ts, sm(sample_Nt));
subplot(4,1,2)

 rectangle = zeros(size(ts_freqs));
 
 filter_n = 100;
 rectangle(1:filter_n)=1;
 rectangle(end-filter_n+1:end)=1;
 
 testf = fftshift(ts_freqs); 
  plot(ts_freqs, abs(fftshift(d_ts*fft(sm(sample_Nt))))); hold on;
  plot(testf(filter_n)*ones(1,10),linspace(0,max(abs(fftshift(d_ts*fft(sm(sample_Nt))))),10));
 subplot(4,1,3)
  plot(ts_freqs,abs(ifftshift(rectangle.*fft(sm(sample_Nt)))));
  
 subplot(4,1,4)
   plot(ts,real(ifft(rectangle.*d_ts.*fft(sm(sample_Nt))))); 
 
 
%  sm_filtered = lowpass(sm(sample_Nt),1E4,1/d_ts,'steepness',.95); 
%  plot(ts, sm_filtered);
%  plot(ts, fftshift(fft(sm_filtered)));


% RVP
figure(31); clf;
subplot(3,1,1)
sRVP = exp(0*1j*pi*(ts_freqs).^2/alpha); %apply RVP filter !!!!!!!!!!!!

plot(ts_freqs,real(sRVP),ts_freqs,imag(sRVP),'--');

subplot(3,1,2)
sbsRVP = sb.*sRVP; %apply RVP filter

plot(ts_freqs, abs(sb),'b',ts_freqs,abs(sbsRVP),'r');

subplot(3,1,3); clf;
plot(ts,real(sm(sample_Nt)),'b'); hold on;
plot(ts,real(ifft(fftshift(sbsRVP))),'r')

%% (FREQUENCY DOMAIN) point source test source 
fmax = f0+BW;
lambdaMax = c/fmax;
kp =  linspace(2*pi*f0/c,2*pi*(f0+BW)/c,300).'; 
Lx = 0.6;
Ly = 1;
% 
xpath = -Lx/2:lambdaMax/4:Lx/2;
ypath = -1.5*ones(size(xpath));

% xpath =.005+[-fliplr(0:lambdaMax/4:Lx/2), lambdaMax/4:lambdaMax/4:Lx/2];
% ypath = -1.5*ones(size(xpath));

x1 = mean(xpath);           %measurement offset x
y1 = mean(ypath);            %measurement offset y   

% expected resolution
% xResEst = max(abs(y1*(lambdaMax)/Lx),lambdaMax/4);
xResEst = (lambdaMax/2);

yResEst = c/(2*BW);

display(['Estimated X (cross-range) Resolution: ', num2str(xResEst)]);
display(['Estimated Y (range) Resolution:       ', num2str(yResEst)]);

%measurement track path

%scene 
% maxRange =  c*ts(end)/2;
display(['Max Range (Determined by Pulse Repetition Time):', num2str(maxRange)]);

nXMeas = 2*ceil(Lx/(xResEst)); %double res
nYMeas = 2*ceil(Ly/(yResEst));

xvec = linspace(-Lx/2,Lx/2,nXMeas);
yvec = linspace(-Ly/2,Ly/2,nYMeas); %!not set to max range
%max range

%cell dimensions 
dx = abs(mean(diff(xvec)));
dy = abs(mean(diff(yvec)));

%scene center
flightAltitude = 0.001; %flight height

%  R=sqrt(flightAltitude^2+y1.^2); %%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!
% R = 0.01;


%imaging regime size



%max range measurable with beat frequency
rmax = c*(1/2)*(1/d_th)*T/(2*BW);

[xgrid,ygrid] = meshgrid(xvec, yvec);

target = targetAsymSmiley(3*floor(xResEst/dx)*dx,3*floor(yResEst/dy)*dy,xgrid,ygrid); %weird asymmetric target at resolution
%make 
figure(4); clf; %target scenes
subplot(2,1,1)
imagesc(xvec,yvec,target)
 axis tight; axis xy; axis equal;
% generate received_signal

sb=zeros(numel(kp),numel(xpath));

kpmtx = permute(kp,[2,3,1]);
     step = 2 %updates at # percent
for ii=1:numel(xpath)
    distance = sqrt((xgrid-xpath(ii)).^2 + (ygrid-ypath(ii)).^2); %from center
    
      %test straight up mixing real components
    sb(:,ii) = squeeze(squeeze(sum(target...
        .*(1./distance).^2.*exp(-1j*kpmtx.*(2*distance))...
        ,[1,2])));   %mix down to baseband

    if ~mod(ii,floor(numel(xpath)*step/100))
    disp(['percent done: ', num2str(ii/numel(xpath)*100)]) 
    end
%     Sb = ifftshift(fft(sb(:,ii))); %take fft to get to frequency regime
%     [~,indx] = max(abs(Sb(Nt/2+1:end))); %find max frequencie position (ignore negative frequencies)
%     image(5); subplot(3,1,1);
%     hold on
%     plot(ts,abs(Sb));
%     
%     subplot(3,1,2)
%     hold on;
%     
%     fb=ts_freqs(Nt/2+indx); %calc baseband frequency
%     r=c*T*fb/(2*BW); %determine range 
%     plot(distance(target(:)==1),r,'-o');
%     axis equal
%     drawnow
end
figure(4), subplot(2,1,2);
imagesc(xpath,ts, angle(sb)), axis xy;
%% (FREQUENCY DOMAIN) backprojection 

%reconstruction zone
%scene 

 [xi,yi] = meshgrid(-Lx/2:dx/2:Lx/2, (-Ly/2:dy/2:Ly/2)); %grid of xi, yi positions
% [xi,yi] = meshgrid(linspace(-Lx/2,Lx/2,nXMeas), linspace(-Ly/2,Ly/2,nYMeas)); %grid of xi, yi positions

%generate H
fest = zeros(numel(xi),1); %initialize estimate signal
% H=zeros(length(k)*length(xpath),numel(xi(:))); 

tic
display('')
display('Reconstruct:')
for ii=1:numel(xi)
Di = sqrt((xi(ii)-xpath).^2 + (yi(ii)-ypath).^2); %distance to all points from xpath

fest(ii) =  sum(((2*Di).*exp(1j*kp.*(2*Di))).*sb,[1,2]); %calc measurment matrix
%  fest(ii) =  sum(((Di).*exp(1j*kp.*(2*Di))).*sb,[1,2]); %calc measurment matrix

   if ~mod(ii,floor(numel(xi)*10/100))
    disp(['percent done: ', num2str((ii/numel(xi))*100)]); 
    end

end
toc
image = reshape(fest, [size(xi,1),size(xi,2)]);
image= abs(image)/max(abs(image(:)));

figure(44)
subplot(2,1,2)
% imagesc(xi(1,:), yi(:,1),image)
imagesc(xvec, yvec,image)
ax =gca;
axis xy; axis tight; axis equal;
 grid on; 
 ax.LineWidth = 0.5;
 ax.GridColor = 'r';
%  labelsXTick = linspace(min(xvec),max(xvec),floor(Lx/xResEst/10));
%  labelsYTick = linspace(min(yvec),max(yvec),floor(Ly/yResEst/4));
%  xticks(labelsXTick);
%  yticks(labelsYTick);
%    yticklabel(1:20:end)
%    xticklabel(1:20:end)
 %solve for image
% Kr = (4*pi/c)*(BW/T)*(fc/(BW/T)-ts);
% 
% Ky = sqrt(Kr.^2-Kx.^2)
%calculate fields at measurement position
%% (FREQUENCY DOMAIN) Matched Filter 

%reconstruction zone
%scene 

[xi,yi] = meshgrid(-Lx/2:dx/2:Lx/2, (-Ly/2:dy/2:Ly/2)); %grid of xi, yi positions
% [xi,yi] = meshgrid(linspace(-Lx/2,Lx/2,nXMeas), linspace(-Ly/2,Ly/2,nYMeas)); %grid of xi, yi positions

%generate H
fest = zeros(numel(xi),1); %initialize estimate signal
% H=zeros(length(k)*length(xpath),numel(xi(:))); 

display('')
display('Reconstruct:')
tic
for ii=1:numel(xpath)
Di = sqrt((xi-xpath(ii)).^2 + (yi-ypath(ii)).^2); %distance to all points from xpath
for jj=1:length(kpmtx)
% Hii  =  ((1./Di).^2.*exp(-1j*kpmtx(jj).*(2*Di))); %calc measurment matrix
%  Hii  =  ((Di).*exp(-1j*kpmtx(jj).*(2*Di))); %calc measurment matrix
  Hii  =  (exp(-1j*kpmtx(jj).*(2*Di))); %calc measurment matrix
fest =  fest+ conj(Hii(:)).*sb(jj, ii);
end

   if ~mod(ii,floor(numel(xpath)*10/100))
    disp(['percent done: ', num2str((ii/numel(xpath))*100)]); 
    end

end
toc

image = reshape(fest, [size(xi,1),size(xi,2)]);
image= abs(image)/max(abs(image(:)));

figure(44)
subplot(2,1,2)
% imagesc(xi(1,:), yi(:,1),image)
imagesc(xvec, yvec,image)
ax =gca;
axis xy; axis tight; axis equal;
 grid on; 
 ax.LineWidth = 0.5;
 ax.GridColor = 'r';
%  labelsXTick = linspace(min(xvec),max(xvec),floor(Lx/xResEst/10));
%  labelsYTick = linspace(min(yvec),max(yvec),floor(Ly/yResEst/4));
%  xticks(labelsXTick);
%  yticks(labelsYTick);
%    yticklabel(1:20:end)
%    xticklabel(1:20:end)
 %solve for image
% Kr = (4*pi/c)*(BW/T)*(fc/(BW/T)-ts);
% 
% Ky = sqrt(Kr.^2-Kx.^2)
%calculate fields at measurement position

%% point source test source
fmax = f0+BW;
lambdaMax = c/fmax;

Lx = 2.5*0.6;
Ly = 1;

xpath = -Lx/2:lambdaMax/4:Lx/2;
ypath = -1.5*ones(size(xpath));

x1 = mean(xpath);           %measurement offset x
y1 = mean(ypath);            %measurement offset y   

% expected resolution
% xResEst = abs(y1*lambdaMax/Lx);
xResEst = lambdaMax/Lx;
yResEst = c/(2*BW);

display(['Estimated X (cross-range) Resolution: ', num2str(xResEst)]);
display(['Estimated Y (range) Resolution:       ', num2str(yResEst)]);

%measurement track path

%scene 
maxRange =  c*ts(end)/2;
display(['Max Range (Determined by Pulse Repetition Time):', num2str(maxRange)]);

nXMeas = 2*ceil(Lx/(xResEst));
nYMeas = 2*ceil(Ly/(yResEst));

xvec = linspace(-Lx/2,Lx/2,nXMeas);
yvec = linspace(-Ly/2,Ly/2,nYMeas); %!not set to max range
%max range

%cell dimensions 
dx = abs(mean(diff(xvec)));
dy = abs(mean(diff(yvec)));

%scene center
flightAltitude = 0.001; %flight height

%  R=sqrt(flightAltitude^2+y1.^2); %%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!
% R = 0.01;


%imaging regime size



%max range measurable with beat frequency
rmax = c*(1/2)*(1/d_th)*T/(2*BW);

[xgrid,ygrid] = meshgrid(xvec, yvec);

target = targetAsymSmiley(xResEst,yResEst,xgrid,ygrid); %weird asymmetric target at resolution
%make 
figure(4); clf; %target scenes
subplot(2,1,1)
imagesc(xvec,yvec,target)
 axis tight; axis xy; axis equal;
% generate received_signal

sb=zeros(numel(ts),numel(xpath));

     step = 2 %updates at # percent
for ii=1:numel(xpath)
    distance = sqrt((xgrid-xpath(ii)).^2 + (ygrid-ypath(ii)).^2); %from center
    tdelay  = 2*distance/c; %time delay to each point
    tmtx = repmat(permute(ts,[1,3,2]),size(xgrid,1), size(xgrid,2)); %make matrix of time points
    
%     %straight up mixing real components
%     sb(:,ii) = squeeze(squeeze(sum(target...
%         .*1./distance.*real(exp(1j*2*pi*ft(tmtx).*(tmtx)))...
%         .*1./distance.*real(exp(1j*2*pi*ft(tmtx-tdelay).*(tmtx-tdelay)))...
%         ,[1,2])));   %mix down to baseband

%     %straight up mixing real components
%     sb(:,ii) = squeeze(squeeze(sum(target...
%         .*1./distance.*real(exp(1j*2*pi*ft(tmtx)))...
%         .*1./distance.*real(exp(1j*2*pi*ft(tmtx-tdelay)))...
%         ,[1,2])));   %mix down to baseband

%       %test straight up mixing real components
%     sb(:,ii) = squeeze(squeeze(sum(target...
%         .*1./distance.*exp(1j*2*pi*ft(tdelay))...
%         ,[1,2])));   %mix down to baseband
%     
%       %test straight up mixing real components
%     sb(:,ii) = squeeze(squeeze(sum(target...
%         .*1./distance.*exp(-1j*2*pi*ft(tdelay))...
%         ,[1,2])));   %mix down to baseband


    %using equation
       %straight up mixing real components
%     sb(:,ii) = squeeze(squeeze(sum(target...
%        .*exp(1j*2*pi*(f0*tdelay+alpha.*tdelay.*tmtx-(1/2)*alpha.*tdelay.^2))...
%         ,[1,2])));   %mix down to baseband
%     
%        sb(:,ii) = squeeze(squeeze(sum(target...
%        .*((1./distance).^2).*exp(1j*2*pi*(f0*tdelay+alpha.*tdelay.*tmtx))...
%         ,[1,2])));   %mix down to baseband
% 

sb(:,ii) = squeeze(squeeze(sum(target...
       .*exp(1j*2*pi*((fc)*tdelay+alpha.*tdelay.*tmtx))...
        ,[1,2])));   %mix down to baseband        

    if ~mod(ii,floor(numel(xpath)*step/100))
    disp(['percent done: ', num2str(ii/numel(xpath)*100)]) 
    end
%     Sb = ifftshift(fft(sb(:,ii))); %take fft to get to frequency regime
%     [~,indx] = max(abs(Sb(Nt/2+1:end))); %find max frequencie position (ignore negative frequencies)
%     image(5); subplot(3,1,1);
%     hold on
%     plot(ts,abs(Sb));
%     
%     subplot(3,1,2)
%     hold on;
%     
%     fb=ts_freqs(Nt/2+indx); %calc baseband frequency
%     r=c*T*fb/(2*BW); %determine range 
%     plot(distance(target(:)==1),r,'-o');
%     axis equal
%     drawnow
end
figure(4), subplot(2,1,2);
imagesc(xpath,ts, angle(sb)), axis xy;



%% Backprojection

%reconstruction zone
%scene 

[xi,yi] = meshgrid(-Lx/2:dx:Lx/2, (-Ly/2:dy:Ly/2)); %grid of xi, yi positions
% [xi,yi] = meshgrid(linspace(-Lx/2,Lx/2,nXMeas), linspace(-Ly/2,Ly/2,nYMeas)); %grid of xi, yi positions

%generate H
fest = zeros(numel(xi),1); %initialize estimate signal
% H=zeros(length(k)*length(xpath),numel(xi(:))); 

display('')
display('Reconstruct:')
for ii=1:numel(xpath)
Di = sqrt((xi(:)-xpath(ii)).^2 + (yi(:)-ypath(ii)).^2).'; %distance to all points from xpath
tdelayi=2*Di/c; %time delay
%measurement matrix; note di is n x 1 and k is 1 x n vector
% H(1+(ii-1)*length(k):ii*length(k),:) = (1./Di.*exp(-1j*k.*Di)...
%                                       .*1./Di.*exp(-1j*k.*Di)).';

% Hjj  = ((1./Di).*exp(-1j*k.*Di)).*((1./Di).*exp(-1j*k.*Di));
% fest = Hjj(:)'*meas
% 
% Hii  =((1./Di).*exp(1j*2*pi*ft(ts.').*(ts.'))).*((1./Di).*exp(1j*2*pi*ft(ts.'-tdelayi).*(ts.'-tdelayi))); %calc measurment matrix
%  fest = fest+ Hii'*sb(:,ii);

% Hii  =((1./Di).*real(exp(1j*2*pi*ft(ts.')))).*((1./Di).*real(exp(1j*2*pi*ft(ts.'-tdelayi)))); %calc measurment matrix
%   fest = fest+ Hii'*sb(:,ii);

% Hii  =((1./Di).*exp(1j*2*pi*ft(tdelayi))); %calc measurment matrix
%   fest = fest+ Hii'*sb(:,ii);


% Hii  =exp(1j*2*pi*(f0*tdelayi+alpha.*tdelayi.*(ts.'))); %calc measurement matrix
% fest = fest + Hii'*sb(:,ii);

Hii  =exp(1j*2*pi*((fc)*tdelayi+alpha.*tdelayi.*(-ts).')); %calc measurment matrix
% fest = fest+ Hii'*sb(:,ii);
% 
for jj = 1:numel(xi)
   fest(jj)= fest(jj) + sum(conv(sb(:,ii),conj(Hii(:,jj)),'same'));
end

   if ~mod(ii,floor(numel(xpath)*10/100))
    disp(['percent done: ', num2str((ii/numel(xpath))*100)]); 
    end



% does not need ranging (c/(2*BW*1/T))... solely based on phase info 
end

image = reshape(fest, [size(xi,1),size(xi,2)]);
% 
% image= ifftshift(fft(fftshift(image,1),[],1),1);
%  image= ifftshift(fft(fftshift(image,2),[],2),2);

image= abs(image)/max(abs(image(:)));


figure(44)
subplot(2,1,2)
% imagesc(xi(1,:), yi(:,1),image)
imagesc(xvec, yvec,image)
ax =gca;
axis xy; axis tight; axis equal;
 grid on; 
 ax.LineWidth = 0.5;
 ax.GridColor = 'r';
%  labelsXTick = linspace(min(xvec),max(xvec),floor(Lx/xResEst/10));
%  labelsYTick = linspace(min(yvec),max(yvec),floor(Ly/yResEst/4));
%  xticks(labelsXTick);
%  yticks(labelsYTick);
%    yticklabel(1:20:end)
%    xticklabel(1:20:end)
 %solve for image
% Kr = (4*pi/c)*(BW/T)*(fc/(BW/T)-ts);
% 
% Ky = sqrt(Kr.^2-Kx.^2)
%calculate fields at measurement position


%% Reconstruct g= Hf
%reconstruction zone
%scene 

[xi,yi] = meshgrid(-Lx/2:xResEst/2:Lx/2, (-Ly/2:yResEst/2:Ly/2)); %grid of xi, yi positions
% [xi,yi] = meshgrid(linspace(-Lx/2,Lx/2,nXMeas), linspace(-Ly/2,Ly/2,nYMeas)); %grid of xi, yi positions

%generate H
fest = zeros(numel(xi),1); %initialize estimate signal
% H=zeros(length(k)*length(xpath),numel(xi(:))); 


display('')
display('Reconstruct:')
for ii=1:numel(xpath)
    
Di = sqrt((xi(:)-xpath(ii)).^2 + (yi(:)-ypath(ii)).^2).'; %distance to all points from xpath
tdelayi=2*Di/c; %time delay
%measurement matrix; note di is n x 1 and k is 1 x n vector
% H(1+(ii-1)*length(k):ii*length(k),:) = (1./Di.*exp(-1j*k.*Di)...
%                                       .*1./Di.*exp(-1j*k.*Di)).';

% Hjj  = ((1./Di).*exp(-1j*k.*Di)).*((1./Di).*exp(-1j*k.*Di));
% fest = Hjj(:)'*meas
% 
% Hii  =((1./Di).*exp(1j*2*pi*ft(ts.').*(ts.'))).*((1./Di).*exp(1j*2*pi*ft(ts.'-tdelayi).*(ts.'-tdelayi))); %calc measurment matrix
%  fest = fest+ Hii'*sb(:,ii);

% Hii  =((1./Di).*real(exp(1j*2*pi*ft(ts.')))).*((1./Di).*real(exp(1j*2*pi*ft(ts.'-tdelayi)))); %calc measurment matrix
%   fest = fest+ Hii'*sb(:,ii);

Hii  =((1./Di).*real(exp(1j*2*pi*ft(-tdelayi)))); %calc measurment matrix
%   fest = fest+ Hii'*sb(:,ii);


% Hii  =exp(1j*2*pi*(f0*tdelayi+alpha.*tdelayi.*(ts.'))); %calc measurement matrix
% fest = fest + Hii'*sb(:,ii);

% Hii  =exp(1j*2*pi*((fc)*tdelayi+alpha.*tdelayi.*(-ts).')); %calc measurment matrix
% % fest = fest+ Hii'*sb(:,ii);
% % 
% sbase==exp(2*pi*f0*tdelasb(:,ii)
% for jj = 1:numel(xi)
%    fest(jj)= fest(jj) + sum(conv(sb(:,ii),conj(Hii(:,jj)),'same'));
% end

   if ~mod(ii,floor(numel(xpath)*10/100))
    disp(['percent done: ', num2str((ii/numel(xpath))*100)]); 
    end



% does not need ranging (c/(2*BW*1/T))... solely based on phase info 
end

image = reshape(fest, [size(xi,1),size(xi,2)]);
% 
% image= ifftshift(fft(fftshift(image,1),[],1),1);
%  image= ifftshift(fft(fftshift(image,2),[],2),2);

image= abs(image)/max(abs(image(:)));


figure(44)
subplot(2,1,2)
% imagesc(xi(1,:), yi(:,1),image)
imagesc(xvec, yvec,image)
ax =gca;
axis xy; axis tight; axis equal;
 grid on; 
 ax.LineWidth = 0.5;
 ax.GridColor = 'r';
%  labelsXTick = linspace(min(xvec),max(xvec),floor(Lx/xResEst/10));
%  labelsYTick = linspace(min(yvec),max(yvec),floor(Ly/yResEst/4));
%  xticks(labelsXTick);
%  yticks(labelsYTick);
%    yticklabel(1:20:end)
%    xticklabel(1:20:end)
 %solve for image
% Kr = (4*pi/c)*(BW/T)*(fc/(BW/T)-ts);
% 
% Ky = sqrt(Kr.^2-Kx.^2)
%calculate fields at measurement position

%% reconstruct with RMA

% [xi,yi] = meshgrid(-Lx/2:dx:Lx/2, (-Ly/2:dy:Ly/2)); %grid of xi, yi positions
dxpath     =mean(abs(diff(xpath)));

%zero pad
pad_Nx = 2.^(nextpow2(size(sb,2))); 
pad_Nr = 2.^(nextpow2(size(sb,1)));

% dxi = Lx/(pad_Nx-1);
% dri = Ly/(pad_Nx-1);
% Kr spatial frequency rage (min, max); from
% exp(1j*2*pi*(f0*tau+alpha*tau*ts) = exp(1j*2*Kr*r)
KrMin = (4*pi/c)*alpha*(f0/alpha+ts(1).');
KrMax = (4*pi/c)*alpha*(f0/alpha+ts(end).');
Krc = (KrMax+KrMin)/2;
dkr   = abs(mean(diff((4*pi/c)*alpha*(f0/alpha+ts.'))));

% dr    = (KrMax-KrMin)/(pad_Nr-1);
%  [Kx, Kr] = meshgrid(-2*pi/(2*dx):2*pi/Lx:2*pi/(2*dx),(4*pi/c)*alpha*(f0/alpha+ts.'));
[Kx, Kr] = meshgrid(-2*pi/(2*dxpath):2*pi/Lx:2*pi/(2*dxpath),(KrMin:dkr:KrMax));

%Calculate associated Ky component
Ky=real(sqrt(Kr.^2-Kx.^2)-Krc);
% Ky=real(Ky);

%plot fields pre-RMA application
figure(5)
subplot(1,3,1)
imagesc(xvec, ts, real(sb));  axis xy; axis tight; 
title('real(sb)')
subplot(1,3,2)
imagesc(xvec, ts, imag(sb));  axis xy; axis tight; 
title('imag(sb)')
subplot(1,3,3)
imagesc(xvec, ts, angle(sb)); axis xy; axis tight;
title('phase(sb)')

%conjugate signal 

sb_kr_xn = ifftshift(fft(sb*exp(-2*pi*fc/c),[],1),1); %!!!!!!!!!!!!!!!!!!!!!!

%slow axis

%fft slow axis

 sb_kr_kxn= ifftshift(fft(fftshift(sb_kr_xn,2),[],2),2);

% plot fields post fft
figure(6)
subplot(4,1,1)
imagesc(Kx(1,:),Kr(:,1), abs(sb_kr_kxn)); axis xy; axis tight;

title('spatial frequenc ydomain (phase)')
%shift
% sb_kr_kxn_mf=sb_kr_kxn.*exp(-1j*R*sqrt(Kr.^2-Kx.^2)-1j*Kx*x1);
sb_kr_kxn_mf=sb_kr_kxn.*exp(-1j*y1*sqrt(Kr.^2-Kx.^2)-1j*Kx*x1);

% add zero padding
% sb_kr_kxn_mf = padarray(sb_kr_kxn_mf,[ceil((pad_Nr-size(sb,1))/2), ceil((pad_Nx-size(sb,2))/2)],0,'pre'); %note, adding paddig to end here Kr
% % sb_kr_kxn = padarray(sb_kr_kxn,[floor((pad_Nr-size(sb,1))), floor((pad_Nx-size(sb,2))/2)],0,'post');
% sb_kr_kxn_mf = padarray(sb_kr_kxn_mf,[floor((pad_Nr-size(sb,1))/2), floor((pad_Nx-size(sb,2))/2)],0,'post');

% Kyi = padarray(Ky,[ceil((pad_Nx-size(sb,1))/2), ceil((pad_Nx-size(sb,2))/2)],0,'pre'); %note, adding paddig to end here Kr
% % sb_kr_kxn = padarray(sb_kr_kxn,[floor((pad_Nr-size(sb,1))), floor((pad_Nx-size(sb,2))/2)],0,'post');
% Kyi = padarray(Kyi,[floor((pad_Nx-size(sb,1))/2), floor((pad_Nx-size(sb,2))/2)],0,'post');

% Kyi = padarray(Ky,[ceil((pad_Nr-size(sb,1))/2), ceil((pad_Nx-size(sb,2))/2)],0,'pre');
% Kyi = padarray(Kyi,[floor((pad_Nr-size(sb,1))/2), floor((pad_Nx-size(sb,2))/2)],0,'post');

numSampling = 2^nextpow2(size(sb_kr_kxn_mf,1)+1000); %opt 2 (sample at next power of 2)
%  numSampling = size(sb_kr_kxn_mf,1); %opt 2 (sample at next power of 2)


% numSampling = 2^nextpow2(length(yvec)); %opt 2 (sample at next power of 2)

% numSampling = 2^(nextpow2(size(sb_kr_kxn,1))); %opt 2 (sample at next power of 2) 
Ky_min = min(min(Ky(Ky~=0)));
Ky_min = min(min(Ky));
Ky_max = max(max(Ky));

dky_resample = ((Ky_max-Ky_min-1/(Ky_max-Ky_min))/numSampling);
Ky_resample = linspace(-Ky_max,Ky_max-1/(Ky_max-Ky_min), numSampling);

figure; 
imagesc(Kx(1,:),Ky_resample(:,1), abs(sb_kr_kxn_mf)); axis xy; axis tight;
title('mag after stolt interpolation')
xlabel('Kx(rad/m)')
ylabel('Ky(rad/m)')

% Sxy=padarray(Sxy,[ceil((pad-size(measurements,1))/2), ceil((pad-size(measurements,2))/2)],0,'pre');
% Sxy=padarray(Sxy,[floor((pad-size(measurements,1))/2), floor((pad-size(measurements,2))/2)],0,'post');


figure(6)
subplot(4,1,2)
imagesc(Kx(1,:),Kr(:,1), angle(sb_kr_kxn_mf)); axis xy; axis tight;
title('mag post matched filter')
xlabel('Kx(rad/m)')
ylabel('Kr(rad/m)')

subplot(4,1,3)
imagesc(Kx(1,:),Kr(:,1), Ky); axis xy; axis tight;
title('abs(Ky)')
xlabel('Kx(rad/m)')
ylabel('Kr(rad/m)')


% dy_rescale=

% dy_resample=2*pi/(Ky_max-Ky_min);
% Ky_resample = (0:dy_resample:yvec(end)).';

sb_kx_ky = zeros(length(Ky_resample),size(sb_kr_kxn_mf,2)); 

for jj = 1:size(sb_kr_kxn_mf,2)
    indx_vec = squeeze(squeeze(real(Ky(:,jj)+Krc)~=0));    
    if ~(sum(indx_vec)==0 || sum(indx_vec)==1) %check that there are enough points to interpolate
        sb_kx_ky(:,jj) = interp1(Ky(indx_vec,jj), sb_kr_kxn_mf(indx_vec,jj), Ky_resample, 'linear',0);
    end
end

    %debug stolt interp step: plot interpolatio for z-slice
            jj=110;
%             indx_vec = squeeze(squeeze(real(Ky(:,jj))~=0));
            %debug
            figure(300); clf; 
            plot(squeeze(squeeze(Ky(:,jj))),real(squeeze(sb_kr_kxn_mf(:,jj))));
            hold on;            
            plot(real(Ky_resample(:)),real(squeeze(squeeze(sb_kx_ky(:,jj)))),'-o')
            drawnow; xlim([-Ky_max, Ky_max])            
%      end 
% end
figure(6)
subplot(4,1,4)
imagesc(Kx(1,:),Ky_resample(:,1).', abs(sb_kx_ky)); axis xy; axis tight;
title('mag after stolt interpolation')
xlabel('Kx(rad/m)')
ylabel('Ky(rad/m)')

% 
% dfy_resample=mean(diff(Ky_resample))/(2*pi);
% 
% (numSampling-1)*1/dfy_resample

% yveci=-2*pi/(2*dky_resample)-2*pi/(2*Ky_max-dky_resample):2*pi/(2*Ky_max+dky_resample):2*pi/(2*dky_resample);
yveci=linspace(-2*pi/(2*dky_resample),2*pi/(2*dky_resample),numSampling);


% xveci=linspace(-Lx/2,Lx/2,pad_Nx);
dKx = abs(mean(diff(Kx(1,:))));
BKx = (max(max(Kx))-min(min(Kx)));
xveci =( -2*pi/(2*dKx):2*pi/(BKx):2*pi/(2*dKx)); 

sb_kx_ky(find(isnan(sb_kx_ky))) = 0;

% 
% % sb_kx_y = (ifft(sb_kx_ky,[],1));
% % sb_x_y  = (ifft(sb_kx_y,[],2));
% 
sb_kx_y = ifftshift(ifft(fftshift(sb_kx_ky,1),[],1),1);
% % sb_kx_y = fftshift(ifft(sb_kx_ky,[],1),1);

% figure(10)
% subplot(2,1,1)
% imagesc(Kx(1,:),yveci,abs(sb_kx_y)); axis xy; axis tight;
% title('mag after fft Y direction')
% xlabel('Kx(rad/m)')
% ylabel('Ky(rad/m)')
% 
% subplot(2,1,2)
% imagesc(Kx(1,:),yveci,abs(fftshift(sb_kx_y,1))); axis xy; axis tight;
% title('mag after fft Y direction')
% xlabel('Kx(rad/m)')
% ylabel('Ky(rad/m)')


figure(7)
subplot(2,1,1)
imagesc(Kx(1,:),yveci,abs(sb_kx_y)); axis xy; axis tight;
title('mag after fft Y direction')
xlabel('Kx(rad/m)')
ylabel('Ky(rad/m)')

sb_x_y  =ifftshift(fft(fftshift(sb_kx_y,2),[],2),2);

subplot(2,1,2)
imagesc(xveci(1:end/20:end),yveci,abs(sb_x_y)); axis xy; axis tight; axis equal;
title('Reconstructed Image (fft X and Y)')
% % ylim([0,max(max(yveci))])
ylim([-4,4])
xlabel('X(m)')
ylabel('Y(m)')


%interpolate






% subplot(4,1,4)
% imagesc(xpath, yvec,abs(Sb_x_y)), axis xy, axis tight; axis equal;



