%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FMCW SAR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data and set parameters
c=2.99792458*10^8;

%load data

%% test data: target at 50 cm

%modulation
f0  = 24E9;            % GHz
BW  = 2E9;             % GHz
T   = 0.125E-3;        % ms Sweep Time !!!!!
Nt   = 256%1024;           % number of measurement points (mixed down)
th  = 0:1/(26E9)/2:T-1/(26E9)/2;  % fast time sampling (for simulator only)
d_th = mean(diff(th)); % time step (for simulator only)
alpha=(BW/T);          % frequency ramp rate

% baseband sampling
ts  =linspace(0,T,Nt); % Baseband time sampling sampling
d_ts = mean(diff(ts)); % Baseband time step 

% Simulation parameters

% Frequency ramp
ft =@(t) (f0 + 1/2*BW*mod(t,T)/T);   %chirp ramp
st =exp(1.0j*2*pi*ft(th).*th);       %chirp functon

% Received signal
% Delay time
tau = 2*.06/c%2E-9 %ns
sr = exp(1.0j*2*pi*ft(th-tau).*(th-tau)); %delayed signal
sm = real(st).*real(sr);                  %mixed signal, note mixed in real domain

% figure
% test_sm = exp(1j*2*pi*(f0*tau+alpha*th*tau-(1/2)*alpha*tau^2));
% plot(th(sample_Nt),real(sm(sample_Nt))); hold on;
% plot(th(sample_Nt), test_sm(sample_Nt),'r');

figure(1); clf; %plot frequency ramp
subplot(3,1,1)
plot(th,ft(th)); hold on;
plot(th,ft(th-tau)); 


subplot(3,1,2) %plot mixed siginal signal 
plot(th,real(sm))

subplot(3,1,3); %plot out signal in frequency domain
th_freqs = -1/(2*d_th):1/(T):1/(2*d_th)-1/T; 
st_FFT=ifftshift(fft(st));
plot(th_freqs,abs(st_FFT)); hold on;

figure(2); clf; %plot mixed signal in frequency domain
subplot(2,1,1)
plot(th_freqs,real(sm));
subplot(2,1,2)
plot(th_freqs,ifftshift(fft(sm)));

sample_Nt=floor(1:length(th)/Nt:length(th)); %range 

figure(3); clf; %plot sampled signal
subplot(2,1,1)
plot(th(sample_Nt),sm(sample_Nt))
subplot(2,1,2)

ts_freqs=-1/(2*d_ts):1/T:1/(2*d_ts);
sb=ifftshift(fft(sm(sample_Nt)))
plot(ts_freqs,abs(sb)); %frequency demand

%calc distance
[~,indx] = max(abs(sb(Nt/2+1:end)));
fb=ts_freqs(Nt/2+indx);

r=c*T*fb/(2*BW)
%beat frequency

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



%% point source test source
fmax = f0+BW;
lambdaMax = c/fmax;

%measurement track path

%scene 
xvec = 0:0.005:1;
yvec = 0:0.0075:.2 %!not set to max range

% yvec = -.1:lambdaMax/2:.1;

%scene center
flightAltitude = 0.01; %flight height
yCenter = mean(yvec) ; %offset
 R=sqrt(flightAltitude^2+yCenter^2); %%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!
% R = 0.01;
xpath = xvec(1:end);


%max range measurable with beat frequency
rmax = c*(1/2)*(1/d_th)*T/(2*BW);

%make 
figure(4); clf; %target scene

[xgrid,ygrid] = meshgrid(xvec, yvec);

target = zeros(size(xgrid(:,:)));
target(floor(4*end/5),floor(4*end/5)) = 1;
target(floor(end/3),floor(end/3)) = 1;

subplot(2,1,1)
imagesc(xvec-x1,yvec-y1,target)
axis equal; axis tight; axis xy;
% generate received_signal


sb=zeros(numel(ts),numel(xpath));

 status = floor(linspace(0,numel(xpath),11));
for ii=1:numel(xpath)
    distance = sqrt(R^2+(xgrid-xpath(ii)).^2 + ygrid.^2); %from center
    tdelay  = 2*distance/c; %time delay to each point
    tmtx = repmat(permute(ts,[1,3,2]),size(xgrid,1), size(xgrid,2)); %make matrix of time points
    
%     %straight up mixing real components
%     sb(:,ii) = squeeze(squeeze(sum(target...
%         .*1./distance.*real(exp(1j*2*pi*ft(tmtx).*(tmtx)))...
%         .*1./distance.*real(exp(1j*2*pi*ft(tmtx-tdelay).*(tmtx-tdelay)))...
%         ,[1,2])));   %mix down to baseband

    %using equation
       %straight up mixing real components
%     sb(:,ii) = squeeze(squeeze(sum(target...
%        .*exp(1j*2*pi*(f0*tdelay+alpha.*tdelay.*tmtx-(1/2)*alpha.*tdelay.^2))...
%         ,[1,2])));   %mix down to baseband
%     
       sb(:,ii) = squeeze(squeeze(sum(target...
       .*((1./distance).^2).*exp(1j*2*pi*(f0*tdelay+alpha.*tdelay.*tmtx))...
        ,[1,2])));   %mix down to baseband
% 
%        sb(:,ii) = squeeze(squeeze(sum(target...
%        .*exp(1j*2*pi*(f0*tdelay+alpha.*tdelay.*tmtx))...
%         ,[1,2])));   %mix down to baseband    


    if sum(ii == status) == 1
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

%% Reconstruct g= Hf
%reconstruction zone

[xi,yi] = meshgrid(xvec, yvec); %grid of xi, yi positions
%generate H
fest = zeros(numel(xi),1); %initialize estimate signal
% H=zeros(length(k)*length(xpath),numel(xi(:))); 

for ii=1:numel(xpath)
    
Di = sqrt(R.^2+(xpath(ii)-xi(:)).^2 + yi(:).^2).'; %distance to all points from xpath
tdelayi=2*Di/c; %time delay
%measurement matrix; note di is n x 1 and k is 1 x n vector
% H(1+(ii-1)*length(k):ii*length(k),:) = (1./Di.*exp(-1j*k.*Di)...
%                                       .*1./Di.*exp(-1j*k.*Di)).';
% Hjj  = ((1./Di).*exp(-1j*k.*Di)).*((1./Di).*exp(-1j*k.*Di));
% fest = Hjj(:)'*meas

% Hii  =((1./Di).*exp(1j*2*pi*ft(ts.').*(ts.'))).*((1./Di).*exp(1j*2*pi*ft(ts.'-tdelayi).*(ts.'-tdelayi))); %calc measurment matrix
% fest = fest+ Hii'*sb(:,ii);

Hii  =exp(1j*2*pi*(f0*tdelayi+alpha.*tdelayi.*ts.')); %calc measurment matrix
fest = fest+ Hii'*sb(:,ii);

% does not need ranging (c/(2*BW*1/T))... solely based on phase info 
end

image = reshape(fest, [size(xi,1),size(xi,2)])
image= abs(image)/max(abs(image(:)));

figure(4)
subplot(2,1,2)
imagesc(xvec-x1, yvec-y1,image)
axis equal; axis xy; axis tight;
%solve for image
% Kr = (4*pi/c)*(BW/T)*(fc/(BW/T)-ts);
% 
% Ky = sqrt(Kr.^2-Kx.^2)
%calculate fields at measurement position

%% reconstruct with RMA
% f  = linspace(f0,f0+BW,Nt); %freq Range
x1 = mean(xpath);           % center of image, x
y1 = mean(yvec);            % center of image, y                

%imaging regime size
Lx = abs(xvec(end)-xvec(1));
Ly = abs(yvec(end)-yvec(1));

%stepsize of measurements   
dx = abs(mean(diff(xpath)));

%zero pad
pad_Nx = 2.^(nextpow2(size(sb,2))); 
pad_Nr = 2.^(nextpow2(size(sb,1)));

dxi = Lx/(pad_Nx-1);

% Kr spatial frequency rage (min, max); from
% exp(1j*2*pi*(f0*tau+alpha*tau*ts) = exp(1j*2*Kr*r)
KrMin = (4*pi/c)*alpha*(f0/alpha+ts(1).');
KrMax = (4*pi/c)*alpha*(f0/alpha+ts(end).');
dkr   = abs(mean(diff((4*pi/c)*alpha*(f0/alpha+ts.'))));

% dr    = (KrMax-KrMin)/(pad_Nr-1);

% [Kx, Kr] = meshgrid(-2*pi/(2*dx):2*pi/(Lx):2*pi/(2*dx),(4*pi/c)*alpha*(f0/alpha+ts.'));
%   [Kx, Kr] = meshgrid(linspace(-2*pi/(2*dx),2*pi/(2*dx)-2*pi/(Lx),size(sb,2)),(4*pi/c)*alpha*(f0/alpha+ts.'));
% [Kx, Kr] = meshgrid(linspace(-2*pi/(2*dxi),2*pi/(2*dxi),pad_Nx),linspace(KrMin,KrMax, pad_Nr));
[Kx, Kr] = meshgrid(linspace(-2*pi/(2*dxi),2*pi/(2*dxi),pad_Nx),0:dkr:KrMax);



%Calcualte associated Ky component
Ky=sqrt(Kr.^2-Kx.^2);
Ky=real(Ky);

%plot fields pre-RMA application
figure(5)
subplot(1,3,1)
imagesc(xvec, ts, real(sb));  axis xy; axis tight; 

subplot(1,3,2)
imagesc(xvec, ts, imag(sb));  axis xy; axis tight; 

subplot(1,3,3)
imagesc(xvec, ts, angle(sb)); axis xy; axis tight;

%conjugate signal 
wc = 2*pi*f';
wc = repmat(wc,1,size(sb,2));

sb_kr_xn = conj(sb.*exp(-j*wc.*ts')); %!!!!!!!!!!!!!!!!!!!!!!

%slow axisimagesc(xpath,Kr(1,:), abs(sb_kr_xn)); axis xy; axis tight;  

%fft slow axis
clear sb_kr_kxn;
sb_kr_kxn= ifftshift(fft(fftshift(sb_kr_xn,2),[],2),2);


% add zero padding
sb_kr_kxn = padarray(sb_kr_kxn,[ceil((length(Kr(:,1))-size(sb,1))), ceil((pad_Nx-size(sb,2))/2)],0,'pre'); %note, adding paddig to end here Kr
% sb_kr_kxn = padarray(sb_kr_kxn,[floor((pad_Nr-size(sb,1))), floor((pad_Nx-size(sb,2))/2)],0,'post');
sb_kr_kxn = padarray(sb_kr_kxn,[0, floor((pad_Nx-size(sb,2))/2)],0,'post');



% Sxy=padarray(Sxy,[ceil((pad-size(measurements,1))/2), ceil((pad-size(measurements,2))/2)],0,'pre');
% Sxy=padarray(Sxy,[floor((pad-size(measurements,1))/2), floor((pad-size(measurements,2))/2)],0,'post');

% plot fields post fft
figure(6)
subplot(4,1,1)
imagesc(Kx(1,:),Kr(:,1), abs(sb_kr_kxn)); axis xy; axis tight;
imagesc(abs(sb_kr_kxn)); axis xy; axis tight;

title('spatial frequenc ydomain (phase)')
%shift
sb_kr_kxn_mf=sb_kr_kxn.*exp(1j*R*sqrt(Kr.^2-Kx.^2)-1j*1/sqrt(2)*Kx*x1);

subplot(4,1,2)
imagesc(Kx(1,:),Kr(:,1), angle(sb_kr_kxn_mf)); axis xy; axis tight;
title('mag post matched filter')
xlabel('Kx(rad/m)')
ylabel('Kr(rad/m)')

% Ky_min = min(min(Ky(Ky~=0)));
% Ky_min =0;
% Ky_max = max(max(Ky));

subplot(4,1,3)
imagesc(Kx(1,:),Kr(:,1), Ky); axis xy; axis tight;
title('abs(Ky)')
xlabel('Kx(rad/m)')
ylabel('Kr(rad/m)')

% Ky = padarray(Ky,[0, ceil((pad_Nx-size(sb,2))/2)],0,'pre');
% Ky = padarray(Ky,[floor((pad_Nr-size(sb,1))), floor((pad_Nx-size(sb,2))/2)],0,'post');

numSampling = 2^nextpow2(size(sb_kr_kxn_mf,1)+10000); %opt 2 (sample at next power of 2)
 numSampling = size(sb_kr_kxn_mf,1); %opt 2 (sample at next power of 2)


% numSampling = 2^nextpow2(length(yvec)); %opt 2 (sample at next power of 2)

% numSampling = 2^(nextpow2(size(sb_kr_kxn,1))); %opt 2 (sample at next power of 2) 


dky_resample = 2*Ky_max/(numSampling-1);
Ky_resample = linspace(-Ky_max,Ky_max, numSampling).';
% dy_rescale=

% dy_resample=2*pi/(Ky_max-Ky_min);
% Ky_resample = (0:dy_resample:yvec(end)).';

sb_kx_ky = zeros(length(Ky_resample),size(sb_kr_kxn_mf,2)); 

% Stolt interpolation
for jj = 1:size(sb_kr_kxn_mf,2)
   indx_vec = squeeze(squeeze(real(Ky(:,jj))~=0));
    if ~(sum(indx_vec)==0 || sum(indx_vec)==1) %check that there are enough points to interpolate
    sb_kx_ky(:,jj) = interp1(Ky(indx_vec,jj), sb_kr_kxn_mf(indx_vec,jj), Ky_resample, 'linear'); 
         end 
end

    %debug stolt interp step: plot interpolatio for z-slice
            jj=156;
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
imagesc(Kx(1,:),Ky_resample(:,1), abs(sb_kx_ky)); axis xy; axis tight;
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
BKx = 2*max(max(Kx));
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
imagesc(xveci,yveci,abs(sb_x_y)); axis xy; axis tight; axis equal;
title('Reconstructed Image (fft X and Y)')
% % ylim([0,max(max(yveci))])
ylim([-.1,.1])
xlabel('X(m)')
ylabel('Y(m)')


%interpolate






% subplot(4,1,4)
% imagesc(xpath, yvec,abs(Sb_x_y)), axis xy, axis tight; axis equal;



