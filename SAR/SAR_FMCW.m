%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FMCW SAR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data and set parameters
c=2.99792458*10^8;

%load data


%% test data: target at 50 cm

%modulation
f0  = 24E9; %GHz
BW  = 2E9;  %GHz
T   = 0.125E-3; % ms Sweep Time
Nt   = 1024;
th  = 0:1/(26E9)/2:T;
d_th = mean(diff(th));

% baseband sampling
ts  =linspace(0,T,Nt);
d_ts = mean(diff(ts));
% delay time
tau =2E-9 %ns
%Frequency ramp
ft =@(t) (f0 + 1/2*BW*mod(t,T)/T);
st =exp(1.0j*2*pi*ft(th).*th);

% received signal
sr = exp(1.0j*2*pi*ft(th-tau).*(th-tau));
sm = real(st).*real(sr);

figure(1); clf;
subplot(3,1,1)
plot(th,ft(th)); hold on;
plot(th,ft(th-tau)); 

subplot(3,1,2)
plot(th,real(sm))

subplot(3,1,3);
th_freqs = -1/(2*d_th):1/(T):1/(2*d_th);
st_FFT=ifftshift(fft(st));
plot(th_freqs,abs(st_FFT)); hold on;

figure(2); clf;
subplot(2,1,1)
plot(th_freqs,real(sm));
subplot(2,1,2)
plot(th_freqs,ifftshift(fft(sm)));

sample_Nt=floor(1:length(th)/Nt:length(th));

figure(3); clf;
subplot(2,1,1)
plot(th(sample_Nt),sm(sample_Nt))
subplot(2,1,2)

ts_freqs=-1/(2*d_ts):1/T:1/(2*d_ts);
sb=ifftshift(fft(sm(sample_Nt)))
plot(ts_freqs,abs(sb));


%calc distance
[~,indx] = max(abs(sb(Nt/2+1:end)));
fb=ts_freqs(Nt/2+indx)

r=c*T*fb/(2*BW)
%beat frequency

%% point source test source
fmax = f0+BW;
lambdaMax = c/fmax;

%measurement track path

%scene 
xvec = 0:.075:1;
yvec = -.1:lambdaMax/2:.1;

xpath = xvec;

[xgrid,ygrid,tau] = meshgrid(xvec, yvec,0:1/(2*26E9):4E-9);
target = zeros(size(xgrid));
target(floor(end/3),floor(end/3),105) = 1;

%calculate fields at measurement position
target.*exp(-1.0j*2*pi*ft(tau)*tau)*exp(-1.0*j*pi*ft(tau
