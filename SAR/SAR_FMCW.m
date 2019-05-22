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

%scene center
flightAltitude = 0.5; %flight height
yCenter =  0.5; %offset
R=sqrt(flightAltitude^2+yCenter^2);

xpath = xvec;

yvec = 0:0.025:2 %!not set  to max range

%max range measurable with beat frequency
rmax = c*(1/2)*(1/d_th)*T/(2*BW);

%make 
figure(4); clf;

[xgrid,ygrid] = meshgrid(xvec, yvec);

target = zeros(size(xgrid(:,:)));
target(floor(end/3),floor(end/3)) = 1;
subplot(2,1,1)
imagesc(xvec,yvec,target)
axis equal; axis tight; axis xy;
% geneate received_signal

sb=zeros(numel(ts),numel(xpath));

figure(5); clf;
for ii=1:numel(xpath)
    distance = sqrt(R^2+(xgrid-xpath(ii)).^2 + ygrid.^2); %from center
    tdelay  = 2*distance/c; %time delay to each point
    tmtx = repmat(permute(ts,[1,3,2]),size(xgrid,1), size(xgrid,2));
    
    sb(:,ii) = squeeze(squeeze(sum(target...
        .*1./distance.*real(exp(1j*2*pi*ft(tmtx).*(tmtx)))...
        .*1./distance.*real(exp(1j*2*pi*ft(tmtx-tdelay).*(tmtx-tdelay)))...
        ,[1,2])));
    
    Sb = ifftshift(fft(sb(:,ii)));
    [~,indx] = max(abs(Sb(Nt/2+1:end)));
    subplot(3,1,1)
    hold on
    plot(ts,abs(Sb));
    
    subplot(3,1,2)
    hold on;
    
    fb=ts_freqs(Nt/2+indx);
    r=c*T*fb/(2*BW);
    plot(distance(target(:)==1),r,'-o');
    axis equal
    drawnow
end

%% Reconstruct g= Hf
%reconstruction zone

[xi,yi] = meshgrid(xvec, yvec);
%generate H
fest = zeros(numel(xi),1);
% H=zeros(length(k)*length(xpath),numel(xi(:))); 

for ii=1:numel(xpath)
    
Di = sqrt(R.^2+(xpath(ii)-xi(:)).^2 + yi(:).^2).';
tdelayi=2*Di/c;
%measurement matrix; note di is n x 1 and k is 1 x n vector
% H(1+(ii-1)*length(k):ii*length(k),:) = (1./Di.*exp(-1j*k.*Di)...
%                                       .*1./Di.*exp(-1j*k.*Di)).';
% Hjj  = ((1./Di).*exp(-1j*k.*Di)).*((1./Di).*exp(-1j*k.*Di));
% fest = Hjj(:)'*meas

Hii  = ((1./Di).*exp(1j*2*pi*ft(ts.').*(ts.'))).*((1./Di).*exp(1j*2*pi*ft(ts.'-tdelayi).*(ts.'-tdelayi)));
fest = fest+ Hii'*sb(:,ii);

end

image=reshape(fest, [size(xi,1),size(xi,2)])
image= abs(image)/max(abs(image(:)));

figure(4)
subplot(2,1,2)
imagesc(xvec, yvec,image)
axis equal; axis xy; axis tight;
%solve for image
% Kr = (4*pi/c)*(BW/T)*(fc/(BW/T)-ts);
% 
% Ky = sqrt(Kr.^2-Kx.^2)
%calculate fields at measurement position

%% reconstruct with RMA
f = linspace(f0,f0+BW,Nt);
x1 = mean(xpath);

dx=mean(diff(xpath));
Lx=abs(xpath(end)-xpath(1));

Ly=abs(yvec(end)-yvec(1));
dy=abs(mean(diff(yvec)));


[Kx, Ky] = meshgrid(-2*pi/(2*dx):2*pi/(Lx):2*pi/(2*dx),-2*pi/(2*dy):2*pi/(Ly):2*pi/(2*dy));
Kr = 4*pi*f/c;


figure(5)
% subplot(3,1,1)
% imagesc(abs(sb(1:114,:))); axis equal; axis tight; 

%fast axis 
sb_f_x=ifftshift(fft(sb,[],1),1);
sb_kr_x=sb_f_x*(2*pi/c);
subplot(4,1,1)
imagesc(xpath,Kr, abs(sb_kr_x)); axis xy; axis tight;  
%slow axis
sb_kr_kx= ifftshift(fft(fftshift(sb_kr_x,2),[],2),2);
subplot(4,1,2)
imagesc(Kx,Kr, abs(sb_kr_kx)); axis xy; axis tight;
sb_kx_ky=sb_kr_kx.*exp(1j*R*sqrt((Kr.').^2-Kx.^2)*y1-1j*Kx*x1);

subplot(4,1,3)
imagesc(Kx, Kr, abs((sb_kx_ky))); axis xy; axis tight;axis equal;


Sb_x_y= ifftshift(ifft(fftshift(ifftshift(ifft(fftshift(sb_kx_ky,1),[],1),1),2),[],2),2);


subplot(4,1,4)
imagesc(xpath, yvec,abs(Sb_x_y)), axis xy, axis tight; axis equal;



