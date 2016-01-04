%-------------------------------------------------------------------------------
%SAR Range Migration Algorithm 
%J. Gollub
%-------------------------------------------------------------------------------

%% load data and set parameters
c=3e8;
load('C:\Users\Jonah Gollub\Downloads\drill_data.mat');
measurements=data.measurements;

X=data.X/1000;
dx=X(1,2)-X(1,1);
Lx=X(1,end)-X(1,1);
[ynum,xnum]=size(X);

Y=data.Y/1000;
dy=Y(2,1)-Y(1,1);
Ly=Y(end,1)-Y(1,1);

f=data.f;
BW=f(end)-f(1);
clear data;
load('C:\Users\Jonah Gollub\Downloads\drill_background.mat');
%probe phase
load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data (1)\Scans\TEMP\probephase.mat');
pc=repmat(permute(probe_phase_meas,[3,2,1]),[size(X),1]);
measurements=(measurements-data.measurements).*exp(-j*0*.65*(pc));
%note: .65 is used to estimate effect of horns from different horn cal ata

%can apply offset if necessary to account for any missed calibration phase
%in cables, etc.

z_offset=0;
%%
%upsample in x and y

pad=2^nextpow2(max(numel(X(1,:)),numel(Y(:,1))));
pad=2^7;

Lx_pad=dx*(pad-1);
Ly_pad=dy*(pad-1);

% calc kx & ky vector information
[kux,kuy,~]=meshgrid(-(2*pi)/(2*dx):2*pi/(Lx_pad):(2*pi)/(2*dx),...
    -(2*pi)/(2*dy):2*pi/(Ly_pad):(2*pi)/(2*dy),...
    1:numel(f));

%calc free space k vector
k=2*pi*f/c;
k=repmat(permute(k,[1,3,2]),[pad,pad,1]);

%calculate plane wave decomposition (FFT). Note we have to be careful about
%shifting in fft because we are only acting along two of the dimensions.
%Therefor it is easier to use circshift (as opposed to ifftshif which acts
%on all dimentions.

Sxy=(1/(Lx_pad*Ly_pad))*circshift(fft2(measurements,pad,pad),[pad/2,pad/2]);

shiftx=dx*(pad/2-floor(xnum/2));
shifty=dy*(pad/2-floor(ynum/2));

Sxy=Sxy.*exp(-1.0j*kux*(shiftx)).*exp(-1.0j*kuy*(shifty)); %we must phase shift such that panel is centered

    %plot decomposition
    figHandle=figure(1);
    scrn = get( groot, 'Screensize'); scrn(1)=2*scrn(3)/3;  scrn(3)=scrn(3)/3;
    set(figHandle,'Position',scrn); clf; subplot(2,2,1);
    imagesc(-Lx_pad/2:dx:Lx_pad/2,-Ly_pad/2:dy:Ly_pad/2,abs(ifft2(fftshift(Sxy(:,:,1)))));
    title('Padded Input Fields f(1)'); axis equal; axis tight; xlabel('x (m)'); ylabel('y (m)')
    subplot(2,2,2);
    imagesc(kux(1,:),kuy(:,1),abs(Sxy(:,:,1)));
    title('FFT of Fields f(1)'); axis equal; axis tight; xlabel('kx (m)'); ylabel('ky (m)')
    drawnow
    
%calculate min max kz wavenumber
Kz=2*k;
kuz=sqrt((2*k).^2-kux.^2-kuy.^2);

%interpolate to evenly spaced grid
Srmg=zeros(size(kux));
for ii=1:size(kux,2)
    for jj=1:size(kuy,1)
        Srmg(jj,ii,:)=interp1(squeeze(kuz(jj,ii,:)),...
            squeeze(Sxy(jj,ii,:)),...
            squeeze(Kz(jj,ii,:)),'linear');
    end
end

Srmg(find(isnan(Srmg))) = 0; %set all Nan values to 0

%apply inverst FFT to get image
fxy = ifftn(fftshift(Srmg.*exp(-1.0j*Kz*z_offset)));
image=(abs(fxy)/max(abs(fxy(:))));
%% plot 

%labeling
xx=linspace(-Lx_pad/2,Lx_pad/2,size(fxy,2));
yy=linspace(-Ly_pad/2,Ly_pad/2,size(fxy,1));
zz=linspace(0,c*numel(f)/(2*BW),numel(f));

    figure(1); subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',image.^2,'xdata',xx,'Ydata',yy,'Zdata',zz);
    axis equal; axis tight; view(3);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');

%% plot subimage
xmin=-0.15;
xmax=0.25;
ymin=-0.25;
ymax=0;
zmin=0.7;
zmax=.85;

subimage=image(yy>ymin & yy<ymax, xx>xmin & xx<xmax, zz>zmin & zz<zmax);
% subimage=subimage/max(subimage(:));
    figure(1); subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',subimage.^2,...
        'xdata',xx(xx>xmin & xx<xmax),...
        'Ydata',yy(yy>ymin & yy<ymax),...
        'Zdata',zz(zz>zmin & zz<zmax)...
        );
    zlabel('downrange (m)'); xlabel(' x crossrange (m)');  ylabel('y crossrange (m)');
    axis equal; axis tight; view(3);

