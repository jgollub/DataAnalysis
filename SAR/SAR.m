%-------------------------------------------------------------------------------
%SAR Range Migration Algorithm 
%Jonah Gollub
%-------------------------------------------------------------------------------


%% load data and set parameters
c=2.99792458*10^8;

% data 401 freq pts
% load('C:\Users\Jonah Gollub\Downloads\drill_data.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\D_38cmPastWall_60X60cm_WallMeasSeries_14-Feb-2019_3_32.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Saw_46cmAway_60X60cm_CorrectAlign_09-Feb-2019_4_59.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\ReferenceFrame_382mmPastWall_60X60cm_WallMeasSeries_19-Feb-2019_3_21.mat');

% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\2By4_InsideWall_60X60cm_WallMeasSeries_21-Feb-2019_2_56.mat');
 load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Baseline_InsideWall_60X60cm_RealWallMeasSeries_23-Feb-2019_1_29.mat');
%   load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Plywood_InsideWall_60X60cm_RealWallMeasSeries_24-Feb-2019_2_56.mat');

% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Drywall_InsideWall_60X60cm_RealWallMeasSeries_25-Feb-2019_19_40.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_26-Feb-2019_6_13.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\RandomObjects_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_26-Feb-2019_19_39.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Fiducial_wDryWall_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_28-Feb-2019_19_20.mat');
  
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\ShiftedRight_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_27-Feb-2019_4_48.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\ShiftedRight_wDryWall_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_27-Feb-2019_21_13.mat');

%101 freq  pts
%  load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Fiducial_wTwosidedDryWall_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_01-Mar-2019_18_9.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\MonsterScan_Fiducial_wTwosidedDryWall_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_04-Mar-2019_2_46.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\TheWall_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_05-Mar-2019_3_22.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\TheWall_HighPowerUnLeveled_CrackedCementBlock_InsideWall_60X60cm_RealWallMeasSeries_07-Mar-2019_0_15.mat');
%   load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\TheWall_HighPowerUnLeveled_CrackedCementBlock_DryWall_60X60cm_RealWallMeasSeries_08-Mar-2019_0_55.mat');


% Measurement positions

 
 
%choose object position; Choose frequency points to use. 
z_offset=.56;

%freq range
freqMin=17.5E9;
freqMax=26.5E9;

pick=logical((data.f>=freqMin).*(data.f<=freqMax)).';
pick(2:2:end)=logical(0);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pick(3:4:end)=logical(0);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
f=data.f(pick);

BW=f(end)-f(1);

X=data.X/1000;
Y=data.Y/1000;

dx=abs(X(1,2)-X(1,1));
%x range
% xMin =-.3%-inf;
% xMax = .3%inf;
xMin =  -inf;
xMax =   inf;


xRange = logical((X(1,:)>=xMin).*(X(1,:)<=xMax)).';

%y range
% yMin = -.3%-inf;
% yMax =  .3%inf;
yMin = -inf;
yMax =  inf;
yRange = logical((Y(:,1)>=yMin).*(Y(:,1)<=yMax));


X=X(yRange,xRange);
Y=Y(yRange,xRange);


%measurements
measurements=data.measurements(yRange,xRange, pick);


%subtract background
%   load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\TheWall_35cmAway_60X60cm_WallMeasSeries_13-Feb-2019_1_1');
 
%correct for probe phase
cal=load('C:\Users\Jonah Gollub\Documents\code\data Analysis\SAR\Calibration Data\extractedHornCalRad2019-2-18.mat');

indxTemp=[];
for ii=1:numel(f)
indxTemp(ii)=find(cal.f==f(ii));
end

measurements=(measurements)./repmat(permute(exp(1j*2*cal.phaseCalData(indxTemp)),[3,2,1]),size(measurements,1),size(measurements,2));

% if f==cal.f(pick)
%  measurements=(measurements)./exp(1j*2*cal.phaseCalData(pick));
%  error('cal frequencies not contained in measurements') 
% end

clear data;
%initilize variables
k =[]; kux = []; kuy =[]; kuz = []; pad = []; Sxy = [];

%% Perform RMA algorithm (functionalized form below)

% upsample in x and y

Lx=X(1,end)-X(1,1);
[ynum,xnum]=size(X);

dy=abs(Y(2,1)-Y(1,1));
Ly=Y(end,1)-Y(1,1);

pad=2^nextpow2(max(numel(X(1,:)),numel(Y(:,1))));
pad=2^8;

dx=Lx/(pad-1);
dy=Ly/(pad-1);

%calc free space k vector
k=2*pi*f/c;
k=repmat(permute(k,[1,3,2]),[pad,pad,1]);

% [kux,kuy,~]=meshgrid(-(2*pi)/(2*dx):2*pi/(Lx):(2*pi)/(2*dx),...
%                      -(2*pi)/(2*dy):2*pi/(Ly):(2*pi)/(2*dy),...
%                      1:numel(f));

[kux,kuy,~]=meshgrid(linspace(-(2*pi)/(2*dx),(2*pi)/(2*dx),pad),...
                     linspace(-(2*pi)/(2*dy),(2*pi)/(2*dy),pad),...
                     1:numel(f));

%calculate plane wave decomposition (FFT). Note we have to be careful about
%shifting in fft because we are only acting along two of the dimensions.
%Therefor it is easier to use circshift (as opposed to ifftshif which acts
%on all dimentions.

Sxy=fft(fftshift(measurements),[],1);
Sxy=ifftshift(fft(Sxy,[],2));
Sxy=padarray(Sxy,[ceil((pad-size(measurements,1))/2), ceil((pad-size(measurements,2))/2)],0,'pre');
Sxy=padarray(Sxy,[floor((pad-size(measurements,1))/2), floor((pad-size(measurements,2))/2)],0,'post');

% Sxy=(1/(Lx_pad*Ly_pad))*circshift(fft2(measurements,pad,pad),[pad/2,pad/2]);
% 
% shiftx=dx*(pad/2-floor(xnum/2));
% shifty=dy*(pad/2-floor(ynum/2));
% 
% Sxy=Sxy.*exp(-1.0j*kux*(shiftx)).*exp(-1.0j*kuy*(shifty)); %we must phase shift such that panel is centered

%     plot decomposition
    figHandle=figure(1);
    scrn = get( groot, 'Screensize');  scrn(1)=2*scrn(3)/3;  scrn(3)=scrn(3)/3;
    set(figHandle,'Position',scrn); clf; subplot(2,2,1);
%     imagesc(-Lx_pad/2:dx:Lx_pad/2,-Ly_pad/2:dy:Ly_pad/2,abs(ifft2(fftshift(Sxy(:,:,1)))));
    imagesc(-Lx/2:dx:Lx/2,-Ly/2:dy:Ly/2,abs(ifftshift(ifft2(fftshift(Sxy(:,:,1))))));
    title('Padded Input Fields f(1)'); axis equal; axis tight; xlabel('x (m)'); ylabel('y (m)')
    subplot(2,2,2);
    imagesc(kux(1,:),kuy(:,1),abs(Sxy(:,:,1)));
    title('FFT of Fields f(1)'); axis equal; axis tight; xlabel('kx (m)'); ylabel('ky (m)')
    drawnow
    
%calculate min max kz wavenumber

kuz=sqrt((2*k).^2-kux.^2-kuy.^2);
kuz=real(kuz); % ignore evanescent fields

kuzMin = min(kuz(:));
kuzMax = max(kuz(:));

%ensure sufficient sampling of Kz 
minSampling = min(nonzeros(diff(kuz,1,3))); %only look end of matrix (high frequency region)
numSampling = max(size(measurements,3),ceil((kuzMax-kuzMin)/minSampling)); %opt 1 (sample at minimum)
numSampling = 2^nextpow2(numSampling); %opt 2 (sample at next power of 2)


Kz=linspace(kuzMin,kuzMax,numSampling);

Sxy=Sxy.*exp(1.0j*(kuz)*z_offset);

%interpolate to evenly spaced grid

Srmg = zeros(size(kuz,1),size(kuz,2),length(Kz));
for ii=1:size(kux,2)
    for jj=1:size(kuy,1)
        indx_vec = squeeze(squeeze(real(kuz(jj,ii,:))~=0));
        if ~(sum(indx_vec)==0 || sum(indx_vec)==1) %check that there are enough points to interpolate
           Srmg(jj,ii,:)=interp1(squeeze(squeeze(kuz(jj,ii,indx_vec))),...
                                 squeeze(Sxy(jj,ii,indx_vec)),...
                                 Kz.','linear');
        end
    end
end
% %debug: plot interpolatio for z-slice
%             ii=44;
%             jj=44;
%             %debug
%             figure(3); cla; 
%             plot(squeeze(squeeze(kuz(jj,ii,indx_vec))),squeeze(Sxy(jj,ii,indx_vec)))
%             hold on;            
%             Srmg(jj,ii,find(isnan(Srmg(jj,ii,:))))=0;
%             plot(real(Kz(:)).',squeeze(squeeze(Srmg(jj,ii,:))),'-o')
%             drawnow; 

Srmg(find(isnan(Srmg))) = 0; %set all Nan values to 0

% k_upsampled=linspace(min(k(:)), max(k(:)),4*size(measurements,3));
% k_upsampled=repmat(permute(k_upsampled,[1,3,2]),[pad,pad,1]);

%   Srmg=Srmg.*exp(1.0*j*(kuz)*z_offset);

%apply inverst FFT to get image
  fxy = fftshift(ifftn(Srmg));
% fxy = ifftn(ifftshift(Srmg));
image=(abs(fxy)/max(abs(fxy(:))));

%% Functional Form of RMA code options variation 1

[image,Lx,Ly,numSampling]=reconstructRF(X,Y,f,measurements,'z_offset',z_offset,'Algorithm','RMA');

%% Functional Form of RMA code options variation 2 (using GPU, but doesn't speed things up currently)

% [image,Lx,Ly,numSampling]=reconstructRF(X,Y,f,measurements,'z_offset',z_offset,'Algorithm','RMA_GPU');

%% Movie of images varying offset parameter
clear F;
figure(101); clf;

z_vec = 0.46:.1:.76

video_filename='z_offset';
path = 'D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\';

for qq=1:length(z_vec)
    z_offset = z_vec(qq)
    
    [image,Lx,Ly,numSampling]=reconstructRF(X,Y,f,measurements,'z_offset',z_offset,'Algorithm','RMA');
    
    
    % plot in loop
    
    %labeling
    xx=linspace(-Lx/2,Lx/2,size(image,2));
    yy=linspace(-Ly/2,Ly/2,size(image,1));
    
    zz=linspace(-c*numel(f)/(2*BW)/2,c*numel(f)/(2*BW)/2,numSampling);
    %  zz=linspace(0,c*numel(f)/(2*BW),numel(f));
    
    % plot subimage
    xmin=min(X(1,:));
    xmax=max(X(1,:));
    ymin=min(Y(:,1));
    ymax=max(Y(:,1));
    zmin=-.2;
    zmax=0.5;
    
    subimage=image(yy>ymin & yy<ymax, xx>xmin & xx<xmax, zz>zmin & zz<zmax);
    subimage=subimage/max(subimage(:));
    figure(101); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',subimage.^1.5,...
        'xdata',xx(xx>xmin & xx<xmax),...
        'Ydata',yy(yy>ymin & yy<ymax),...
        'Zdata',zz(zz>zmin & zz<zmax)...
        );
    zlabel('downrange (m)'); xlabel(' x crossrange (m)');  ylabel('y crossrange (m)');
    axis equal; axis tight; 
    
    %view angle
    view(180,-90);% view(3);
    
    title(['Offset distance = ',num2str(z_offset,'%4.4f')])
    
    F(qq) = getframe(gcf);
    drawnow;
    
end
writerObj = VideoWriter([path,video_filename]);
writerObj.FrameRate = 4;

% open the video writer
open(writerObj);

% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

%% Perform RMA algorithm

% upsample in x and y

Lx=X(1,end)-X(1,1);
[ynum,xnum]=size(X);

dy=abs(Y(2,1)-Y(1,1));
Ly=Y(end,1)-Y(1,1);

pad=2^nextpow2(max(numel(X(1,:)),numel(Y(:,1))));
pad=2^8;

dx=Lx/(pad-1);
dy=Ly/(pad-1);

%calc free space k vector
k=2*pi*f/c;
k=repmat(permute(k,[1,3,2]),[pad,pad,1]);

% [kux,kuy,~]=meshgrid(-(2*pi)/(2*dx):2*pi/(Lx):(2*pi)/(2*dx),...
%                      -(2*pi)/(2*dy):2*pi/(Ly):(2*pi)/(2*dy),...
%                      1:numel(f));

[kux,kuy,~]=meshgrid(linspace(-(2*pi)/(2*dx),(2*pi)/(2*dx),pad),...
                     linspace(-(2*pi)/(2*dy),(2*pi)/(2*dy),pad),...
                     1:numel(f));

%calculate plane wave decomposition (FFT). Note we have to be careful about
%shifting in fft because we are only acting along two of the dimensions.
%Therefor it is easier to use circshift (as opposed to ifftshif which acts
%on all dimentions.

Sxy=fft(fftshift(measurements),[],1);
Sxy=ifftshift(fft(Sxy,[],2));
Sxy=padarray(Sxy,[ceil((pad-size(measurements,1))/2), ceil((pad-size(measurements,2))/2)],0,'pre');
Sxy=padarray(Sxy,[floor((pad-size(measurements,1))/2), floor((pad-size(measurements,2))/2)],0,'post');

% Sxy=(1/(Lx_pad*Ly_pad))*circshift(fft2(measurements,pad,pad),[pad/2,pad/2]);
% 
% shiftx=dx*(pad/2-floor(xnum/2));
% shifty=dy*(pad/2-floor(ynum/2));
% 
% Sxy=Sxy.*exp(-1.0j*kux*(shiftx)).*exp(-1.0j*kuy*(shifty)); %we must phase shift such that panel is centered

%     plot decomposition
    figHandle=figure(1);
    scrn = get( groot, 'Screensize');  scrn(1)=2*scrn(3)/3;  scrn(3)=scrn(3)/3;
    set(figHandle,'Position',scrn); clf; subplot(2,2,1);
%     imagesc(-Lx_pad/2:dx:Lx_pad/2,-Ly_pad/2:dy:Ly_pad/2,abs(ifft2(fftshift(Sxy(:,:,1)))));
    imagesc(-Lx/2:dx:Lx/2,-Ly/2:dy:Ly/2,abs(ifftshift(ifft2(fftshift(Sxy(:,:,1))))));
    title('Padded Input Fields f(1)'); axis equal; axis tight; xlabel('x (m)'); ylabel('y (m)')
    subplot(2,2,2);
    imagesc(kux(1,:),kuy(:,1),abs(Sxy(:,:,1)));
    title('FFT of Fields f(1)'); axis equal; axis tight; xlabel('kx (m)'); ylabel('ky (m)')
    drawnow
    
%calculate min max kz wavenumber

kuz=sqrt((2*k).^2-kux.^2-kuy.^2);
kuz=real(kuz); % ignore evanescent fields

kuzMin = min(kuz(:));
kuzMax = max(kuz(:));

%ensure sufficient sampling of Kz 
minSampling = min(nonzeros(diff(kuz,1,3))); %only look end of matrix (high frequency region)
numSampling = max(size(measurements,3),ceil((kuzMax-kuzMin)/minSampling)); %opt 1 (sample at minimum)
numSampling = 2^nextpow2(numSampling); %opt 2 (sample at next power of 2)


Kz=linspace(kuzMin,kuzMax,numSampling);

Sxy=Sxy.*exp(1.0j*(kuz)*z_offset);

%interpolate to evenly spaced grid

Srmg = zeros(size(kuz,1),size(kuz,2),length(Kz));
for ii=1:size(kux,2)
    for jj=1:size(kuy,1)
        indx_vec = squeeze(squeeze(real(kuz(jj,ii,:))~=0));
        if ~(sum(indx_vec)==0 || sum(indx_vec)==1) %check that there are enough points to interpolate
           Srmg(jj,ii,:)=interp1(squeeze(squeeze(kuz(jj,ii,indx_vec))),...
                                 squeeze(Sxy(jj,ii,indx_vec)),...
                                 Kz.','linear');
        end
    end
end
% %debug: plot interpolatio for z-slice
%             ii=44;
%             jj=44;
%             %debug
%             figure(3); cla; 
%             plot(squeeze(squeeze(kuz(jj,ii,indx_vec))),squeeze(Sxy(jj,ii,indx_vec)))
%             hold on;            
%             Srmg(jj,ii,find(isnan(Srmg(jj,ii,:))))=0;
%             plot(real(Kz(:)).',squeeze(squeeze(Srmg(jj,ii,:))),'-o')
%             drawnow; 

Srmg(find(isnan(Srmg))) = 0; %set all Nan values to 0

% k_upsampled=linspace(min(k(:)), max(k(:)),4*size(measurements,3));
% k_upsampled=repmat(permute(k_upsampled,[1,3,2]),[pad,pad,1]);

%   Srmg=Srmg.*exp(1.0*j*(kuz)*z_offset);

%apply inverst FFT to get image
  fxy = fftshift(ifftn(Srmg));
% fxy = ifftn(ifftshift(Srmg));
image=(abs(fxy)/max(abs(fxy(:))));
 
%% plot 

%labeling
xx=linspace(-Lx/2,Lx/2,size(image,2));
yy=linspace(-Ly/2,Ly/2,size(image,1));

 zz=linspace(-c*numel(f)/(2*BW)/2,c*numel(f)/(2*BW)/2,numSampling);
%  zz=linspace(0,c*numel(f)/(2*BW),numel(f));

    figure(1); cla; subplot(2,2,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',image.^1.5,'xdata',xx,'Ydata',yy,'Zdata',zz);
    axis equal; axis tight; view(3);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');

     drawnow;
 
%% plot subimage

xmin=min(X(1,:));
xmax=max(X(1,:));
ymin=min(Y(:,1));
ymax=max(Y(:,1));
zmin=-0.15;
zmax=.15;

subimage=image(yy>ymin & yy<ymax, xx>xmin & xx<xmax, zz>zmin & zz<zmax);
 subimage=subimage/max(subimage(:));
      figure(102); %subplot(2,1,2); 
    cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',subimage.^1.5,...
        'xdata',xx(xx>xmin & xx<xmax),...
        'Ydata',yy(yy>ymin & yy<ymax),...
        'Zdata',zz(zz>zmin & zz<zmax)...
        );
    zlabel('downrange (m)'); xlabel(' x crossrange (m)');  ylabel('y crossrange (m)');
    axis equal; axis tight; view(180,-90); %view(180,0)

%% Reconstruct g=Hf backprojection
%reconstruction zone

xvec = -0.3:0.0075:0.3;
yvec = -0.3:0.0075:0.3;
zvec = 0.45:0.015:0.7;

[xi,yi,zi] = meshgrid(xvec, yvec, zvec);

kr=(2*pi*f/c).';
%generate H !!! too large keep in memory
% H=zeros(length(kr)*numel(X),numel(xi(:)), 'single'); 

% for ii=1:numel(X)
% 
%         Di = sqrt((X(ii)-xi(:)).^2+(Y(ii)-yi(:)).^2 + zi(:).^2).';
%         
%         %measurement matrix; note di is n x 1 and k is 1 x n vector
%         H(1+(ii-1)*length(kr):ii*length(kr),:) =...
%           single(((1./Di).*exp(-1j*kr.*Di)).*((1./Di).*exp(-1j*kr.*Di)));
% 
%       if mod(ii,1)==0;
%           ii
%       end
% end 


fest=zeros(numel(xi),1);
size(fest)

for ii=1:numel(xi)

        Di = sqrt((X(:)-xi(ii)).^2+(Y(:)-yi(ii)).^2 + zi(ii).^2).';
        
        %measurement matrix; note di is 1 x n and k is n x 1 vector
        %calulate fest one row at a time (to conserve memory)
         fii= single(((1./Di).*exp(-1j*kr.*Di)).*((1./Di).*exp(-1j*kr.*Di)));

          fest(ii) = fii(:)'*reshape(permute(measurements,[3, 1, 2]),[],1);
      if mod(ii,100)==0;
          ii
      end
end 
clear f_image;
f_image=reshape(fest,length(yvec),length(xvec),length(zvec));
f_image=abs(f_image)/max(abs(f_image(:)));

%% plot g=HF

%labeling
xx=xvec;
yy=yvec;
zz=zvec;
%  zz=linspace(0,c*numel(f)/(2*BW),numel(f));

    figure(1); cla; subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',f_image.^1.5,'Xdata',xx,'Ydata',yy,'Zdata',zz);
    axis equal; axis tight; view(180,-90);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');

     drawnow;
