%-------------------------------------------------------------------------------
%SAR Range Migration Algorithm 
%J. Gollub
%-------------------------------------------------------------------------------


%% load data and set parameters
c=2.99792458*10^8;

% data
% load('C:\Users\Jonah Gollub\Downloads\drill_data.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\D_38cmPastWall_60X60cm_WallMeasSeries_14-Feb-2019_3_32.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Saw_46cmAway_60X60cm_CorrectAlign_09-Feb-2019_4_59.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\ReferenceFrame_382mmPastWall_60X60cm_WallMeasSeries_19-Feb-2019_3_21.mat');

% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\2By4_InsideWall_60X60cm_WallMeasSeries_21-Feb-2019_2_56.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Baseline_InsideWall_60X60cm_RealWallMeasSeries_23-Feb-2019_1_29.mat');
  load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Plywood_InsideWall_60X60cm_RealWallMeasSeries_24-Feb-2019_2_56.mat');

% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Drywall_InsideWall_60X60cm_RealWallMeasSeries_25-Feb-2019_19_40.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_26-Feb-2019_6_13.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\RandomObjects_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_26-Feb-2019_19_39.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\Fiducial_wDryWall_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_28-Feb-2019_19_20.mat');
  
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\ShiftedRight_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_27-Feb-2019_4_48.mat');
% load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\ShiftedRight_wDryWall_HighPowerUnLeveled_Drywall_InsideWall_60X60cm_RealWallMeasSeries_27-Feb-2019_21_13.mat');

%choose object position; Choose frequency points to use. 
z_offset=.56;
pick=logical((data.f>17E9).*(data.f<27E9));

measurements=data.measurements(:,:, pick);
f=data.f(pick);
BW=f(end)-f(1);

% Measurement positions
X=data.X/1000; 
dx=abs(X(1,2)-X(1,1));
Lx=X(1,end)-X(1,1);
[ynum,xnum]=size(X);

Y=data.Y/1000;
dy=abs(Y(2,1)-Y(1,1));
Ly=Y(end,1)-Y(1,1);

clear data;
%subtract background
%   load('D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\TheWall_35cmAway_60X60cm_WallMeasSeries_13-Feb-2019_1_1');
 
%correct for probe phase
cal=load('C:\Users\Jonah Gollub\Documents\code\data Analysis\SAR\Calibration Data\extractedHornCalRad2019-2-18.mat');
if f==cal.f(pick)
 measurements=(measurements)./exp(1j*2*cal.phaseCalData(pick));
 error('cal frequencies not contained in measurements') 
end

%initilize variables
k =[]; kux = []; kuy =[]; kuz = []; pad = []; Sxy = [];
%% Perform RMA algorithm

% clear F;
% figure(101); clf;
% 
%  z_vec = -.6:.05:.3
% 
% video_filename='z_offset';
% path = 'D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\';
% 
% for qq=1:length(z_vec) 
% z_offset = z_vec(qq)    
% 
% kmeas=2*pi*f/c;
% kmeas=repmat(permute(kmeas,[1,3,2]),[size(measurements,1),size(measurements,2),1]);
% 
% measurementsShifted=measurements./exp(-1.0j*(kmeas/.7).*cal_offset);
% measurementsShifted=measurements./exp(-1.0j*(kmeas/.7).*(-0.45));

% upsample in x and y

pad=2^nextpow2(max(numel(X(1,:)),numel(Y(:,1))));
pad=2^7;

Lx_pad=dx*(pad-1); 
Ly_pad=dy*(pad-1);

%calc free space k vector
k=2*pi*f/c;
k=repmat(permute(k,[1,3,2]),[pad,pad,1]);

[kux,kuy,~]=meshgrid(-(2*pi)/(2*dx):2*pi/(Lx_pad):(2*pi)/(2*dx),...
                     -(2*pi)/(2*dy):2*pi/(Ly_pad):(2*pi)/(2*dy),...
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
    imagesc(-Lx_pad/2:dx:Lx_pad/2,-Ly_pad/2:dy:Ly_pad/2,abs(ifftshift(ifft2(fftshift(Sxy(:,:,1))))));
    title('Padded Input Fields f(1)'); axis equal; axis tight; xlabel('x (m)'); ylabel('y (m)')
    subplot(2,2,2);
    imagesc(kux(1,:),kuy(:,1),abs(Sxy(:,:,1)));
    title('FFT of Fields f(1)'); axis equal; axis tight; xlabel('kx (m)'); ylabel('ky (m)')
    drawnow
    
%calculate min max kz wavenumber

kuz=sqrt((2*k).^2-kux.^2-kuy.^2);
kuz=real(kuz); %!!!!!!!!!!!!!!!!!!! ignore evanescent fields

Kz=linspace(min(kuz(:)),max(kuz(:)),4*size(measurements,3));

% dt=c/(BW);
% Kz=min(kuz(:)):4*pi/(c*dt):max(kuz(:));


Sxy=Sxy.*exp(1.0j*(kuz)*z_offset);

%interpolate to evenly spaced grid

Srmg = zeros(size(kuz,1),size(kuz,2),length(Kz));
for ii=1:size(kux,2)
    for jj=1:size(kuy,1)
     %           kuz(jj,ii,real(kuz(jj,ii,:))==0)=0;
        indx_vec = squeeze(squeeze(real(kuz(jj,ii,:))~=0));
        if ~sum(indx_vec)==0
 
           Srmg(jj,ii,:)=interp1(squeeze(squeeze(kuz(jj,ii,indx_vec))),...
                squeeze(Sxy(jj,ii,indx_vec)),...
                Kz.','linear');
            
%             %debug
%             figure(3); cla; 
%             plot(squeeze(squeeze(kuz(jj,ii,indx_vec))),squeeze(Sxy(jj,ii,indx_vec)))
%             hold on;            
%             Srmg(jj,ii,find(isnan(Srmg(jj,ii,:))))=0;
%             plot(real(Kz(:)).',squeeze(squeeze(Srmg(jj,ii,:))),'-o')
%             drawnow; pause;
        end
    end
end

Srmg(find(isnan(Srmg))) = 0; %set all Nan values to 0

% k_upsampled=linspace(min(k(:)), max(k(:)),4*size(measurements,3));
% k_upsampled=repmat(permute(k_upsampled,[1,3,2]),[pad,pad,1]);

%   Srmg=Srmg.*exp(1.0*j*(kuz)*z_offset);

%apply inverst FFT to get image
  fxy = fftshift(ifftn(Srmg));
% fxy = ifftn(ifftshift(Srmg));
image=(abs(fxy)/max(abs(fxy(:))));

% plot in loop
% 
% %labeling
% xx=linspace(-Lx_pad/2,Lx_pad/2,size(fxy,2));
% yy=linspace(-Ly_pad/2,Ly_pad/2,size(fxy,1));
% 
%  zz=linspace(-c*numel(f)/(2*BW)/2,c*numel(f)/(2*BW)/2,numel(f));
% % zz=linspace(0,c*numel(f)/(BW),numel(f));
% 
% % plot subimage
% xmin=min(X(1,:));
% xmax=max(X(1,:));
% ymin=min(Y(:,1));
% ymax=max(Y(:,1));
% zmin=-3.5;
% zmax=3.5;
% 
% subimage=image(yy>ymin & yy<ymax, xx>xmin & xx<xmax, zz>zmin & zz<zmax);
% subimage=subimage/max(subimage(:));
%    figure(101); cla;
%     set(gcf,'color','white'); colormap('parula');
%     vol3d('Cdata',subimage.^2,...
%         'xdata',xx(xx>xmin & xx<xmax),...
%         'Ydata',yy(yy>ymin & yy<ymax),...
%         'Zdata',zz(zz>zmin & zz<zmax)...
%         );
%     zlabel('downrange (m)'); xlabel(' x crossrange (m)');  ylabel('y crossrange (m)');
%     axis equal; axis tight; view(3);
%      
%      title(['Offset distance = ',num2str(z_offset,'%4.4f')])
% 
%     F(qq) = getframe(gcf);
%      drawnow;
% 
% end
%   writerObj = VideoWriter([path,video_filename]);
%   writerObj.FrameRate = 4;
% 
%   % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);
  

%% plot 

%labeling
xx=linspace(-Lx_pad/2,Lx_pad/2,size(fxy,2));
yy=linspace(-Ly_pad/2,Ly_pad/2,size(fxy,1));

 zz=linspace(-c*numel(f)/(2*BW)/2,c*numel(f)/(2*BW)/2,numel(f)*4);
%  zz=linspace(0,c*numel(f)/(2*BW),numel(f));

    figure(1); cla; subplot(2,2,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',image.^1,'xdata',xx,'Ydata',yy,'Zdata',zz);
    axis equal; axis tight; view(3);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');

     drawnow;
 
%% plot subimage
xmin=min(X(1,:));
xmax=max(X(1,:));
ymin=min(Y(:,1));
ymax=max(Y(:,1));
zmin=0.05;
zmax=.25;

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

%% move along path plot 
% clear F;
% figure(101); clf;
% 
% video_filename='move_through_wall';
% path = 'D:\Dropbox (Duke Electric & Comp)\Through Wall\9GHz-401pts\';
% 
% z_vec=linspace(0,0.15,24)
% 
% for qq=1:length(z_vec) 
%  z_vec(qq) 
% 
% 
% xmin=min(X(1,:));
% xmax=max(X(1,:));
% ymin=min(Y(:,1));
% ymax=max(Y(:,1));
% zmin=z_vec(qq) + 0;
% zmax=z_vec(qq) + 0.25;
% 
% subimage=image(yy>ymin & yy<ymax, xx>xmin & xx<xmax, zz>zmin & zz<zmax);
%  subimage=subimage/max(subimage(:));
%  
%      figure(101); 
%     cla;
%     set(gcf,'color','white'); colormap('parula');
%     vol3d('Cdata',subimage.^1.5,...
%         'xdata',xx(xx>xmin & xx<xmax),...
%         'Ydata',yy(yy>ymin & yy<ymax),...
%         'Zdata',zz(zz>zmin & zz<zmax)...
%         );
%     zlabel('downrange (m)'); xlabel(' x crossrange (m)');  ylabel('y crossrange (m)');
%     axis equal; axis tight; view(180,-90);
% 
%      title(['Offset distance = ',num2str(z_vec(qq),'%4.4f')])
% 
%     F(qq) = getframe(gcf);
%      drawnow;
% 
% end
%   writerObj = VideoWriter([path,video_filename]);
%   writerObj.FrameRate = 4;
% 
%   % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
l% end
% % close the writer object
% close(writerObj);

