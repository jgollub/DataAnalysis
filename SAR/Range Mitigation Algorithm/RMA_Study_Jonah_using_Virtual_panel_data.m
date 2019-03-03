%% parameters for panel data taken from NFS
c=3e8;
% 
load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Near Field Scans\RevA4\6_feeds_NFS-29-Jan-2015.mat')
A1data=data;
clear data;
freqnum=101;
f=A1data.f;

y_apt=-A1data.X*1e-3;  %neg sign for import
z_apt=A1data.Y*1e-3;

yn_apt=numel(y_apt(1,:));
zn_apt=numel(z_apt(:,1));

ymin=y_apt(1,1);
ymax=y_apt(1,end);
zmin=z_apt(1,1);
zmax=z_apt(end,1);

dy=abs(y_apt(1,2)-y_apt(1,1));
dz=(z_apt(2,1)-z_apt(1,1));

SizeY=.56;
SizeZ=.4;
Q=100;
feeds=6;

% %Target Sphere
% x_offset=1;
% y_offset=0;
% z_offset=0;
 probe_offset_y=0;

%% target block
x_offset=1;

y_offset=0;
y_range=.05;

z_offset=0;
z_range=.05;

width=.02;
target_res=.01;
%reconstruction Plane
x_offset_recon=x_offset;



%% parameter data for virtualizer panel and probes
SizeY=.2;
SizeZ=.2;
Q=100;
feeds=1;
%f= linspace(17.5e9, 26.5e9, 101);
f=A1data.f;
dy=.005;
dz=.005;


panel = create_panel('fsweep',f,'SizeY',SizeY,'SizeZ',SizeZ,'ElementSizeY', ...
    dy, 'ElementSizeZ', dz,'Q',Q);
panel.u=[0; 1; 0]; 
panel.v=[0; 0; 1];
panel =panel_feed(panel,'type', 'circular','feedLocs',feeds);

probe = create_panel('type', 'probe', 'SizeY', .012, 'SizeZ', .012,...
    'ElementSizeY', .002, 'ElementSizeZ', .002, 'fsweep', f);


y_apt=panel.y; 
z_apt=panel.z;

yn_apt=numel(y_apt(1,:));
zn_apt=numel(z_apt(:,1));

ymin=y_apt(1,1);
ymax=y_apt(1,end);

zmin=z_apt(1,1);
zmax=z_apt(end,1);

dky=2*pi/(ymax-ymin);
dkz=2*pi/(zmax-zmin);

%% Resolution Target

target=BlockTarget(x_offset, 0, width, width, target_res/6, width);
% target.xmat(~any(target.shape,2),:)=[];
% target.xmat(:,~any(target.shape,1))=[];

[locs, coord_grid] = test_space2(x_offset_recon,width*6, y_offset,...
    6*width,z_offset, 6*width, target_res);

sigma = ones(size(target.locs(:,1)));

try figure(scene_plot)
catch
    scene_plot=figure();
    scrsize=get(0,'ScreenSize');
    figsize=scene_plot.Position(3:4);
    scene_plot.Position(1:2)=[scrsize(3)-figsize(1),scrsize(4)-figsize(2)-80];      
   
end
%plot region of reconstruction
clf
panel_plot(scene_plot,probe,panel);

hold on;
%scatter3(locs.locs(:,1),locs.locs(:,2),locs.locs(:,3),1);

%plot target structure
hold on;
scatter3(target.locs(:,1),target.locs(:,2),target.locs(:,3),10,'red');


%% First calculate fields on the aperture directly

 %Directly calculated fourier componentS at each position on the aperture for 
 %all scattering target points

 k=2*pi*f/c;
 
 R_nx=@(x,y,z,xt,yt,zt) sqrt((xt-x).^2+(yt-y).^2+(zt-z).^2);
 
 S_xyk=zeros([size(y_apt),numel(k)]);
 S_kxkyk=0*S_xyk;
 tracker=0;
 total_steps=yn_apt*zn_apt*numel(k);

 %probe
 R_Rx=R_nx(0, 0, 0, target.xmat,target.ymat,target.zmat);
 tic
 for iy=1:yn_apt
     for iz=1:zn_apt
         for jk=1:numel(k)
             % for each point of the panel aperture, calculate the fields due to
             % all of the target points
             R_Tx=R_nx(0, y_apt(1,iy), z_apt(iz,1), target.xmat,target.ymat,target.zmat);
             
             S_xyk(iz,iy,jk)=sum(sum(sum((1/(R_Tx.*R_Rx)).*exp(-1j*k(jk)*(R_Tx+R_Rx)))));
             
             tracker=tracker+1;
             if mod(tracker,floor(total_steps/10))==0
                 fprintf('Calculating S_xyk is %3.1f %% done.\n', (tracker/total_steps)*100)
             end
         end
     end
     
 end
 
 tracker=0;
 for jk=1:numel(k)
     S_kxkyk(:,:,jk)=ifftshift(fft2(fftshift(squeeze(S_xyk(:,:,jk)))));
     tracker=tracker+1;
     if mod(tracker,floor(numel(k)/5))==0
         fprintf('Calculating FFT of S_kxkyk is %3.1f %% done.\n', (tracker/numel(k))*100)
     end
 end
 fprintf('Calc time is %f s',toc) 
 figure(98)
imagesc(abs(squeeze(S_xyk(:,:,1))));
 
 figure(99)
imagesc(abs(squeeze(S_kxkyk(:,:,1))))
axis equal; axis tight;
 
 %% accepting fields on the aperture
 
 
 %% K space decomposition



% transform NFS data to wavenumber space
FTpanelA1=[];
PSFspace=[];
PSFfreq=[];
for ifeed=1:feeds
for ifreq=1:numel(f)
%FTpanelA1(:,:,ifreq,ifeed)=ifftshift(fft2(fftshift(A1data.measurements(:,:,ifreq,ifeed))))*dy*dz;
%fix here need to find fields from panels
FTpanelA1(:,:,ifreq,ifeed)=ifftshift(fft2(fftshift(A1data.measurements(:,:,ifreq,ifeed))))*dy*dz;

% %apply phase conjugation
%  PSFspace(:,:,ifreq,ifeed)=conj(FTpanelA1(:,:,ifreq,ifeed)).*FTpanelA1(:,:,ifreq,ifeed);
% %PSFspace(:,:,ifreq,ifeed)=FTpanelA1(:,:,ifreq,ifeed).';


%PSFfreq(:,:,ifreq,ifeed)=ifftshift(ifft2(fftshift(PSFspace(:,:,ifreq,ifeed))));

end
end
figure(100)
imagesc(abs(FTpanelA1(:,:,51,1)))
axis equal; axis tight;

% figure(101)
% imagesc(abs(PSFfreq(:,:,51,1)))
% axis equal; axis tight;
 
 %% then the field received at the aperture is 
 E_k_apt=FTpanelA1(:,:,:,1).*S_kxkyk;
 
 figure(102)
imagesc(abs(E_k_apt(:,:,51)))
axis equal; axis tight;
 

%% interpolation 
dky=2*pi/abs(ymax-ymin);
dkz=2*pi/abs(zmax-zmin);
% [Ky, Kz, K]=meshgrid(-2*pi/(2*dy):dky:2*pi/(2*dy),-2*pi/(2*dz):dkz:2*pi/(2*dz),k);
ky=linspace(-2*pi/(2*dy),2*pi/(2*dy),yn_apt);
kz=linspace(-2*pi/(2*dz),2*pi/(2*dz),zn_apt);
[Kz,Ky,K] = ndgrid(kz,ky,k);

Kx=sqrt(K.^2-Ky.^2-Kz.^2)+K;

% evanescent fields decay and do not need to be included. Hence if imag(Kx)
% must be real.

% [Zreal,Yreal,Kreal]=find(imag(Kx)==0);

Kx_flat = squeeze(sum(Kx,3));
[Zr,Yr]=find(imag(Kx_flat)==0);

% figure()
% imagesc(imag(Kx_flat)==0)
% vol3d('cdata',imag(Kx_flat)==0);
% axis equal

%Interpolation of k's into new k-space
Kx(imag(Kx)~=0)=0;
Kx_min=min(min(min((Ky(Ky~=0)))));
Kx_max=max(max(max(Ky)));


% kx_ref = squeeze(Kx(round(numel(zn_apt)/2),round(numel(yn_apt)/2),:));
kx_ref = linspace(Kx_min,Kx_max,numel(k));

x0 = 1; 

%central depth... add in noise here if need be
E_comp=E_k_apt.*conj(squeeze(FTpanelA1(:,:,:,1)));%.* exp(1j*real(Kx).*x0);
E_interp=E_comp.*0;

for m = 1:numel(Zr)
    E_interp(Zr(m),Yr(m),:) = interp1(...
        squeeze(Kx(Zr(m),Yr(m),:)).',abs(squeeze(E_comp(Zr(m),Yr(m),:))).',kx_ref,'nearest');
    
    display(['Interpolation ... ',num2str(round(m/numel(Zr)*100)),'%'])
end

% for m=1:numel(Zreal)
%     E_interp(Zreal(m),Yreal(m),Kreal)=interp1(...
%         squeeze(Kx(Zreal(m),Yreal(m),Kreal)).',squeeze(E_comp(Zreal(m),Yreal(m),Kreal)).',kx_ref,'nearest');
%     
%     E_interp=interp3(Zreal,Yreal,Kreal, ...
%         squeeze(Ky(Zreal,Yreal,Kreal),squeeze(E_comp(Zreal(m),Yreal(m),Kreal(m))).',ky_ref,'nearest');
%     
%     
%     display(['Interpolation ... ',num2str(round(m/numel(Xreal)*100)),'%'])
% end
%% calculate scene using inverse FFT

E_interp(isnan(E_interp)) = 0;

E = fftshift(fftshift(fftshift(...
        ifftn(ifftshift(ifftshift(ifftshift(...
            E_interp...
        ,1),2),3))...
    ,1),2),3);

display(['Decoding and imaging time : ',num2str(toc),' s'])

dkx = mean(diff(kx_ref));
x = linspace(-pi/dkx,pi/dkx,numel(kx_ref)) + x0;
figure(scene_plot)
Elog = 20*log10(abs(E) ./ max(max(max(abs(E)))));
 [X,Y,Z] = meshgrid(x,y_apt(1,:),z_apt(:,1));
p = patch(isosurface(X,Y,Z,permute(Elog,[2 3 1]),-3));
set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',1);
p = patch(isosurface(X,Y,Z,permute(Elog,[2 3 1]),-10));
set(p,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
camlight left
lighting gouraud
box on
grid on
 axis([0, 1.25, -y_apt(1,1), -y_apt(1,end), z_apt(1,1), z_apt(end,1)])
view(3); 
daspect([1 1 1]);


%% calculate fields
g=[];
%take measurement
g = forward_model(probe, panel, sigma, target.locs);

g=permute(g, [ 3 1 2]);
g=g(:);

    % Reconstruction using Matched Filter
    %   f_est=image_recon(probe, panel, g, locs.locs, 'matched_filter');
    %Calc expected fields
    fields_probe=dipoles_to_fieldsEXP3(probe,locs.locs);
    fields_panel=dipoles_to_fieldsEXP3(panel,locs.locs);
    
    H=makeH_faceted(fields_probe,fields_panel);
    f_est=H'*g;
    f_est_plot=reshape(f_est,[locs.ny,locs.nx,locs.nz])
% 3D plotting
figure(1)    
cla
    upsamp_3D_plot=[];
    upsample_3D=1;
    [f_est_plot_3D] = upsample_image_3D(abs(f_est_plot), upsample_3D)/max(max(max(abs(f_est_plot))));
    
    f_est_norm=(f_est_plot_3D-min(f_est_plot_3D(:)))./(max(f_est_plot_3D(:))-min(f_est_plot_3D(:)));
    
    vol3d('cdata',f_est_plot_3D,...
        'XData',locs.xrange,...
        'YData',locs.yrange,...
        'ZData',locs.zrange,...
        'Alpha',f_est_norm.^3 ...
        );
    
    axis equal;
    axis tight;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(-48,26)
    title(['Reconstruction from Data']);
    
    hold on;
    %plot target pts
    scatter3(target.locs(:,1),target.locs(:,2),target.locs(:,3),10,'red');
axis equal
    

%%



s_rec=[];
%apply g

for ifeed=1:feeds
    for ifreq=1:numel(f)
        s_rec(:,:,ifreq,ifeed)=PSFspace(:,:,ifreq,ifeed)*g(ifreq+(ifeed-1)*numel(f));
    end
end

k_yR=[];

for ifreq=1:numel(f)
k_yR(:,:,ifreq)=sqrt(2*pi*f(ifreq)/c-ky_R.^2-kz_R.^2);
end


yval=1;

E_kspace=[];
for ifeed=1:feeds
for ifreq=1:numel(f)
E_kspace(:,:,ifreq,ifeed)=s_rec(:,:,ifreq,ifeed).*exp(1j*(k_yR(:,:,ifreq)+2*pi*f(ifreq)/c)*yval);
end
end
E_kspace=sum(E_kspace,3);
E_kspace=sum(E_kspace,4);
E_xyz=ifftshift(ifft2(fftshift(E_kspace)));

figure(1)
imagesc(ky_R(1,:),kz_R(:,1), abs(E_kspace))
axis equal 
axis tight

figure(2)
imagesc(ky_R(1,:),kz_R(:,1), abs(E_xyz))
axis equal 
axis tight


figure()
imagesc((A1data.X(1,:)*1e-3),(A1data.Y(:,1)*1e-3), abs(A1data.measurements(:,:,freqnum,1)))
axis equal 
axis tight

%% apply measurement, i.e. MW(k)*g(k)

%% Add phase influece of scene

%% FT-1 to get scene by slice.

