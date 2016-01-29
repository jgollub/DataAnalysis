data=load('D:\Dropbox (Duke Electric & Comp)\WPT (1) (1)\Printed Cavity\Hughes_Latest\Tx_Panel_101_both_pol.mat')
figure(1); clf;
f_num=1;
c=3e8;
measurements=flip(flip(data.measurements,1),2);
X=flip(data.X,2);
Y=flip(data.Y,1);
imagesc(X(1,:),Y(:,1),abs(measurements(:,:,f_num,1)));
axis xy; axis equal; axis tight;

pad=2^nextpow2(max(numel(X(1,:)),numel(Y(:,1))));

dx=(X(1,2)-X(1,1))/1e3;
dy=(Y(2,1)-Y(1,1))/1e3;
Lx=(X(1,end)-X(1,1))/1e3; Ly=(Y(end,1)-Y(1,1))/1e3;
dFx=dx*(pad)/Lx;
dFy=dy*(pad)/Ly;

kxVec=-2*pi/(2*dx):2*pi/(dFx*Lx):(2*pi/(2*dx)-2*pi/(dFy*Lx));
kyVec=-2*pi/(2*dy):2*pi/(dFy*Ly):(2*pi/(2*dy)-2*pi/(dFy*Ly));
[kx, ky]=meshgrid(kxVec,kyVec);

k0=2*pi*data.f(f_num)/c;
kz=sqrt(k0.^2-kx.^2-ky.^2);
% end
%propogate fields
figure(2); clf;
 for z=-.1:.001:0;
% z=-0.02;
subplot(1,2,1);
 
fE=fft2(measurements(:,:,f_num,1),pad,pad);
fE=fE.*exp(-1.0j*kx*(Lx/2)).*exp(-1.0j*ky*(Ly/2));
imagesc(abs(ifftshift(fE)));
% imagesc(ifftshift(abs(fE))); imagesc(abs(ifft2(fE)));
% imagesc(abs(ifft2(fft2(measurements(:,:,jj,ii),pad,pad))))
axis xy; axis equal; axis tight;
fEz=fftshift(ifftshift(fE).*exp(-1j*(kz*z)));
% imagesc(real(kz));
% imagesc(abs(fEz));
E=ifft2(fEz);
subplot(1,2,2);
imagesc(X(1,:),Y(:,1),abs(E));
axis xy; axis equal; axis tight;
 drawnow

pause
 end
