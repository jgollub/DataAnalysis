data=load('D:\Dropbox (Duke Electric & Comp)\WPT (1) (1)\Printed Cavity\Hughes_Latest\Tx_Panel_101_both_pol.mat')
figure(1); clf;
f_num=1;
c=3e8;
measurements=flip(flip(data.measurements,1),2);
X=flip(data.X,2);
Y=flip(data.Y,1);
imagesc(X(1,:),Y(:,1),abs(measurements(:,:,f_num,1)));
axis xy; axis equal; axis tight;

dx=(X(1,2)-X(1,1))/1e3;
dy=(Y(2,1)-Y(1,1))/1e3;
Lx=(X(1,end)-X(1,1))/1e3; Ly=(Y(end,1)-Y(1,1))/1e3;
kxVec=-2*pi/(2*dx):2*pi/(Lx):2*pi/(2*dx);
kyVec=-2*pi/(2*dy):2*pi/(Ly):2*pi/(2*dy);
[kx, ky]=meshgrid(kxVec,kyVec);

for ii=1:1
k0=2*pi*data.f(ii)/c;
kz=sqrt(k0.^2-kx.^2-ky.^2);
end
%propogate fields
figure(2); clf;
% for z=-.07:.001:0;
z=-.06
subplot(1,2,1);
for ii=1:1
fE=fft2(fftshift(measurements(:,:,1,ii)));
imagesc(kxVec,kyVec,abs(ifftshift(fE)));
axis xy; axis equal; axis tight;
fE=fftshift(ifftshift(fE).*exp(-1j*kz*z));
E=ifftshift(ifft2(fE));
end
subplot(1,2,2);
imagesc(kxVec,kyVec,abs(E));
axis xy; axis equal; axis tight;

