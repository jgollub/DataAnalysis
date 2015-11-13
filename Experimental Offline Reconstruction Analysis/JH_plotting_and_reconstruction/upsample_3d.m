
function [Aup] = upsample_3d(A,upsample_rate)

Nx = size(A,2);
Ny = size(A,1);
Nz = size(A,3);

XI = linspace(1,Nx,round(upsample_rate*Nx));
YI = linspace(1,Ny,round(upsample_rate*Ny));
ZI = linspace(1,Nz,round(upsample_rate*Nz));
[XI YI ZI] = meshgrid(XI,YI,ZI);

Aup = interp3(A,XI,YI,ZI,'linear');
