
function [Aup] = upsample_image(A,upsample_rate)

Nx = size(A,2);
Ny = size(A,1);

[x y] = meshgrid(linspace(0,Nx,Nx),linspace(0,Ny,Ny));
[xup yup] = meshgrid(linspace(0,Nx,round(upsample_rate*Nx)),linspace(0,Ny,round(upsample_rate*Ny)));

Aup = interp2(x,y,A,xup,yup,'*spline');