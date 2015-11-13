function [ Aup ] = upsample_image_3D(A, upsample_rate)

Nx = size(A,2);
Ny = size(A,1);
Nz = size(A,3);

[x y z] = meshgrid(linspace(0,Nx,Nx),linspace(0,Ny,Ny),linspace(0,Nz,Nz));
[xup yup zup] = meshgrid(linspace(0,Nx,round(upsample_rate*Nx)),linspace(0,Ny,round(upsample_rate*Ny)),linspace(0,Nz,round(upsample_rate*Nz)));

Aup=interp3(x,y,z,A,xup,yup,zup,'linear');
end

