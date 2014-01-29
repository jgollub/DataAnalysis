function voxelimage_projections(S,Az,El,Z,threshold)
% voxelimage_projections.m, John Hunt, 04.03.13
%draws a four-panel figure with one orthoscopic and three planar projections
%of the data in the scattering density matrix S

%S  - 3D scattering density matrix 
%Az - 2D azimuth plaid matrix of the type produced by 'meshgrid'
%El - 2D elevation plaid matrix of the type produced by 'meshgrid'
%Z  - 1D vector of ranges to each range plane
%threshold - [0,1] threshold the image (optional, default is zero)

if (nargin==5)
    t = threshold;
else
     t = 0;
end
            
subplot(2,2,1)
axis equal
view([45,45])
voxelimage(S,Az,El,Z,t)

subplot(2,2,2)
axis equal
voxelimage(S,Az,El,Z,t)
view([90,0])

subplot(2,2,3)
axis equal
voxelimage(S,Az,El,Z,t)
view([0,90])

subplot(2,2,4)
axis equal
view([0,0])
voxelimage(S,Az,El,Z,t)


end
