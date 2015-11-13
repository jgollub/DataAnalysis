function [] = voxelimage(S,Az,El,Z,threshold)
% voxelimage.m, John Hunt, 04.03.13
%S  - 3D scattering density matrix 
%Az - 2D azimuth plaid matrix of the type produced by 'meshgrid'
%El - 2D elevation plaid matrix of the type produced by 'meshgrid'
%Z  - 1D vector of ranges to each range plane
%threshold - [0,1] threshold the image (optional, default is 0) 

daz = Az(1,2) - Az(1,1);
del = El(2,1) - El(1,1);
dz = Z(2) - Z(1);

Az = repmat(Az,[1,1,length(Z)]);
El = repmat(El,[1,1,length(Z)]);
Z1(1,1,:) = Z;
Z = repmat(Z1,[size(Az,1),size(Az,2),1]);
X = Z.*tan(Az);
Y = Z.*tan(El)./cos(Az);

a = gca;
%cla

S = abs(S);
S = S./maxall(S);
S = S.^2; %reduce the intensity of low-scattering voxels so they don't wash out the image
 
for nx=1:size(S,2) 
    for ny=1:size(S,1)
        for nz=1:size(S,3)
            if (nargin==5)
                 t = threshold;
            else
                 t = 0;
            end
            if S(ny,nx,nz)>t
                 [x,y,z] = voxelpatch(X(ny,nx,nz),Y(ny,nx,nz),Z(ny,nx,nz),daz,del,dz);
%                  [x,y,z] = voxelpatch(X(nx,ny,nz),Y(nx,ny,nz),Z(nx,ny,nz),daz,del,dz);
                p = patch(x,z,y,'w');
                 set(p,'FaceAlpha',S(ny,nx,nz));
%                  set(p,'FaceAlpha',S(nx,ny,nz));
                set(p,'Edgecolor','none');
            end
        end
    end
end

xlabel('X')
ylabel('Z')
zlabel('Y')
axis([minall(X) maxall(X) minall(Z)-dz/2 maxall(Z)+dz/2 minall(Y) maxall(Y)])
set(a,'Color',[0,0,0])

end

%%
function [X Y Z] = voxelpatch(xc,yc,zc,daz,del,dz)
% voxelpatch.m, John Hunt, 04.03.13
%
%[X Y Z] = voxelpatch(xc,yc,zc,daz,del,dz)
%
%Inputs:
%daz         - voxel azimuth (angle in xz-plane) length
%del         - voxel elevation (angle above xz-plane) length
%dz         - voxel z side length
%(xc,yc,zc) - voxel center location

%% find vectors along each of the four sides of the voxel that are not in the
%range planes

%rotation matrix that rotates the vector at the center of the voxel in the
%azimuth direction
Raz = rot('y',daz/2);

r0 = [xc;yc;zc];
r01 = Raz'*r0;
r02 = Raz'*r0;
r03 = Raz*r0;
r04 = Raz*r0;

%find the vectors around which to rotate the r0x vectors into the elevation direction
e01 = cross(r01,[0;1;0]);
e02 = cross(r02,[0;1;0]);
e03 = cross(r03,[0;1;0]);
e04 = cross(r04,[0;1;0]);

%rotate r0x vectors in the elevation direction (angle away from the xz-plane)
r01 = rot(e01,-del/2)*r01;
r02 = rot(e02,del/2)*r02;
r03 = rot(e03,del/2)*r03;
r04 = rot(e04,-del/2)*r04;

%% find the intersections of these edge vectors with the range planes
r1 = (zc-dz/2)/r01(3)*r01;
r2 = (zc-dz/2)/r02(3)*r02;
r3 = (zc-dz/2)/r03(3)*r03;
r4 = (zc-dz/2)/r04(3)*r04;
r5 = (zc+dz/2)/r01(3)*r01;
r6 = (zc+dz/2)/r02(3)*r02;
r7 = (zc+dz/2)/r03(3)*r03;
r8 = (zc+dz/2)/r04(3)*r04;

%%
     %top                     %bottom                  %front                   %back                    %left                    %right
X = [r1(1) r5(1) r8(1) r4(1); r2(1) r3(1) r7(1) r6(1); r1(1) r4(1) r3(1) r2(1); r5(1) r6(1) r7(1) r8(1); r3(1) r4(1) r8(1) r7(1); r1(1) r2(1) r6(1) r5(1)]';
Y = [r1(2) r5(2) r8(2) r4(2); r2(2) r3(2) r7(2) r6(2); r1(2) r4(2) r3(2) r2(2); r5(2) r6(2) r7(2) r8(2); r3(2) r4(2) r8(2) r7(2); r1(2) r2(2) r6(2) r5(2)]';
Z = [r1(3) r5(3) r8(3) r4(3); r2(3) r3(3) r7(3) r6(3); r1(3) r4(3) r3(3) r2(3); r5(3) r6(3) r7(3) r8(3); r3(3) r4(3) r8(3) r7(3); r1(3) r2(3) r6(3) r5(3)]';

end

%%
function R = rot(axis,angle)

if ischar(axis)
    switch axis
        case 'x'
            R = [1 0 0;0 cos(angle) -sin(angle);0 sin(angle) cos(angle)];
        case 'y'
            R = [cos(angle) 0 sin(angle);0 1 0; -sin(angle) 0 cos(angle)];
        case 'z'
            R = [cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0; 0 0 1];
    end
else
        u = axis./sqrt(sum(axis.^2));
        R = [cos(angle)+u(1)^2*(1-cos(angle)),         u(1)*u(2)*(1-cos(angle))-u(3)*sin(angle), u(1)*u(3)*(1-cos(angle))+u(2)*sin(angle);...
             u(2)*u(1)*(1-cos(angle))+u(3)*sin(angle), cos(angle)+u(2)^2*(1-cos(angle)),         u(2)*u(3)*(1-cos(angle))-u(1)*sin(angle);...
             u(3)*u(1)*(1-cos(angle))-u(2)*sin(angle), u(3)*u(2)*(1-cos(angle))+u(1)*sin(angle), cos(angle)+u(3)*(1-cos(angle))];
end
           
end

%%
function [mina] = minall(A)

dim = length(size(A));
for n=1:dim
     A = squeeze(min(A));
end
    mina = A;
end

%%
function [maxa] = maxall(A)

dim = length(size(A));
for n=1:dim
    A = squeeze(max(A));
end
    maxa = A;
end
