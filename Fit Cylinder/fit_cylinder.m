%% dummy data for testing
xoff=-3;
yoff=-2;
radius=1;
zoff=1;
[xx,yy,zz]=meshgrid(-10:.1:10,-10:.1:10,-10:.1:10);

Cp=[1 0 0 -xoff;
    0 1 0 -yoff;
    0 0 0 0;
    0 0 0  -(radius^2)];

%rot z
t1=[ cos(pi/4) -sin(pi/4)  0  0;
     sin(pi/4)  cos(pi/4)  0  0;
     0          0          1  0;
     0          0          0  1];

%rot y
t2=[ cos(pi/4)     0      sin(pi/4)  0;
     0             1      0          0;
     -sin(pi/4)    0     cos(pi/4)   0;
     0             0     0           1];

 %rot x
t3=[ 1           0        0           0;
     0         cos(pi/4)  -sin(pi/4)  0;
     0         sin(pi/4)  cos(pi/4)   0;
     0           0         0          1]; 
 
 t4=[1           0        0          3
     0           1        0          3;
     0           0        1          0;
     0           0        0          1]; 
 
 M=t1*t2*t3*t4';
 
Cdummy=inv(M)*Cp*inv(M).';


dummy_grid=[xx(:),yy(:),zz(:),1*ones(size(xx(:)))].';
dummy_quadric=sum(dummy_grid.*(Cdummy*dummy_grid),1);
dummy_selectpts=(abs(dummy_quadric))<.05 & zz(:).'<0;
 xdata=xx(dummy_selectpts).';
 ydata=yy(dummy_selectpts).';
 zdata=zz(dummy_selectpts).';

figure(1); clf; hold on;
scatter3(xdata,ydata,zdata,10,'filled','y'); 
 axis equal; axis tight; xlabel('x');ylabel('y'); zlabel('z');


scatter3(0,0,0,50,'blue','filled')
% pts=(abs((zz-zoff).^2+(yy-yoff).^2-radius^2)<.1 & zz<zoff);
% 

% scatter3(xdata,ydata,zdata,10,'filled','k'); 
data=[xdata,ydata,zdata,ones(size(xdata))].';
%% Experimental Data
Optical_Scan='D:\MetaImager\CREAform\Cylinder Experiment\dataforpassivecalibration\creaform.txt';
Optical_Data=dlmread(Optical_Scan,'\t',0,0);
opt_fiducial = Optical_Data(Optical_Data(:,9)==0,1:3)/1000; %also convert m
target_pts=opt_fiducial(Optical_Data(:,9)<1 ,[1 2 3]); %any coded fiducial less than 7 is coord sys or bar 

%isolate cylinder
target_pts=target_pts(target_pts(:,3)>.8 & target_pts(:,2)>-1.6,:);

figure(2); clf;
scatter3(target_pts(:,1),target_pts(:,2),target_pts(:,3),5,'filled','k');
axis equal; axis tight; xlabel('x');ylabel('y'); zlabel('z');

data=[target_pts,ones(size(target_pts(:,1)))].';
delta=.005;
[xx, yy, zz]=meshgrid(min(target_pts(:,1)):delta:max(target_pts(:,1)),...
                      min(target_pts(:,2)):delta:max(target_pts(:,2)),...
                      (-.8+min(target_pts(:,3))):delta:(+.8+max(target_pts(:,3))));
hold on;
%  scatter3(xx(:),yy(:),zz(:),1,'filled','g');
%% solve for best cylindrical fit
figure(2)
%function for vectorizing lower triangle of matrix
vech=@(M) M(tril(true(size(M))));

%populate A matrix
A=zeros(size(data,2),10);
for ii=1:size(data,2);
    M=data(:,ii)*(data(:,ii).');
    A(ii,:)=vech(M); 
end

%SVD
[U,S,V]=svd(A,0);
% extract most degenerate fit (least square using SVD)
C=V(:, end);
C=  [ C(1) .5*C(2) .5*C(3) .5*C(4); ...
      .5*C(2) C(5) .5*C(6) .5*C(7); ...
      .5*C(3) .5*C(6) C(8) .5*C(9); ...
      .5*C(4) .5*C(7) .5*C(9) C(10) ];

grid=[xx(:),yy(:),zz(:),1*ones(size(xx(:)))].';
quadric=sum(grid.*(C*grid),1);
selectpts=abs(quadric)<.005^2;

hold on;
scatter3(xx(selectpts),yy(selectpts),zz(selectpts),5,'r');
 axis equal; axis tight; xlabel('x');ylabel('y'); zlabel('z');
