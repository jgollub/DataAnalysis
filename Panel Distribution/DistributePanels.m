figure(201); clf;

%constants
cm=1e-2;

%system size
panelSize=13*cm;
stepy=13*cm;
stepz=14*cm;
systemWidth=1.4;
systemHeight=1.75;

[yPos,zPos]=meshgrid(...
    linspace(-systemWidth/2,systemWidth/2,floor(systemWidth/stepy)),...
    linspace(-systemHeight/2,systemHeight/2,floor(systemHeight/stepz))....
      );
xPos=zeros(size(yPos));
  
if abs((zPos(1,1)-zPos(2,1)))<panelSize | abs((yPos(1,1)-yPos(1,2)))<panelSize
    error('aperture area is not large enough to fit panels')
end

  scatter3(xPos(:),yPos(:),zPos(:),10,'k');
  hold on;


  grid=zeros(size(yPos));
  
  %set Tx positions
  TxGrid=logical(grid);
  TxGrid([3,10],2:9)=logical(1);
%     TxGrid([5],4:2:6)=logical(1);
  TxGrid([6],3:2:9)=logical(1);
  TxGrid([7],2:2:9)=logical(1);
    scatter3(xPos(TxGrid(:)),yPos(TxGrid(:)),zPos(TxGrid(:)),10,'r','filled');

    axis equal;
    view(90,0);
%Rx sub-panel

  %set Tx positions
  RxGrid=logical(grid);
  RxGrid(:,[1,end])=logical(1);
%      TxGrid([5],4:2:6)=logical(1);
  RxGrid([1 2 4 5],1:2:9)=logical(1);
   RxGrid([8 9],2:2:9)=logical(1);
  RxGrid([5 6],[2 4 6 8])=logical(1);
  RxGrid([end end-1],1:2:9)=logical(1);
 
    scatter3(xPos(RxGrid(:)),yPos(RxGrid(:)),zPos(RxGrid(:)),10,'blue','filled');

    axis equal;
    view(90,0);
%Tx panels

% %panel Positions
% [TX_x,TX_y,TX_z]=meshgrid(0,linspace(-0.6,0.6,8),linspace(-.7,.7,3));
% TX_Positions=horzcat(TX_x(:),TX_y(:),TX_z(:));
% 
% [RX_x,RX_y,RX_z]=meshgrid(0,linspace(-.8,.8,5),linspace(-1,1,14));
% RX_Positions=horzcat(RX_x(:),RX_y(:),RX_z(:));
% 
% figure(101); clf; hold on;
% scatter3(RX_x(:),RX_y(:),RX_z(:),10,'b','filled');
% scatter3(TX_x(:),TX_y(:),TX_z(:),10,'r','filled');



%%%%%%%%%%%%%%

% yr=[3*15/2, 3*15/2, 3*15/2, 3*15/2, 3*15/2, 3*15/2];
% zr=[3*15/2, 5*15/2, 11*15/2, 13*15/2,  19*15/2, 21*15/2];
% ty=[yr,yr+6*15/2, yr+12*15/2, yr+18*15/2]/100;
% tz=[zr,zr,zr,zr]/100;
% 
% yt=[15/2, 15/2, 15/2, 15/2, 15/2, 15/2];
% zt=[3*15/2, 5*15/2, 11*15/2, 13*15/2,  19*15/2, 21*15/2];
% yt=[yt,yt+4*15/2, yt+6*15/2,yt+10*15/2,yt+12*15/2,yt+16*15/2,yt+18*15/2,yt+22*15/2];
% zt=[zt,zt,zt,zt,zt,zt,zt,zt];
% yt2=[3*15/2, 3*15/2, 3*15/2, 3*15/2, 3*15/2, 3*15/2];
% zt2=[1*15/2, 7*15/2, 9*15/2,  15*15/2, 17*15/2,  23*15/2];
% yt2=[yt2, yt2+6*15/2, yt2+12*15/2, yt2+18*15/2];
% zt2=[zt2,zt2,zt2,zt2];
% ry=[yt,yt2]/100;
% rz=[zt,zt2]/100;



%%%%%%%%%%%%%%%%


yr1=[15/2, 15/2, 15/2*3, 15/2*3, 15/2*5, 15/2*5];
zr1=[3*15/2, 5*15/2, 1*15/2, 7*15/2,  3*15/2, 5*15/2];
yt1=[15/2*3, 15/2*3];
zt1=[3*15/2, 5*15/2];
yr1=[yr1,yr1+3*15, yr1+6*15, yr1+9*15, yr1-1*15, yr1+10*15];
zr1=[zr1,zr1,zr1,zr1, zr1+7*15, zr1+7*15];
yt1=[yt1,yt1+3*15, yt1+6*15, yt1+9*15, yt1-1*15, yt1+10*15];
zt1=[zt1,zt1,zt1,zt1, zt1+7*15, zt1+7*15];
yr2=[15/2, 15/2*3, 15/2*3, 15/2*5, 15/2*5, 15/2*7];
zr2=[3*15/2, 5*15/2, 1*15/2, 5*15/2,  1*15/2, 3*15/2];
yt2=[15/2*3, 15/2*5];
zt2=[3*15/2, 3*15/2];
yr2=[yr2,yr2+4*15,yr2+8*15,yr2+2*15,yr2+6*15,yr2+4*15];
zr2=[zr2+4*15, zr2+4*15, zr2+4*15, zr2+7*15, zr2+7*15, zr2+10*15];
yt2=[yt2,yt2+4*15,yt2+8*15,yt2+2*15,yt2+6*15,yt2+4*15];
zt2=[zt2+4*15, zt2+4*15, zt2+4*15, zt2+7*15, zt2+7*15, zt2+10*15];
ty=[yt1,yt2]/100;tz=[zt1,zt2]/100;ry=[yr1,yr2]/100;rz=[zr1,zr2]/100;

%% %%%%%%%%%%%%%%%%%%%%


A=18;

yr2=[A/2, A/2, A/2*3, A/2*3, A/2*5, A/2*5];
zr2=[3*A/2, 5*A/2, 1*A/2, 7*A/2,  3*A/2, 5*A/2];
yt2=[A/2*3, A/2*3];
zt2=[3*A/2, 5*A/2];

yr1=[A/2, A/2*3, A/2*3, A/2*5, A/2*5, A/2*7];
zr1=[3*A/2, 5*A/2, 1*A/2, 5*A/2,  1*A/2, 3*A/2];
yt1=[A/2*3, A/2*5];
zt1=[3*A/2, 3*A/2];

step_diag_a=[2,-2]; 
step_diag_b=[2,-2]; 
diag1_pos=[0 0];
diag2_pos=[-3.5 -.5];
diag3_pos=[-6 -2];
diag4_pos=[-5.5 -6.5];
diag1_num=[0:1];
diag2_num=[0:3];
diag3_num=[0:3];
diag4_num=[0:1];

diag1=[0*diag1_num', (diag1_num.*step_diag_a(1)+diag1_pos(1))', (diag1_num*step_diag_a(2)+diag1_pos(2))'];
diag2=[0*diag2_num', (diag2_num.*step_diag_b(1)+diag2_pos(1))', (diag2_num*step_diag_b(2)+diag2_pos(2))'];
diag3=[0*diag3_num', (diag3_num.*step_diag_a(1)+diag3_pos(1))', (diag3_num*step_diag_a(2)+diag3_pos(2))'];
diag4=[0*diag4_num', (diag4_num.*step_diag_b(1)+diag4_pos(1))', (diag4_num*step_diag_b(2)+diag4_pos(2))'];

temp=[diag1; diag3; diag2; diag4]*20;
rx=zeros(72,1);
ry=rx; rz=rx;
tx=zeros(24,1);
ty=tx; tz=tx;
for ii=1:6
rx((ii-1)*6+1:ii*6)=temp(ii,1)+0.*(1:6)';
ry((ii-1)*6+1:ii*6)=temp(ii,2)+yr1;
rz((ii-1)*6+1:ii*6)=temp(ii,3)+zr1;

tx((ii-1)*2+1:ii*2)=temp(ii,1)+0.*(1:2)';
ty((ii-1)*2+1:ii*2)=temp(ii,2)+yt1;
tz((ii-1)*2+1:ii*2)=temp(ii,3)+zt1;

end

for ii=7:12
rx((ii-1)*6+1:ii*6)=temp(ii,1);
ry((ii-1)*6+1:ii*6)=temp(ii,2)+yr2;
rz((ii-1)*6+1:ii*6)=temp(ii,3)+zr2;

tx((ii-1)*2+1:ii*2)=temp(ii,1)+0.*(1:2)';
ty((ii-1)*2+1:ii*2)=temp(ii,2)+yt2;
tz((ii-1)*2+1:ii*2)=temp(ii,3)+zt2;
end
rx=rx/100;
ry=ry/100;
rz=rz/100;

tx=tx/100;
ty=ty/100;
tz=tz/100;

max(ry)-min(ry)
max(rz)-min(rz)

clf
scatter3(rx,ry,rz,10,'g','filled'); hold on;
scatter3(tx,ty,tz,10,'r','filled');
    axis equal;
    view(90,0);


%% %%%%%%%%%%%%%%%%%%%%


A=18;%18

yr2=[A/2, A/2, A/2*3, A/2*3, A/2*5, A/2*5];
zr2=[3*A/2, 5*A/2, 1*A/2, 7*A/2,  3*A/2, 5*A/2];
yt2=[A/2*3, A/2*3];
zt2=[3*A/2, 5*A/2];

yr1=[A/2, A/2*3, A/2*3, A/2*5, A/2*5, A/2*7];
zr1=[3*A/2, 5*A/2, 1*A/2, 5*A/2,  1*A/2, 3*A/2];
yt1=[A/2*3, A/2*5];
zt1=[3*A/2, 3*A/2];

step_diag_a=[2,-2]; 
step_diag_b=[2,-2]; 
diag1_pos=[0 0];
diag2_pos=[-3.5 -.5];
diag3_pos=[-6 -2];
diag4_pos=[-5.5 -6.5];
diag1_num=[0:1];
diag2_num=[0:3];
diag3_num=[0:3];
diag4_num=[0:1];

r_amp=0;
diag1=[0*diag1_num', (diag1_num.*step_diag_a(1)+r_amp*rand([1,2])+diag1_pos(1))', (diag1_num*step_diag_a(2)+r_amp*rand([1,2])+diag1_pos(2))'];
diag2=[0*diag2_num', (diag2_num.*step_diag_b(1)+r_amp*rand([1,4])+diag2_pos(1))', (diag2_num*step_diag_b(2)+r_amp*rand([1,4])+diag2_pos(2))'];
diag3=[0*diag3_num', (diag3_num.*step_diag_a(1)+r_amp*rand([1,4])+diag3_pos(1))', (diag3_num*step_diag_a(2)+r_amp*rand([1,4])+diag3_pos(2))'];
diag4=[0*diag4_num', (diag4_num.*step_diag_b(1)+r_amp*rand([1,2])+diag4_pos(1))', (diag4_num*step_diag_b(2)+r_amp*rand([1,2])+diag4_pos(2))'];

temp=[diag1; diag3; diag2; diag4]*(A+2); %2
rx=zeros(72,1);
ry=rx; rz=rx;
tx=zeros(24,1);
ty=tx; tz=tx;
for ii=1:6
rx((ii-1)*6+1:ii*6)=temp(ii,1)+0.*(1:6)';
ry((ii-1)*6+1:ii*6)=temp(ii,2)+yr1;
rz((ii-1)*6+1:ii*6)=temp(ii,3)+zr1;

tx((ii-1)*2+1:ii*2)=temp(ii,1)+0.*(1:2)';
ty((ii-1)*2+1:ii*2)=temp(ii,2)+yt1;
tz((ii-1)*2+1:ii*2)=temp(ii,3)+zt1;

end

for ii=7:12
rx((ii-1)*6+1:ii*6)=temp(ii,1);
ry((ii-1)*6+1:ii*6)=temp(ii,2)+yr2;
rz((ii-1)*6+1:ii*6)=temp(ii,3)+zr2;

tx((ii-1)*2+1:ii*2)=temp(ii,1)+0.*(1:2)';
ty((ii-1)*2+1:ii*2)=temp(ii,2)+yt2;
tz((ii-1)*2+1:ii*2)=temp(ii,3)+zt2;
end
rx=(rx-mean(rx))/100;
ry=(ry-mean(ry))/100;
rz=(rz-mean(rz))/100;

tx=(tx-mean(tx))/100;
ty=(ty-mean(ty))/100;
tz=(tz-mean(tz))/100;

max(ry)-min(ry)
max(rz)-min(rz)

figure(999)
clf

scatter3(rx,ry,rz,10,'g','filled'); hold on;
scatter3(tx,ty,tz,10,'r','filled');
    axis equal;
    view(90,0);
    ylabel('y axis')
    zlabel('z axis')
ylim([-1.25 1.25])
zlim([-1.25 1.25])



%%  Sub-Assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% orientation 
%layout
cm=1e-2;
spacing_y=17*cm;
spacing_z=17*cm;

[Ys,Zs]=meshgrid(-spacing_y:spacing_y:spacing_y,spacing_z:-spacing_z:-spacing_z);
keep=logical(ones(1,9));
keep(5)=logical(0); 
Xs=zeros(size(Ys));
Ys=Ys(keep);
Zs=Zs(keep);

Or1=logical([0 1 0 1 0 0 0 0]); 
Or2=logical([0 0 0 1 0 0 1 0]);
Or3=logical([0 1 0 0 1 0 0 0]); 
Or4=logical([0 0 0 0 1 0 1 0]); 
Or={Or4,Or4,Or2,Or2,...
    Or2,Or1,Or1,Or1,...
    Or3,Or3,Or1, Or1};

% figure(200); clf; hold on; 
% 
% for ii=1:4
%     subplot(2,2,ii); hold on;
% scatter3(Xs(~OrTest{ii}),Ys(~OrTest{ii}),Zs(~OrTest{ii}),10,'g')
% scatter3(Xs(OrTest{ii}),Ys(OrTest{ii}),Zs(OrTest{ii}),10,'r')
% view(90,0);
% title(['Or ',num2str(ii)]);
% view(90,0);
% end
y_range=1-17*1.5*cm;
z_range=1.08-17*1.5*cm;
[Ya, Za]=meshgrid(linspace(-y_range,y_range,3),linspace(z_range,-z_range,4));
% 
% shift_Ya=[+.2;+.1;0.05;+.2;...
%             0;0.05;-.05;0;...
%             -.2;-0.05;-.1;-.2];           

shift_Ya=[+.2;+.12;0.08;+.2;...
            0;0.05;-.05;0;...
            -.2;-0.07;-.13;-.2];         

 Ya=Ya(:)+shift_Ya;
 
 rng(1);
shift_Za=0.02*(rand(12,1)-.5)
Za=Za(:)+shift_Za;

ty=zeros(24,1);
tz=ty;
tx=ty;
ry=zeros(72,1);
rz=ry;
rx=ry;

for ii=1:12
ty(2*(ii-1)+1:2*ii)=Ys(Or{ii})+Ya(ii);
tz(2*(ii-1)+1:2*ii)=Zs(Or{ii})+Za(ii);
    
ry(6*(ii-1)+1:6*ii)=Ys(~Or{ii})+Ya(ii);
rz(6*(ii-1)+1:6*ii)=Zs(~Or{ii})+Za(ii);
end

figure(201); clf;
scatter3(tx,ty,tz,10,'g'); hold on;
scatter3(rx,ry,rz,10,'r')
    view(90,0)
    axis equal; axis tight;