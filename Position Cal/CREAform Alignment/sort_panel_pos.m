close all
data_panel=load('C:\Users\Jonah Gollub\Downloads\duke_system_fiducials.mat')


datapos =[9 9 5 5 1 2 2 6 6 10 ...
     1 10 7 7 3 8 4 11 8 4 ...
    12 12 3 11];

sortpos=zeros(12,2,3);

for i=1:24
if sortpos(datapos(i),1,1)==0;
    
sortpos(datapos(i),1,1)=data_panel.x0(i);
sortpos(datapos(i),1,2)=data_panel.y0(i);
sortpos(datapos(i),1,3)=data_panel.z0(i);

else
sortpos(datapos(i),2,1)=data_panel.x0(i);
sortpos(datapos(i),2,2)=data_panel.y0(i);
sortpos(datapos(i),2,3)=data_panel.z0(i);
end
end

panelcenter=[];
panelcenter_I=[];


for i=1:12
[~, panelcenter_I]=min(sortpos(i,:,1));
panelcenter(i,:)=sortpos(i,panelcenter_I,:);
end


%find B3 position for ref
panel_origin=panelcenter(7,:)
%subtract bow of circuit board
bevelthickness=1;

%panelcenter(1,3)=panelcenter((1,3)
%panelcenter(2,3)=panelcenter((2,3)
panelcenter(3,3)=panelcenter(3,3)-2.5+bevelthickness;
panelcenter(4,3)=panelcenter(4,3)-1.5+bevelthickness;
panelcenter(5,3)=panelcenter(5,3)-3+bevelthickness;
%panelcenter(6,3)=panelcenter(6,3);
panelcenter(7,3)=panelcenter(7,3)-1.5+bevelthickness;
panelcenter(8,3)=panelcenter(8,3)-2.5+bevelthickness;
panelcenter(9,3)=panelcenter(9,3)-1.5+bevelthickness;
panelcenter(10,3)=panelcenter(10,3)-3.5+bevelthickness;
panelcenter(11,3)=panelcenter(11,3)-1.5+bevelthickness;
%panelcenter(12,3)=panelcenter((12,3);

figure; 
view([0, 90])
xlim([0,2000])
ylim([0,2000])

hold on
for d=1:12; 
%scatter3(sortpos(d,1,1),sortpos(d,1,2),sortpos(d,1,3))
scatter3(panelcenter(d,1),panelcenter(d,2),panelcenter(d,3))
xlabel('x')
ylabel('y')
zlabel('z')
pause(.1)
end

%xpos

max(panelcenter(:,1))-min(panelcenter(:,1))


panelcenter(1:4,1)-mean(panelcenter(1:4,1))
panelcenter(5:8,1)-mean(panelcenter(5:8,1))
panelcenter(9:12,1)-mean(panelcenter(9:12,1))

%ypos

panelcenter([1 5 9],2)-mean(panelcenter([1 5 9],2))
panelcenter([2 6 10],2)-mean(panelcenter([2 6 10],2))
panelcenter([3 7 11],2)-mean(panelcenter([3 7 11],2))
panelcenter([4 8 12],2)-mean(panelcenter([4 8 12],2))


%zpos

panelcenter(:,3)-mean(panelcenter(:,3))

panelcenter(1:4,3)-mean(panelcenter(1:4,3))
panelcenter(5:8,3)-mean(panelcenter(5:8,3))
panelcenter(9:12,3)-mean(panelcenter(9:12,3))

%translate to center of B3
%translate x
panelcenter(:,1)=panelcenter(:,1)-panelcenter(7,1);
%translate y
panelcenter(:,2)=panelcenter(:,2)-panelcenter(7,2);


figure(2); 
view([0, 90])
xlim([-1500,1000])
ylim([-1500,1000])
hold on;
for d=1:12; 
%scatter3(sortpos(d,1,1),sortpos(d,1,2),sortpos(d,1,3))
scatter3(panelcenter(d,1),panelcenter(d,2),panelcenter(d,3))
xlabel('x')
ylabel('y')
zlabel('z')
pause(.1)
end


%translate to virtualizer coordinates
panelcenter=[panelcenter(:,3), panelcenter(:,1), panelcenter(:,2)];

figure(3);
cla
view([-90, 0])
zlim([-1500,1000])
ylim([-1500,1000])
hold on;
for d=1:12; 
%scatter3(sortpos(d,1,1),sortpos(d,1,2),sortpos(d,1,3))
scatter3(panelcenter(d,1),panelcenter(d,2),panelcenter(d,3))
xlabel('x')
ylabel('y')
zlabel('z')
pause(.1)
end

%convert to meters
panelcenter=panelcenter/1000;

rot=fitNormal(panelcenter,1);

%panelcenter_rot=[atand(rot(2)/rot(3)), atand(rot(3)/rot(2)), atand(rot(2)/rot(1))]


%% find antenna positions

close all

data_antenna=load('C:\Users\Jonah Gollub\Downloads\duke_system_w_covers_fiducials (1).mat')

figure(1)
cla

subplot(1,2,1)
scatter3(data_antenna.x0,data_antenna.y0,data_antenna.z0)
view([0, 90])
xlim([0,2000])
ylim([0,2000])
axis equal


subplot(1,2,2)

scatter3(data_antenna.x0,data_antenna.y0,data_antenna.z0)
max(data_antenna.z0)
min(data_antenna.z0)
view([0, 0])
xlim([0,2000])
zlim([-115,-70])
axis auto

datapos=[1 4 5 6 7 8 ...
         11 13 16 17 18 ...
         19 21 22 23 24 25 ...
         27 28 29 30 31 33 35 ...
         ];
     
     antenna_ordering=[16 17 13 7 20 10 6 3 18 21 14 23 11 8 4 1 5 24 22 19 15 12 9 2];
     
%%
d=datapos;
subplot(1,2,1)
hold on
scatter3(data_antenna.x0(d),data_antenna.y0(d),data_antenna.z0(d),'filled');
subplot(1,2,2)
hold on
scatter3(data_antenna.x0(d),data_antenna.y0(d),data_antenna.z0(d),'filled');

% 
%      datapanel=[2 3 9 10 12 15 20 26 32 34]
     
     antenna_pos=[]
     antenna_pos.x0=data_antenna.x0(datapos);
     antenna_pos.y0=data_antenna.y0(datapos);
     antenna_pos.z0=data_antenna.z0(datapos);

for i=1:24
new_anntenna_pos.x0(i)=antenna_pos.x0(find(antenna_ordering==i));
new_anntenna_pos.y0(i)=antenna_pos.y0(find(antenna_ordering==i));
new_anntenna_pos.z0(i)=antenna_pos.z0(find(antenna_ordering==i));
end

figure(3)
cla
for i=1:24
    hold on
scatter3(new_anntenna_pos.x0(i),new_anntenna_pos.y0(i),...
    new_anntenna_pos.z0(i))
view([0, 90])
xlim([0,2000])
ylim([0,2000])
axis equal
pause(.25)
end

%account for fiducials offset
new_anntenna_pos.x0=new_anntenna_pos.x0+8;
new_anntenna_pos.y0=new_anntenna_pos.y0-8;

%translate to NFS plane
new_anntenna_pos.x0=new_anntenna_pos.x0-panel_origin(1);
new_anntenna_pos.y0=new_anntenna_pos.y0-panel_origin(2);
new_anntenna_pos.z0=new_anntenna_pos.z0-panel_origin(3)+55;

figure(2); 
cla
for d=1:24; 
%scatter3(sortpos(d,1,1),sortpos(d,1,2),sortpos(d,1,3))
scatter3(new_anntenna_pos.x0(d),new_anntenna_pos.y0(d),new_anntenna_pos.z0(d))
view([0, 90])
xlim([-1500,1000])
ylim([-1500,1000])
hold on;
xlabel('x')
ylabel('y')
zlabel('z')
pause(.1)
end

%translate to virtualizer coordinates
final_anntenna_pos=[new_anntenna_pos.z0(:),...
                    new_anntenna_pos.x0(:), ...
                    new_anntenna_pos.y0(:)...
                    ];

%convert to meters
final_anntenna_pos=final_anntenna_pos/1000;

figure(5)
for d=1:24; 
hold on
scatter3(final_anntenna_pos(d,1),...
    final_anntenna_pos(d,2),...
    final_anntenna_pos(d,3)...
    )
view([-90, 0])
zlim([-1.500,1.000])
ylim([-1.500,1.000])

xlabel('x')
ylabel('y')
zlabel('z')
pause(.1)
end



