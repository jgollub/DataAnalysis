close all

data=load('C:\Users\Jonah Gollub\Downloads\duke_system_w_covers_fiducials (1).mat')

figure(1)
cla

subplot(1,2,1)
scatter3(data.x0,data.y0,data.z0)
view([0, 90])
xlim([0,2000])
ylim([0,2000])
axis equal


subplot(1,2,2)

scatter3(data.x0,data.y0,data.z0)
max(data.z0)
min(data.z0)
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
scatter3(data.x0(d),data.y0(d),data.z0(d),'filled');
subplot(1,2,2)
hold on
scatter3(data.x0(d),data.y0(d),data.z0(d),'filled');

% 
%      datapanel=[2 3 9 10 12 15 20 26 32 34]
     
     antenna_pos=[]
     antenna_pos.x0=data.x0(datapos);
     antenna_pos.y0=data.y0(datapos);
     antenna_pos.z0=data.z0(datapos);

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


