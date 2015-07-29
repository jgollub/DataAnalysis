xplot_range=[-.1 .2];
yplot_range=[.2 .5];
zplot_range=[-1.2,-.9];

Center_Pos=[.07 .32 -1.07];
% plot Sphere calculated positions
%MapSpecularyReflect_sphere only
probe_group_color=[];
Legend_fill=[];
for j_probe=1:6

    figure(Sphere_Pos_fig)
    hold all;
    probe_group_color(j_probe,:)=[1-j_probe/Num_Probes,j_probe/Num_Probes,0];
    Legend_fill{j_probe}=['\color[rgb]{',num2str(probe_group_color(j_probe,1)),',',num2str(probe_group_color(j_probe,2)),',',num2str(probe_group_color(j_probe,3)),'}Probe ',num2str(j_probe)]

subplot(2,2,1)
hold all;
plot(sphere_center(:,1,j_probe),sphere_center(:,2,j_probe),'Marker', '.','Color',probe_group_color(j_probe,:) ,'LineStyle','none')
xlabel('X (m)');
ylabel('Y (m)');
xlim(xplot_range);
ylim(yplot_range);
axis equal;
title('XY Alignment')
drawnow;

subplot(2,2,2)
hold all;
plot(sphere_center(:,2,j_probe),sphere_center(:,3,j_probe),'Marker', '.','Color',probe_group_color(j_probe,:), 'LineStyle','none')
xlabel('Y (m)');
ylabel('Z (m)');
xlim(yplot_range);
ylim(zplot_range);
title('YZ Alignment')
axis equal;
drawnow;

subplot(2,2,3)
hold all;
plot(sphere_center(:,3,j_probe),sphere_center(:,1,j_probe),'Marker', '.','Color',probe_group_color(j_probe,:) ,'LineStyle','none')
xlabel('Z (m)');
ylabel('X (m)');
xlim(zplot_range);
ylim(xplot_range);
title('ZX Alignment');
axis equal;
drawnow;
end
%


%%%% Average all panels each probe
% Legend_fill=[];
% for j_probe=1:6
%     
%     figure(Sphere_Pos_fig)
%     
%     probe_group_color(j_probe,:)=[1-j_probe/Num_Probes,j_probe/Num_Probes,0];
%     Legend_fill{j_probe}=['\color[rgb]{',num2str(probe_group_color(j_probe,1)),',',num2str(probe_group_color(j_probe,2)),',',num2str(probe_group_color(j_probe,3)),'}Probe ',num2str(j_probe)]
%     
%     subplot(2,2,1)
%     hold all;
% 
%     plot(mean(sphere_center(:,1,j_probe)),mean(sphere_center(:,2,j_probe)),'Marker', '.','Color',probe_group_color(j_probe,:) ,'LineStyle','none')
%     
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     xlim(xplot_range);
%     ylim(yplot_range);
%     axis equal;
%     title('XY Alignment')
%     drawnow;
%     
%     subplot(2,2,2)
%     hold all;
% 
%     plot(mean(sphere_center(:,2,j_probe)),mean(sphere_center(:,3,j_probe)),'Marker', '.','Color',probe_group_color(j_probe,:), 'LineStyle','none')
%     
%     xlabel('Y (m)');
%     ylabel('Z (m)');
%     xlim(yplot_range);
%     ylim(zplot_range);
%     title('YZ Alignment')
%     axis equal;
%     drawnow;
%     
%     subplot(2,2,3)
%     hold all;
% 
%     plot(mean(sphere_center(:,3,j_probe)),mean(sphere_center(:,1,j_probe)),'Marker', '.','Color',probe_group_color(j_probe,:) ,'LineStyle','none')
%     
%     xlabel('Z (m)');
%     ylabel('X (m)');
%     xlim(zplot_range);
%     ylim(xplot_range);
%     title('ZX Alignment');
%     drawnow;
%     axis equal;
%     
% end


%%%%by panel

% probe_group_color=[];
% for j_panel=1:12
% 
%     figure(Sphere_Pos_fig)
%     hold all;
%     probe_group_color(j_panel,:)=[1-j_panel/Num_Panels,j_panel/Num_Panels,0];
%     Legend_fill{j_panel}=['\color[rgb]{',num2str(probe_group_color(j_panel,1)),',',num2str(probe_group_color(j_panel,2)),',',num2str(probe_group_color(j_panel,3)),'}Panel ',num2str(j_panel)]
% 
%     subplot(2,2,1)
%     hold all;
%     plot(mean(sphere_center(j_panel,1,:)),mean(sphere_center(j_panel,2,:)),'Marker', '.','Color',probe_group_color(j_panel,:) ,'LineStyle','none')
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     xlim(xplot_range);
%     ylim(yplot_range);
%     axis equal;
% 
%     title('XY Alignment')
%     drawnow;
% 
%     subplot(2,2,2)
%     hold all;
%     plot(mean(sphere_center(j_panel,2,:)),mean(sphere_center(j_panel,3,:)),'Marker', '.','Color',probe_group_color(j_panel,:), 'LineStyle','none')
%     xlabel('Y (m)');
%     ylabel('Z (m)');
%     axis equal;
% 
%     xlim(yplot_range);
%     ylim(zplot_range);
%     title('YZ Alignment')
%     drawnow;
% 
%     subplot(2,2,3)
%     hold all;
%     plot(mean(sphere_center(j_panel,3,:)),mean(sphere_center(j_panel,1,:)),'Marker', '.','Color',probe_group_color(j_panel,:) ,'LineStyle','none')
%     xlabel('Z (m)');
%     ylabel('X (m)');
%     axis equal;
% 
%     xlim(zplot_range);
%     ylim(xplot_range);
%     title('ZX Alignment');
%     drawnow;
% end



LegendBar=legend(Legend_fill(:))

set(LegendBar, 'Position', [.7 .15 .15 .3])

figure(Sphere_Pos_fig)
subplot(2,2,1)

hold all;
plot(Center_Pos(1),Center_Pos(2),'Marker','+','Color',[0 0 0],'MarkerSize',20);
axis equal

figure(Sphere_Pos_fig)
subplot(2,2,2)
hold all;
plot(Center_Pos(2),Center_Pos(3),'Marker','+','Color',[0 0 0],'MarkerSize',20)
axis equal

figure(Sphere_Pos_fig)
subplot(2,2,3)
hold all;
plot(Center_Pos(3),Center_Pos(1),'Marker','+','Color',[0 0 0],'MarkerSize',20)
axis equal