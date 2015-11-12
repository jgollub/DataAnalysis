function [reduced_region_fields, reduced_region_positions, fid_physical_position, panel_type] = processSlottedMCPanel(fields, xx, yy, zz)

% Finding the Board
%%%%% Idea: start from the corners, pick the fiducials and zero pad the rest
board_size=0.11; % meters

%%% Left Bottom Corner
%%%% Find the indices for the board edges
value_y=-board_size/2;
value_z=-board_size/2;

[~,index_y_lb]=min(abs(yy(:)-value_y));
[~,index_z_lb]=min(abs(zz(:)-value_z));

% scatter3(0,yy(index_y_lb),zz(index_z_lb),100,'green','square','LineWidth',2)


%%% Right Bottom Corner
%%%% Find the indices for the board edges
value_y=board_size/2; 
value_z=-board_size/2;

[~,index_y_rb]=min(abs(yy(:)-value_y));
[~,index_z_rb]=min(abs(zz(:)-value_z));


% scatter3(0,yy(index_y_rb),zz(index_z_rb),100,'green','square','LineWidth',2)


%%% Left Top Corner
%%%% Find the indices for the board edges
value_y=-board_size/2;
value_z=board_size/2;

[~,index_y_lt]=min(abs(yy(:)-value_y));
[~,index_z_lt]=min(abs(zz(:)-value_z));


% scatter3(0,yy(index_y_lt),zz(index_z_lt),100,'green','square','LineWidth',2)

%%% Right Top Corner
%%%% Find the indices for the board edges
value_y=board_size/2;
value_z=board_size/2;

[~,index_y_rt]=min(abs(yy(:)-value_y));
[~,index_z_rt]=min(abs(zz(:)-value_z));


% scatter3(0,yy(index_y_rt),zz(index_z_rt),100,'green','square','LineWidth',2)

% Define coordinates for wrt board edges
if abs(index_y_lb)>abs(index_y_lt)
y_left_limit=index_y_lb;
else
y_left_limit=index_y_lt;
end

if abs(index_y_rb)>abs(index_y_rt)
y_right_limit=index_y_rb;
else
y_right_limit=index_y_rt;
end

if abs(index_z_lb)>abs(index_z_rb)
z_bottom_limit=index_z_rb;
else
z_bottom_limit=index_z_lb;
end

if abs(index_z_lt)>abs(index_z_rt)
z_top_limit=index_z_lt;
else
z_top_limit=index_z_rt;
end

% %plot outline
% line([0, 0 ,0, 0, 0],...
%     [yy(y_left_limit), yy(y_right_limit), yy(y_right_limit), yy(y_left_limit),yy(y_left_limit)], ...
%     [zz(z_top_limit),zz(z_top_limit),zz(z_bottom_limit),zz(z_bottom_limit),zz(z_top_limit)],'color','blue','LineWidth',1)

% Assign New Coordinates Enclosing the Panel
reduced_region_positions.y=yy((yy>yy(y_left_limit) & yy<yy(y_right_limit)) & (zz<zz(z_top_limit)& zz>zz(z_bottom_limit)));
reduced_region_positions.z=zz((yy>yy(y_left_limit) & yy<yy(y_right_limit)) & (zz<zz(z_top_limit)& zz>zz(z_bottom_limit)));
reduced_region_positions.x=xx((yy>yy(y_left_limit) & yy<yy(y_right_limit)) & (zz<zz(z_top_limit)& zz>zz(z_bottom_limit)));

reduced_region_fields=fields(yy>yy(y_left_limit) & yy<yy(y_right_limit) & (zz<zz(z_top_limit)& zz>zz(z_bottom_limit)));
 

hold on;  
scatter3(reduced_region_positions.x(:),reduced_region_positions.y(:),reduced_region_positions.z(:),40,20*log10(abs(reduced_region_fields(:))),'filled')
        axis normal; colormap('hot');set(gcf,'color','w');
        xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
        view(90,0); hold on;
 
% scatter3(fiducial.x(:),fiducial.y(:),fiducial.z(:),20,20*log10(abs(reduced_region)),'filled')
 % At this point we need to decide whether this is Tx or Rx panel
 %%%Adding vertically
 y_limit=.03;
 z_limit=.03;

 MC_vertical_test=sum(abs(reduced_region_fields(abs(reduced_region_positions.y)>y_limit)));
 
%%%Adding horizontally

 MC_horizontal_test=sum(abs(reduced_region_fields(abs(reduced_region_positions.z)>z_limit)));
 
 %%%% Choose = Tx or Rx

 mm=0.001;

 if MC_vertical_test>MC_horizontal_test
     panel_type='Rx';
         
     %fiducial left top
     fid_physical_position.lt.y= -18*mm;
     fid_physical_position.lt.z= 39.49*mm;
     %fiducial right top
     fid_physical_position.rt.y= 18*mm;
     fid_physical_position.rt.z= 39.49*mm;
     %fiducial left bottom
     fid_physical_position.lb.y= -10*mm;
     fid_physical_position.lb.z= -39.54*mm;
     %fiducial right bottom
     fid_physical_position.rb.y= 18*mm;
     fid_physical_position.rb.z= -39.54*mm;
 else
     panel_type='Tx';
     
      %fiducial left top
     fid_physical_position.lt.y= -39.49*mm;
     fid_physical_position.lt.z= 18*mm;
     %fiducial right top
     fid_physical_position.rt.y= 39.54*mm;
     fid_physical_position.rt.z= 18*mm;
     %fiducial left bottom
     fid_physical_position.lb.y= -39.49*mm;
     fid_physical_position.lb.z= -18*mm;
     %fiducial right bottom
     fid_physical_position.rb.y= 39.54*mm;
     fid_physical_position.rb.z= -10*mm;

 end