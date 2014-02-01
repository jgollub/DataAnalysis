%reconstruct image from each panel of 1-12 panel array given measurement
%Target Radius
R_sph=.1262; %m
Sphere_Pos_fig=figure(99);

xplot_range=[-.1 .2];
yplot_range=[.2 .5];
zplot_range=[-1.4,-.8];

sphere_center=[];
%parse data. each chunck of H matrix-> 101 freqs X 6 switches= 606 measurement modes per panel
freqs=101;
RF_pathPerPanel=6;
Num_Probes=6;

% rf data
H = scene_data.H;
g = scene_data.obj_saved(1).measurement;
Z = scene_data.Z;
Az = scene_data.Az;
El = scene_data.El;

%specific slice to image
upsample=1;
Panel_Probe_combo=1:1;
slices=1:6;

%Load Panel Positions
load('Panel_Positions.mat');

%load probe positions
load('Probe_Positions.mat');

%size of H matrix
Num_Panels=(size(scene_data.H,1)/(freqs*RF_pathPerPanel*Num_Probes))
Num_figs=Num_Panels/2

%kinect scene data
if isfield(scene_data,'rgb')
    sn = 1;
    rgb = scene_data(sn).rgb;
    xyz = scene_data(sn).xyz;
    azimuth_min = scene_data(sn).Az_extent(1);
    azimuth_max = scene_data(sn).Az_extent(2);
    elevation_min = -scene_data(sn).El_extent(1);
    elevation_max = -scene_data(sn).El_extent(2);
    z_min = scene_data(sn).Z_extent(1);
    z_max = scene_data(sn).Z_extent(2);
    objs = scene_data(sn).objs;
else
    fprintf('%s\n%s\n','Sorry, there''s no Kinect data in this dataset :(.','We can still look at the RF data though.')
end


%figure to plot subfigures into

Obj_funct_Plot=figure(36)
for j_probe=1:6;
    for i_panel=1:12;
        
        figure(semilogy)
        hold on;
        [obj obj_debias objfunc]=singlePanel_recon(H,g,Az,El,Z,Num_Panels,freqs,j_probe,i_panel);
        
        semilogy(abs(objfunc(2:end)-objfunc(1:(end-1))))
        
        xlabel('twist iteration')
        ylabel('magnitude of objective function change')
        drawnow
        
        %% reshape reconstructed image vector into 3D scene, plot
        obj = abs(obj);
        nz = length(Z);
        nel = size(El,1);
        naz = size(Az,2);
        nae = nel*naz;
        obj3D = zeros(nel,naz,nz);
        for zn=1:nz
            objd = obj((1:nae)+nae*(zn-1));
            obj3D(:,:,zn) = reshape(objd,nel,naz);
        end
        
        %%smooth object before testing for max region
        %Smoothed_Obj=conv2(single(obj3D(:,:,slices)),fspecial('gaussian',7, 1),'same');
        %[maxplot_val max_index]=max(max(max(obj3D(:,:,slices))))
        
        [maxplot_val max_index]=max(max(max(obj3D(:,:,slices))))
        recon_fig=figure(35);
        for i_slice=slices
            subplot(ceil(length(slices)/2),2,i_slice);
            set(gca,'YDir','reverse')
            % subplot(2,Num_figs,plotpos);
            Az_range=tan(Az(1,:))*Z(i_slice);
            El_range=tan(El(:,1))*Z(i_slice);
            Recon_obj=upsample_image(obj3D(:,:,i_slice),upsample);
            
            imagesc(Az_range,El_range,Recon_obj,[0 maxplot_val]);
            axis xy
            axis equal
            axis tight
            
            
            title(['slice num: ',num2str(i_slice),'; Z_pos', num2str(Z(i_slice))]);
            xlabel('Az (m)');
            ylabel('El (m)');
            
            %solve for max reflection
            if i_slice==max_index
                Spot_loc=LocateReflPeak(Az_range, El_range,Recon_obj,recon_fig);
                if length(Spot_loc(:,1))>1
                    error('more than one max found')
                end
                reflect_pos(i_slice,:)=Spot_loc.';
                
                %specular reflection of sphere
                Panel_vec=[PanelPositions(i_panel), PanelPositions(i_panel+Num_Panels), -PanelPositions(i_panel+2*Num_Panels)]; % note: -z to correct for left-handed coordinate system
                Probe_vec=[ProbePos(j_probe), ProbePos(j_probe+Num_Probes), -ProbePos(j_probe+2*Num_Probes)]; %
                
                Linc_vec=[Spot_loc -Z(i_slice)]-Panel_vec; %(-Z(slice)) COR for LH System
                Lref_vec=[Spot_loc -Z(i_slice)]-Probe_vec;% (-Z(slice)) COR for LH System
                
                norm_vec=-(Linc_vec+Lref_vec)./norm(Linc_vec+Lref_vec);
                
                %center of sphere
                S_vec=[Spot_loc -Z(i_slice)]-R_sph*norm_vec; %(-Z(slice)) COR for LH System
            end
            
        end
        
        %Create colorbar
        
        B=colorbar;
        caxis([0,maxplot_val]);
        set(B, 'Position', [.8314 .11 .0581 .8150])
        for i=1:length(slices)
            pos=get(subplot(ceil(length(slices)/2),2,i), 'Position');
            axes(subplot(ceil(length(slices)/2),2,i))
            set(subplot(ceil(length(slices)/2),2,i), 'Position', [pos(1)-.1 pos(2) pos(3) pos(4)])
        end
        
        %save calculated sphere positions
        sphere_center(i_panel,1:3,j_probe)=S_vec;
    end
    
    
    % plot Sphere calculated positions
    
    figure(Sphere_Pos_fig)
    hold all;
    probe_group_color=[1-j_probe/Num_Probes,1-j_probe/Num_Probes,0];
    
    subplot(2,2,1)
    hold all;
    plot(sphere_center(:,1,j_probe),sphere_center(:,2,j_probe),'Marker', '.','Color',probe_group_color ,'LineStyle','none')
    xlabel('X (m)');
    ylabel('Y (m)');
    xlim(xplot_range);
    ylim(yplot_range);
    title('XY Alignment')
    drawnow;
    
    subplot(2,2,2)
    hold all;
    plot(sphere_center(:,2,j_probe),sphere_center(:,3,j_probe),'Marker', '.','Color',probe_group_color, 'LineStyle','none')
    xlabel('Y (m)');
    ylabel('Z (m)');
    xlim(yplot_range);
    ylim(zplot_range);
    title('YZ Alignment')
    drawnow;
    
    subplot(2,2,3)
    hold all;
    plot(sphere_center(:,3,j_probe),sphere_center(:,1,j_probe),'Marker', '.','Color',probe_group_color ,'LineStyle','none')
    xlabel('Z (m)');
    ylabel('X (m)');
    xlim(zplot_range);
    ylim(xplot_range);
    title('ZX Alignment');
    drawnow;
    
    %suptitle('Reconstructed Center of Sphere')
    
    % set(gcf,'NextPlot','add');
    % axes;
    % htitle=title('Individual Panel Reconstructions');
    % set(gca,'Visible', 'off');
    % set(htitle,'Visible', 'on');
    %plot_range_slices(obj3D,Az,El,Z,3)
    % figure(34)
    % set(gcf,'NextPlot','add')
    % hold on;
    % if i==1
    %     sumPlot=upsample_image(obj3D(:,:,slice),upsample)./Num_Panels;
    % elseif i~1
    %     sumPlot=sumPlot+upsample_image(obj3D(:,:,slice),upsample)./Num_Panels;
    
    % imagesc(tan(Az(1,:))*Z(slice),tan(El(:,1))*Z(slice),sumPlot);
    % axis xy
    % title('Sum of Individual Panel Reconstructions')
    % xlabel('X (m)');
    % ylabel('Y (m)');
end
