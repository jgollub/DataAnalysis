%% Reconstruct Full Scene

%panel positions
%At startup: initialize panels
debugVal=1
%%Choose parameters which panels to generate lookup panels for
sysConfig.numFeeds    =6; 
sysConfig.nProbes     =24;
sysConfig.nPanels     =12;
sysConfig.nFeeds      =6;
sysConfig.nFrequencies=96;
sysConfig.fIndices    =1:96;

% choosePanels = [A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4]
  choosePanels = logical([ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]);
  chooseProbes = logical(ones(24,1).');
  sysConfig.nPanels            = numel(find(choosePanels));
  

  %Measurement data
  Exp_data=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\alignment\depth_resolution_A_23-Mar-2015.mat');
  Exp_data.g =reshape(Exp_data.g,101,sysConfig.nFeeds,12,24);
  Exp_data.g =Exp_data.g(6:101,:,choosePanels,chooseProbes); %remove bad measurement data
  Exp_data.g =Exp_data.g(:);
  
  % Near Field Scan (NFS) data location for panels
filepathTx='C:\Users\Jonah Gollub\Documents\MetaImagerData\SemiRigid_NFS_scans\';
filepathRx='C:\Users\Jonah Gollub\Documents\MetaImagerData\Dipole_Antenna_Scans\';
folder='D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Panel Lookup Tables\Saved Data';

fprintf('NFS scan data: \n\t %s \nLookup table will be saved: \n\t %s\n',filepathTx, folder)


%% load panels


if debugVal
try figure(NFS_plot); catch NFS_plot=figure();
    scrsize=get(0,'ScreenSize');
    figsize=NFS_plot.Position(3:4);
    NFS_plot.Position(1:2)=[scrsize(3)-figsize(1),scrsize(4)-figsize(2)-80];
end
end
clf;
%load and set panel
fprintf('Creating panels...'); mark=clock;


counter=0;
for en=find(choosePanels)
    counter=counter+1;

    % determine file name
    [ pan_index,pan_letter]=ind2sub([4 3],en);
    filename=[char(64+pan_letter),num2str(pan_index)];
    FullPath=dir([filepathTx,filename,'.mat']);
    fprintf(['Loading Tx file: ',FullPath.name,'\n']);
    
    % load NFS data
    NFS_data=load([filepathTx,FullPath.name]);

    % Create dipole equivalent panel
    importedPanel = import_panel(NFS_data.data);

        %rotation
     importedPanel=panel_rotate(...p
        importedPanel,...
        Exp_data.panelRotation(en,1),...
        Exp_data.panelRotation(en,2),...
        Exp_data.panelRotation(en,3) ...
        );
    
        importedPanel=panel_offset(...
        importedPanel,...
        Exp_data.panelPosition(en,1),...
        Exp_data.panelPosition(en,2),...
        Exp_data.panelPosition(en,3) ...
        );
    

    %move panels but shift offset to experimental function
    for je=1:sysConfig.numFeeds
        fprintf('calculating for feed %d \n',je)    
         panels{je,en}=subpanel(importedPanel, je);
        panels(je,en)=frequency_clip(panels(je,en),6:101); %clip to number divisible by 4 and elim bad pts
        fprintf('Convert loaded data to dipoles data: %.3f mins\n', toc/60);
        
        % plot effective dipoles of panels for frequency 1
        if debugVal
            figure(NFS_plot)
            subplot(numel(find(choosePanels)),6,(en-1)*sysConfig.numFeeds+je)
            scatter(panels{je,counter}.y(:),panels{je,counter}.z(:),1,abs(panels{je,counter}.dipoles.y(51,:)));
            hold on
            feedLocation=locate(panels{je,counter});
            plot(feedLocation(2),feedLocation(3),'*r')
            text(feedLocation(2),feedLocation(3),[' ', num2str(je)],'Color','r')
            xlabel('Y'); ylabel('Z')
            axis xy; axis equal; axis tight;
            set(gca,'xdir','reverse')
            title(['Panel ',...
                char(64+pan_letter),num2str(pan_index),'_',num2str(je)])
            drawnow
        end
    end
end

sysConfig.freqs=panels{1,1}.f;

figure(10)
clf
hold on

% scatter3(subPanelPosition(:,1),subPanelPosition(:,2),subPanelPosition(:,3),30,'g');
panel_plot(figure(10),panels);
axis equal
drawnow

%% load probes

probes=cell(24,1);

NFS_plot_probes=figure(11);

%filenames

fileRx=cell(24,1);
for i=1:24
    fileRx{i}=sprintf('R%.2i',i);
end

%load data
counter=0;
for en=find(chooseProbes)
    counter=counter+1;
    
    FullPath=dir([filepathRx,fileRx{en},'.mat']);
    fprintf(['Loading Rx file: ',FullPath.name,'\n']);
    
    % load NFS data
    NFS_data=load([filepathRx,FullPath.name]);
    
    % Create dipole equivalent panel
    importedProbe = import_panel(NFS_data);
    importedProbe.type='probe';
    
       %rotation
    importedProbe=panel_rotate(...
        importedProbe,...
        Exp_data.probeRotation(en,1),...
        Exp_data.probeRotation(en,2),...
        Exp_data.probeRotation(en,3) ...
        );
    
    importedProbe=panel_offset(...
        importedProbe,...
       Exp_data.probePosition(en,1),...
       Exp_data.probePosition(en,2),...
       Exp_data.probePosition(en,3) ...
        );
    
 
    
    %move panels but shift offset to experimental function
    
    probes(en)=frequency_clip({importedProbe},sysConfig.fIndices); %clip to number divisible by 4 and elim bad pts
    fprintf('Convert loaded data to dipoles data: %.3f mins\n', toc/60);
    
    % plot effective dipoles of probes for frequency 1
    if debugVal 
        figure(NFS_plot_probes)
        subplot(6,4,counter)
        scatter(probes{counter}.y(:),probes{counter}.z(:),1,abs(probes{counter}.dipoles.y(51,:)));
        hold on
        feedLocation=locate(probes{counter});
        plot(feedLocation(2),feedLocation(3),'*r')
        text(feedLocation(2),feedLocation(3),[' ', num2str(en)],'Color','r')
        xlabel('Y'); ylabel('Z')
        axis xy; axis equal; axis tight;
        set(gca,'xdir','reverse')
        title(['Probe ', fileRx{i}])
        drawnow
    end
end

sysConfig.freqs=probes{1,1}.f;

figure(10)
hold on
panel_plot(figure(10), probes{:});
axis equal
drawnow


%% save import panel and probe data
save('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\PanelProbeVirtualizerData.mat',...
    '-v7.3','panels','probes')

%% fake target scene

% define the resolution for the the g 
voxel_size = .0025; % ~voxel size

depth_offset=mean(Exp_data.extent(1,:));
depth=Exp_data.extent(1,2)-Exp_data.extent(1,1);
width=Exp_data.extent(2,2)-Exp_data.extent(2,1);
height=Exp_data.extent(3,2)-Exp_data.extent(3,1);

%calculate imaging domain
[imgDomain_fg, ~] = test_space(depth_offset, depth, width, height, voxel_size);

%shift panels randomly to test algorithm [1 cm std]
rng(2);
% shift=0*randn(sysConfig.nPanels,3); 
shift=0.02*ones(sysConfig.nPanels,3); 
rng(2);
rotation=0*pi/180*randn(sysConfig.nPanels,3); 
shiftedPanels=panels;

for i=1:12
    for j=1:6
        shiftedPanels(j,i)=panel_offset(...
            panels(j,i),...
            shift(i,1),...
            shift(i,2),...
            shift(i,3) ...
            );
         shiftedPanels(j,i)=panel_rotate(...
            shiftedPanels(j,i),...
            rotation(i,1),...
            rotation(i,2),...
            rotation(i,3),...
             locate(shiftedPanels(1,i)));
%         panel_plot(figure(22),shift_panels_test(j,i));
%         drawnow
    end
end

% start with 10 points
numScatterers=10;
imgDomain_sigma=[];
minSeparation=.1;
rng(2);
% scattererIndex=randi(length(imgDomain),1,numScatteres);
for i=1:numScatterers
    x_range=-(min(imgDomain_fg(:,1))-max(imgDomain_fg(:,1))).*rand(1,1)+min(imgDomain_fg(:,1));
    y_range=-(min(imgDomain_fg(:,2))-max(imgDomain_fg(:,2))).*rand(1,1)+min(imgDomain_fg(:,2));
    z_range=-(min(imgDomain_fg(:,3))-max(imgDomain_fg(:,3))).*rand(1,1)+min(imgDomain_fg(:,3));
    scattererPositions(i,:)=[x_range, y_range,z_range];
    
    if i>1
        sct_distance=sum(bsxfun(@minus,scattererPositions(1:i-1,:),scattererPositions(i,:)).^2,2).^(1/2);
        %values for f of small region (log)
        
        while any(sct_distance<minSeparation);
            x_range=-(min(imgDomain_fg(:,1))-max(imgDomain_fg(:,1))).*rand(1,1)+min(imgDomain_fg(:,1));
            y_range=-(min(imgDomain_fg(:,2))-max(imgDomain_fg(:,2))).*rand(1,1)+min(imgDomain_fg(:,2));
            z_range=-(min(imgDomain_fg(:,3))-max(imgDomain_fg(:,3))).*rand(1,1)+min(imgDomain_fg(:,3));
            scattererPositions(i,:)=[x_range, y_range,z_range]
            sct_distance=sum(bsxfun(@minus,scattererPositions(1:i-1,:),scattererPositions(i,:)).^2,2).^(1/2)
        end
        
    end
    
end

%in order to maintain normalization between experiment and virtualized g
imgDomain_virtualizer=zeros(2*numScatterers,3);

imgDomain_dummy(:,:,1)= scattererPositions;
imgDomain_dummy(:,:,2)= bsxfun(@plus,scattererPositions,[R_l, R_l,R_l]);
imgDomain_virtualizer(1:2:2*numScatterers,:)=imgDomain_dummy(:,:,1);
imgDomain_virtualizer(2:2:2*numScatterers,:)=imgDomain_dummy(:,:,2);

imgDomain_sigma=zeros(length(imgDomain_virtualizer),1);
imgDomain_sigma(1:2:2*numScatterers)=1;

g=forward_model(probes, shiftedPanels(:), imgDomain_sigma, imgDomain_virtualizer);
g=permute(g,[3 2 1]);
g=g(:);
panelLocations=locate(shiftedPanels(1,:));
PerceivedPanelLocations=locate(panels(1,:));

%save backup g
save('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\LastSimulatedg.mat',...
    'g',...
    'panelLocations',...
    'shift',...
    'rotation',...
    'PerceivedPanelLocations',...
    'scattererPositions')

%% %%% Reconstruct over all panel/probes ====================================


%% Setup the imaging domain
  
% define the imaging domain
R_l = .007; % ~voxel size

depth_offset=mean(Exp_data.extent(1,:));
depth=Exp_data.extent(1,2)-Exp_data.extent(1,1);
width=Exp_data.extent(2,2)-Exp_data.extent(2,1);
height=Exp_data.extent(3,2)-Exp_data.extent(3,1);

%calculate imaging domain
[imgDomain, domainGrid] = test_space(depth_offset, depth, width, height, R_l);

%% Generate lookup tables
% do this at system startup

%Resolutions Sampling Values for panel
  dR_panel     = 0.2; % meters
  dAz_panel    = 0.03; %radians
  dEl_panel    = 0.025; %radians
   NFS_distance = 0.086; %meters

panelTableRes = [dR_panel, dAz_panel, dEl_panel]; %Radial/ Azimuth/ Elevation 
panelTables      = makeLookUpTablesExperimental(panels(:,choosePanels), NFS_distance,imgDomain, panelTableRes);

probeTableRes = [.2, .1, .1]; %can be reduced for probes
probeTables   = makeLookUpTablesExperimental(probes, NFS_distance,imgDomain, probeTableRes);

% save lookup tables
save('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\LastBuiltLookupTables.mat','panelTables','probeTables');

%% load previous data panels, g, lookup tables etc

% load panel probe data from file (preprocessed near field scan data)
load('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\PanelProbeVirtualizerData.mat');
sysConfig.freqs=panels{1,1}.f;
% load last simulated g
load('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\LastSimulatedg.mat');

% load last lookup tables
load('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\LastBuiltLookupTables.mat');

%% Scene Pre-Calculations =================================================
% choosePanels = [A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4]
  choosePanels = logical([ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]);
  chooseProbes = logical(ones(24,1).');
  sysConfig.nPanels            = numel(find(choosePanels));
  


nLevels =4; % nearest neighbor subdomains used in interpolation, can be 1, 2, 3, or 4
delta   =4; % error factor: small delta -> more accurate, big delta -> faster
% get appropriate subdomain and frequency bin sizes, can override if
% desired 
[partitioning.R_k, ~, partitioning.S, partitioning.T] =getDomainSizes(probes, panels, imgDomain, delta);

R_k=0.04;
% Restructure imaging domain for recon, divvy the fine imaging grid by subdomain
fprintf('Breaking imaging domain into subdomains...'); mark = clock();
partitioning.scene_cfg0     =prepScene(imgDomain, R_k, R_l, nLevels);
partitioning.large_k        =extract_subs(partitioning.scene_cfg0.locs_kl(1,:));
fprintf('finished in %0.2f seconds\n', etime(clock(),mark));

%% recon full region
mark=clock();
[f_mf,indices, locs, timing] = reconstructionSelectPanels(sysConfig,...
                                            g,... %Exp_data.g,...
                                            imgDomain,...
                                            probes, probeTables,chooseProbes,...
                                            panels, panelTables,choosePanels,...
                                            partitioning,...
                                            R_k...
                                            );
                                        
partitioning.indices=indices;         
partitioning.locs=locs;
fprintf('[Done] %f(s)\n', etime(clock(),mark));
%% Plot Results ===========================================================
fprintf('Plotting results... '); 
fg =figure(29);

fg.Position=[37    90   561   840];
subplot(2,1,1)
cla
plotScatters(f_mf, indices,imgDomain, domainGrid, 'Type','linear');
% f_plot=volume_plot2(fg, imgDomain, f_mf, indices, domainGrid, 'log',2); 
title('look-up table')
view(-33,16)
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
alpha('color');

hold on;
scatter3(scattererPositions(:,1),scattererPositions(:,2),scattererPositions(:,3),4,'r')

% Filter and Threshold ===================================================
figure(fg);
subplot(2,1,2)
cla

f_plot=plotScatters(f_mf, indices,imgDomain, domainGrid, 'Type','log');
hold on
scatter3(scattererPositions(:,1),scattererPositions(:,2),scattererPositions(:,3),4,'r')

% Locate scatters & calculate reduced region to reconstruct over
figure(fg);
subplot(2,1,2)

num_points=10;
subImgDomain=cell(length(num_points));
lookWithinRadius=0.03;

f_plot_dummy=f_plot;

solvedGlobalPosition=zeros(num_points,3);

%estimate position of target scatteres
[solvedGlobalPosition, subImgDomain, subIndices]=filterScatterers(f_mf,...
               indices, ...
               imgDomain, ...
               domainGrid,...
               15,...
               num_points,...
               lookWithinRadius ...
               ); 
           
 partitioning.subIndices=subIndices;
 
%plot positions of scatterers that have been solved for globally
 scatter3(solvedGlobalPosition(:,1),solvedGlobalPosition(:,2),solvedGlobalPosition(:,3),4,'black')
for i=1:10
    subplot(2,1,2)
 smallReconstructionRegionPlot(figure(fg), solvedGlobalPosition(i,:),lookWithinRadius,num2str(i))
end


 %sort positions of the scatterers such that they can be compared agains
 %the real positions (for fake_g)
 A1=repmat(scattererPositions,1,1, ...
     length(solvedGlobalPosition));
 for i=1:length(scattererPositions)
     B1=A1(:,:,i)-repmat(solvedGlobalPosition(i,:),length(solvedGlobalPosition),1);
     [~, sortIndex(i)]=min(arrayfun(@(idx) norm(B1(idx,:)), 1:size(B1,1)));
 end
diff=scattererPositions(sortIndex,:)-solvedGlobalPosition;
diff=diff(~isnan(solvedGlobalPosition(:,1)),:); 
MSE_pts=sum(sum((diff).^2,2),1)/length(scattererPositions)

 %% recon with makeH directly calculated H
 
 collected_small_imgDomain=cat(1,subImgDomain{:});
 subIndices=[];
 counter=0;
 for i=1:numel(subImgDomain)
 subIndices{i}=(1+counter):(counter+length(subImgDomain{i}))
 counter=counter+length(subImgDomain{i})
 end
 
 panelFields = dipoles_to_fieldsEXP3(panels, collected_small_imgDomain);
 
 probeFields = dipoles_to_fieldsEXP3(probes, collected_small_imgDomain);

 H=makeH_faceted(probeFields,panelFields);
 
 f_abs_est=cgs(H'*H,H'*g,1e-3,30);
 figure(37)
 scatter3(collected_small_imgDomain(:,1),...
          collected_small_imgDomain(:,2),...
          collected_small_imgDomain(:,3),abs(f_abs_est).^2./max(abs(f_abs_est).^2))
hold on

scatter3(scattererPositions(:,1),scattererPositions(:,2),scattererPositions(:,3),4,'r')

for i=1:num_points
     %solve for max position using centroid
      x_centroid=sum(subImgDomain{i}(:,1).*abs(f_abs_est(subIndices{i})).^2)/sum(abs(f_abs_est(subIndices{i})).^2);
      y_centroid=sum(subImgDomain{i}(:,2).*abs(f_abs_est(subIndices{i})).^2)/sum(abs(f_abs_est(subIndices{i})).^2);
      z_centroid=sum(subImgDomain{i}(:,3).*abs(f_abs_est(subIndices{i})).^2)/sum(abs(f_abs_est(subIndices{i})).^2);
      solvedGlobalPosition_absH(i,:)=[x_centroid,y_centroid,z_centroid];
end   
hold on
  scatter3(solvedGlobalPosition_absH(:,1),solvedGlobalPosition_absH(:,2),solvedGlobalPosition_absH(:,3),4,'blue')

%% Reconstruct for single Tx/Rx and all Rx/Tx; find position of scatters; cycle
%  through all and adjust positions
figure(30)
clf


results=struct([]);

% solved_panel_position=Exp_data.panelPosition;
spottedPosition=zeros(12,3,num_points);
panelsPrimed=panels;


% for i=1:2
%     for j=1:6
%         panelsPrimed(j,i)=panel_offset(...
%             panelsPrimed(j,i),...
%             .03,...
%             .03,...
%             .03 ...
%             );
%          shiftedPanels(j,i)=panel_rotate(...
%             shiftedPanels(j,i),...
%             rotation(i,1),...
%             rotation(i,2),...
%             rotation(i,3),...
%              locate(shiftedPanels(1,i)));
% %         panel_plot(figure(22),shift_panels_test(j,i));
% %         drawnow
%     end
% end

MSE_start=sum(sum((PerceivedPanelLocations-panelLocations).^2,2),1)/length(panelLocations);
MSE_solved=[];
for iterate=1:10%3
    
    Mean_error(iterate,:)= mean(abs(locate(panelsPrimed(1,:))-locate(shiftedPanels(1,:))));
    for i_panel=1:2%12
        fprintf('panel %i\n',i_panel)
        %select panels
        choosePanels=logical(zeros(12,1).');
        choosePanels(i_panel)=logical(1);
        chooseSubpanels=cat(1,choosePanels,choosePanels,choosePanels,choosePanels,choosePanels,choosePanels);
        chooseSubpanels=chooseSubpanels(:); %choosing one panel (six feeds at time)
        for j_region=1:num_points
 
% p=ones(10,1);
% p(j_region)=0;
% p=logical(p);
% clear dummy
% dummy(:,:,1)= scattererPositions(sortIndex(p),:);
% dummy(:,:,2)= bsxfun(@plus,dummy(:,:,1),[R_l, R_l,R_l]);
% imgDomain_virtualizer(1:2:2*(numScatterers-1),:)=dummy(:,:,1);
% imgDomain_virtualizer(2:2:2*(numScatterers-1),:)=dummy(:,:,2);
% imgDomain_sigma=zeros(length(imgDomain_virtualizer),1);
% imgDomain_sigma(1:2:2*numScatterers)=1;
% 
% t=ones(12,1);
% t(i_panel)=0;
% t=logical(t);
% g_subtract=forward_model(probes, shiftedPanels(:), imgDomain_sigma, imgDomain_virtualizer);
% g_subtract=permute(g_subtract,[3 2 1]);
% gprime=g-g_subtract(:);

            [f_mf, subIndices{j_region},subImgDomain{j_region}, timing] = reconstructionSelectPanels(sysConfig,...
                g,... Exp_data.g,...
                subImgDomain{j_region},...
                probes, probeTables,chooseProbes,...
                panelsPrimed, panelTables(chooseSubpanels),choosePanels,...
                partitioning,...
                R_k,...
                'Reconstruction','Direct','inputIndices',subIndices{j_region});
%             subIndices{j_region}
            % find best guess at target position
           % subplot(num_points,sysConfig.nPanels,i_panel+sysConfig.nPanels*(j_region-1))
            [spottedPosition(i_panel,:,j_region), ~] = filterScatterers(f_mf,...
                                                            subIndices{j_region},... subIndices,...
                                                            imgDomain,...
                                                            domainGrid,...
                                                            15, ...
                                                            1, ...
                                                            lookWithinRadius ...
                                                            );                                                    
        figure(30)
        cla
        hold on
              scatter3(spottedPosition(i_panel,1,j_region),...
                    spottedPosition(i_panel,2,j_region),...
                    spottedPosition(i_panel,3,j_region),...
                    30,'green')

                  scatter3(solvedGlobalPosition(j_region,1),...
                    solvedGlobalPosition(j_region,2),...
                    solvedGlobalPosition(j_region,3),...
                    30,'red')
                    smallReconstructionRegionPlot(figure(30),...
                        solvedGlobalPosition(j_region,:),...
                        lookWithinRadius,num2str(j_region))
                    plotScatters(f_mf,subIndices{j_region}, ...
                        imgDomain, domainGrid, 'Type','log');
                    
                scatter3(scattererPositions(:,1),scattererPositions(:,2),scattererPositions(:,3),4,'k')
                    
       xlim([min(imgDomain(subIndices{j_region},1)),max(imgDomain(subIndices{j_region},1))]);
       ylim([min(imgDomain(subIndices{j_region},2)),max(imgDomain(subIndices{j_region},2))]);
       zlim([min(imgDomain(subIndices{j_region},3)),max(imgDomain(subIndices{j_region},3))]);

            axis equal
                  drawnow
%                     scatter3(imgDomain(subIndices{j_region},1),...
%                         imgDomain(subIndices{j_region},2),...
%                         imgDomain(subIndices{j_region},3),3,'green')
%                     scatter3(subLocs(:,1),subLocs(:,2),subLocs(:,3),3,'blue')
%                      figure(28)
%   
%         [f_plot] = plotScatters(f_mf,subIndices, imgDomain, domainGrid, 'Type','log');
%         hold on
%         scatter3(spottedPosition(i_panel,1,j_region),...
%                     spottedPosition(i_panel,2,j_region),...
%                     spottedPosition(i_panel,3,j_region),...
%                     30,'green')
        smallReconstructionRegionPlot(figure(29), solvedGlobalPosition(j_region,:),lookWithinRadius,num2str(j_region))

            %check that the spotted position is in the expected range or set to NaN
        if  norm(spottedPosition(i_panel,:,j_region)-solvedGlobalPosition(j_region,:))>lookWithinRadius
            %pause
            recordAllSpottedPosition(i_panel,:,j_region)=spottedPosition(i_panel,:,j_region) %
            spottedPosition(i_panel,:,j_region)=NaN;
        end
        end
    end
    
    
    spottedPosition
    
    figure(31)
    hold on;
    for i=1:num_points
        marker_color=(1:12)*3;
        marker_size=500*ones(12,1);
        scatter3(spottedPosition(:,1,i),...
            spottedPosition(:,2,i),...
            spottedPosition(:,3,i),...
            marker_size,...
            'k',...
            'x'...
            );
        
        %calculate variation of located positions for each point
        temp_region=spottedPosition(:,1,i);
        x_std=std(temp_region(~isnan(temp_region)));
        temp_region=spottedPosition(:,2,i);
        y_std=std(temp_region(~isnan(temp_region)));
        temp_region=spottedPosition(:,3,i);
        z_std=std(temp_region(~isnan(temp_region)));
        fprintf('for Point %i ..., \nX variance is: %.5f\nY variance is: %.3f\nZ variance is: %.3f\n'...
            ,i, x_std, y_std,z_std)
    end
    scatter3(spottedPosition(:,1),spottedPosition(:,2),spottedPosition(:,3),4,'r') %!!!!!
    %% Calculate shift of panels and move ======================================
    %determine guess for how panels should be moved to improve image
    
    K_U=[];
    K_r=[];
    lrms=[];
     
    pts_global_exist=logical(all(~isnan(solvedGlobalPosition(:,1)),2)); 
    for i=1:sysConfig.nPanels
        %remove NaN valued points
        pts_local_exist=logical(squeeze(all(~isnan(spottedPosition(i,:,:)))));
        use_pts=logical(pts_local_exist.*pts_global_exist);
        %    solve for trransformation using Kabsch alg and also set LS_weighting----haven't yet
        [K_U(:,:,i) K_r(:,i) lrms(i)]=Kabsch(...
            permute(spottedPosition(i,:,use_pts),[3 2 1]).',...
            solvedGlobalPosition(use_pts,:).'...
            );
        
          %check to see that nothing is too far out of bounds, and throw it out if it's s
       if (norm(K_r(:,i))> .06) | (trace(K_U(:,:,i))<2.5)
         K_r(:,i)=0;
            fprintf('movement out of bounds for panel %i, set to zero\n',i)
            K_U(:,:,i)=eye(3);
            fprintf('rotation out of bounds for panel %i, set to identity\n',i)
        end
    end
 results(iterate).K_U=K_U;
 results(iterate).K_r=K_r;
    % move copy of panels (move_panels) to new positions
   
    %move panels
    for i=1:sysConfig.nPanels
        for j=1:sysConfig.nFeeds
            panelsPrimed(j,i)=transformPanelArbitraryUr(panelsPrimed(j,i),K_U(:,:,i),K_r(:,i));
        end
    end
    solvedPanelLocations=locate(panelsPrimed(1,:))
    %plot MSE
    disp('MSE (position only)')
    results(iterate).MSE_solved=sum(sum((solvedPanelLocations-panelLocations).^2,2),1)/length(panelLocations)
    figure(26)
    plot(1:numel(results(:).MSE_solved),[results(:).MSE_solved])
    title('MSE Vs. Iteration (Position Only)')
    xlabel('Iteration')
    ylabel('MSE')

figure(29)
subplot(2,1,2)
hold on
     scatter3(panelLocations(:,1),panelLocations(:,2),panelLocations(:,3),4,'black')
     scatter3(solvedPanelLocations(:,1),solvedPanelLocations(:,2),solvedPanelLocations(:,3),4,'red')
     scatter3(PerceivedPanelLocations(:,1),PerceivedPanelLocations(:,2),PerceivedPanelLocations(:,3),4,'blue')
axis equal;
     %     panel_plot(figure(11),panelsPrimed);
%     axis equal
%     drawnow
%     hold on
%     panel_plot(figure(11),panels);
    
    %% Test reconstruction of new image =========================================
    choosePanels = logical([ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]);
    chooseProbes = logical(ones(24,1).');
    
    % Scene Pre-Calculations =================================================
    nLevels =4; % nearest neighbor subdomains used in interpolation, can be 1, 2, 3, or 4
    delta   =4; % error factor: small delta -> more accurate, big delta -> faster
    % get appropriate subdomain and frequency bin sizes, can override if
    % desired
    [partitioning.R_k, ~, partitioning.S, partitioning.T] =getDomainSizes(probes, panelsPrimed, imgDomain, delta);
    
    R_k=0.04;
    % Restructure imaging domain for recon, divvy the fine imaging grid by subdomain
    fprintf('Breaking imaging domain into subdomains...'); mark = clock();
    partitioning.scene_cfg0     =prepScene(imgDomain, R_k, R_l, nLevels);
    partitioning.large_k        =extract_subs(partitioning.scene_cfg0.locs_kl(1,:));
    fprintf('finished in %0.2f seconds\n', etime(clock(),mark));
    
    [f_mf, Indices, locs, timing] = reconstructionSelectPanels(sysConfig,...
        g,... Exp_data.g,...
        imgDomain,...
        probes, probeTables,chooseProbes,...
        panelsPrimed, panelTables,choosePanels,...
        partitioning,...
        R_l...
        );
     figure(35)
     cla
    [f_plot] = plotScatters(f_mf,indices,imgDomain, domainGrid, 'Type','linear');
    
     hold on
   % scatter3(scattererIndex(:,1),scattererIndex(:,2),scattererIndex(:,3),4,'r')
    drawnow
    
    %% Locate scatters & calculate reduced region to reconstruct over
    f_plot_dummy=f_plot;
 
solvedGlobalPosition=zeros(num_points,3);

[solvedGlobalPosition, subImgDomain, subIndices] = filterScatterers(f_mf,...
                                        indices,...
                                        imgDomain,...
                                        domainGrid,...
                                        15, ...
                                        num_points,...
                                        lookWithinRadius ...
                                        );
figure(29)
 scatter3(solvedGlobalPosition(:,1),solvedGlobalPosition(:,2),solvedGlobalPosition(:,3),10,'black')
end