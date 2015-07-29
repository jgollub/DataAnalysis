% Reconstruct Full Scene

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
  Exp_data=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Measured Scenes\hayrettin\g_mannequin_clean.mat');
  Exp_data.g =reshape(Exp_data.g,101,sysConfig.nFeeds,12,24);
  Exp_data.g =Exp_data.g(6:101,:,choosePanels,chooseProbes); %remove bad measurement data
  Exp_data.g =Exp_data.g(:);
  
  % Near Field Scan (NFS) data location for panels
filepathTx='C:\Users\Jonah Gollub\Documents\MetaImagerData\SemiRigid_NFS_scans\';
filepathRx='C:\Users\Jonah Gollub\Documents\MetaImagerData\Dipole_Antenna_Scans\';
folder='D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Panel Lookup Tables\Saved Data';

fprintf('NFS scan data: \n\t %s \nLookup table will be saved: \n\t %s\n',filepathTx, folder)

%% load previous data panels, g, lookup tables etc

% load panel probe data from file (preprocessed near field scan data)
load('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\PanelProbeVirtualizerData.mat');
sysConfig.freqs=panels{1,1}.f;
% load last simulated g
load('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\LastSimulated_Mannequin.mat');

% load last lookup tables
load('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\LastBuiltLookupTables.mat');


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
     importedPanel=panel_rotate(...
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

%% Build Targets ==========================================================

x=Exp_data.xGrid(Exp_data.roi(:,1));
y=Exp_data.yGrid(Exp_data.roi(:,2));
z=Exp_data.zGrid(Exp_data.roi(:,3));
extent=[min(x),max(x);min(y),max(y);min(z),max(z)]; 
offset=[mean(extent(1,:)),mean(extent(2,:)),mean(extent(3,:))];
% define the resolution for the the g 
voxel_size = .0025; % ~voxel size

% Parameters for the target
rescale   =1;
scene_res =0.002;
off.x     =1;
off.y     =0;
off.z     =offset(3)-0.05;

% Parameters for the target
rescale   =1;
scene_res =0.002;
off.x     =1;
off.y     =0;
off.z     =offset(3)-0.05;
% Create mannequin target
testMan          = objectCreator('Man_001',0.9,rescale);
% testMan1.vert     = testMan1.vert;
% testMan           = testMan1;
% testMan.vert(:,1) = -1*testMan1.vert(:,1); 
% testMan.vert(:,2) = testMan1.vert(:,3); 
% testMan.vert(:,3) = -1*testMan1.vert(:,2); 
testMan           = objectMover(testMan,off);
% % Rotation  matrix
% theta = -7;
% RX    = testMan.vert(:,1)*cosd(theta)-testMan.vert(:,3)*sind(theta);
% RZ    = testMan.vert(:,1)*sind(theta)+testMan.vert(:,3)*cosd(theta);
% testMan.vert(:,1) = RX;
% testMan.vert(:,3) = RZ;

targetMan  = ZBuffer(testMan, scene_res);
targetMan0  = ZBuffer(testMan,  0.005);

target.locs   = [targetMan.locs];
target.sigma  = [(targetMan.sigma)'];
target.sigma  = target.sigma';

target0.locs   = [targetMan0.locs];
target0.sigma  = [(targetMan0.sigma)'];
target0.sigma  = target0.sigma';

% Collects objecs in the scene
objects = testMan;

% Experimental extent
x=target0.locs(:,1);
y=target0.locs(:,2);
z=target0.locs(:,3);
extent=[min(x),max(x);min(y),max(y);min(z),max(z)];
range=[extent(1,2)-extent(1,1),extent(2,2)-extent(2,1),extent(3,2)-extent(3,1)];
offset=[mean(extent(1,:)),mean(extent(2,:)),mean(extent(3,:))];

% figure();
% scatter3(target0.locs(:,1),target0.locs(:,2),target0.locs(:,3),5)
%% Define ROI =============================================================
fprintf('Selecting ROI... '); mark=clock();
R_l    =0.007;
% R_r    =0.04;
padding=R_r/sqrt(2);
S=2;

% Fine
[imgDomain,domainGrid]=test_space2(offset(1),range(1),offset(2),range(2),offset(3),range(3),R_l);
imgDomain=imgDomain.locs;

% scatter3(imgDomain(:,1),imgDomain(:,2),imgDomain(:,3),5)
% panel_plot(figure(),probes)
g=forward_model({probes{:}}, {panels{:}}, target.sigma, target.locs );
g0 = vec(permute(g, [3 2 1]));
g=g0;
%save backup g
save('D:\MetaImager\Virtualizer Data\Imported Panels and Probes\LastSimulated_Mannequin.mat',...
    'g');


%% %%% Reconstruct over all panel/probes ====================================




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

R_k=0.04
% Restructure imaging domain for recon, divvy the fine imaging grid by subdomain
fprintf('Breaking imaging domain into subdomains...'); mark = clock();
partitioning.scene_cfg0     =prepScene(imgDomain, R_k, R_l, nLevels);
partitioning.large_k        =extract_subs(partitioning.scene_cfg0.locs_kl(1,:));
fprintf('finished in %0.2f seconds\n', etime(clock(),mark));

%% recon full region
mark=clock();
[f_mf,indices, locs, timing] = reconstructionSelectPanels(sysConfig,...
                                            g,... %Exp_data.g,...
                                            target.locs,... imgDomain,...
                                            probes, probeTables,chooseProbes,...
                                            panels, panelTables,choosePanels,...
                                            partitioning,...
                                            R_k...
                                            );
                                        
partitioning.indices=indices;         
partitioning.locs=locs;
fprintf('[Done] %f(s)\n', etime(clock(),mark));

recon_fig1=figure();
volume_plot(recon_fig1,imgDomain, f_mf, indices, domainGrid, 10); title('look-up table')
plotScatters(f_mf, indices,imgDomain, domainGrid, 'Type','log');


%% ==============================
sysConfig.numFeeds    =6; 
sysConfig.nProbes     =24;
sysConfig.nPanels     =12;
sysConfig.nFeeds      =6;
gtest=reshape(g,sysConfig.nFrequencies,sysConfig.nPanels*sysConfig.nFeeds,sysConfig.nProbes);   

gmean=squeeze(mean(abs(gtest).^2,1))
max_gmean=max(max(gmean))
dB_gmean=10*log10(gmean/max_gmean);
threshold=-5;
dB_gmean(dB_gmean>threshold);

[I,J]=ind2sub(size(dB_gmean),find(dB_gmean>threshold));


g_threshold=[];
for i=1:length(I)
g_threshold(1+(i-1)*sysConfig.nFrequencies:i*sysConfig.nFrequencies)=gtest(:,I(i),J(i));
end


chooseProbes_Threshold=zeros(24,1);
chooseProbes_Threshold(unique(J))=1;
chooseProbes_Threshold=logical(chooseProbes_Threshold);

choosepanels_Threshold=zeros(72,1);
choosepanels_Threshold(unique(I))=1;
choosepanels_Threshold=logical(choosepanels_Threshold);
%% recon full region
mark=clock();
[f_mf,indices, locs, timing] = reconstructionSelectPanels2(sysConfig,...
                                            g,... %Exp_data.g,...
                                            target.locs,... imgDomain,...
                                            probes, probeTables,chooseProbes_Threshold,...
                                            panels, panelTables,choosepanels_Threshold,...
                                            partitioning,...
                                            R_k...
                                            );
                                        
partitioning.indices=indices;         
partitioning.locs=locs;
fprintf('[Done] %f(s)\n', etime(clock(),mark));

recon_fig2=figure();
volume_plot(recon_fig1,imgDomain, f_mf, indices, domainGrid, 10); title('look-up table')
plotScatters(f_mf, indices,imgDomain, domainGrid, 'Type','log');
