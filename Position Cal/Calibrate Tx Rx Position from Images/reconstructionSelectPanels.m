%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
function [f_mf, indices,locs, timing] = reconstructionSelectPanels(sysConfig,...
    full_g,...
    Domain,...
    probes,probeTables,ChooseRx,...
    panels,panelTables,ChooseTx,...
    partitioning,...
    padding,...
    varargin...
    )

parser=inputParser;
defaultRecon='Partitioning';
defaultInputIndices=[];
addParameter(parser,'Reconstruction',defaultRecon);
addParameter(parser,'inputIndices',defaultInputIndices);
parse(parser,varargin{:});

%% Measurement ============================================================

g =reshape(full_g,sysConfig.nFrequencies,sysConfig.nFeeds,sysConfig.nPanels,sysConfig.nProbes);
g =g(sysConfig.fIndices,:,ChooseTx,ChooseRx);
g =g(:);

if strcmp(parser.Results.Reconstruction,'Partitioning')
    % Restructure the panels and probes for reconstruction
    % partition the frequency sweep
    freqs_st     =reshape(sysConfig.freqs, partitioning.T, partitioning.S);
    % unpack the individual transmit and receive antennas, used for recon
    primedPanels =unpack_duke(panels(:,ChooseTx), freqs_st);
    primedProbes =unpack_duke(probes(ChooseRx), freqs_st);
    % fprintf('[Done] %f(s)\n',etime(clock(),mark));
    
    %kinect_data   =cat(2,Exp_data.xGrid(:),Exp_data.yGrid(:),Exp_data.zGrid(:));
    %kinect_offset = mean(kinect_data);
    %% Pre-recon ==============================================================
    % do for every measurement
    fprintf('Pre-recon calculation... '); mark=clock();
    % re-order measurement vector according to the partitioning rules
    g_prep    =reorderG(g, primedProbes, primedPanels, freqs_st);
    % identify coarse bins that are sufficiently close to the kinect
    %padding   =R_k; % 'buffer' around kinect 'best-guess'
    % keep only subdomains close to the kinect return
    % roiIdx    =surf_constrain(partition.large_k, kinect_data, padding);
    roiIdx    =surf_constrain(partitioning.large_k, Domain, padding);
    
    scene_cfg =levelCleaning(partitioning.scene_cfg0, roiIdx);
    % prune beam patterns before calling GPU functions
    locs_k    =extract_subs(scene_cfg.locs_kl(1,:));
    [panelTablesPr, panelMem] =pruneBeamPattern(panelTables, primedPanels, locs_k);
    [probeTablesPr, probeMem] =pruneBeamPattern(probeTables, primedProbes, locs_k);
    % grab indices
    indices   =cell2mat(scene_cfg.indices_kl(1,:)');
    locs      =cell2mat(scene_cfg.locs_kl(1,:))';
    fprintf('[Done] %f(s)\n', etime(clock(),mark));
    fprintf('Lookup table memory footprint = %0.0f MB\n', (panelMem + probeMem)/1e6);
    
    %% Matched filter APB =====================================================
    % apply the adjoint of the measurement matrix to the measurement vector
    [f_mf, ~, timing] =multAPB_duke(g_prep, ...
        primedProbes, ...
        probeTablesPr, ...
        primedPanels, ...
        panelTablesPr, ...
        freqs_st, ...
        scene_cfg, ...
        true...
        );
    fprintf('matched filter completed in %0.3f seconds\n', timing.A+timing.B);
    
    %%========================================================================
elseif strcmp(parser.Results.Reconstruction,'Direct')
    
    indices=parser.Results.inputIndices; % must have been process from a previous run
    locs=  Domain;
    timing=NaN;
%     primedPanels =unpack_duke(panels(:,ChooseTx), freqs_st);
%     primedProbes =unpack_duke(probes(ChooseRx), freqs_st);
    panelFields = dipoles_to_fieldsEXP3(panels(:,ChooseTx), Domain);
    probeFields = dipoles_to_fieldsEXP3(probes(ChooseRx), Domain);
    
    H=makeH_faceted(probeFields,panelFields);
    
    f_mf=cgs(H'*H,H'*g,1e-3,30);
     
end
end
