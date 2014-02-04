freqpt=401;
scalingbar=[0 .01];
figure(1)
% 
% subplot(2,3,1)
% 
% imagesc(twoPanel.X(1,:),twoPanel.Y(:,1),abs(twoPanel.measurements(:,:,freqpt)),scalingbar)
% set(gca,'YDir','normal')
% axis equal
% axis tight
% 
% subplot(2,3,4)
% imagesc(onePanel.X(1,:),onePanel.Y(:,1),abs(onePanel.measurements(:,:,freqpt)),scalingbar)
% set(gca,'YDir','normal')
% axis equal
% axis tight
% 
% subplot(1,3,2:3)
% imagesc(onePanel.X(1,:),onePanel.Y(:,1),(abs(onePanel.measurements(:,:,freqpt))-abs(twoPanel.measurements(:,:,freqpt)))./abs(onePanel.measurements(:,:,freqpt)),[0,1])
% set(gca,'YDir','normal')
% axis equal
% axis tight
% colorbar()

subplot(2,3,1)
imagesc(twoPanel.X(1,:),twoPanel.Y(:,1),angle(twoPanel.measurements(:,:,freqpt)))
set(gca,'YDir','normal')
axis equal
axis tight

subplot(2,3,4)
imagesc(onePanel.X(1,:),onePanel.Y(:,1),angle(onePanel.measurements(:,:,freqpt)))
set(gca,'YDir','normal')
axis equal
axis tight



Phase_Dif=angle(onePanel.measurements(:,:,freqpt))-angle(twoPanel.measurements(:,:,freqpt));

IM_mask=ones(size(Phase_Dif));                     %Mask (if applicable)
%%
IM_mag=abs(Phase_Dif);                             %Magnitude image
IM_phase=angle(Phase_Dif);                         %Phase image

%%  Set parameters
max_box_radius=4;                           %Maximum search box radius (pixels)
threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image

%% Unwrap
residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
[IM_unwrapped, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping

%% Display results
figure; imagesc(residue_charge), colormap(gray), axis square, axis off, title('Phase residues (charged)');
figure; imagesc(branch_cuts), colormap(gray), axis square, axis off, title('Branch cuts');
figure; imagesc(immultiply(IM_phase,IM_mask)), colormap(gray), axis square, axis off, title('Wrapped phase');
tempmin=min(min(IM_unwrapped));          %This bit is just done to create a pleasing display when a mask is used
temp=(IM_unwrapped==0);
temp_IM=IM_unwrapped;
temp_IM(temp)=tempmin;
figure; imagesc(temp_IM), colormap(gray), axis square, axis off, title('Unwrapped phase');



subplot(1,3,2:3)
imagesc(onePanel.X(1,:),onePanel.Y(:,1), Phase_Dif,[.1 ,1])
set(gca,'YDir','normal')
axis equal
axis tight
colorbar()