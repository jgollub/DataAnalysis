freqpt=401;
scalingbar=[0 .01];
figure(1)

subplot(2,3,1)

imagesc(twoPanel.X(1,:),twoPanel.Y(:,1),abs(twoPanel.measurements(:,:,freqpt)),scalingbar)
set(gca,'YDir','normal')
axis equal
axis tight

subplot(2,3,4)
imagesc(onePanel.X(1,:),onePanel.Y(:,1),abs(onePanel.measurements(:,:,freqpt)),scalingbar)
set(gca,'YDir','normal')
axis equal
axis tight

subplot(1,3,2:3)
imagesc(onePanel.X(1,:),onePanel.Y(:,1),(abs(onePanel.measurements(:,:,freqpt))-abs(twoPanel.measurements(:,:,freqpt)))./abs(onePanel.measurements(:,:,freqpt)),[0,1])
set(gca,'YDir','normal')
axis equal
axis tight
colorbar()

% 
% imagesc(twoPanel.X(1,:),twoPanel.Y(:,1),angle(twoPanel.measurements(:,:,freqpt)))
% set(gca,'YDir','normal')
% axis equal
% axis tight
% 
% subplot(2,3,4)
% imagesc(onePanel.X(1,:),onePanel.Y(:,1),angle(onePanel.measurements(:,:,freqpt)))
% set(gca,'YDir','normal')
% axis equal
% axis tight
% 
% subplot(1,3,2:3)
% imagesc(onePanel.X(1,:),onePanel.Y(:,1),angle(onePanel.measurements(:,:,freqpt))-angle(twoPanel.measurements(:,:,freqpt)))
% set(gca,'YDir','normal')
% axis equal
% axis tight
% colorbar()