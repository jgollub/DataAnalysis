panels_primed=panels;
  MSE_panel=zeros(length(solvedGlobalPosition),3);
for iterate=1:10
      MSE_panel(iterate,:)= mean(abs(locate(panels_primed(1,:))-locate(shiftedPanels(1,:))));
      
    %% recon full region
mark=clock();
[f_mf,indices, locs, timing] = reconstructionSelectPanels(sysConfig,...
                                            g,... %Exp_data.g,...
                                            imgDomain,...
                                            probes, probeTables,chooseProbes,...
                                            panels_primed, panelTables,choosePanels,...
                                            partitioning,...
                                            R_k...
                                            );
                                        
partitioning.indices=indices;         
partitioning.locs=locs;
fprintf('[Done] %f(s)\n', etime(clock(),mark));

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
MSE_pts(iterate)=sum(sum((diff).^2,2),1)/length(scattererPositions)

      [K_U K_r]=Kabsch(...
          solvedGlobalPosition.',  ...        
          scattererPositions(sortIndex,:).' ...
            );

    for i=1:sysConfig.nPanels
        for j=1:sysConfig.nFeeds
            panels_primed(j,i)=transformPanelArbitraryUr(panels_primed(j,i),K_U,K_r);
        end
    end
    
end