%load data
close all
%expdata=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Measured Scenes\Anechoic Chamber Upgrade\newTiles.mat')

expdata=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Measured Scenes\inSceneCalibrationTarget\bigEggCrate.mat')
%expdata=load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Measured Scenes\Semi-Rigid&Phase_stable\Mannequin Bust Measurement\Mannequin_with_Gun_Kinect.mat')

%add dummy panel rotation matrix intilized to zero---> will be included in

for
%%panel=      [A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4]
% choosePanels=[  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1];
choosePanels=ones(12,1);
choosePanel(panel_set)
%%probe=     [P1, P2, P3, p4, p5, p6];
chooseProbes=[ 1, 1,  1,  1,  1,  1];

[f_est,f_est_plot,locs]=Offline_Recon(...
    expdata,...
    'Probes', chooseProbes,...
    'Panels', choosePanels,...
    'ProbePositions', probePosition,...
    'PanelPositions', panelPosition,...
    'PanelRotations', panelRotation,...
    'Recon_Type','partition',...
    'ROI',1 ...
   ... 'Recon_Type','TwIST_TV',...
   ... 'GridSpacing', 0.01...
   ...'XRange',[1.3 1.37]...
   ... 'YRange',[-0.4 -0.2],...
   ... 'ZRange',[0.37 0.47]...
    );
    
%%Plotting

 %3D plotting
 if isempty(f_est_plot)
     
            try figure(image_plot)
            catch
                image_plot=figure();
                scrsize=get(0,'ScreenSize');
                figsize=image_plot.Position(3:4);
                image_plot.Position(1:2)=[scrsize(3)-figsize(1),scrsize(4)/2-figsize(2)-80];
            end
            clf
            subplot(1,3,1); volume_plot(...
                image_plot, ...
                locs.locs, ...
                f_est.f_mf, ...
                cell2mat(locs.indices_kl), ...
                [locs.ny locs.nx locs.nz], ...
                'linear' ...
            );
            
            title('matched filter');
            subplot(1,3,2); volume_plot(... 
                image_plot, ...
                locs.locs, ...
                f_est.f_gmres, ...
                cell2mat(locs.indices_kl), ...
                [locs.ny locs.nx locs.nz], ...
                'linear' ...
            );
            title('least square');

 else

upsample_3D=2;
upsample_2D=4;
 
upsamp_2D_plot=[];
f_est_plot_3D=[];


%3D plot
[f_est_plot_3D] = upsample_image_3D(f_est_plot, upsample_3D);

 figure(105);
      clf;
     
      vol3d('cdata',abs(f_est_plot_3D),...
        'XData',[min(locs.xrange) max(locs.xrange)],...
        'YData',[min(locs.yrange) max(locs.yrange)],...
        'ZData',[min(locs.zrange) max(locs.zrange)], ...
        'alpha',(abs(f_est_plot_3D)/max(max(max(abs(f_est_plot_3D))))) ...
        );
    axis equal;
    axis tight;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(-48,14)
    title(['Reconstruction from Experimental Data']);
    fprintf('Scene reconstructed in %f minutes\n', toc/60);
end    
%% Plot all slices
% 
% for i=1:nx
%      upsamp_2D_plot(:,:,i)=upsample_image(permute(f_est_plot(:,i,:),[3 1 2]),upsample_2D);
%       %  upsamp_2D_plot(:,:,i)=permute(f_est_plot(:,i,:),[3 1 2]);
% end
% %find max in slice for normalization
% [maxplot_val max_index]=max(max(max(abs(upsamp_2D_plot(:,:,nx)))));
% 
% for i_slice=1:nx
%     %     subplot(ceil(length(slices)/2),2,mod(i_slice,6)+1);
%     
%     % subplot(2,Num_figs,plotpos);
%     
%     %plot image
%     figure()
%     images(i)=imagesc(y_extent, z_extent,abs(upsamp_2D_plot(:,:,i_slice)),[0 maxplot_val+maxplot_val/3]);
%     colorbar
%     
%     slice_pos=unique(x);
%     title(['slice num: ',num2str(i_slice),'; X Pos:', num2str(slice_pos(i_slice))]);
%     xlabel('Y (m)');
%     ylabel('Z (m)');
%     
%     axis xy
%     set(gca,'XDir','reverse')
%     axis equal
%     axis tight
%     
% end
% 
% %plot sum of 3 slices with highest "energy." at 5mm resolution this
% %consists with the range resolution of 
% figure(107)
% clf
% 
% for i=1:nx
% slice_energy(i)=max(max(permute(f_est_plot(:,i,:),[3 1 2])))
% end
% [~,max_slices]=sort(slice_energy);
% slices_to_plot=max_slices(end-1:end);
% [maxplot_val max_index]=max(max(sum(permute(f_est_plot(:,slices_to_plot,:),[3 1 2]),3)./3));
% 
% sum_plot=upsample_image(sum(permute(f_est_plot(:,slices_to_plot,:),[3 1 2]),3)./3,upsample_2D);
% imagesc(y_extent, z_extent,abs(sum_plot),[0 maxplot_val]);
% colorbar
% 
% title(['Sum Plot']);
% xlabel('Y (m)');
% ylabel('Z (m)');
% 
% axis xy
% set(gca,'XDir','reverse')
% axis equal
% axis tight

% % if debugVal
% %     figure(102);
% %     clf;
% %     vol3d('cdata',reshape(abs(H(1, :)),...
% %         max(roi(:,2))-min(roi(:,2))+1,max(roi(:,1))-min(roi(:,1))+1, max(roi(:,3))-min(roi(:,3))+1),...
% %         'XData',[min(x) max(x)],...
% %         'YData',[min(y) max(y)],...
% %         'ZData',[min(z) max(z)] ...
% %     );
% %
% %     axis equal;
% %     axis tight;
% %     xlabel('X');
% %     ylabel('Y');
% %     zlabel('Z');
% %     view([-40,40]);
% %     title('H Measurement matrix (Feed 1, Freq 51)')
% % end
% 
% % try figure(Obj_funct_Plot); catch Obj_funct_Plot=figure(); end
% %
% % semilogy(abs(objfunc(2:end)-objfunc(1:(end-1))))
% % xlabel('twist iteration')
% % ylabel('magnitude of objective function change')
% % drawnow
% % 
% % try figure(image_plot); catch image_plot=figure(); end
% % clf
% % 
% % %plot 3D and slices
% % %for plotting with imagesc
% % xVals_plot=unique(locs.locs(:,1));
% % yVals_plot=unique(locs.locs(:,2));
% % zVals_plot=unique(locs.locs(:,3));
% % 
% % x_reshape=length(unique(locs.locs(:,1)));
% % y_reshape=length(unique(locs.locs(:,2)));
% % z_reshape=length(unique(locs.locs(:,3)));
% 
%  
%    
%save(['C:\Users\Jonah Gollub\Desktop\Resolution Targets Solved Data\solvedData',datestr(now,'mmmm-dd-yyyy--HH-MM-AM'),'.mat'],'f_est', 'f_est_plot','nx','ny','nz','x','y','z','x_extent','y_extent','z_extent')    