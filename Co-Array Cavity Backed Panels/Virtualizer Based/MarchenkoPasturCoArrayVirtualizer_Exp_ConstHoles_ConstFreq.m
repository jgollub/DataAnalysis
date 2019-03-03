%% Investigate changing Q with number of frequency pts oversampled, and number of holes set
%% import CoArrays
disp(sprintf(['\n+===================================================+\n',...
               '+===================================================+\n']))

figure(11); clf;
SVD_Plot_data_Norm=[];
SVD_Plot_data=[];
counter=0;
fbegin=18e9;
fend=26.5e9;
num_a1=20;
num_a2=20;


f0=(fbegin+fend)/2;
B=(fend-fbegin);
Q_threshold=floor(num_a2*num_a1*f0/(B));
% Q_cavity_values=[250 500 1000 2000 4000 8000 16000 Q_threshold];
Q_cavity_values=[250 500 1000 2000 4000 8000 Q_threshold];

 nfreqpts=floor(Q_cavity_values(end-1)*B/(f0)); %2x last frequency threshold
%  nfreqpts=floor(Q_threshold*B/(f0)); %2x last frequency threshold

for Q_cavity=Q_cavity_values
    
    counter=counter+1;
    
    disp(sprintf(['\nNUMBER OF FREQUENCY POINTS:\n'...
        'Useful Frequency samples %g.\n'...
        'Selected number of frequency points %g.\n',...
        'Q of cavity %g.\n',...
        'Tau of cavity %g.\n'], ...
        floor(B*Q_cavity/f0),...
        nfreqpts,...
        Q_cavity,...
        Q_cavity/B))
    
    f=linspace(fbegin, fend, floor(nfreqpts));
    
    disp(sprintf(['NUMBER OF HOLES\n',...
        'Useful number of holes = %g. \n',...
        'Selected number of holes P1 = %g.\n'...
        'Selected number of holes P2 = %g.\n', ...
        'P1 X P2 = %g.\n',...
        'Useful Q (given set P1,P2) = %g.\n'], ...
        floor(B*Q_cavity/f0),...
        num_a1,...
        num_a2,...
        num_a1*num_a2,...
        num_a1*num_a2*(f0/B)));
    
    % Setup aperture size in mill's cross pattern
    aperture_size=(3e8/fend)*num_a1/2;
    disp(sprintf(['APPERTURE:\n',...
        'Spacing between points=%g m'],...
        aperture_size/num_a1));
    disp(sprintf('wavelength=%g m\n',3e8/fend));
    
    %=======Mills Cross==========================
    
    % caluculate element spacing
    element_spacing=aperture_size/(num_a1-1);
    
    % Mills cross panel 1
    MCr_array1=zeros(num_a1);
    MCr_array1(2:2:end,1)=1;
    MCr_array1(1:2:end,end)=1;
    
    % MCr_array1=ones(num_a1);
    
    %Mills cross panel 2
    MCr_array2=zeros(num_a2);
    MCr_array2(1,1:2:end)=1;
    MCr_array2(end,2:2:end)=1;
    
    %create panel 1
    panel1_MillsC=create_panel2('type','CoArray','CoArrayElements',MCr_array1,...
        'fsweep',f,...
        'sizeY',aperture_size,...
        'sizeZ',aperture_size,...
        'ElementSizeY',element_spacing,...
        'ElementsizeZ', element_spacing ...
        );
    %feed panel from randomized cavity 1
    panel1_MillsC=panel_feed(panel1_MillsC, 'type','random_cavity_phase_amplitude',...
        'Q_cavity',Q_cavity,...
        'renormalize','self');
    %
    %create panel 2
    panel2_MillsC=create_panel2('type','CoArray','CoArrayElements',MCr_array2,...
        'fsweep',f,...
        'sizeY',aperture_size,...
        'sizeZ',aperture_size,...
        'ElementSizeY',element_spacing,...
        'ElementsizeZ',element_spacing ...
        );
    
    % %feed panel from randomized cavity 2
    panel2_MillsC=panel_feed(panel2_MillsC, 'type','random_cavity_phase_amplitude',...
        'Q_cavity',Q_cavity,...
        'renormalize','self');
    
      %panel2_MillsC=create_panel('type', 'dipole','fsweep',f);
    %  probe=panel_offset(probe, [0,aperture_size/2,0]);
    
    figure(12); cla; scatter3(panel1_MillsC.x,panel1_MillsC.y,panel1_MillsC.z,10,'k');
    hold on;   scatter3(panel2_MillsC.x,panel2_MillsC.y,panel2_MillsC.z,10,'r');
    %   scatter3(panel2_MillsC.x,panel2_MillsC.y,panel2_MillsC.z,10,'blue')
    
    %% Define imaging domain (slice)
    
    % define the absolute imaging domain
    
    x_range=1;
    
    % minimum wavelength
    lambda_min=3E8/max(f);
    
    %calculate the nyquist/SBP sampling rate at the scene =
    %range*(sampling rate)/aperture
    
    R_l=0.5*x_range*(lambda_min/2)/(aperture_size); %Oversampled x2
    
    scene_size=x_range*lambda_min/(element_spacing);
    imgDomain_offset = [x_range 0 0];
    imgDomain_range  = [0  scene_size scene_size];
    
    [imgDomain, domainGrid] =test_space3(imgDomain_offset,imgDomain_range,R_l);
    
    disp(sprintf(['SCENE:\n'...
        'domain y, z size = %g m\n',...
        'sampling dy, dz = %g m\n',...
        'number of voxels = %g \n'], ...
        scene_size, R_l, size(imgDomain.locs,1)));
    
    figure(12); hold on; scatter3(imgDomain.locs(1:2:end,1),imgDomain.locs(1:2:end,2),imgDomain.locs(1:2:end,3),5,'r');
    axis equal
    
    
    %% Define imaging target
    %% Target 3 =================================================================
    % res=(3e8/fend)/aperture_size;
    % target=struct();
    %
    % theta=[];
    % theta(1)=.5;
    % Arc=R_l;
    % a=.25/(2*pi);
    % % b=a+res/(2*pi);
    %
    % r_check=scene_size/2;
    % i=1;
    % while r_check > a*theta(end)
    % theta(i+1)=theta(i)+Arc/(a*theta(i));
    % i=i+1;
    % end
    %
    % yy1=a.*theta.*cos(theta);
    % zz1=a.*theta.*sin(theta);
    % xx1=ones(1,numel(theta));
    %
    %  xx2=b.*theta.*cos(theta);
    %  yy2=b.*theta.*sin(theta);
    %  zz2=ones(1,numel(theta));
    % figure(24); cla; scatter3(xx1,yy1,zz1,10,'r')
    %   hold on; axis equal; scatter3(xx2,yy2,zz2,10,'r')
    %
    % xx=imgDomain.locs(:,1);
    % yy=imgDomain.locs(:,2);
    % zz=imgDomain.locs(:,3);
    % spiral=zeros(numel(xx),1);
    %
    % for i=1:numel(xx1)
    % spiral=logical(sqrt((xx-xx1(i)).^2+(yy-yy1(i)).^2+(zz-zz1(i)).^2)<4*R_l)+spiral;
    % % spiral=logical(sqrt((xx-xx2(i)).^2+(yy-yy1(i)).^2+(zz-zz1(i)).^2)<2*R_l)+spiral;
    % end
    % spiral=logical(spiral);
    % target.locs=imgDomain.locs(spiral,:);
    % target.sigma=ones(size(target.locs));
    %
    %
    % test=figure(12);
    %  hold on; scatter3(target.locs(:,1),target.locs(:,2),target.locs(:,3),'b');
    %
    % %
    % % barwidth=(3e8/fend)/aperture_size;
    % %
    % % disp(sprintf('bar width=%g',barwidth));
    % % object  = objectCreator('Spiral');
    % % object  = objectCreator('ResTarget',1,target_res,barwidth);
    % % num_y=floor(scene_sampling/(6*barwidth));
    % % num_z=floor(scene_sampling/(12*barwidth));
    % %         target_res   = R_l/3;
    % %         targetOffs=[];
    % %
    % %
    % %
    % % target=cell(num_z,num_y);
    % %
    % % for i=1:num_z
    % %     for j=1:num_y
    % %     object  = objectCreator('ResTarget',1,target_res,barwidth);
    % %     targetOffs.x = 1;
    % %     targetOffs.y = scene_sampling/2-6*barwidth*j;
    % %     targetOffs.z = scene_sampling/2-12*barwidth*i;
    % %
    % %         % Move
    % %         object = objectMover(object,targetOffs);
    % %
    % % % Zbuffer
    % % target{i,j}  = ZBuffer(object,target_res);
    % % % target0 = ZBuffer(object,2*target_res);
    % %
    % % test=figure(12);
    % %  hold on; scatter3(target{i,j}.locs(:,1),target{i,j}.locs(:,2),target{i,j}.locs(:,3),'b');
    % axis equal
    % view(-90,0);
    % xlabel('xaxis')
    % ylabel('yaxis')
    % zlabel('zaxis')
    
    %% generate g
    % g=forward_model(panel2_MillsC, panel1_MillsC, target.sigma.*(randn(size(target.sigma))+1i*randn(size(target.sigma))), target.locs);
    %
    %
    % % gij = forward_model(panel2_MillsC, panel1_MillsC, target{i,j}.sigma.*(randn(size(target{i,j}.sigma))+1i*randn(size(target{i,j}.sigma))), target{i,j}.locs);
    % % if i==1 & j==1
    % % g=zeros(size(gij(:)));
    % % end
    % % g=g+gij(:);
    % %
    % %     end
    % % end
    
    %% ================== calculating the fields ========================
    
    panel1_MCr_fields = dipoles_to_fieldsEXP3(panel1_MillsC, imgDomain.locs);
    panel2_MCr_fields = dipoles_to_fieldsEXP3(panel2_MillsC, imgDomain.locs);
    
    %% ================generate H ======================================
    
    
    H_MCr=makeH_faceted(panel1_MCr_fields,panel2_MCr_fields);
    SVD_Plot_data(:,counter)=svd(H_MCr,'econ');
    SVD_Plot_data_Norm(:,counter)=SVD_Plot_data(:,counter)/max(SVD_Plot_data(:,counter));
    
    title('SVD (Cross-Range Slice)')
    
    legend('CoArray')% ,'CoArray (random)') % ,'Conventional, Q=1', 'Conventional, Q=100', 'Conventional, Q=10000')
end

%%
clear SVD_plot_Q
SVD_plot_Q=figure();

% Create axes
axes1 = axes('Parent',SVD_plot_Q);
box(axes1,'on');
hold(axes1,'off');

legend_text=['\bf Q='];

color_palette=[0.8 0.8 1
               0.6 0.6 1
               0.4 0.4 1
               0.2 0.2 1
               0 0 1
               0 0 0.8
               0 0 0.6
               0 0 0.4
               0 0 0.2
               0 0 0];

legend_text=['\bf','Q = '];           
% Create multiple lines using matrix input to plot
plot1 = semilogy(SVD_Plot_data_Norm,'LineWidth',3,'Parent',axes1);
% plot1 = semilogy(SVD_Plot_data/max(max(SVD_Plot_data)),'LineWidth',3,'Parent',axes1);
set(plot1(1),'DisplayName',[legend_text,num2str(Q_cavity_values(1))],...
    'Color',color_palette(1,:));
set(plot1(2),'DisplayName',[legend_text,num2str(Q_cavity_values(2))],...
    'Color',color_palette(2,:));
set(plot1(3),'DisplayName',[legend_text,num2str(Q_cavity_values(3))],...
    'Color',color_palette(3,:));
set(plot1(4),'DisplayName',[legend_text,num2str(Q_cavity_values(4))],...
    'Color',color_palette(4,:));
set(plot1(5),'DisplayName',[legend_text,num2str(Q_cavity_values(5))],'Color',color_palette(5,:));
set(plot1(6),'DisplayName',[legend_text,num2str(Q_cavity_values(6))],'Color',color_palette(6,:));
%set(plot1(7),'DisplayName',[legend_text,num2str(Q_cavity_values(7))],'Color',color_palette(7,:));
 set(plot1(7),'DisplayName',[legend_text,num2str(Q_cavity_values(7)),' (NMf_0/B) '],'LineStyle','--','Color',[1 0 0]);

% Create title

title(['Normalized SVD Vs. Q (# Freq. Pts = ', num2str(nfreqpts), ', N & M = ', num2str(num_a1),' Irises)']);
xlabel('\bf Measurement Mode Index')
ylabel('\bf Singular Value')
% Create legend
legend(axes1,'show');
set(gca,'fontsize',16)
set(gca,
%%
clear SVD_plot_Q
SVD_plot_Q=figure();

% Create axes
axes1 = axes('Parent',SVD_plot_Q);
box(axes1,'on');
hold(axes1,'off');

legend_text=['\bf Q='];

color_palette=jet(7);

legend_text=['\bf','Q = '];           
% Create multiple lines using matrix input to plot
 plot1 = semilogy(SVD_Plot_data_Norm,'LineWidth',3,'Parent',axes1);
% plot1 = semilogy(SVD_Plot_data/max(max(SVD_Plot_data)),'LineWidth',3,'Parent',axes1);
set(plot1(1),'DisplayName',[legend_text,num2str(Q_cavity_values(1))],...
    'Color',color_palette(1,:));
set(plot1(2),'DisplayName',[legend_text,num2str(Q_cavity_values(2))],...
    'Color',color_palette(2,:));
set(plot1(3),'DisplayName',[legend_text,num2str(Q_cavity_values(3))],...
    'Color',color_palette(3,:));
set(plot1(4),'DisplayName',[legend_text,num2str(Q_cavity_values(4))],...
    'Color',color_palette(4,:));
set(plot1(5),'DisplayName',[legend_text,num2str(Q_cavity_values(5))],'Color',color_palette(5,:));
set(plot1(6),'DisplayName',[legend_text,num2str(Q_cavity_values(6))],'Color',color_palette(6,:));
%set(plot1(7),'DisplayName',[legend_text,num2str(Q_cavity_values(7))],'Color',color_palette(7,:));
 set(plot1(7),'DisplayName',[legend_text,num2str(Q_cavity_values(7)),' (NMf_0/B) '],'LineStyle','--','Color',[1 0 0]);

% Create title

title(['Normalized SVD Vs. Q (# Freq. Pts = ', num2str(nfreqpts), ', N & M = ', num2str(num_a1),' Irises)']);
xlabel('\bf Measurement Mode Index')
ylabel('\bf Singular Value')
% Create legend
legend(axes1,'show');
set(gca,'fontsize',16)


%%
sampconst=1;
% [edges,N]=histcounts(S_probe,10);
figure(15); clf;
for i=1:size(SVD_Plot_data,2)
subplot(4,ceil(size(SVD_Plot_data,2)/4),i)
MarchenkoPasturOverlayHist(SVD_Plot_data(:,i),50,sampconst);
% title(['Q =',num2str(Q_cavity_values(i))]);
end
% figure(16); clf;
% MarchenkoPasturOverlayHist(SVD_Plot_data_Norm(:,1),50,.3*sampconst);

%% reconstruct image
% g=g(:);
% f_est_MF=H_MCr'*g;
% f_est_CGS = cgs(H_MCr'*H_MCr, H_MCr'*g, 1e-3, 20);
% figure(27); subplot (1,2,1);
% cla; scatter3(imgDomain.locs(:,1),...
%                  imgDomain.locs(:,2), ...
%                  imgDomain.locs(:,3), ...
%                  30, ...
%                  (abs(f_est_MF).^2),...
%                  'filled');
% title('MF');
% view(-90,0);
% axis equal
% axis tight
% figure(27); subplot (1,2,2);
% cla; scatter3(imgDomain.locs(:,1),...
%                  imgDomain.locs(:,2), ...
%                  imgDomain.locs(:,3), ...
%                  30, ...
%                  (abs(f_est_CGS).^2),...
%                  'filled');
% title('CGS');
% view(-90,0);
% axis equal
% axis tight
%
% %%---------------------
%
% figure(28); subplot (1,2,1);
% f_est_MF_grid=reshape(f_est_MF,[domainGrid(1) domainGrid(3)]);
%
% cla; imagesc(imgDomain.locs(:,2),...
%                  imgDomain.locs(:,3), ...
%                  (abs(f_est_MF_grid).^2)...
%                   );
%               axis equal
%               axis tight
% title('MF');
% axis equal
% axis tight
% figure(28); subplot (1,2,2);
% f_est_CGS_grid=reshape(f_est_CGS,[domainGrid(1) domainGrid(3)]);
%
% cla; imagesc(imgDomain.locs(:,2),...
%                  imgDomain.locs(:,3), ...
%                  (abs(f_est_CGS_grid).^2)...
%                  );
% title('CGS');
% axis equal
% axis tight