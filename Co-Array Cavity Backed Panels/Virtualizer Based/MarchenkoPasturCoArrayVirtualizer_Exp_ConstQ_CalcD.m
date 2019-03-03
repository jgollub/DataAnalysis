% Run coarray distribution in virtualizer 
 

%%need to check
%orientation of coarray feeding into virtualizer
%check normalization currently set to self

% %% Choose coarray
% [array1,array2,bestcvl]=gencoarray([10 10],[10 10],0,1)
% l1_z=size(array1,1)/2;
% l1_y=size(array1,2)/2;
% l2_z=size(array2,1)/2;
% l2_y=size(array2,2)/2;
% 
% array1=array1(1:l1_z,1:l1_y);
% array2=array2(1:l2_z,1:l2_y);
% 
% figure(9); cla; subplot(1,2,1); cla;
% imagesc(array1); axis equal; axis tight;
% subplot(1,2,2); hold on; 
% imagesc(array2); axis equal; axis tight;
% save('C:\Users\Jonah Gollub\Documents\code\data Analysis\Co-Array\last_CoArray.mat','array1','array2');
% % 
% close all
%% import CoArrays
figure(12); clf; 
Q_cavity=1047; %10000

fbegin=18e9;
fend=26.5e9;
B=(fend-fbegin);
f0=(fbegin+fend)/2;
% f=linspace(18e9, 26e9, floor(Q_cavity/2.5));
 nfreqpts_match=floor(B*Q_cavity/(f0));
%nfreqpts_match=floor(B*Q_cavity/(f0));


freq_value=[100 200 400 800 1600 3200 nfreqpts_match]
% freq_value=[800 1200 1600 2000 2400 2800 nfreqpts_match]

SVD_Plot_data=NaN(max(freq_value),numel(freq_value));
SVD_Plot_data_Norm=NaN(max(freq_value),numel(freq_value));
counter=0;
disp(sprintf(['\n\n+===================================================+\n\n']))

for nfreqpts=freq_value
counter=counter+1;
% select number of frequency points
disp(sprintf(['\nNUMBER OF FREQUENCY POINTS:\n'...
              'Useful Frequency samples %g.\n'...
              'Selected number of frequency points %g.\n',...
              'Q of cavity %g.\n',...
              'Tau of cavity %g.\n'], ...
              floor(B*Q_cavity/f0),...
              nfreqpts,...
              Q_cavity,...
              2*Q_cavity/B));

sampconst=1;
f=linspace(fbegin, fend, floor(nfreqpts)*sampconst);

%  number of holes
%  num_a1=2*floor(sqrt( floor(2*(fend-fbegin)*Q_cavity/(fbegin+fend))));
%  num_a1=6;
%  num_a2=6;

num_a1=round(sqrt(B*Q_cavity/(f0)));
num_a2=num_a1;



disp(sprintf(['NUMBER OF HOLES\n',...
              'Useful number of holes = %g. \n',...
              'Selected number of holes P1 = %g.\n'...
              'Selected number of holes P2 = %g.\n', ...
              'P1 X P2 = %g.\n', ...
              'Useful Q (given set P1,P2) = %g.\n'], ...
              floor(pi*B*Q_cavity/f0),...
              num_a1,...
              num_a2,...
              num_a1*num_a2,...
              num_a1*num_a2*(f0/(pi*B))...
              ));

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
panel1_MillsC=panel_feed(panel1_MillsC, 'type','random_cavity_phase_amplitude', ...
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
% 
% %feed panel from randomized cavity 2    
panel2_MillsC=panel_feed(panel2_MillsC, 'type','random_cavity_phase_amplitude',...
    'Q_cavity',Q_cavity,...
    'renormalize','self');

%  panel2_MillsC=create_panel('type', 'dipole','fsweep',f);
%  probe=panel_offset(probe, [0,aperture_size/2,0]);

figure(12); cla; scatter3(panel1_MillsC.x,panel1_MillsC.y,panel1_MillsC.z,10,'k');
%       hold on;   scatter3(panel2_MillsC.x,panel2_MillsC.y,panel2_MillsC.z,10,'r');
   scatter3(panel2_MillsC.x,panel2_MillsC.y,panel2_MillsC.z,10,'blue')

%% Define imaging domain (slice)

% define the absolute imaging domain
x_range=1;
imgDomain_offset=[x_range 0 0];

% minimum wavelength
lambda_min=3E8/max(f);

%calculate the nyquist/SBP sampling rate at the scene =
%range*(sampling rate)/aperture

R_l=.5*x_range*(lambda_min/2)/(aperture_size); %oversample x 2??

scene_size=x_range*lambda_min/(element_spacing);
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
figure(11);

H_MCr=makeH_faceted(panel1_MCr_fields,panel2_MCr_fields);

if size(H_MCr,2)>size(H_MCr,1)
    SVD_Plot_data(1:size(H_MCr,1),counter)=svd(H_MCr,'econ');
else
    SVD_Plot_data(1:size(H_MCr,2),counter)=svd(H_MCr,'econ');
end

SVD_Plot_data_Norm(1:numel(SVD_Plot_data(:,counter)),counter)=...
                    SVD_Plot_data(:,counter)/max(SVD_Plot_data(:,counter));

drawnow;

end
%% plot SVD Versus number Frequencies
clear SVD_Plot
SVD_Plot=figure(); clf;
% Create axes

axes1 = axes('Parent',SVD_Plot);
box(axes1,'on');
hold(axes1,'off');

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
 color_palette=flip(color_palette(1:6,:),1)

legend_text=['\bf','# Freq. Pts. = '];
% Create multiple lines using matrix input to plot
resort_plot_data=[flip(SVD_Plot_data_Norm(:,1:end-1),2), SVD_Plot_data_Norm(:,end)];
%resort_plot_data=[flip(SVD_Plot_data(:,1:end-1),2)/max(max(SVD_Plot_data)), SVD_Plot_data(:,end)/max(max(SVD_Plot_data))];
plot1 = plot(resort_plot_data(:,2:end),'LineWidth',3,'Parent',axes1);
set(plot1(1),'DisplayName',[legend_text,num2str(freq_value(5))],...
    'Color',color_palette(1,:));
set(plot1(2),'DisplayName',[legend_text,num2str(freq_value(4))],...
    'Color',color_palette(2,:));
set(plot1(3),'DisplayName',[legend_text,num2str(freq_value(3))],...
    'Color',color_palette(3,:));
set(plot1(4),'DisplayName',[legend_text,num2str(freq_value(2))],...
    'Color',color_palette(4,:));
 set(plot1(5),'DisplayName',[legend_text,num2str(freq_value(1))],'Color',color_palette(5,:));
% set(plot1(6),'DisplayName',[legend_text,num2str(freq_value(1))],'Color',color_palette(6,:));
%    set(plot1(7),'DisplayName',[legend_text,num2str(freq_value(7))],'Color',[0 0 0.600000011920929]);
set(plot1(6),'DisplayName',[ legend_text,num2str(freq_value(7)), ' (BQ/( f0))'],'LineStyle','--','Color',[1 0 0]);
 
% Create title
title(['Normalized SVD Vs. Number of Freq. Pts. (Q = ', num2str(Q_cavity), ', N & M = ', num2str(num_a1),' Irises)']);

% Create legend
legend(axes1,'show');
%% plot SVD Versus number Frequencies
clear SVD_Plot
SVD_Plot=figure(); clf;
% Create axes

axes1 = axes('Parent',SVD_Plot);
box(axes1,'on');
hold(axes1,'off');

color_palette=jet(7)
 color_palette=flip(color_palette(1:6,:),1)

legend_text=['\bf','# Freq. Pts. = '];
% Create multiple lines using matrix input to plot
resort_plot_data=[flip(SVD_Plot_data_Norm(:,1:end-1),2), SVD_Plot_data_Norm(:,end)];
%resort_plot_data=[flip(SVD_Plot_data(:,1:end-1),2)/max(max(SVD_Plot_data)), SVD_Plot_data(:,end)/max(max(SVD_Plot_data))];
plot1 = plot(resort_plot_data(:,2:end),'LineWidth',3,'Parent',axes1);
set(plot1(1),'DisplayName',[legend_text,num2str(freq_value(5))],...
    'Color',color_palette(1,:));
set(plot1(2),'DisplayName',[legend_text,num2str(freq_value(4))],...
    'Color',color_palette(2,:));
set(plot1(3),'DisplayName',[legend_text,num2str(freq_value(3))],...
    'Color',color_palette(3,:));
set(plot1(4),'DisplayName',[legend_text,num2str(freq_value(2))],...
    'Color',color_palette(4,:));
 set(plot1(5),'DisplayName',[legend_text,num2str(freq_value(1))],'Color',color_palette(5,:));
% set(plot1(6),'DisplayName',[legend_text,num2str(freq_value(1))],'Color',color_palette(6,:));
set(plot1(7),'DisplayName',[legend_text,num2str(freq_value(7))],'Color',[0 0 0.600000011920929]);
set(plot1(6),'DisplayName',[ legend_text,num2str(freq_value(7)), ' (BQ/( f0))'],'LineStyle','--','Color',[1 0 0]);
 
% Create title
title(['Normalized SVD Vs. Number of Freq. Pts. (Q = ', num2str(Q_cavity), ', N & M = ', num2str(num_a1),' Irises)']);

% Create legend
legend(axes1,'show');
set(gca,'fontsize',16)
 
 
 %%
figure(15); clf;
for i=1:size(SVD_Plot_data,2)
subplot(4,2,i)
d=freq_value(i)/nfreqpts_match/(4);
MarchenkoPasturOverlayHist(SVD_Plot_data(~isnan(SVD_Plot_data(:,i)),i),50,d);
 title(['# Freq =',num2str(freq_value(i)),' \lambda= ',num2str(d)]);
xlim([0 5])

end
 
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