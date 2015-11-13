% Run coarray distribution in virtualizer 

%% import CoArrays
figure(12); clf; 
Q_cavity=2000; %10000

fbegin=18e9;
fend=26.5e9;
B=(fend-fbegin);
f0=(fbegin+fend)/2;

% minimum wavelength
lambda_min=3E8/fend;

% f=linspace(18e9, 26e9, floor(Q_cavity/2.5));
nfreqpts_match=floor(B*Q_cavity/f0);


%freq_value=[ 100 200 400 800 1600 3200 nfreqpts_match]
freq_value=[200 400 800 1600 3200 nfreqpts_match]

%% aperture setup ============================================================
%  number of holes
num_a1=floor(sqrt(B*Q_cavity/f0));
num_a2=num_a1;

% Setup aperture size in mill's cross pattern
aperture_size=lambda_min*(num_a1/2-1);

% calculate element spacing
element_spacing=aperture_size/(num_a1/2-1);
panel_positions_spacing=aperture_size/(num_a1-1);

disp(sprintf(['\n\n+===================================================+\n\n']))
disp(sprintf(['APPERTURE:\n',...
              'Spacing between points=%g m'],...
              element_spacing));    

% Mills cross panel 1
MCr_array1=zeros(num_a1);
MCr_array1(2:2:end,1)=1;
MCr_array1(1:2:end,end)=1;

% MCr_array1=ones(num_a1);

%Mills cross panel 2
MCr_array2=zeros(num_a2);
MCr_array2(1,2:2:end)=1;
MCr_array2(end,1:2:end)=1;

  
          
%% Define imaging domain (slice) ========================================

% define the absolute imaging domain
x_range=1;
imgDomain_offset=[x_range 0 0];

%calculate the nyquist/SBP sampling rate at the scene =
%range*(sampling rate)/aperture

R_l=x_range*(lambda_min/4)/(aperture_size); %sample at nyquist for imaging

%cross range extent of unaliased region
scene_size=x_range*lambda_min/(element_spacing/2); % points are offset
imgDomain_range  = [0 scene_size scene_size];

[imgDomain, domainGrid] =test_space3(imgDomain_offset,imgDomain_range,R_l);


disp(sprintf(['SCENE:\n'...
              'domain y, z size = %g m\n',...
              'sampling dy, dz = %g m\n',...
              'number of voxels = %g \n'], ...
              scene_size, R_l, size(imgDomain.locs,1)));

figure(12); hold on; scatter3(imgDomain.locs(1:2:end,1),imgDomain.locs(1:2:end,2),imgDomain.locs(1:2:end,3),5,'r'); 
axis equal                 
%% Define imaging target ====================================================
figure(12); clf;
res=x_range*(3e8/fend)/aperture_size;
target=struct();
Arc=R_l;
[imgDomain_target] =test_space3(imgDomain_offset,imgDomain_range,R_l/2); %oversample to ensure accurate g

%creates spiral at resolution limit
theta=[];
theta(1)=1;
a=2*res/(2*pi);

r_check=scene_size/2;
i=1;
while r_check > a*theta(end)
theta(i+1)=theta(i)+Arc/(a*theta(i));
i=i+1;
end

yy1=a.*theta.*cos(theta);
zz1=a.*theta.*sin(theta);
xx1=ones(1,numel(theta));

figure(24); cla; 
scatter3(xx1,yy1,zz1,10,'r')

xx=imgDomain_target.locs(:,1);
yy=imgDomain_target.locs(:,2);
zz=imgDomain_target.locs(:,3);
spiral=zeros(numel(xx),1);

for i=1:numel(xx1)
% spiral=logical(sqrt((xx-xx1(i)).^2+(yy-yy1(i)).^2+(zz-zz1(i)).^2)<res/2)+spiral;
spiral=logical(sqrt((xx-xx1(i)).^2+(yy-yy1(i)).^2+(zz-zz1(i)).^2)<res/5)+spiral;
end
spiral=logical(spiral);
target.locs=imgDomain_target.locs(spiral,:);
target.sigma=ones(size(target.locs));

target.sigma=target.sigma.*(randn(size(target.sigma))+1i*randn(size(target.sigma)));


test=figure(12);
 hold on; scatter3(target.locs(:,1),target.locs(:,2),target.locs(:,3),'b'); 
 
axis equal
view(-90,0);
xlabel('xaxis')
ylabel('yaxis')
zlabel('zaxis')
%% sweep through different frequency sweeps          
disp(sprintf('wavelength=%g m\n',3e8/fend));
SVD_Plot_data=NaN(max(freq_value),numel(freq_value));
SVD_Plot_data_Norm=NaN(max(freq_value),numel(freq_value));

f_est_MF=zeros(length(imgDomain.locs), numel(freq_value));
f_est_CGS=zeros(length(imgDomain.locs), numel(freq_value));
counter=0;

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
              Q_cavity/(pi*f0)));

f=linspace(fbegin, fend, floor(nfreqpts));

disp(sprintf(['NUMBER OF HOLES\n',...
              'Useful number of holes = %g. \n',...
              'Selected number of holes P1 = %g.\n'...
              'Selected number of holes P2 = %g.\n', ...
              'P1 X P2 = %g.\n', ...
              'Useful Q (given set P1,P2) = %g.\n'], ...
              floor(B*Q_cavity/f0),...
              num_a1,...
              num_a2,...
              num_a1*num_a2,...
              num_a1*num_a2*(f0/(pi*B))...
              ));

%======= generate panels==========================

%create panel 1
panel1_MillsC=create_panel2('type','CoArray','CoArrayElements',MCr_array1,...
                     'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',panel_positions_spacing,...
                     'ElementsizeZ', panel_positions_spacing ...
                      );
%feed panel from randomized cavity 1  
rng(1) 
panel1_MillsC=panel_feed(panel1_MillsC, 'type','random_cavity_phase_amplitude', 'Q_cavity',Q_cavity,'renormalize','self');
% 
%create panel 2
panel2_MillsC=create_panel2('type','CoArray','CoArrayElements',MCr_array2,...
                     'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',panel_positions_spacing,...
                     'ElementsizeZ',panel_positions_spacing ...
                      );
% 
% %feed panel from randomized cavity 2    
rng(2) 
panel2_MillsC=panel_feed(panel2_MillsC, 'type','random_cavity_phase_amplitude', 'Q_cavity',Q_cavity,'renormalize','self');

%  panel2_MillsC=create_panel('type', 'dipole','fsweep',f);
%  probe=panel_offset(probe, [0,aperture_size/2,0]);

figure(12); cla; scatter3(panel1_MillsC.x,panel1_MillsC.y,panel1_MillsC.z,10,'k');
%       hold on;   scatter3(panel2_MillsC.x,panel2_MillsC.y,panel2_MillsC.z,10,'r');
   scatter3(panel2_MillsC.x,panel2_MillsC.y,panel2_MillsC.z,10,'blue')

%% generate g
g=forward_model(panel2_MillsC, panel1_MillsC,target.sigma, target.locs);


% gij = forward_model(panel2_MillsC, panel1_MillsC, target{i,j}.sigma.*(randn(size(target{i,j}.sigma))+1i*randn(size(target{i,j}.sigma))), target{i,j}.locs);
% if i==1 & j==1
% g=zeros(size(gij(:)));
% end
% g=g+gij(:);
% 
%     end
% end

%% ================== calculating the fields ========================

panel1_MCr_fields = dipoles_to_fieldsEXP3(panel1_MillsC, imgDomain.locs);
panel2_MCr_fields = dipoles_to_fieldsEXP3(panel2_MillsC, imgDomain.locs);

%% ================generate H ======================================
figure(11);

H_MCr=makeH_faceted(panel1_MCr_fields,panel2_MCr_fields);

if length(imgDomain.locs)> nfreqpts
    SVD_Plot_data(1:size(H_MCr,1),counter)=svd(H_MCr,'econ');
    SVD_Plot_data_Norm(1:numel(SVD_Plot_data(:,counter)),counter)=...
        SVD_Plot_data(:,counter)/max(SVD_Plot_data(:,counter));
else
    SVD_Plot_data(1:length(imgDomain.locs),counter)=svd(H_MCr,'econ');
    SVD_Plot_data_Norm(:,counter)=...
        SVD_Plot_data(:,counter)/max(SVD_Plot_data(:,counter));
end

%% reconstruct image
g=g(:);
SNR=15; %dB
g_noise=mean(abs(g))*10^(-SNR/20)*(1/sqrt(2))*(randn(size(g))+1.0i*randn(size(g)));
g=g%+g_noise; %!!!!!!!!!!!!
f_est_MF(:,counter)=H_MCr'*g;
mat=H_MCr'*H_MCr;
mat=mat+(norm(mat)*10^(-SNR/10))*eye(size(mat));
% f_est_CGS(:,counter) = cgs(mat, f_est_MF(:,counter), 1e-3, 20);
f_est_CGS(:,counter) = pcg(mat, f_est_MF(:,counter), 1e-3, 20);
mat=0;

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
plot1 = plot(resort_plot_data,'LineWidth',3,'Parent',axes1);
set(plot1(1),'DisplayName',[legend_text,num2str(freq_value(5))],...
    'Color',color_palette(1,:));
set(plot1(2),'DisplayName',[legend_text,num2str(freq_value(4))],...
    'Color',color_palette(2,:));
set(plot1(3),'DisplayName',[legend_text,num2str(freq_value(3))],...
    'Color',color_palette(3,:));
set(plot1(4),'DisplayName',[legend_text,num2str(freq_value(2))],...
    'Color',color_palette(4,:));
 set(plot1(5),'DisplayName',[legend_text,num2str(freq_value(1))],'Color',color_palette(5,:));
%   set(plot1(6),'DisplayName',[legend_text,num2str(freq_value(2))],'Color',color_palette(6,:));
%  set(plot1(7),'DisplayName',[legend_text,num2str(freq_value(7))],'Color',[0 0 0.600000011920929]);
set(plot1(6),'DisplayName',[ legend_text,num2str(freq_value(6)), ' (BQ/( f0))'],'LineStyle','--','Color',[1 0 0]);
 
% Create title
title(['Normalized SVD Vs. Number of Freq. Pts. (Q = ', num2str(Q_cavity), ', N & M = ', num2str(num_a1),' Irises)']);

% Create legend
legend(axes1,'show');
%%
% hline = findobj(gcf, 'type', 'line');
% color_palette=bone(6);
% set(hline(2),'LineStyle','-','Color',color_palette(1,:));
% set(hline(3),'LineStyle','-','Color',color_palette(2,:));
% set(hline(4),'LineStyle','-','Color',color_palette(3,:));
% set(hline(5),'LineStyle','-','Color',color_palette(4,:));
%  set(hline(6),'LineStyle','-','Color',color_palette(5,:));
% %  set(hline(7),'LineStyle','-','Color',color_palette(6,:));
% %  set(plot1(7),'DisplayName',[legend_text,num2str(freq_value(7))],'Color',[0 0 0.600000011920929]);
% set(hline(1),'LineStyle','--','Color',[1 0 0]);
%  
% % Create title
% title(['Normalized SVD Vs. Number of Freq. Pts. (Q = ', num2str(Q_cavity), ', N & M = ', num2str(num_a1),' Irises)']);
% 
% % Create legend
% legend(axes1,'show');

%% plot Marchenko-pastur distribution
figure(15); clf;
for i=1:size(SVD_Plot_data,2)
subplot(4,2,i)
samp_const=.25;
d=freq_value(i)/nfreqpts_match;
MarchenkoPasturOverlayHist(SVD_Plot_data(~isnan(SVD_Plot_data(:,i)),i),50,d); %1/4
 title(['# Freq =',num2str(freq_value(i)),' D= ',num2str(d)]);
 xlim([0 3])
end
%% plot reconstruction
figure(28); hold on; clf;
colormap gray
for i=1:size(f_est_MF,2)
 subplot(2,size(f_est_MF,2),size(f_est_MF,2)+i);

% scatter3(imgDomain.locs(:,1),...
%                  imgDomain.locs(:,2), ...
%                  imgDomain.locs(:,3), ...
%                  30, ...
%                  (abs(f_est_MF).^2),...
%                  'filled'); 
%          title('MF');
%          view(-90,0);
%          axis equal
%          axis tight
imagesc(unique(imgDomain.locs(:,2)),...
                 unique(imgDomain.locs(:,3)), ...
                 squeeze(reshape((abs(f_est_MF(:,i)).^2), domainGrid(1), domainGrid(3), domainGrid(2)))...
                 );
              title('MF');
              set(gca, 'XDir','reverse')
              axis equal
              axis tight

figure(28); 
subplot(2,size(f_est_MF,2),i);
cla; 
% scatter3(imgDomain.locs(:,1),...
%                  imgDomain.locs(:,2), ...
%                  imgDomain.locs(:,3), ...
%                  30, ...
%                  (abs(f_est_CGS).^2),...
%                  'filled'); 
% title('CGS'); view(-90,0); axis equal; axis tight;
             
imagesc(unique(imgDomain.locs(:,2)),...
        unique(imgDomain.locs(:,3)), ...
                 squeeze(reshape((abs(f_est_CGS(:,i)/max(abs(f_est_CGS(:,i))))), domainGrid(1), domainGrid(3), domainGrid(2)))...
                 );
             set(gca, 'XDir','reverse'); axis equal; axis tight
             title('CGS');

end
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