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
Q_cavity=20000; %10000
% f=linspace(18e9, 26e9, floor(Q_cavity/2.5));
fbegin=17e9;
fend=26.5e9;
nfreqpts=floor(2*(fend-fbegin)*Q_cavity/(fbegin+fend));
disp(sprintf('number of measurement modes=%g\n',nfreqpts));
sampconst=0.5;
f=linspace(fbegin, fend, floor(nfreqpts)*sampconst);
% f=linspace(fbegin, fend,101);

num_a1=floor(sqrt(nfreqpts));

num_a2=num_a1;
num_a1*num_a2

aperture_size=(3e8/fend)*num_a1/2;
disp(sprintf('spacing between points=%g',aperture_size/num_a1));
disp(sprintf('wavelength=%g',3e8/fend));

% element_size=aperture_size/(length(array1)-1);
% 
% panel1=create_panel2('type','CoArray','CoArrayElements', array1,...
%                      'fsweep',f,...
%                      'sizeY',aperture_size,...
%                      'sizeZ',aperture_size,...
%                      'ElementSizeY',element_size,...
%                      'ElementsizeZ',element_size ...
%                       );
%  panel1=panel_feed(panel1,'type','random_cavity_phase_amplitude', 'Q_cavity',Q_cavity,'renormalize','self');
% 
% panel2=create_panel2('type','CoArray','CoArrayElements',array2,...
%                      'fsweep',f,...
%                      'sizeY',aperture_size,...
%                      'sizeZ',aperture_size,...
%                      'ElementSizeY',element_size,...
%                      'ElementsizeZ',element_size ...
%                       );
% panel2=panel_feed(panel2, 'type','random_cavity_phase_amplitude','Q_cavity',Q_cavity,'renormalize','self');
% 
% figure(12); scatter3(panel1.x,panel1.y,panel1.z,10,'k'); 
% hold on;    scatter3(panel2.x,panel2.y,panel2.z,20,'r');
% legend('panel 1', 'panel 2')


% panel1B=create_panel2('type','CoArray','CoArrayElements', array1,...
%                      'fsweep',f,...
%                      'sizeY',.2,...
%                      'sizeZ',.2,...
%                      'ElementSizeY',.2/19,...
%                      'ElementsizeZ',.2/19 ...
%                       );
%  panel1B=panel_feed(panel1B,'type','random_cavity_phase', 'renormalize','self');
% 
% panel2B=create_panel2('type','CoArray','CoArrayElements',array2,...
%                      'fsweep',f,...
%                      'sizeY',.2,...
%                      'sizeZ',.2,...
%                      'ElementSizeY',.2/19,...
%                      'ElementsizeZ',.2/19 ...
%                       );
% panel2B=panel_feed(panel2B, 'type','random_cavity_phase','renormalize','self');

% num_a1=length(find(array1))
% num_a2=length(find(array2))

% num_a1=30;

%==========Random distribution
% 
% 
% array1=zeros(floor(aperture_size/element_size));
% array2=array1;
% num_a1=floor(sqrt((f(end)-f(1))*Q_cavity/mean(f)));
% num_a2=num_a1;
% 
% pos1=randperm(numel(array1));
% pos1=pos1(1:num_a1);
% 
% pos2=randperm(numel(array2));
% pos2=pos2(1:num_a2);
% 
% rnd_array1=zeros(size(array1));
% rnd_array1(pos1)=1;
% 
% rnd_array2=zeros(size(array2));
% rnd_array2(pos2)=1;
% 
% element_size=aperture_size/(num_a1-1);
% 
% panel1_random=create_panel2('type','CoArray','CoArrayElements',rnd_array1,...
%                      'fsweep',f,...
%                      'sizeY',aperture_size,...
%                      'sizeZ',aperture_size,...
%                      'ElementSizeY',element_size,...
%                      'ElementsizeZ',element_size ...
%                       );
% panel1_random=panel_feed(panel1_random, 'type','random_cavity_phase_amplitude', 'Q_cavity',Q_cavity,'renormalize','self');
% 
% panel2_random=create_panel2('type','CoArray','CoArrayElements',rnd_array2,...
%                      'fsweep',f,...
%                      'sizeY',aperture_size,...
%                      'sizeZ',aperture_size,...
%                      'ElementSizeY',element_size,...
%                      'ElementsizeZ',element_size ...
%                       );
% panel2_random=panel_feed(panel2_random, 'type','random_cavity_phase_amplitude','Q_cavity',Q_cavity, 'renormalize','self');
% 
% 
% figure(10); cla; scatter3(panel1_random.x,panel1_random.y,panel1_random.z,10,'k');
%       hold on;   scatter3(panel2_random.x,panel2_random.y,panel2_random.z,10,'r');
% 
%  probe=create_panel('type', 'dipole','fsweep',f);
% probe=panel_offset(probe, [0,aperture_size/2,0]);
% 
% figure(20); cla; scatter3(panel1_random.x,panel1_random.y,panel1_random.z,10,'k');
%       hold on;   scatter3(probe.x,probe.y,probe.z,50,'r');
%==========Experimental ===============================


% 
%   % Near Field Scan (NFS) data location for panels
% filepathTx='C:\Users\Jonah Gollub\Documents\MetaImagerData\SemiRigid_NFS_scans\';
% filepathRx='C:\Users\Jonah Gollub\Documents\MetaImagerData\Dipole_Antenna_Scans\';
% folder='D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Panel Lookup Tables\Saved Data';
% 
% fprintf('NFS scan data: \n\t %s \nLookup table will be saved: \n\t %s\n',filepathTx, folder)
% 
% 
%     % determine file name
%     [ pan_index,pan_letter]=ind2sub([4 3],1);
%     filename=[char(64+pan_letter),num2str(pan_index)];
%     FullPath=dir([filepathTx,filename,'.mat']);
%     fprintf(['Loading Tx file: ',FullPath.name,'\n']);
%     
%     % load NFS data
%     NFS_data=load([filepathTx,FullPath.name]);
%     
%     % Create dipole equivalent panel
%     importedPanel = import_panel(NFS_data.data);
% 
% % for je=1:6
%     importedPanel=subpanel(importedPanel, 1);
%     importedPanel=frequency_clip(importedPanel)
% % end
    %=======Mills Cross==========================
% num_a1=10;

element_size=aperture_size/(num_a1-1);

MCr_array1=zeros(num_a1);
MCr_array1(2:2:end,1)=1;
MCr_array1(1:2:end,end)=1;

MCr_array2=zeros(num_a2);
MCr_array2(1,2:2:end)=1;
MCr_array2(end,1:2:end)=1;

panel1_MillsC=create_panel2('type','CoArray','CoArrayElements',MCr_array1,...
                     'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',element_size,...
                     'ElementsizeZ', element_size ...
                      );
                  
panel1_MillsC=panel_feed(panel1_MillsC, 'type','random_cavity_phase_amplitude', 'Q_cavity',Q_cavity,'renormalize','self');

panel2_MillsC=create_panel2('type','CoArray','CoArrayElements',MCr_array2,...
                     'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',element_size,...
                     'ElementsizeZ',element_size ...
                      );
                  
panel2_MillsC=panel_feed(panel2_MillsC, 'type','random_cavity_phase_amplitude', 'Q_cavity',Q_cavity,'renormalize','self');
          
figure(12); cla; scatter3(panel1_MillsC.x,panel1_MillsC.y,panel1_MillsC.z,10,'k');
      hold on;   scatter3(panel2_MillsC.x,panel2_MillsC.y,panel2_MillsC.z,10,'r');

%============================================
% 
% % create standard panel (for comparison)
% Conv_PanelQ001=create_panel('type','random',...
%                         'Q',.001,...
%                         'fsweep',f,...
%                      'sizeY',aperture_size,...
%                      'sizeZ',aperture_size,...
%                      'ElementSizeY',element_size,...
%                      'ElementsizeZ',element_size ...
%                          );
%  Conv_Panel001=panel_feed(Conv_PanelQ001);    
%  
%  Conv_Panel100=create_panel('type','random',...
%                         'Q',100,...
%                         'fsweep',f,...
%                      'sizeY',aperture_size,...
%                      'sizeZ',aperture_size,...
%                      'ElementSizeY',element_size,...
%                      'ElementsizeZ',element_size ...
%                          );
%  Conv_Panel100=panel_feed(Conv_Panel100);      
%  
%   Conv_Panel10000=create_panel('type','random',...
%                         'Q',10000,...
%                         'fsweep',f,...
%                      'sizeY',aperture_size,...
%                      'sizeZ',aperture_size,...
%                      'ElementSizeY',element_size,...
%                      'ElementsizeZ',element_size ...
%                          );
%  Conv_Panel10000=panel_feed(Conv_Panel10000);     
%  
%% Define imaging domain (slice)

% define the absolute imaging domain
x_ext=1;
lambda_min=3E8/max(f);
R_l=.5*x_ext*lambda_min/(aperture_size);
% R_l =0.005;

scene_sampling=x_ext*lambda_min/(element_size);
imgDomain_offset = [x_ext 0 0];
imgDomain_range  = [0  scene_sampling scene_sampling];
                      
[imgDomain, domainGrid] =test_space3(imgDomain_offset,imgDomain_range,R_l);

figure(12); hold on; scatter3(imgDomain.locs(1:2:end,1),imgDomain.locs(1:2:end,2),imgDomain.locs(1:2:end,3),5,'r'); 
axis equal


%% Define imaging target
%% Target 3 =================================================================
res=(3e8/fend)/aperture_size;
target=struct();

theta=[];
theta(1)=.5;
Arc=R_l;
a=.25/(2*pi);
% b=a+res/(2*pi);

r_check=scene_sampling/2;
i=1;
while r_check > a*theta(end)
theta(i+1)=theta(i)+Arc/(a*theta(i));
i=i+1;
end

yy1=a.*theta.*cos(theta);
zz1=a.*theta.*sin(theta);
xx1=ones(1,numel(theta));

 xx2=b.*theta.*cos(theta);
 yy2=b.*theta.*sin(theta);
 zz2=ones(1,numel(theta));
figure(24); cla; scatter3(xx1,yy1,zz1,10,'r')
  hold on; axis equal; scatter3(xx2,yy2,zz2,10,'r')

xx=imgDomain.locs(:,1);
yy=imgDomain.locs(:,2);
zz=imgDomain.locs(:,3);
spiral=zeros(numel(xx),1);

for i=1:numel(xx1)
spiral=logical(sqrt((xx-xx1(i)).^2+(yy-yy1(i)).^2+(zz-zz1(i)).^2)<4*R_l)+spiral;
% spiral=logical(sqrt((xx-xx2(i)).^2+(yy-yy1(i)).^2+(zz-zz1(i)).^2)<2*R_l)+spiral;
end
spiral=logical(spiral);
target.locs=imgDomain.locs(spiral,:);
target.sigma=ones(size(target.locs));


test=figure(12);
 hold on; scatter3(target.locs(:,1),target.locs(:,2),target.locs(:,3),'b'); 
 
% 
% barwidth=(3e8/fend)/aperture_size;
% 
% disp(sprintf('bar width=%g',barwidth));
% object  = objectCreator('Spiral');
% object  = objectCreator('ResTarget',1,target_res,barwidth);
% num_y=floor(scene_sampling/(6*barwidth));
% num_z=floor(scene_sampling/(12*barwidth));
%         target_res   = R_l/3;
%         targetOffs=[];
% 
%         
%         
% target=cell(num_z,num_y);
% 
% for i=1:num_z
%     for j=1:num_y
%     object  = objectCreator('ResTarget',1,target_res,barwidth);
%     targetOffs.x = 1;
%     targetOffs.y = scene_sampling/2-6*barwidth*j;
%     targetOffs.z = scene_sampling/2-12*barwidth*i;   
%        
%         % Move
%         object = objectMover(object,targetOffs);
%  
% % Zbuffer
% target{i,j}  = ZBuffer(object,target_res);
% % target0 = ZBuffer(object,2*target_res);
% 
% test=figure(12);
%  hold on; scatter3(target{i,j}.locs(:,1),target{i,j}.locs(:,2),target{i,j}.locs(:,3),'b'); 
axis equal
view(-90,0);
xlabel('xaxis')
ylabel('yaxis')
zlabel('zaxis')


%% generate g
g=forward_model(panel2_MillsC, panel1_MillsC, target.sigma.*(randn(size(target.sigma))+1i*randn(size(target.sigma))), target.locs);


% gij = forward_model(panel2_MillsC, panel1_MillsC, target{i,j}.sigma.*(randn(size(target{i,j}.sigma))+1i*randn(size(target{i,j}.sigma))), target{i,j}.locs);
% if i==1 & j==1
% g=zeros(size(gij(:)));
% end
% g=g+gij(:);
% 
%     end
% end

%% ================== calculating the fields ========================
% panel1fields = dipoles_to_fieldsEXP3(panel1, imgDomain.locs);
% panel2fields = dipoles_to_fieldsEXP3(panel2, imgDomain.locs);

% panel1Bfields = dipoles_to_fieldsEXP3(panel1B, imgDomain.locs);
% panel2Bfields = dipoles_to_fieldsEXP3(panel2B, imgDomain.locs);

% panel1_rnd_fields = dipoles_to_fieldsEXP3(panel1_random, imgDomain.locs);
% panel2_rnd_fields = dipoles_to_fieldsEXP3(panel2_random, imgDomain.locs);
%  probe_fields= dipoles_to_fieldsEXP3(probe, imgDomain.locs);

panel1_MCr_fields = dipoles_to_fieldsEXP3(panel1_MillsC, imgDomain.locs);
panel2_MCr_fields = dipoles_to_fieldsEXP3(panel2_MillsC, imgDomain.locs);

% convPanelFields001 = dipoles_to_fieldsEXP3(Conv_Panel001, imgDomain.locs);
% convPanelFields100 = dipoles_to_fieldsEXP3(Conv_Panel100, imgDomain.locs);
% convPanelFields10000 = dipoles_to_fieldsEXP3(Conv_Panel10000, imgDomain.locs);

% panel_exp_fields = dipoles_to_fieldsEXP3(importedPanel, imgDomain.locs);

%% ================generate H ======================================
figure(11); 

% H=makeH_faceted(panel1fields,panel2fields);
% S=svd(H) ;
% semilogy(S,'r')

% HB=makeH_faceted(panel1Bfields,panel2Bfields);
% SB=svd(HB) ;
% semilogy(SB,'m')

% H_rnd=makeH_faceted(panel1_rnd_fields,panel2_rnd_fields);
% S_rnd=svd(H_rnd);
% semilogy(S_rnd,'green'); hold on;

% H_probe=makeH_faceted(panel1_rnd_fields,probe_fields);
% S_probe=svd(H_probe);
% semilogy(S_probe,'green'); hold on;

H_MCr=makeH_faceted(panel1_MCr_fields,panel2_MCr_fields);
S_Mc=svd(H_MCr,'econ');
semilogy(S_Mc,'green'); hold on;
% 
% H_exp=makeH_faceted(panel_exp_fields,probe_fields);
% S_exp=svd(H_exp,'econ');
% semilogy(S_exp,'green'); hold on;
% H_Conv001=makeH_faceted(convPanelFields001,convPanelFields001); 
% S_Conv001=svd(H_Conv001);
% hold on; semilogy(S_Conv001,'k')
% 
% H_Conv100=makeH_faceted(convPanelFields100,convPanelFields100); 
% S_Conv100=svd(H_Conv100);
% hold on; semilogy(S_Conv100,'k')
% 
% H_Conv10000=makeH_faceted(convPanelFields10000,convPanelFields10000); 
% S_Conv10000=svd(H_Conv10000);
% hold on; semilogy(S_Conv10000,'k')


% Hpanel=makeH_faceted(panelfields,panelfields2);

title('SVD (Cross-Range Slice)')

legend('CoArray')% ,'CoArray (random)') % ,'Conventional, Q=1', 'Conventional, Q=100', 'Conventional, Q=10000')


% [edges,N]=histcounts(S_probe,10);
figure(15); cla;
MarchenkoPasturOverlayHist(S_Mc,50,.5*sampconst);

%% reconstruct image
g=g(:);
f_est_MF=H_MCr'*g;
f_est_CGS = cgs(H_MCr'*H_MCr, H_MCr'*g, 1e-3, 20);
figure(27); subplot (1,2,1);
cla; scatter3(imgDomain.locs(:,1),...
                 imgDomain.locs(:,2), ...
                 imgDomain.locs(:,3), ...
                 30, ...
                 (abs(f_est_MF).^2),...
                 'filled'); 
title('MF');
view(-90,0);
axis equal
axis tight
figure(27); subplot (1,2,2);
cla; scatter3(imgDomain.locs(:,1),...
                 imgDomain.locs(:,2), ...
                 imgDomain.locs(:,3), ...
                 30, ...
                 (abs(f_est_CGS).^2),...
                 'filled'); 
title('CGS');
view(-90,0);
axis equal
axis tight

%%---------------------

figure(28); subplot (1,2,1);
f_est_MF_grid=reshape(f_est_MF,[domainGrid(1) domainGrid(3)]);

cla; imagesc(imgDomain.locs(:,2),...
                 imgDomain.locs(:,3), ...
                 (abs(f_est_MF_grid).^2)...
                  ); 
              axis equal
              axis tight
title('MF');
axis equal
axis tight
figure(28); subplot (1,2,2);
f_est_CGS_grid=reshape(f_est_CGS,[domainGrid(1) domainGrid(3)]);

cla; imagesc(imgDomain.locs(:,2),...
                 imgDomain.locs(:,3), ...
                 (abs(f_est_CGS_grid).^2)...
                 ); 
title('CGS');
axis equal
axis tight