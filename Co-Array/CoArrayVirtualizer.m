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
f=linspace(18e9, 26e9, 801);
Q_cavity=100;
aperture_size=.2;
element_size=aperture_size/(length(array1)-1);
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

num_a1=6;
num_a2=num_a1;

pos1=randperm(numel(array1));
pos1=pos1(1:num_a1);

pos2=randperm(numel(array2));
pos2=pos2(1:num_a2);

rnd_array1=zeros(size(array1));
rnd_array1(pos1)=1;

rnd_array2=zeros(size(array2));
rnd_array2(pos2)=1;

panel1_random=create_panel2('type','CoArray','CoArrayElements',rnd_array1,...
                     'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',element_size,...
                     'ElementsizeZ',element_size ...
                      );
panel1_random=panel_feed(panel1_random, 'type','random_cavity_phase_amplitude', 'Q_cavity',Q_cavity,'renormalize','self');

panel2_random=create_panel2('type','CoArray','CoArrayElements',rnd_array2,...
                     'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',element_size,...
                     'ElementsizeZ',element_size ...
                      );
panel2_random=panel_feed(panel2_random, 'type','random_cavity_phase_amplitude','Q_cavity',Q_cavity, 'renormalize','self');


figure(10); cla; scatter3(panel1_random.x,panel1_random.y,panel1_random.z,10,'k');
      hold on;   scatter3(panel2_random.x,panel2_random.y,panel2_random.z,10,'r');

probe=create_panel('type', 'dipole','fsweep',f);
probe=panel_offset(probe, [0,aperture_size/2,0]);

figure(20); cla; scatter3(panel1_random.x,panel1_random.y,panel1_random.z,10,'k');
      hold on;   scatter3(probe.x,probe.y,probe.z,50,'r');


% create standard panel (for comparison)
Conv_PanelQ001=create_panel('type','random',...
                        'Q',.001,...
                        'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',element_size,...
                     'ElementsizeZ',element_size ...
                         );
 Conv_Panel001=panel_feed(Conv_PanelQ001);    
 
 Conv_Panel100=create_panel('type','random',...
                        'Q',100,...
                        'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',element_size,...
                     'ElementsizeZ',element_size ...
                         );
 Conv_Panel100=panel_feed(Conv_Panel100);      
 
  Conv_Panel10000=create_panel('type','random',...
                        'Q',10000,...
                        'fsweep',f,...
                     'sizeY',aperture_size,...
                     'sizeZ',aperture_size,...
                     'ElementSizeY',element_size,...
                     'ElementsizeZ',element_size ...
                         );
 Conv_Panel10000=panel_feed(Conv_Panel10000);     
 
%% Define imaging domain (slice)

% define the absolute imaging domain
R_l =0.02;

imgDomain_offset = [1 0 0];
imgDomain_range  = [0 1 1];
                      
[imgDomain, domainGrid] =test_space3(imgDomain_offset,imgDomain_range,R_l);

figure(12); hold on; scatter3(imgDomain.locs(1:2:end,1),imgDomain.locs(1:2:end,2),imgDomain.locs(1:2:end,3),5,'r'); 
axis equal


%% ================== calculating the fields ========================
% panel1fields = dipoles_to_fieldsEXP3(panel1, imgDomain.locs);
% panel2fields = dipoles_to_fieldsEXP3(panel2, imgDomain.locs);

% panel1Bfields = dipoles_to_fieldsEXP3(panel1B, imgDomain.locs);
% panel2Bfields = dipoles_to_fieldsEXP3(panel2B, imgDomain.locs);

panel1_rnd_fields = dipoles_to_fieldsEXP3(panel1_random, imgDomain.locs);
panel2_rnd_fields = dipoles_to_fieldsEXP3(panel2_random, imgDomain.locs);
probe_fields= dipoles_to_fieldsEXP3(probe, imgDomain.locs);
% 
% convPanelFields001 = dipoles_to_fieldsEXP3(Conv_Panel001, imgDomain.locs);
% convPanelFields100 = dipoles_to_fieldsEXP3(Conv_Panel100, imgDomain.locs);
% convPanelFields10000 = dipoles_to_fieldsEXP3(Conv_Panel10000, imgDomain.locs);



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

H_probe=makeH_faceted(panel1_rnd_fields,probe_fields);
S_probe=svd(H_probe);
semilogy(S_probe,'green'); hold on;

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