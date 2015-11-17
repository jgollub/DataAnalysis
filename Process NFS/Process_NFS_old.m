%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Near Field Scan Data Processing                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Processing Parameters Parameters
%debug plotting?
debug_on=1;

%Output frequencies
f_num=101;
f=linspace(18e9,26.5e9,f_num);

%Number of polarizations
numPol=2;

%constants
c=299792458;    % [m/s] Speed of light

%Folder containing NFS data
Raw_Data_Folder='C:\Users\lab\Documents\MidImager_Data\NFS_data_5';
files=dir([Raw_Data_Folder,'\*.csv']);

%% clean NFS  data (.CSV files) by removing spaces/comments/etc 
mkdir([Raw_Data_Folder,'\CLEANED_CSV']);
mkdir([Raw_Data_Folder,'\CLEANED_CSV\CAL\']);

F_NSI_last=[];
for i=1:length(files)
    file_in=[Raw_Data_Folder,'\',files(i).name];
    file_out=[Raw_Data_Folder,'\CLEANED_CSV\CLEANED_',files(i).name];
    f_NSI=cleanData(file_in,file_out);
    if F_NSI_last~=F_NSI_last
    error('Near field scans do not have the same parameters')
    end
    F_NSI_last=f_NSI;
    fprintf(['FINISHED CSV CLEAN FOR PANEL: ', files(i).name,'\n'])
end

save([Raw_Data_Folder,'\CLEANED_CSV\CAL\NSI_FREQS.mat'],'f_NSI');

%% Probe phase response (NSI system specific)
debug_on=1;


if ~exist('f_NSI')
load([Raw_Data_Folder,'\CLEANED_CSV\CAL\NSI_FREQS.mat'])
disp('WARNING:using loaded f_NSI data')
end

%choose probe correction
use_case=3;


NSI_probe_response=probePhase(f,use_case);

if debug_on
% check phase vs analytic
% probe length
p_length=7.5*2.54/100;
a=0.010668;
m=1; n_indx=0;
k_probe=(2*pi*f/c).*(1-(c./(2*a*f)).^2).^.5;
NSI_probe_response_analytic=exp(-1.0j*k_probe*p_length);

debug_probe=figure(1);
subplot(2,1,1);
plot(f,NSI_probe_response,f,NSI_probe_response_analytic,'--r');
xlabel('frequency'); ylabel('amplitude'); legend('measured','analytic');

subplot(2,1,2)
plot(f,phase(NSI_probe_response),f,phase(NSI_probe_response_analytic),'--r');
xlabel('frequency'); ylabel('phase'); legend('measured','analytic');
end

%% NSI cable phase response (taken from 4 pt measurement)
debug_on=1;
 %file='C:\Users\lab\Documents\MidImager_Data\NFS_OLD_PANELS\Probe_cal\cable_response_101.csv';
file='C:\Users\lab\Documents\MidImager_Data\NFS_data\NSI_CABLE_MEASURED\cable.csv';

[directory,name,ext]=fileparts(file);
file_in=[directory,'\',name,ext];
file_out=[Raw_Data_Folder,'\CLEANED_CSV\CAL\CLEANED_',name,ext];
%clean raw csv file and load
f_cables=cleanData(file_in,file_out);
NSI_cable_measurement=load(file_out);

%Collect appropriate phase info for cables [assuming measuring 2x2 pts]
NSI_cables_response=10.^(NSI_cable_measurement(1:4:numel(f_cables)*4,4)/20).*exp(1.0j*NSI_cable_measurement(1:4:numel(f_cables)*4,5)*pi/180);  

%choose and extract desired frequency values
[test,f_cable_position]=ismember(f,f_cables);
if ~all(test)
   error('frequency mismatch between available cable frequencies and requested frequencies')
end
NSI_cables_response=NSI_cables_response(f_cable_position,:);

if debug_on
   % probe length %note due to the large length of the cable, phase
   % unwrapping is likely to be inconclusive
   cable_length=((20+15+6+6+(15+5)/12)*12*2.54)/100; %47ft of cable (rx 20ft + tx 15ft+6ft +6ft +15in 5in(connectors))
 
   n_indx=1/(.83);
   k_cable=(n_indx*2*pi*f/c);
   NSI_cable_response_analytic=exp(-1.0j*k_cable*cable_length);
   figure(2); hold on;
   title('Cable Response')
   subplot(2,1,1);
   plot(f/1e9, real(NSI_cables_response/max(abs(NSI_cables_response))),'k',f/1e9,real(NSI_cable_response_analytic),'--r')
   xlabel('Frequency (GHz)');
   ylabel('real value');
   subplot(2,1,2);
   plot(f/1e9,phase(NSI_cables_response),f/1e9,phase(NSI_cable_response_analytic),'--r');
   xlabel('Frequency (GHz)')
   ylabel('phase (rad)')
end

%% conector phase 
% Remove calkit
dt=115.881e-12; % [s] Per calkit specsheet
dx=c*dt;        % [m] Calkit pathlength
connector=exp(-1.0j * dx * 2*pi*f/c);
    
if debug_on
debug_connector=figure(3);
subplot(2,1,1);
plot(f,real(connector).','--r');
xlabel('frequency'); ylabel('amplitude'); legend('datasheet');
end

%% Process NSI files; remove excess phase; save as .MAT files
mkdir([Raw_Data_Folder,'\PANEL_FILES_RAW\']);
files=dir([Raw_Data_Folder,'\CLEANED_CSV\*.CSV']);

%make sure to set number of polarizations

for i=1:length(files)
    file_in=[Raw_Data_Folder,'\CLEANED_CSV\',files(i).name];
    file_out=[Raw_Data_Folder,'\PANEL_FILES_RAW\',files(i).name(9:end-4),'.mat'];
    
    %Process each panel and remove excess phase (6+ entries in function)   
    NSI2Panel(file_out,file_in,f_NSI,f,numPol,NSI_probe_response,NSI_cables_response,connector);
   % NSI2Panel(file_out,file_in,f_NSI,f,2,NSI_probe_response, NSI_cables_response);
    fprintf(['FINISHED FILE: ',files(i).name(9:end-4),'.mat \n']);
end

%% Apply NFS rotation to align with physical orientation

p_lt=cell(2,1);
p_rt=cell(2,1);
p_lb=cell(2,1);
p_rb=cell(2,1);
ey=cell(2,1); 
clear Rx Tx X Y measurements f
mkdir(Raw_Data_Folder,'\PANEL_FILES_ALIGNED\')
files=dir([Raw_Data_Folder,'\PANEL_FILES_RAW\*.MAT']);
for i=1:length(files)
    file_in=[Raw_Data_Folder,'\PANEL_FILES_RAW\',files(i).name];
    load(file_in);  %% Load the near-field scan

sampdist=(abs(X(1,1)-X(1,2)))/1000; %NFS sampling period
for loop=1:1:2
   
        figure(4); 
%       range=[-0.0475, -0.0525]; %%back-propagation distance for plane 1 (in meters) - this one is physically measured for the scan
        range=[-0.0595, -0.0555];
       % range=[-0.0595, -0.0645]; %%back-propagation distance for plane 1 (in meters) - this one is physically measured for the scan
           
        ey{loop}=bp(measurements(:,:,:,1),X,Y,range(loop)); %% ey is the summed up back-propagated field - all frequencies and all polarizations
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        %!!! use single polarization to find fiducials
        
        
        %yy=X/1000; 
        %zz=Y/1000;
        
        yy=(X-mean(mean(X,1),2))/1000; 
        zz=(Y-mean(mean(Y,1),2))/1000;
        
        xx=range(loop)*ones(size(yy));
        %plot
        if loop==1
            figure(5);clf;
            figure(4);clf;
        end
        figure(4); hold on;
        subplot(2,2,loop); hold on; cla;
        scatter3(xx(:),yy(:),zz(:),6,20*log10(abs(ey{loop}(:))),'filled')
        axis image; colormap('hot');set(gcf,'color','w');
        title(['Back-propogate ', num2str(range(loop)),' cm']);
        xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
        view(90,0); hold on;

        % process specific panel
           figure(5); 
           
%          [search_fields, search_pos, fid_exact, panel_type] = processMCPanel(ey{loop},xx, yy, zz);
           [search_fields, search_pos, fid_exact, panel_type] = processSlottedMCPanel(ey{loop},xx, yy, zz);
        
%%%%%%%%%%%%%%%%%%%%%%%% For RX PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fiducial_range=0.015; % fiducial guessed to be within this range

search_region = fiducial_region(search_fields,search_pos, fid_exact,fiducial_range);

figure(4);
subplot(2,2,2+loop);  cla; hold on;
scatter3(search_region.lt.x(:),search_region.lt.y(:),search_region.lt.z(:),20,20*log10(abs(search_region.lt.fields)),'filled')
scatter3(search_region.rt.x(:),search_region.rt.y(:),search_region.rt.z(:),20,20*log10(abs(search_region.rt.fields)),'filled')
scatter3(search_region.lb.x(:),search_region.lb.y(:),search_region.lb.z(:),20,20*log10(abs(search_region.lb.fields)),'filled')
scatter3(search_region.rb.x(:),search_region.rb.y(:),search_region.rb.z(:),20,20*log10(abs(search_region.rb.fields)),'filled')


scatter3([search_pos.x(1),search_pos.x(1),search_pos.x(1),search_pos.x(1)],...
    [fid_exact.lt.y, fid_exact.rt.y, fid_exact.lb.y,fid_exact.rb.y],...
    [fid_exact.lt.z,fid_exact.rt.z,fid_exact.lb.z,fid_exact.rb.z], ...
    20,'k')

axis image; colormap('hot');set(gcf,'color','w');
xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
view(90,0); hold on;

 title(['Back-propogate ', num2str(range(loop)),' cm']);

% Localize Maxima (Using Dan's findsubmax function)
%%%%%%%%%%%%%%%% FIRST PLANE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos(loop)= findFiducial( search_region,sampdist);

figure(4); subplot(2,2,2+loop);  hold on; 
scatter3(search_region.lt.x(1),pos(loop).left_b(1),pos(loop).left_b(2),20,'g','*')
scatter3(search_region.lt.x(1),pos(loop).right_b(1),pos(loop).right_b(2),20,'g','*')
scatter3(search_region.lt.x(1),pos(loop).left_t(1),pos(loop).left_t(2),20,'g','*')
scatter3(search_region.lt.x(1),pos(loop).right_t(1),pos(loop).right_t(2),20,'g','*')
axis image; colormap('hot');set(gcf,'color','w');
xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
view(90,0); hold on;

figure(5); hold on;
scatter3(search_region.lt.x(1),pos(loop).left_b(1),pos(loop).left_b(2),80,'g','*')
scatter3(search_region.lt.x(1),pos(loop).right_b(1),pos(loop).right_b(2),80,'g','*')
scatter3(search_region.lt.x(1),pos(loop).left_t(1),pos(loop).left_t(2),80,'g','*')
scatter3(search_region.lt.x(1),pos(loop).right_t(1),pos(loop).right_t(2),80,'g','*')
axis image; colormap('hot');set(gcf,'color','w');
xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
view(90,0); hold on;
% ALIGNMENT CORRECTION (3D)

%%%%%%%%%%%%%%%%%%%%%%%First Plane%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_lt{loop}=[range(loop),pos(loop).left_t(1),pos(loop).left_t(2)];
p_rt{loop}=[range(loop),pos(loop).right_t(1),pos(loop).right_t(2)];
p_lb{loop}=[range(loop),pos(loop).left_b(1),pos(loop).left_b(2)];
p_rb{loop}=[range(loop),pos(loop).right_b(1),pos(loop).right_b(2)];

end

%calculate range vectors
r1=(p_lt{2}-p_lt{1});
r1=r1/norm(r1);

r2=(p_rt{2}-p_rt{1});
r2=r2/norm(r2);

r3=(p_lb{2}-p_lb{1});
r3=r3/norm(r3);

r4=(p_rb{2}-p_rb{1});
r4=r4/norm(r4);

%calculate cross range vectors between points
d1=(p_lt{loop}-p_rt{loop});
d1=d1/norm(d1);

d2=(p_lb{loop}-p_rb{loop});
d2=d2/norm(d2);

d3=(p_lt{loop}-p_lb{loop});   
d3=d3/norm(d3);

d4=(p_rt{loop}-p_rb{loop});
d4=d4/norm(d4);
%set up matrix relation A=TB

v_plot=[      p_lt{1},   r1
              p_rt{1},   r2
              p_lb{1},   r3
              p_rb{1},   r4
              p_rt{1},   d1
              p_rb{1},   d2
              p_lb{1},   d3
              p_rb{1},   d4];
figure(5); hold on;             
quiver3(v_plot(:,1),v_plot(:,2),v_plot(:,3),v_plot(:,4),v_plot(:,5),v_plot(:,6),'LineWidth',2,'Color','g');

if panel_type=='Tx'
            
    % cross range and range vectors to transform to
    d1prime=[0, -1, 0];

    d2prime=[0, (fid_exact.lb.y-fid_exact.rb.y),(fid_exact.lb.z-fid_exact.rb.z)];
    d2prime=d2prime/norm(d2prime);
    
    d3prime=[0,  0, 1];
    d4prime=[0,  0, 1];
    
    r1prime=[-1,  0, 0];
    r2prime=[-1,  0, 0];
    r3prime=[-1,  0, 0];
    r4prime=[-1,  0, 0];
    
    %A & B matrix
    A=[d1prime; d2prime; d3prime; d4prime; r1prime; r2prime; r3prime; r4prime];
    B=[d1;d2;d3;d4;r1;r2;r3;r4];
    
elseif panel_type=='Rx'
    % cross range and range vectors to transform to
    d1prime=[0 -1 0];
    d2prime=[0 -1 0];
    d3prime=[0,  (fid_exact.lt.y-fid_exact.lb.y), (fid_exact.lt.z-fid_exact.lb.z)];
    d3prime=d3prime/norm(d3prime);
    d4prime=[0  0 1];
    
    r1prime=[-1  0 0];
    r2prime=[-1  0 0];
    r3prime=[-1  0 0];
    r4prime=[-1  0 0];
    
    %A & B matrix
    A=[d1prime; d2prime; d3prime; d4prime; r1prime; r2prime; r3prime; r4prime];
    B=[d1;d2;d3;d4;r1;r2;r3;r4];
    
end

vprime_plot=[range(1) fid_exact.lt.y fid_exact.lt.z r1prime
             range(1) fid_exact.rt.y fid_exact.rt.z r2prime
             range(1) fid_exact.lb.y fid_exact.lb.z r3prime
             range(1) fid_exact.rb.y fid_exact.rb.z r4prime
             range(1) fid_exact.rt.y fid_exact.rt.z d1prime
             range(1) fid_exact.rb.y fid_exact.rb.z d2prime
             range(1) fid_exact.lb.y fid_exact.lb.z d3prime
             range(1) fid_exact.rb.y fid_exact.rb.z d4prime];
         
quiver3(vprime_plot(:,1),vprime_plot(:,2),vprime_plot(:,3),vprime_plot(:,4),vprime_plot(:,5),vprime_plot(:,6),'LineWidth',2,'Color','b');

B=B';
A=A';
[U,S,V]=svd(A*B'*inv(B*B'));
T=V*U';

% apply the rotation ========================================================
fmax=max(f);
fmin=min(f);

[ry,rx,nf,p]=size(measurements);
magaddsq=zeros(ry,rx);
new_measurements=zeros(size(measurements));

fieldsPrimed=[];
for frno=1:nf
wavl=c/f(frno);
fieldsPrimed(:,:,1)=transfield(measurements(:,:,frno,1),[0,0,range(2)/wavl], wavl/sampdist,wavl/sampdist); %!!! sign of range(2) is implied (-) may need change later
% fieldsPrimed(:,:,2)=transfield(measurements(:,:,frno,2),[0,0,range(2)/wavl], wavl/sampdist,wavl/sampdist);
fieldsPrimed(:,:,2)=0; %!!!!!!!!!!!!!!!!!! note should take in both polarizations
fieldsPrimed(:,:,3)=0;

fieldsPrimed=rotfield(fieldsPrimed,T,wavl/sampdist,wavl/sampdist);

new_measurements(:,:,frno,1)=fieldsPrimed(:,:,1);
new_measurements(:,:,frno,2)=fieldsPrimed(:,:,2);

magaddsq=magaddsq+sqrt(abs(fieldsPrimed(:,:,1)).^2+abs(fieldsPrimed(:,:,2)).^2);

end

figure(7); clf;
subplot(2,2,1) 
scatter3(xx(:),yy(:),zz(:),10,20*log10(reshape(magaddsq(:,:,1,1),[],1)),'filled');

axis image; colormap('hot');set(gcf,'color','w');
xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
view(90,0); hold on;

%translate to zero yy, zz to zero center
%calculate the translation
% solve for translation ===================================================
subplot(2,2,3)
[search_fields, samePosition, ~, ~] = processMCPanel(magaddsq,xx, yy, zz);   

search_region=[]; yyy=[];zzz=[];E_lb=[]; E_rb=[]; E_lt=[];E_rt=[]; clear pos;
[ search_region ] = fiducial_region(search_fields,samePosition, fid_exact,fiducial_range);

scatter3([search_pos.x(1),search_pos.x(1),search_pos.x(1),search_pos.x(1)],...
    [fid_exact.lt.y, fid_exact.rt.y, fid_exact.lb.y,fid_exact.rb.y],...
    [fid_exact.lt.z, fid_exact.rt.z, fid_exact.lb.z,fid_exact.rb.z], ...
    20,'k')

subplot(1,2,2); hold on;
scatter3(search_region.lt.x(:),search_region.lt.y(:),search_region.lt.z(:),20,20*log10(abs(search_region.lt.fields)),'filled')
scatter3(search_region.rt.x(:),search_region.rt.y(:),search_region.rt.z(:),20,20*log10(abs(search_region.rt.fields)),'filled')
scatter3(search_region.lb.x(:),search_region.lb.y(:),search_region.lb.z(:),20,20*log10(abs(search_region.lb.fields)),'filled')
scatter3(search_region.rb.x(:),search_region.rb.y(:),search_region.rb.z(:),20,20*log10(abs(search_region.rb.fields)),'filled')

scatter3([search_pos.x(1),search_pos.x(1),search_pos.x(1),search_pos.x(1)],...
    [fid_exact.lt.y, fid_exact.rt.y, fid_exact.lb.y,fid_exact.rb.y],...
    [fid_exact.lt.z, fid_exact.rt.z, fid_exact.lb.z,fid_exact.rb.z], ...
    20,'k')

axis image; colormap('hot');set(gcf,'color','w');
xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
view(90,0); hold on;

 title(['Back-propogate after rotation cm']);

 %translate to center
 %determine offset
 
posPrime= findFiducial( search_region,sampdist);
  
hold on; 
scatter3(search_region.lt.x(1),posPrime.left_b(1),posPrime.left_b(2),20,'g','*')
scatter3(search_region.lt.x(1),posPrime.right_b(1),posPrime.right_b(2),20,'g','*')
scatter3(search_region.lt.x(1),posPrime.left_t(1),posPrime.left_t(2),20,'g','*')
scatter3(search_region.lt.x(1),posPrime.right_t(1),posPrime.right_t(2),20,'g','*')
axis image; colormap('hot');set(gcf,'color','w');
xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
view(90,0); hold on;

% ALIGNMENT CORRECTION (3D)

%%%%%%%%%%%%%%%%%%%%%%%First Plane%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcluate average offset
fid=zeros(4,2);
fid(1,:)=-(posPrime.left_b  -  [fid_exact.lb.y, fid_exact.lb.z]);
fid(2,:)=-(posPrime.right_b -  [fid_exact.rb.y,fid_exact.rb.z]);
fid(3,:)=-(posPrime.left_t  -  [fid_exact.lt.y, fid_exact.lt.z]);
fid(4,:)=-(posPrime.right_t -  [fid_exact.rt.y,fid_exact.rt.z]);

transvec=[mean(fid,1),0]; %note bp does X Y translation

magaddsq_final=zeros(ry,rx);
for frno=1:nf
wavl=c/f(frno);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHANGED POLARIZATION
fieldsPrimed(:,:,1)=transfield(new_measurements(:,:,frno,2), transvec/wavl, wavl/sampdist,wavl/sampdist); %!!! sign of range(2) is implied (-) may need change later
fieldsPrimed(:,:,2)=transfield(new_measurements(:,:,frno,1), transvec/wavl, wavl/sampdist,wavl/sampdist);
fieldsPrimed(:,:,3)=0;

new_measurements(:,:,frno,1)=fieldsPrimed(:,:,1);
new_measurements(:,:,frno,2)=fieldsPrimed(:,:,2);
magaddsq_final=magaddsq_final+sqrt(abs(fieldsPrimed(:,:,1)).^2+abs(fieldsPrimed(:,:,2)).^2);

end
fprintf(['Finished translating ', files(i).name ,' SHIFT = [', num2str(transvec),']\n']);


figure(9); clf; search_region=[]; yyy=[];zzz=[];E_lb=[]; E_rb=[]; E_lt=[];E_rt=[]; posLast=posPrime;posPrime=[];%%%%%!!!!!!!!!!!!!!!It was figure(8)

[search_fields, samePosition, ~, ~] = processMCPanel(magaddsq_final,xx, yy, zz);   
search_region = fiducial_region(search_fields,samePosition, fid_exact,fiducial_range);

clf;
scatter3(xx(:),yy(:),zz(:),40,20*log10(magaddsq_final(:)),'filled')
        axis normal; colormap('hot');set(gcf,'color','w');
        xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
        view(90,0); hold on;
 hold on;
scatter3([search_pos.x(1),search_pos.x(1),search_pos.x(1),search_pos.x(1)],...
    [fid_exact.lt.y, fid_exact.rt.y, fid_exact.lb.y,fid_exact.rb.y],...
    [fid_exact.lt.z, fid_exact.rt.z, fid_exact.lb.z,fid_exact.rb.z], ...
    20,'k')

posPrime= findFiducial( search_region,sampdist);

hold on; 
scatter3(search_region.lt.x(1),posPrime.left_b(1),posPrime.left_b(2),20,'g','*')
scatter3(search_region.lt.x(1),posPrime.right_b(1),posPrime.right_b(2),20,'g','*')
scatter3(search_region.lt.x(1),posPrime.left_t(1),posPrime.left_t(2),20,'g','*')
scatter3(search_region.lt.x(1),posPrime.right_t(1),posPrime.right_t(2),20,'g','*')


fid=zeros(4,2);

fid(1,:)=-(posPrime.left_b  -  [fid_exact.lb.y, fid_exact.lb.z]);
fid(2,:)=-(posPrime.right_b -  [fid_exact.rb.y,fid_exact.rb.z]);
fid(3,:)=-(posPrime.left_t  -  [fid_exact.lt.y, fid_exact.lt.z]);
fid(4,:)=-(posPrime.right_t -  [fid_exact.rt.y,fid_exact.rt.z]);
mean(fid,1);
% pause();
drawnow

X=yy*1000;
Y=zz*1000;

measurements=new_measurements;
save([Raw_Data_Folder,'\PANEL_FILES_ALIGNED\',files(i).name], 'measurements', 'X', 'Y', 'f');
 fprintf(['FINISHED ROTATING FILE: ',files(i).name(16:end-4),'.mat \n']);

end

%% Optical Scanning and panel placement
Optical_Scan='C:\Users\lab\Documents\MidImager_Data\Optical_Scanner\MidImager_Positions.txt';
Optical_Data=dlmread(Optical_Scan,'\t',0,0);
opt_fiducial = Optical_Data(Optical_Data(:,9)==0,1:3)/1000; %also convert m
opt_fiducial_coded=Optical_Data(Optical_Data(:,9)>7,[1 2 3 9]); %any coded fiducial less than 7 is coord sys or bar 
opt_fiducial_coded(:,1:3)=opt_fiducial_coded(:,1:3)/1000; %also convert m

figure(10); clf;
hold on
% scatter3(opt_fiducial(:,1),opt_fiducial(:,2),opt_fiducial(:,3), 20,'k');
scatter3(opt_fiducial_coded(:,1),opt_fiducial_coded(:,2),opt_fiducial_coded(:,3), 20,'blue');
view(0,90)

fiducial_coord=zeros(size(opt_fiducial,1),3,size(opt_fiducial_coded,1));
for i=1:size(opt_fiducial_coded,1)
fiducial_coord(:,:,i)=bsxfun(@minus, opt_fiducial, opt_fiducial_coded(i,1:3));
end

%left quadrant

clear lt_fiducial rt_quadrant lb_quadrant rb_quadrant rt_fiducial lt_quadrant 
clear lb_fiducial rb_fiducial center_panel center_panel_global center_panel_global
clear lt_fiducial_global rt_fiducial_global lb_fiducial_global rb_fiducial_global
clear panel_type M T
for i=1:size(opt_fiducial_coded,1)
    %upper left
lt_quadrant=fiducial_coord(fiducial_coord(:,1,i)<0 & fiducial_coord(:,2,i)>0, :,i); 
[~,n_indx]=min(sqrt(sum(lt_quadrant.^2,2)));
lt_fiducial(i,:)=lt_quadrant(n_indx,:)+opt_fiducial_coded(i,1:3);

    %upper righ
rt_quadrant=fiducial_coord(fiducial_coord(:,1,i)>0 & fiducial_coord(:,2,i)>0, :,i); 
[~,n_indx]=min(sqrt(sum(rt_quadrant.^2,2)));
rt_fiducial(i,:)=rt_quadrant(n_indx,:)+opt_fiducial_coded(i,1:3);

   %left bottom
lb_quadrant=fiducial_coord(fiducial_coord(:,1,i)<0 & fiducial_coord(:,2,i)<0, :,i); 
[~,n_indx]=min(sqrt(sum(lb_quadrant.^2,2)));
lb_fiducial(i,:)=lb_quadrant(n_indx,:)+opt_fiducial_coded(i,1:3);
  %right bottom
rb_quadrant=fiducial_coord(fiducial_coord(:,1,i)>0 & fiducial_coord(:,2,i)<0, :,i); 
[~,n_indx]=min(sqrt(sum(rb_quadrant.^2,2)));
rb_fiducial(i,:)=rb_quadrant(n_indx,:)+opt_fiducial_coded(i,1:3);

opt_d1=norm(lt_fiducial(i,:)-rt_fiducial(i,:));
opt_d2=norm(lb_fiducial(i,:)-rb_fiducial(i,:));
opt_d3=norm(lt_fiducial(i,:)-lb_fiducial(i,:));
opt_d4=norm(rt_fiducial(i,:)-rb_fiducial(i,:));
opt_d5=norm(lt_fiducial(i,:)-rb_fiducial(i,:));
opt_d6=norm(rt_fiducial(i,:)-lb_fiducial(i,:));

if (opt_d1>opt_d2) & ((opt_d3+opt_d4)/2 >(opt_d1+opt_d2)/2)
    panel_type{i}='Rx';
    v1(i,:)=(rt_fiducial(i,:)-lt_fiducial(i,:));
    v2(i,:)=(rt_fiducial(i,:)-rb_fiducial(i,:));       
    v3(i,:)=cross(v1(i,:),v2(i,:));
    
    center_panel(i,:)=lt_fiducial(i,:)+v1(i,:)/2-v2(i,:)/2;
    
    v1(i,:)=v1(i,:)/norm(v1(i,:));
    v2(i,:)=v2(i,:)/norm(v2(i,:));
    v3(i,:)=v3(i,:)/norm(v3(i,:));

elseif (opt_d3>opt_d4) & ((opt_d2+opt_d1)/2 >(opt_d3+opt_d4)/2)
    panel_type{i}='Tx';
    v1(i,:)=(rt_fiducial(i,:)-lt_fiducial(i,:));
    v2(i,:)=(lt_fiducial(i,:)-lb_fiducial(i,:));
    v3(i,:)=cross(v1(i,:),v2(i,:)); 
    
    center_panel(i,:)=lb_fiducial(i,:)+v1(i,:)/2+v2(i,:)/2;
    
    v1(i,:)=v1(i,:)/norm(v1(i,:));
    v2(i,:)=v2(i,:)/norm(v2(i,:));
    v3(i,:)=v3(i,:)/norm(v3(i,:));
else
    panel_type{i}='';
    center_panel(i,:)=nan(1,3);
    v1(i,:)=nan(1,3);
    v2(i,:)=nan(1,3);
    v3(i,:)=nan(1,3);
end

end

hold on;
for i=1:size(opt_fiducial_coded,1)
    scatter3(lt_fiducial(i,1),lt_fiducial(i,2),lt_fiducial(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
    scatter3(rt_fiducial(i,1),rt_fiducial(i,2),rt_fiducial(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
    scatter3(lb_fiducial(i,1),lb_fiducial(i,2),lb_fiducial(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
    scatter3(rb_fiducial(i,1),rb_fiducial(i,2),rb_fiducial(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
    scatter3(center_panel(i,1),center_panel(i,2),center_panel(i,3), 20,[1-i*1/12 0 i*1/12],'*');
end

%set origin; choose panel here; build transform
cPan=9;
u=v1(cPan,:);
v=v2(cPan,:);
w=v3(cPan,:);

Q=[[[1 0 0].' , [0 1 0 ].', [0 0 1].' ].',-center_panel(cPan,:).';[0 0 0 1]]; %shift to center
 M=[[u.', v.', w.' ].',[0 0 0].';[0 0 0 1]]; %rotate into center panel position
% M=[[u.', v.', w.' ].',-center_panel(cPan,:).';[0 0 0 1]];
R=[0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1];
T=R*M*Q;
% N=inv(T)';
% N=N(1:3,1:3);

center_panel_global=zeros(size(center_panel));
lt_fiducial_global=zeros(size(center_panel));
rt_fiducial_global=zeros(size(center_panel));
lb_fiducial_global=zeros(size(center_panel));
rb_fiducial_global=zeros(size(center_panel));

for i=1:size(opt_fiducial_coded,1)
    
tmp1=T*[center_panel(i,:),1].';
center_panel_global(i,1:3)=tmp1(1:3);

tmp2=T*[lt_fiducial(i,:),1].';
lt_fiducial_global(i,1:3)=tmp2(1:3);

tmp3=T*[rt_fiducial(i,:),1].';
rt_fiducial_global(i,1:3)=tmp3(1:3);

tmp4=T*[lb_fiducial(i,:),1].';
lb_fiducial_global(i,1:3)=tmp4(1:3);

tmp5=T*[rb_fiducial(i,:),1].';
rb_fiducial_global(i,1:3)=tmp5(1:3);
end

figure(11); clf;
hold on
for i=1:size(opt_fiducial_coded,1)
    scatter3(lt_fiducial_global(i,1),lt_fiducial_global(i,2),lt_fiducial_global(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
    scatter3(rt_fiducial_global(i,1),rt_fiducial_global(i,2),rt_fiducial_global(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
    scatter3(lb_fiducial_global(i,1),lb_fiducial_global(i,2),lb_fiducial_global(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
    scatter3(rb_fiducial_global(i,1),rb_fiducial_global(i,2),rb_fiducial_global(i,3), 20,[1-i*1/12 0 i*1/12],'filled');
      scatter3(center_panel_global(i,1),center_panel_global(i,2),center_panel_global(i,3), 20,[1-i*1/12 0 i*1/12],'*');
end

%determine individual panel transformation

for i=1:size(opt_fiducial_coded,1)

    
    %Rx panel
if (opt_d1>opt_d2) & ((opt_d3+opt_d4)/2 >(opt_d1+opt_d2)/2) 
    panel_type{i}='Rx';
    v2(i,:)=(rt_fiducial_global(i,:)-lt_fiducial_global(i,:));
    v3(i,:)=(rt_fiducial_global(i,:)-rb_fiducial_global(i,:));
    v1(i,:)=cross(v2(i,:),v3(i,:));
    
    v1(i,:)=v1(i,:)/norm(v1(i,:));
    v2(i,:)=v2(i,:)/norm(v2(i,:));
    v3(i,:)=v3(i,:)/norm(v3(i,:));

    
    %Tx Panel
elseif (opt_d3>opt_d4) & ((opt_d2+opt_d1)/2 >(opt_d3+opt_d4)/2)
    panel_type{i}='Tx';
    v2(i,:)=(rt_fiducial_global(i,:)-lt_fiducial_global(i,:));
    v3(i,:)=(lt_fiducial_global(i,:)-lb_fiducial_global(i,:));
    v1(i,:)=cross(v2(i,:),v3(i,:));   
    
    v1(i,:)=v1(i,:)/norm(v1(i,:));
    v2(i,:)=v2(i,:)/norm(v2(i,:));
    v3(i,:)=v3(i,:)/norm(v3(i,:));
else
    panel_type{i}='';
    v1(i,:)=nan(1,3);
    v2(i,:)=nan(1,3);
    v3(i,:)=nan(1,3);
end
u=v1(i,:);
v=v2(i,:);
w=v3(i,:);

Qg=[[[1 0 0].' , [0 1 0 ].', [0 0 1].' ].',center_panel_global(i,:).';[0 0 0 1]]; %shift to panel
Mg=inv([[u.', v.', w.' ].',[0 0 0].';[0 0 0 1]]); %inv because we are rotating panels to their positions (opposite of transformation)
T_Panel(:,:,i)=Qg*Mg; %opposite of translation to coordinate system 
end

mkdir([Raw_Data_Folder,'\OPTICAL_SCAN_POSITIONS\']);
save([Raw_Data_Folder,'\OPTICAL_SCAN_POSITIONS\optical_scan_position.mat'], 'lt_fiducial_global', 'rt_fiducial_global', 'lb_fiducial_global', 'rb_fiducial_global','center_panel_global', 'T_Panel');
 fprintf(['FINISHED PROCESSING OPTICAL DATA. OUTPUT: optical_scan_position.mat \n']);

