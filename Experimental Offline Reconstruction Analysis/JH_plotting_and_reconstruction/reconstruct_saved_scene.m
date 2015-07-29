function [obj3D Az El Z objfunc] = reconstruct_saved_scene(scene_data,dataset_num,regularization,tolerance)
% reconstruct and plot saved scene data
%Inputs:
%Outputs:
%obj3D - reconstructed scene

%kinect scene data
if isfield(scene_data,'rgb')
    sn = 1;
    rgb = scene_data(sn).rgb;
    xyz = scene_data(sn).xyz;
    azimuth_min = scene_data(sn).Az_extent(1);
    azimuth_max = scene_data(sn).Az_extent(2);
    elevation_min = -scene_data(sn).El_extent(1);
    elevation_max = -scene_data(sn).El_extent(2);
    z_min = scene_data(sn).Z_extent(1);
    z_max = scene_data(sn).Z_extent(2);
    objs = scene_data(sn).objs;
else
    fprintf('%s\n%s\n','Sorry, there''s no Kinect data in this dataset :(.','We can still look at the RF data though.')
end

% rf data
H = scene_data.H;
g = scene_data.obj_saved(dataset_num).measurement;
Z = scene_data.Z;
Az = scene_data.Az;
El = scene_data.El;


%% scene reconstruction
S = svd(H,'econ');  
Hp = H./(max(S));
g = g./(max(S));

lam = min(S)/max(S);
Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));

%l-1 norm minimization
%[obj obj_debias objfunc] = TwIST(g,Hp,regularization,'lambda',lam,'Initialization',0,'MaxiterA',2000,'StopCriterion',1,'ToleranceA',tolerance,'Verbose',0);
%TV minimization
[obj obj_debias objfunc] = TwIST(g,Hp,regularization,'lambda',lam,'ToleranceA',tolerance,'Verbose',0,'Phi',Phi,'Psi',Psi,'Monotone',0,'StopCriterion',1);

figure
semilogy(abs(objfunc(2:end)-objfunc(1:(end-1))))
xlabel('twist iteration')
ylabel('magnitude of objective function change')
drawnow

%% reshape reconstructed image vector into 3D scene, plot
obj = abs(obj);
nz = length(Z);
nel = size(El,1);
naz = size(Az,2);
nae = nel*naz;
obj3D = zeros(nel,naz,nz);
for zn=1:nz
    objd = obj((1:nae)+nae*(zn-1));
    obj3D(:,:,zn) = reshape(objd,nel,naz);
end

figure

if isfield(scene_data,'rgb')
    %% plot 3D-rgb kinect scene and overlay RF reconstruction
    axes('pos', [0,0.5,1,0.5]);
    view(-27,34)
    drawnow
    set(gca,'CameraViewAngleMode','manual')

    draw_Kinect_object_scene(xyz,rgb,objs,1:size(objs,3))
    draw_extent_box([azimuth_min azimuth_max]*(pi/180),[elevation_min elevation_max]*(pi/180),[z_min z_max],1,[0 0 1],0)
    plot_recon_surface(obj3D,Az,El,Z,0,4);
    
    %% plot RF scene alone on same figure
    axes('pos', [0,0,1,0.5])
    image(zeros(1,1,3));axis off
    axes('pos', [0,0.05,1,0.4]);
end

%% plot RF scene alone
view(-27,34)
drawnow
set(gca,'CameraViewAngleMode','manual')

plot_recon_surface(obj3D,Az,El,Z,0,4);
view(-27,34);
set(gca,'CameraViewAngleMode','manual');
set(gca,'Color',[0 0 0],'XColor',[1 1 1], 'YColor',[1 1 1], 'ZColor',[1 1 1],'DataAspectRatio',[1 1 1]);
axis equal;axis tight;box on;colormap jet

figure
plot_range_slices(obj3D,Az,El,Z,3)