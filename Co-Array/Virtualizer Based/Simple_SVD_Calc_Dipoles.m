%% Virtualizer measurement simulation
% initialize panels
x_range=0;
y_range=.5;
z_range=.5;
range_res=1e-2;
crossrange_res=.05;
offset=1;

[imgDomain, ~]  = test_space(offset, x_range, y_range, z_range, range_res,crossrange_res);

panel = create_panel('fsweep', linspace(18e9, 26e9, 101));
panel= panel_feed(panel, 'feedLocs', 1);
probe1 = create_panel('type', 'dipole', 'fsweep', linspace(18e9, 26e9, 101));
probe1.type
probe2 = create_panel('type', 'dipole', 'fsweep', linspace(18e9, 26e9, 101));


panelfields = dipoles_to_fieldsEXP3(panel, imgDomain);
probeFields1 = dipoles_to_fieldsEXP3(probe1, imgDomain);
probeFields2 = dipoles_to_fieldsEXP3(probe2, imgDomain);

H=makeH_faceted(probeFields1,probeFields2);
Hp=makeH_faceted(panelfields,probeFields2); 

title('SVD')
S=svd(H) ;
Sp=svd(Hp) ;


% cla
figure();

legend('dipole', 'panel')
semilogy(S,'b')
hold on 
semilogy(Sp,'r')
ylim([10^-7,10^0])
%
% % initialize the target
% object = objectCreator('ResTarget', -1, .002, .02);
% off.x  = 1;
% object = objectMover(object, off);
% target = ZBuffer(object, .002);

% visualize the sensor array and target
figure(301);
cla
hold on;
scatter3(probe1.x,probe1.y,probe1.z,'r');
scatter3(probe2.x,probe2.y,probe2.z,'r');
scatter3(imgDomain(:,1), imgDomain(:,2), imgDomain(:,3), 1, 'filled'); axis equal

% run the forward model simulation
% g = forward_model(probe1, probe2, imgDomain.sigma, imgDomain.locs);


%% Image reconstruction
% generate the field-of-view grid and perform least squares image reconstruction

% [locs, grid]  = test_space(1, 0, .25, .25, .005);
% opts.max_its  = 5;
% f_est         = image_recon(probes, panels, g, locs, 'least_squares', opts);
% figure;
% scatter(locs(:,2), locs(:,3), 50, (abs(f_est).^2), 's', 'filled'); axis equal; colorbar;

