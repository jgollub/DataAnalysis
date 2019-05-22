%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                         %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Bistatic reconstruction code for        %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      3rd floor SAR setup                %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%          (using VNA)                     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      Aaron Diebold, April 2019           %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    Requires lots of RAM, use server!     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

c = physconst('lightspeed');

%% Load data & coordinates
load('C:\Users\avd13\Desktop\Aaron_Imaging_Files\Data\Bistatic\TwoSphere_RaisedFoamTable_09-Apr-2019_12_20.mat');
g = data.measurements; 
f = data.f;
Yr = data.Y1/100;
Yt = data.Y2/100;
Zt = data.Z2/100;
clear data;
load('C:\Users\avd13\Desktop\Aaron_Imaging_Files\Data\Bistatic\Background_RaisedFoamTable_08-Apr-2019_20_3.mat');
bg = data.measurements; clear data;
load('C:\Users\avd13\Documents\MIMO Imaging\two_horns.mat')
load('C:\Users\avd13\Documents\MIMO Imaging\Calibration\cob1.mat')
load('C:\Users\avd13\Documents\MIMO Imaging\Calibration\cob2.mat')

f_indx = 1:4:201;
nf = numel(f_indx);
f = f(f_indx);
lam = c./f;
k = 2*pi./lam;
nm = numel(Yt);

horn_phase = interp1(two_horns.f,two_horns.phase,f);
% horn_offset = 0.0381;

g = g(:,f_indx);
bg = bg(:,f_indx);

g = g./repmat(horn_phase,[size(g,1) 1]);                                   % Compensate for horn phase in signal
bg = bg./repmat(horn_phase,[size(bg,1) 1]);

S = g - bg;
S = reshape(S, [nm*nf 1]);

coords_R = [zeros(1,numel(Yr)); Yr.'; zeros(1,numel(Yr)); ones(1,numel(Yr))];
coords_T = [zeros(1,numel(Yt)); Yt.'; Zt.'; ones(1,numel(Yt))];
coords_R = cob1*coords_R;                                                  % Change of basis from Creaform data
coords_T = cob2*coords_T;
% coords_R(1,:) = coords_R(1,:) + horn_offset;                               % Compensate for horn length in coordinates
% coords_T(1,:) = coords_T(1,:) + horn_offset;

figure
scatter3(coords_R(1,:),coords_R(2,:),coords_R(3,:));
hold on
scatter3(coords_T(1,:),coords_T(2,:),coords_T(3,:));
hold on
scatter3(0,0,0,'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');

h = @(x,y,z,kk) exp(-1j.*kk.*sqrt(x.^2 + y.^2 + z.^2))./(sqrt(x.^2 + y.^2 + z.^2));

%% Define target locations (or region center for experiment) 
x_tar = 1.1;
y_tar = mean(coords_R(2,:));
z_tar = coords_T(3,end);

%% Simulate data?
% S = zeros(nm,nf);
% 
% for ii = 1:nm
%     for jj = 1:nf
%         for kk = 1:numel(x_tar)
%             S(ii,jj) = S(ii,jj) + h(x_tar(kk)-coords_T(1,ii),y_tar(kk)-coords_T(2,ii),z_tar(kk)-coords_T(3,ii),k(jj)).*h(x_tar(kk)-coords_R(1,ii),y_tar(kk)-coords_R(2,ii),z_tar(kk)-coords_R(3,ii),k(jj));
%         end
%     end
% end
% S = reshape(S,[nm*nf 1]);

%% Define scene and H matrix
x_offset = mean(x_tar);
y_offset = mean(y_tar);
z_offset = mean(z_tar);
x_range = -.2:2*mean(lam):.2;
y_range = -.3:2*mean(lam):.3;
z_range = -.3:2*mean(lam):.3;

x_list = x_range + x_offset;
y_list = y_range + y_offset;
z_list = z_range + z_offset;

n_voxels = numel(x_range)*numel(y_range)*numel(z_range);

fprintf('Number of elements in H matrix = %d \n',nm*nf*n_voxels)

[Xr1,~,~,~] = ndgrid(coords_R(1,:),x_list,y_list,z_list);
[Yr1,~,~,~] = ndgrid(coords_R(2,:),x_list,y_list,z_list);
[Zr1,~,~,~] = ndgrid(coords_R(3,:),x_list,y_list,z_list);

[Xt1,Xs,Ys,Zs] = ndgrid(coords_T(1,:),x_list,y_list,z_list);
[Yt1,~,~,~] = ndgrid(coords_T(2,:),x_list,y_list,z_list);
[Zt1,~,~,~] = ndgrid(coords_T(3,:),x_list,y_list,z_list);

tic
for ii = 1:nf
    H_bi(ii,:,:,:,:) = h(Xs-Xt1,Ys-Yt1,Zs-Zt1,k(ii))...
            .*h(Xs-Xr1,Ys-Yr1,Zs-Zr1,k(ii));
end
toc
H_bi = permute(H_bi,[2 1 3 4 5]);

H = reshape(H_bi,[nm*nf n_voxels]);

%% Load new data with same H matrix?
% clear g bg S
% load('C:\Users\avd13\Desktop\Aaron_Imaging_Files\Data\Bistatic\TwoSphere_RaisedFoamTable_09-Apr-2019_12_20.mat');
% g = data.measurements; clear data;
% load('C:\Users\avd13\Desktop\Aaron_Imaging_Files\Data\Bistatic\Background_RaisedFoamTable_08-Apr-2019_20_3.mat');
% bg = data.measurements; clear data;
% 
% g = g(:,f_indx);
% bg = bg(:,f_indx);
% 
% g = g./repmat(horn_phase,[size(g,1) 1]);
% bg = bg./repmat(horn_phase,[size(bg,1) 1]);
% 
% S = g - bg;
% S = reshape(S, [nm*nf 1]);

%% Reconstruct
f_est = H'*S;                   % Matched filter
% f_est = pinv(H)*S;            % SVD truncation pseudoinverse
% f_est = tikRecon_L(H,S);      % Tikhonov-regularized pseudoinverse
% f_est = inv(H'*H)*H'*S;       % Least squares
% f_est = gmres(ctranspose(H)*H,ctranspose(H)*S,1,1e-03,5);    % GMRES

f_est_im = reshape(f_est,[numel(x_range) numel(y_range) numel(z_range)]);

% im1 = squeeze(mean(abs(f_est_im),3));
im1 = squeeze(abs(f_est_im(:,:,floor(size(f_est_im,3)/2))));
[y_im,x_im,im] = gridrefine(y_list,x_list,im1./max(im1(:)),10);     % Interpolate 2D image by a factor

figure
pcolor(y_im,x_im,mag2db(abs(im))); caxis([-15 0])
shading flat; colorbar;
xlabel('y (m)'); ylabel('x (m)');
ax = gca; ax.FontName = 'Times'; ax.FontSize = 14;

%% 3D Plot
if numel(x_list) > 1
    
plotting.style = 'iso';
    
f_est_im = (f_est_im-min(f_est_im(:)))/...
    (max(f_est_im(:))-min(f_est_im(:)));
f_est_im = single(f_est_im);
    
[XpreInt,YpreInt,ZpreInt] = ...
    meshgrid(x_list,y_list,z_list);
    
xInt = min(XpreInt(:)):0.004:max(XpreInt(:));
yInt = min(YpreInt(:)):0.004:max(YpreInt(:));
zInt = min(ZpreInt(:)):0.004:max(ZpreInt(:));
[XInt,YInt,ZInt] = meshgrid(xInt,yInt,zInt);

imagePreInt = permute(abs(f_est_im)./max(abs(f_est_im(:))),[2 1 3]);
imageInt = interp3(XpreInt,YpreInt,ZpreInt,...
    imagePreInt,...
    XInt,YInt,ZInt,...
    'spline');
    
plotting.image_lin = permute(imageInt,[1 3 2]);
plotting.x = zInt;
plotting.y = yInt;
plotting.z = xInt;

if strcmp(plotting.style,'vol3d')
plotting.fax=figure;
    plotting.vol = vol3d('CData',permute(imageInt.^2,[3 1 2]),...
        'XData',[min(YInt(:)),max(YInt(:))],...
        'YData',[min(ZInt(:)),max(ZInt(:))],...
        'ZData',[min(XInt(:)),max(XInt(:))],...
        'Alpha',permute(imageInt.^3,[3 1 2]));
    axis image
    grid on, box on
    view([46,26])
    xlabel('y'), ylabel('z'), zlabel('x')
    colormap hot; colorbar
    
elseif strcmp(plotting.style,'iso')


plotting.levels = [-3.0,-6.0,-9.0,-26.0];
plotting.alphas = [ 1.0, 0.5, 0.3, 0.001];
plotting.colors = {[0.3,0.2,0.7],[0.3,0.4,0.7],[0.3,0.6,0.7],[0.3,0.8,0.7]};
plotting.fax = plotIsoSurf(plotting);

view(29, 9)
xlabel('z'), ylabel('y'), zlabel('x')

end

end









