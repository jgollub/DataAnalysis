
home

reg = 1E-9;
tol = 1E-9;

load('C:\Users\John\Dropbox\MetaImager Data\MetaImager Scenes\2PanelTx_hornRx_Scene_Data\2PaneTx_hornRx_Cross.mat')
H=saveddata.H;
g=saveddata.g;

Z = linspace(1,1.15,8);
[Az El] = meshgrid(linspace(-23,23,32).*(pi/180),linspace(-23,23,32).*(pi/180));

%% TV minimization
S = svd(H);  %move up and out
Hp = H./max(S);
gp = g./max(S);
lam = min(S)/max(S);
Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));
[obj obj_debias objfunc] = TwIST(gp,Hp,reg,'lambda',lam,'ToleranceA',tol,'Verbose',0,'Phi',Phi,'Psi',Psi);

for zn=1:length(Z)
    obj3D(:,:,zn) = reshape(obj((1:32^2)+(zn-1)*32^2),32,32);
end
%% plot
%voxelimage_projections(robj(:,:,:,3),Az,El,Z)
for n=1:8
    subplot(4,2,n)
    imagesc(abs(obj3D(:,:,n)))
    axis equal;axis tight;colormap gray
end
