Num_Panels=12;
numFeeds=6;
probe=1;
panel=6;


H = scene_data(1).H;
g = scene_data(1).obj_saved.measurement;
Z = scene_data(1).Z;
Az = scene_data(1).Az;
El = scene_data(1).El;
freqs_num=length(scene_data(1).F);
tolerance=scene_data(1).obj_saved.tolerance;
%regularization=scene_data(1).obj_saved.regularization;
regularization=.0005

%% select H and g for individual panel reconstruction
Hpanel_i=H(1+(probe-1)*Num_Panels*numFeeds*freqs_num+numFeeds*freqs_num*(panel-1):...
    (probe-1)*Num_Panels*numFeeds*freqs_num+freqs_num*numFeeds*panel,:);
gpanel_i=g(1+(probe-1)*Num_Panels*numFeeds*freqs_num+freqs_num*numFeeds*(panel-1):(probe-1)*Num_Panels*numFeeds*freqs_num+freqs_num*numFeeds*panel);

%TwIST
S=svd(Hpanel_i,'econ');
Hp = Hpanel_i./(max(S));
gpanel_i = gpanel_i./(max(S));
tolerance=5e-5;
lam=5e-5;

[f_est,~,objfunc,~,~,~,~] = TwIST(gpanel_i,Hp,regularization,...
    'lambda',lam,'ToleranceA',tolerance,'StopCriterion',3,'MaxiterA',2e2,'Verbose',0,'ToleranceA',1e-15);

figure(Obj_funct_Plot)
semilogy(abs(objfunc(2:end)-objfunc(1:(end-1))))
xlabel('twist iteration')
ylabel('magnitude of objective function change')
drawnow

% 
% S = svd(Hpanel_i,'econ');
% Hp = Hpanel_i./(max(S));
% gpanel_i = gpanel_i./(max(S));
% 
% lam = min(S)/max(S);
% Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
% Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));

%l-1 norm minimization
%[obj obj_debias objfunc] = TwIST(g,Hp,regularization,'lambda',lam,'Initialization',0,'MaxiterA',2000,'StopCriterion',1,'ToleranceA',tolerance,'Verbose',0);

%TV minimization
%[obj obj_debias objfunc] = TwIST(gpanel_i,Hp,regularization,'lambda',lam,'ToleranceA',tolerance,'Verbose',0,'Phi',Phi,'Psi',Psi,'Monotone',0,'StopCriterion',1);

obj_plot = abs(obj);
nz = length(Z);
nel = size(El,1);
naz = size(Az,2);
nae = nel*naz;
obj3D = zeros(nel,naz,nz);
for zn=1:nz
    objd = obj_plot((1:nae)+nae*(zn-1));
    obj3D(:,:,zn) = reshape(objd,nel,naz);
end

figure()
imagesc(Az(1,:),El(:,1),obj3D(:,:,4))
colorbar
axis xy
