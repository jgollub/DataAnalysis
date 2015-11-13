function [obj obj_debias objfunc]=singlePanel_recon(H,g,Az,El,Z,Num_Panels, freqs,probe,panel)

%parameters
regularization=.0003;
tolerance=1e-5;

%% select H and g for individual panel reconstruction
Hpanel_i=H(1+(probe-1)*Num_Panels*freqs+freqs*(panel-1):(probe-1)*Num_Panels*freqs+freqs*panel,:);
gpanel_i=g(1+(probe-1)*Num_Panels*freqs+freqs*(panel-1):(probe-1)*Num_Panels*freqs+freqs*panel,:);

%TwIST
S=svd(Hpanel_i,'econ');
Hp = Hpanel_i./(max(S));
gpanel_i = gpanel_i./(max(S));
tolerance=5e-5;
lam=5e-5;

[f_est,~,objfunc,~,~,~,~] = TwIST(gpanel_i,Hp,tolerance,...
    'lambda',lam,'StopCriterion',3,'MaxiterA',2e2,'Verbose',0,'ToleranceA',1e-15);

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

%TV minimization
[obj obj_debias objfunc] = TwIST(gpanel_i,Hp,regularization,'lambda',lam,'ToleranceA',tolerance,'Verbose',0,'Phi',Phi,'Psi',Psi,'Monotone',0,'StopCriterion',1);
end
