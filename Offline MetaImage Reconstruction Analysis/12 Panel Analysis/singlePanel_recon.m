function [obj obj_debias objfunc]=singlePanel_recon(H,g,Az,El,Z,Num_Panels, freqs,probe,panel)

%parameters
regularization=.0003;
tolerance=1e-5;

%% select H and g for individual panel reconstruction
Hpanel_i=H(1+probe*Num_Panels*freqs*(panel-1):probe*Num_Panels*freqs*panel,:);
gpanel_i=g(1+probe*Num_Panels*freqs*(panel-1):probe*Num_Panels*freqs*panel,:);

S = svd(Hpanel_i,'econ');
Hp = Hpanel_i./(max(S));
gpanel_i = gpanel_i./(max(S));

lam = min(S)/max(S);
Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));

%l-1 norm minimization
%[obj obj_debias objfunc] = TwIST(g,Hp,regularization,'lambda',lam,'Initialization',0,'MaxiterA',2000,'StopCriterion',1,'ToleranceA',tolerance,'Verbose',0);

%TV minimization
[obj obj_debias objfunc] = TwIST(gpanel_i,Hp,regularization,'lambda',lam,'ToleranceA',tolerance,'Verbose',0,'Phi',Phi,'Psi',Psi,'Monotone',0,'StopCriterion',1);
end
