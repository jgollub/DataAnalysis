    Hbasis =H;
    %% TV minimization
    S = svd(Hbasis);  %move up and out
    Hbasis = Hbasis./max(S);
    g = g./max(S);
    lam = min(S)/max(S);
    Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
    Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z));
    [obj obj_debias objfunc] = TwIST(g,Hbasis,reg,'lambda',lam,'ToleranceA',tol,'Verbose',0,'Phi',Phi,'Psi',Psi);
   