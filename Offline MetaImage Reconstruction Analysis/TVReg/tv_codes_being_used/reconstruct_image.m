function [obj3D] = reconstruct_image(g, reg, tol, basis, H, Az, El, Z)

maxv=0;
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%     gA = g(1:middle);
%     
%     HA= H(1:mid,:);

    %% psuedo-inverse
    % obj = Hinv*g;

    %% ISTA+1DTV
    %[obj objective] = ista_1DTV(G,H,2.8E-4,5E-1,alpha,2,10);
    % obj = obj.*norm_const'; %uncomment this line if the energy in the measurement matrix columns was normalized

    %% ISTA+point motion
    %[obj objective] = ista_point_motion(G,H,length(distances),length(phi2),2.8E-4,20E-1,alpha,2,10);
    % obj = obj.*norm_const'; %uncomment this line if the energy in the measurement matrix columns was normalized

    %% Twist
    fprintf('%s','TwIST reconstruction...')
 
%     Hbasis = H*(~wave)+HW*wave;
    Hbasis =H;
    %% TV minimization
    S = svd(Hbasis);  %move up and out
    Hbasis = Hbasis./(max(S)*(1+1e-5));
    g = g./(max(S)*(1+1e-5));
    lam = min(S)/max(S);
    Phi = @(x) TVnorm_RI(x,size(El,1),size(Az,2),length(Z));
    Psi = @(x,tau) mycalltoTVnewMatrix_RI(x,tau,10,size(El,1),size(Az,2),length(Z)); % 10 is number of TV iterations per TwIST iterations
    [obj obj_debias objfunc] = TwIST(g,Hbasis,reg,'lambda',lam,'ToleranceA',tol,'Verbose',0,'Phi',Phi,'Psi',Psi,'Monotone',0,'StopCriterion',1);
   
   % [obj obj_debias objfunc] = TwIST(g,Hbasis,reg,'lambda',lam,'ToleranceA',5E-9,'Verbose',0,'Phi',Phi,'Psi',Psi);
    
    %% l1
%     [obj obj_debias objfunc] = TwIST(g,Hbasis,reg,'lambda',1e-4,'ToleranceA',5E-12,'Verbose',0);
    %[obj obj_debias objfunc] = TwIST(g,Hbasis,reg,'lambda',1e-4,'StopCriterion',1,'MaxiterA',2000,'ToleranceA',tol,'Verbose',0);
    %obj = W*obj;
    fprintf('%s\n','done.')
    
    nmaxv = max(0,max(max(max(abs(obj)))));
    if nmaxv>maxv
        maxv = nmaxv;
    end
    
    %% Reshape and plot
    figure(8)
    plot(objfunc)
    figure(7)
    nae = size(El,1)*size(Az,2);
    for nz=1:length(Z)
        
        objd = obj((1:nae)+nae*(nz-1));
%         if wave  
%            objd = W*objd; 
%         end
        obj3D(:,:,nz) = reshape(objd,size(El,1),size(Az,2));
        %subplot(ceil(sqrt(length(Z))),ceil(sqrt(length(Z))),nz)
        subplot(ceil(length(Z)/2),2,nz)
        imagesc(abs(obj3D(:,:,nz)),'Xdata',[Az(1,1) Az(1,end)].*180/pi,'YData',[El(end,1) El(1,1)].*180/pi) %unflipped !!!      
        colormap jet
        %         colormap gray
        %axis equal
        caxis([0 nmaxv]); 
        axis xy; axis equal; axis tight
        xlabel('Azimuth (degrees)')
        ylabel('Elevation (degrees)')
        title(['range: ', num2str(Z(nz))])
      
    end
    figure; imagesc(abs(sum(obj3D,3)));
    %fps = 1/toc
   
    
%     sv = input(['Press ''s'' to save.'], 's');
%     if strcmp(sv,'s')
%         obj_saved(:,:,:,ns) = obj3D;
%         ns = ns+1;
%     end
    
%     changeparams = input('Enter 1 to change parameters: ');
%     while changeparams
%         %% change params
%         regn = input(['Current regularization is ' num2str(reg) '. Enter a number to change: ']);
%         toln = input(['Current tolerance is ' num2str(tol) '. Enter a number to change: ']);
%         waven = input(['Current basis is ' basis '. Enter 1 to change: ']);
%         rerecon = input('Enter 1 to reconstruct image with new parameters: ');
        
%         if wave
%             basis = 'wavelet';
%         else
%             basis = 'Haar';
%         end
% 
%         if ~isnan(toln)
%            tol = toln;
%         end
% 
%         if ~isnan(regn)
%            reg = regn;
%         end
% 
%         if waven
%            wave = ~wave;
%         end
% 
%         if rerecon
%             %%
%                 fprintf('%s','Repeating TwIST reconstruction...')
%                 [obj obj_debias objfunc] = TwIST(g,Hbasis,reg,'lambda',1e-4,'StopCriterion',1,'MaxiterA',2000,'ToleranceA',tol,'Verbose',0);
%                 fprintf('%s\n','done.')
%                 nmaxv = max(0,max(max(max(abs(obj)))));
%                 if nmaxv>maxv
%                     maxv = nmaxv;
%                 end
%                 %Reshape and plot
%                 figure(8)
%                 plot(objfunc)
% 
%                 figure(7)
%                 nae = size(El,1)*size(Az,2);
%                 for nz=1:length(Z)
% 
%                     objd = obj((1:nae)+nae*(nz-1));
%                     if wave  
%                        objd = W*objd; 
%                     end
%                     obj3D(:,:,nz) = reshape(objd,size(El,1),size(Az,2));
%                     %subplot(ceil(sqrt(length(Z))),ceil(sqrt(length(Z))),nz)
%                     subplot(ceil(length(Z)/2),2,nz)
%                     imagesc(abs(obj3D(:,:,nz)),'Xdata',[Az(1,1) Az(1,end)].*180/pi,'YData',[El(end,1) El(1,1)].*180/pi) %unflipped !!!
%                     colormap gray
%                     %axis equal
%                     caxis([0 nmaxv]); 
%                     axis xy; axis equal; axis tight
%                     xlabel('Azimuth (degrees)')
%                     ylabel('Elevation (degrees)')
%                     title(['range: ', num2str(Z(nz))])
%                 end
%         end
    
  %  changeparams = input('Enter 1 to change parameters: ');
        
    end
    


