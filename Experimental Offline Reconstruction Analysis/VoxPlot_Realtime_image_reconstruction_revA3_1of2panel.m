function [obj3D,Az,El,Z,F,g,reg]=VoxPlot_Realtime_image_reconstruction_revA3_1of2panel(frm,obj_saved)
%frm=4;
tic
home
%close all

%% define imaging volume
rmin = 0.98;
rmax = 1.1;
nfreqs = 101;
reg = obj_saved.regularization(frm);
take_bkg = 0;
wave = 0;

% %% import measurement matrix
% %filename = 'C:\Users\Lab\Documents\Metaimaging\MeasurementMatricies\FullVectorField_1-2m';
% [H, F, Az, El, Z] = measurement_matrix_build_3D1of2PanelVertical(rmin,rmax,nfreqs,-22*pi/180,22*pi/180,-26*pi/180,26*pi/180,data1p);%filename);
% %H rows are: all frequencies, switch 1; all frequencies switch 2 ...
% %H columns are: all Elevation (y cosine) angles for Azimuth 1, z 1; all Elevation (y cosine) angles for Azimuth 2, z 1...

%% convert H to wavelet basis
H=obj_saved.H;
Z=obj_saved.Z;
F=obj_saved.F;
Az=obj_saved.Az;
El=obj_saved.El;

%%%%

HW = zeros(size(H));
for n=1:length(Z)
[Fw, W] = HaarWaveletTransform(ones(32,32));
W=inv(W);
HW(:,(1:32^2)+(n-1)*32^2) = H(:,(1:32^2)+(n-1)*32^2)*W;
end

%alpha = 1.67e-04;%10*eigs(H'*H,1);

fsamples = length(F);
fstart = F(1)/1E9;
fstop = F(end)/1E9;

% %% initialize instraments
% delete(instrfind) %delete any existing instrament objects 
% vobj_switch = agilent_11713C_switchdriver_startVISAObject; %open switch communications
% vobj_vna    = agilent_E8364B_NA_startVISAObject;           %open vna communications
% [buffersize, f] = Initialize_8364B(vobj_vna, fsamples, fstart, fstop); % setup VNA scan
% agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,1,[1:6])

% some tricks
%SVD truncation
% k = 80;
% g = S(1:k,1:k)\U(:,1:k)'*g;
% H = V(:,1:k)';

%normalize energy in columns
% for n=1:size(H,2)
%     norm_const(n) = sum(abs(H(:,n)).^2)^(1/2);
%     H(:,n) = H(:,n)./norm_const(n);
% end
    

%% take background
% if take_bkg
%     fprintf('%s','Collecting background data...')
%     background = zeros(length(F),6);
%     avg_factor = 5;
%     for sn=1:6
%         agilent_11713C_switchdriver_openChannelNumbers(vobj_switch,1,sn)
%         for an=1:avg_factor
%             %collect data from VNA
%             background(:,sn) = background(:,sn) + transpose(Read_8364B(vobj_vna,buffersize));
%         end
%         background(:,sn) = background(:,sn)./avg_factor;
%         agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,1,sn)
% 
%         figure(1)
%         subplot(3,2,sn)
%         xlabel('frequency (GHz)')
%         ylabel(['S_{7' num2str(sn) '} dB'])
%         plot(F/1E9,20.*log10(abs(background(:,sn))));
%         drawnow
%     end
%     fprintf('%s\n','done.')
% end
% 
% %% collect image 
% figure(3)
% colormap gray
% 
 tic
 maxv=0;
% ns = 1;
% for n=1:frames
%     
%     G = zeros(length(F),6);
%     g = zeros(length(F),6);
%     fprintf('%s','Collecting data...')
%     for sn=1:6
%         agilent_11713C_switchdriver_openChannelNumbers(vobj_switch,1,sn)
%         %collect data from VNA
%         G(:,sn) = Read_8364B(vobj_vna,buffersize);
%         agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,1,sn)
%         g(:,sn) = G(:,sn) - background(:,sn);
%       
%         figure(1);
%         subplot(3,2,sn);
%          hold all
%         xlabel('frequency (GHz)')
%         ylabel(['S_{7' num2str(sn) '} dB'])
%         %plot(F,20.*log10(abs(G(:,sn)))); hold on
%         plot(F./1E9,20.*log10(abs(g(:,sn))),'r'); 
%         ylim([-90 -40])
%         drawnow
%     end
%     fprintf('%s\n','done.')
%     
%    
%     g = g(:); %now g is a column vector with all fequencies for switch 1, then all freq for switch 2 etc., same as the rows of H 
g=obj_saved.measurement(:,frm);
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
    Hbasis = H*(~wave)+HW*wave;
    obj = TwIST(g,Hbasis,reg,'lambda',1e-4,'ToleranceA',5E-8,'Verbose',0);
    %obj = W*obj;
    fprintf('%s\n','done.')
    
     nmaxv = max(0,max(max(max(abs(obj)))));
    if nmaxv>maxv
        maxv = nmaxv;
    end
    
    %% Reshape and plot
    figure(frm)
    nae = size(El,1)*size(Az,2);
    for nz=1:length(Z)
        
        objd = obj((1:nae)+nae*(nz-1));
        if wave  
           objd = W*objd; 
        end
        obj3D(:,:,nz) = reshape(objd,size(El,1),size(Az,2));
        %subplot(ceil(sqrt(length(Z))),ceil(sqrt(length(Z))),nz)
        subplot(5,2,nz)
        imagesc(flipud(abs(obj3D(:,:,nz))),'Xdata',[Az(1,1) Az(1,end)].*180/pi,'YData',[El(end,1) El(1,1)].*180/pi)
        colormap gray
        %axis equal
        caxis([0 nmaxv]); 
        axis xy
        xlabel('Azimuth (degrees)')
        ylabel('Elevation (degrees)')
        title(['range: ', num2str(Z(nz))])
      
    end
    fps = 1/toc;
   
   % sv = input(['Press ''s'' to save.'], 's');
   % regn = input(['Current regularization is ' num2str(reg) '. Enter a number to change: ']);
   % if wave
   %    basis = 'wavelet';
   % else
   %    basis = 'Haar';
   % end
   % waven = input(['Current basis is ' basis '. Enter 1 to change: ']);
    
   % if strcmp(sv,'s')
   %     if ns==1;
   %         obj_saved.F = F;
   %         obj_saved.Az = Az;
   %         obj_saved.El = El;
   %         obj_saved.Z = Z;
   %         obj_saved.H = H;
   %     end
   %     obj_saved.measurement(:,ns) = g;
   %     obj_saved.regularization(ns) = reg;
   %     obj_saved.reconstructed(:,:,:,ns) = obj3D;
   %     ns = ns+1;
   % end
   % if ~isnan(regn)
   %    reg = regn;
   % end
   % if waven
   %    wave = ~wave;
   % end
    
    %home
   

    

