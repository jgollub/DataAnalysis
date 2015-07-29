%% define imaging volume
rmin = 0.98;
rmax = 1.1;
nfreqs = 401;
reg = 8E-8;
take_bkg = 0;
wave = 0;

%% import measurement matrix
%filename = 'C:\Users\Lab\Documents\Metaimaging\MeasurementMatricies\FullVectorField_1-2m';
[H, F, Az, El, Z] = measurement_matrix_build_3D1of2PanelVerticalVaryFrqPts(rmin,rmax,nfreqs,-22*pi/180,22*pi/180,-26*pi/180,26*pi/180,data);%filename);
%H rows are: all frequencies, switch 1; all frequencies switch 2 ...
%H columns are: all Elevation (y cosine) angles for Azimuth 1, z 1; all Elevation (y cosine) angles for Azimuth 2, z 1...

%% convert H to wavelet basis
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

%% initialize instruments
delete(instrfind) %delete any existing instrament objects 
vobj_switch = agilent_11713C_switchdriver_startVISAObject; %open switch communications
vobj_vna    = agilent_E8364B_NA_startVISAObject;           %open vna communications
[buffersize, f] = Initialize_8364B(vobj_vna, fsamples, fstart, fstop); % setup VNA scan
agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,1,[1:6])

%% take background401
if take_bkg
    fprintf('%s','Collecting background401 data...')
    background401 = zeros(length(F),6);
    avg_factor = 5;
    for sn=1:6
        agilent_11713C_switchdriver_openChannelNumbers(vobj_switch,1,sn)
        for an=1:avg_factor
            %collect data from VNA
            background401(:,sn) = background401(:,sn) + transpose(Read_8364B(vobj_vna,buffersize));
        end
        background401(:,sn) = background401(:,sn)./avg_factor;
        agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,1,sn)

        figure(1)
        subplot(3,2,sn)
        xlabel('frequency (GHz)')
        ylabel(['S_{7' num2str(sn) '} dB'])
        plot(F/1E9,20.*log10(abs(background401(:,sn))));
        drawnow
    end
    fprintf('%s\n','done.')
end

%initialize measurement matrix to full 401 pts (later will be reduced to
%less points depending on query)
data1p=data;

% Set up the movie.
writerObj = VideoWriter('C:\Users\Lab\Dropbox\MetaImager Project\MetaImager Data\MetaImager Scenes\Vary Num Freq Modes for  Reconstructed Image\Sphere1-401FreqptsA'); % Name it.
writerObj.FrameRate = 1; % How many frames per second.
open(writerObj); 

%vary number of points used to reconstruct image
%for rangeval=[4]
for rangeval=[40 20 8 4 2 1];
nfreqs=400/rangeval +1

%pick out background values to normalize with
background=background401(1:rangeval:end,:);

%pick out matrix field values
data1p.Ey_intp_origin_panel_center=data.Ey_intp_origin_panel_center(:,:,:,1:rangeval:end,:);
%pick out frequency points
data1p.F=data.F(1:rangeval:end);
%run image reconstruct
Realtime_image_reconstruction_revA3_1of2panelVaryFrqPts(rmin,rmax,nfreqs,data1p,background)

   %make movie frame
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end

close(writerObj);