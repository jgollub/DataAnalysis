%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Near Field Scan Image Reconstruction Algorithm
%Jonah Gollub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [image,Lx,Ly,numSampling] = reconstructRF(X, Y, f, measurements, varargin)

% constants
c=2.99792458*10^8;

%parse inputs

defaultPad=2^nextpow2(max(numel(X(1,:)),numel(Y(:,1))));
ReconstructionTypes = {'RMA', 'gHF', 'SUM'};

p = inputParser;
addOptional(p, 'Algorithm', 'RMA', @(x) find(strcmp(findReconstructionTypes,x)));
addOptional(p, 'Pad', defaultPad, @(x) rem(x,1)==0);
addOptional(p, 'z_offset', 0, @(x) isnumeric(x));

parse(p,varargin{:});

pad      = p.Results.Pad;
alg      = p.Results.Algorithm;
z_offset = p.Results.z_offset;

switch alg
    case 'RMA'
        %% Perform RMA algorithm
               
        Lx=X(1,end)-X(1,1);
        [ynum,xnum]=size(X);
        
        dy=abs(Y(2,1)-Y(1,1));
        Ly=Y(end,1)-Y(1,1);
        
        dx=Lx/(pad-1);
        dy=Ly/(pad-1);

        %calc free space k vector
        k=2*pi*f/c;
        k=repmat(permute(k,[1,3,2]),[pad,pad,1]);
        
        [kux,kuy,~]=meshgrid(linspace(-(2*pi)/(2*dx),(2*pi)/(2*dx),pad),...
            linspace(-(2*pi)/(2*dy),(2*pi)/(2*dy),pad),...
            1:numel(f));
        
        %calculate plane wave decomposition (FFT). Note we have to be careful about
        %shifting in fft because we are only acting along two of the dimensions.
        %Therefor it is easier to use circshift (as opposed to ifftshif which acts
        %on all dimentions.
        
        Sxy=fft(fftshift(measurements),[],1);
        Sxy=ifftshift(fft(Sxy,[],2));
        Sxy=padarray(Sxy,[ceil((pad-size(measurements,1))/2), ceil((pad-size(measurements,2))/2)],0,'pre');
        Sxy=padarray(Sxy,[floor((pad-size(measurements,1))/2), floor((pad-size(measurements,2))/2)],0,'post');
        
        % %     plot decomposition
        %     figHandle=figure(1);
        %     scrn = get( groot, 'Screensize');  scrn(1)=2*scrn(3)/3;  scrn(3)=scrn(3)/3;
        %     set(figHandle,'Position',scrn); clf; subplot(2,2,1);
        % %     imagesc(-Lx_pad/2:dx:Lx_pad/2,-Ly_pad/2:dy:Ly_pad/2,abs(ifft2(fftshift(Sxy(:,:,1)))));
        %     imagesc(-Lx/2:dx:Lx/2,-Ly/2:dy:Ly/2,abs(ifftshift(ifft2(fftshift(Sxy(:,:,1))))));
        %     title('Padded Input Fields f(1)'); axis equal; axis tight; xlabel('x (m)'); ylabel('y (m)')
        %     subplot(2,2,2);
        %     imagesc(kux(1,:),kuy(:,1),abs(Sxy(:,:,1)));
        %     title('FFT of Fields f(1)'); axis equal; axis tight; xlabel('kx (m)'); ylabel('ky (m)')
        %     drawnow
        %
        %calculate min max kz wavenumber
        
        kuz=sqrt((2*k).^2-kux.^2-kuy.^2);
        kuz=real(kuz); % ignore evanescent fields
        
        kuzMin = min(kuz(:));
        kuzMax = max(kuz(:));
        
        %ensure sufficient sampling of Kz
        minSampling = min(nonzeros(diff(kuz,1,3))); %only look end of matrix (high frequency region)
        numSampling = max(size(measurements,3),ceil((kuzMax-kuzMin)/minSampling)); %opt 1 (sample at minimum)
        numSampling = 2^nextpow2(numSampling); %opt 2 (sample at next power of 2)
        
        Kz=linspace(kuzMin,kuzMax,numSampling);
        
        Sxy=Sxy.*exp(1.0j*(kuz)*z_offset);
        
        %interpolate to evenly spaced grid
        
        Srmg = zeros(size(kuz,1),size(kuz,2),length(Kz));
        for ii=1:size(kux,2)
            for jj=1:size(kuy,1)
                indx_vec = squeeze(squeeze(real(kuz(jj,ii,:))~=0));
                if ~(sum(indx_vec)==0 || sum(indx_vec)==1) %check that there are enough points to interpolate
                    Srmg(jj,ii,:)=interp1(squeeze(squeeze(kuz(jj,ii,indx_vec))),...
                        squeeze(Sxy(jj,ii,indx_vec)),...
                        Kz.','linear');
                end
            end
        end
        % %debug: plot interpolatio for z-slice
        %             ii=44;
        %             jj=44;
        %             %debug
        %             figure(3); cla;
        %             plot(squeeze(squeeze(kuz(jj,ii,indx_vec))),squeeze(Sxy(jj,ii,indx_vec)))
        %             hold on;
        %             Srmg(jj,ii,find(isnan(Srmg(jj,ii,:))))=0;
        %             plot(real(Kz(:)).',squeeze(squeeze(Srmg(jj,ii,:))),'-o')
        %             drawnow;
        
        Srmg(find(isnan(Srmg))) = 0; %set all Nan values to 0
        
        %apply inverst FFT to get image
        fxy = fftshift(ifftn(Srmg));
        
        image=(abs(fxy)/max(abs(fxy(:))));
        
end
end

