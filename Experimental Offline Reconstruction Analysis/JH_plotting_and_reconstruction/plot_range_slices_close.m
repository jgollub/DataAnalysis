function plot_range_slices(obj3D,Az,El,Z,upsample,range_indexes,grid)
%plot a 3D scene in multipule range slices
%Inputs:
%obj3D - a 3D matrix where the first dimension is elevation, the second is azimuth, and the third is range
%Az,El - azimuth and elevation, 2D plaid matricies of the type produced by meshdrid()
%Z - 1D array of ranges
%upsample - image upsampling factor
%range_indexes - 1D array of range indicies to plot

obj3D = abs(obj3D);

%plot a 3D image as 2D range slices
if nargin>=6
    ranges = range_indexes;
else
    ranges = 1:size(obj3D,3);
end

Nz = length(ranges);
cmax = maxall(obj3D(:,:,ranges));

if nargin>=5
    for zn=1:length(ranges)
        if length(ranges)~=1
            subplot(ceil(Nz/2),2,zn);
        end
        
        
        %% cart coords
        if nargin==7
            z = 115;
            l = 1;
            s = 1;
            imagesc(tan(Az(1,:))*z,tan(El(:,1))*z,upsample_image(obj3D(:,:,ranges(zn)),upsample))
            xlabel('cm')
            ylabel('cm')
            if strcmp(grid,'h')
                % horizontal grid!!!!!!!!!!!!
                line([0 0],[-l l],'Color','r')
                xn = 0;
                small = 1;
                while xn<abs(tan(Az(1,1))*z+3)
                    xn = xn+s;
                    Azn = atan(xn/z);
                    if small
                        line([xn xn],[-l l]./2,'Color','r')
                        line(-[xn xn],[-l l]./2,'Color','r')
                        small = 0;
                    else
                        line([xn xn],[-l l],'Color','r')
                        line(-[xn xn],[-l l],'Color','r')
                        small = 1;
                    end
                end
                
            elseif strcmp(grid,'v')
                % vertical grid!!!!!!!!!!!!
                line([-l l],[0 0],'Color','r')
                yn = 0;
                small = 1;
                while yn<abs(tan(El(1,1))*z+3)
                    yn = yn+s;
                    if small
                        line([-l l]./2,[yn yn],'Color','r')
                        line([-l l]./2,-[yn yn],'Color','r')
                        small = 0;
                    else
                        line([-l l],[yn yn],'Color','r')
                        line([-l l],-[yn yn],'Color','r')
                        small = 1;
                    end
                end
            end
            
        else
            % angular coords
            imagesc(Az(1,:)*180/pi,El(:,1)*180/pi,upsample_image(obj3D(:,:,ranges(zn)),upsample)) 
            xlabel('Az')
            ylabel('El')
        end
        %%
        title(['z = ', num2str(Z(zn))])
        axis equal;axis tight;
        caxis([0 cmax]);
    end
    
elseif nargin==4
    for zn=ranges
        subplot(ceil(Nz/2),2,zn)
        imagesc(Az(1,:)*180/pi,El(:,1)*180/pi,obj3D(:,:,zn))
        title(['z = ', num2str(Z(zn))])
        axis equal;axis tight;
        caxis([0 cmax]);
    end
    
elseif nargin==2
    for zn=ranges
        subplot(ceil(Nz/2),2,zn)
        imagesc(upsample_image(obj3D(:,:,zn),upsample))
        title(['z = ', num2str(Z(zn))])
        axis equal;axis tight;
        caxis([0 cmax]);
    end
    
elseif nargin==1
    for zn=ranges
        subplot(ceil(Nz/2),2,zn)
        imagesc(obj3D(:,:,zn))
        title(['z = ', num2str(Z(zn))])
        axis equal;axis tight;
        caxis([0 cmax]);
    end
end
