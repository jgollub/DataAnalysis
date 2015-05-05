%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


function [f_plot, maxVal_position] = filterScatters(f_raw, indices, imgDomain, domainGrid, R_l,threshold,R_fine)
        maxVal_position=NaN;
        f = zeros(size(imgDomain, 1), 1);
        f(indices) = abs(f_raw).^2;
        f_plot=reshape(f,domainGrid);
       
        x_centroid=sum(imgDomain(indices,1).*abs(f_raw(:)).^2)/sum(abs(f_raw).^2);
        y_centroid=sum(imgDomain(indices,2).*abs(f_raw(:)).^2)/sum(abs(f_raw).^2);
        z_centroid=sum(imgDomain(indices,3).*abs(f_raw(:)).^2)/sum(abs(f_raw).^2);
        
        maxVal_position=[x_centroid, y_centroid,z_centroid];
        
        f_db       = 20*log10(abs(f_plot));
        thresh     = max(max(max(f_db))) -threshold;
%         f          = ones(size(f_plot))*thresh;
%         f_plot = 20*log10(abs(f_plot));
        f_db  (f_db  <thresh)= thresh;
            
%         %upsample to determine max position
%         if ~isnan(R_fine)
%             X_subImgDomain=[min(imgDomain(indices,1)):R_fine:max(imgDomain(indices,1))];
%             Y_subImgDomain=[min(imgDomain(indices,2)):R_fine:max(imgDomain(indices,2))];
%             Z_subImgDomain=[min(imgDomain(indices,3)):R_fine:max(imgDomain(indices,3))];
%             
%             
%             %convert to log scale and correct indexing
%             f2_db       = 20*log10(abs(f_raw));
%             thresh     = max(max(max(f2_db))) -threshold;
%             %         f2          = ones(size(f_raw))*thresh;
%             %         db = 20*log10(abs(f_raw));
%             f2_db(f2_db<thresh)= thresh;
%             
%             %define fine grid for interpolation
%             [Xf, Yf, Zf]=meshgrid(X_subImgDomain,Y_subImgDomain,Z_subImgDomain);
%             XI=[Xf(:),Yf(:),Zf(:)];
%             XYZf=[imgDomain(indices,1),imgDomain(indices,2),imgDomain(indices,3)];
%             f_up=griddatan(XYZf, f2_db, XI);
%             
%             [f_max, f_I]=max(f_up);
%             maxVal_position=[Xf(f_I),Yf(f_I),Zf(f_I)];
%         else
%             [max_value,index_max]=max(f_plot(:));
%             maxVal_position=imgDomain(index_max,:);
%         end
end