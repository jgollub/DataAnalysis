function [peakPos peak_Pixels]=LocateReflPeak(x_range, y_range,Recon_obj,recon_fig)


%set function parameters and find max

thres = (max([min(max(Recon_obj,[],1))  min(max(Recon_obj,[],2))]));
max_val=max(max(Recon_obj));

filt = (fspecial('gaussian',2,1));
%filt = (fspecial('average', 1));
edge =1;
res=2;
%use centroid fitting to find max position (between pixels)
%%find peak position

%scale function for plotting

[ysize xsize]=size(Recon_obj);

%y-coordinate
yaxis_slope=(y_range(end)-y_range(1))/(ysize-1);
y_c=y_range(1)-yaxis_slope;
y_pos=@(y) yaxis_slope*y+y_c;

%x-coordinate
xaxis_slope=(x_range(end)-x_range(1))/(xsize-1);
x_c=x_range(1)-xaxis_slope;
x_pos=@(x) xaxis_slope*x+x_c;

%find peak position in pixels
figure(36)
%peak_Pixels=FastPeakFind(Recon_obj, thres, filt, edge, res).';
[~, peak_Pixels(1)]=max(max(Recon_obj,[],1))
[~, peak_Pixels(2)]=max(max(Recon_obj,[],2))

figure(recon_fig)
for i=1:length(peak_Pixels)/2
    
    peakPos(i,:)=[x_pos(peak_Pixels(i,1)) y_pos(peak_Pixels(i,2))]; %x & y pos
    if nargin>3
        figure(recon_fig);
        hold on;
        plot(peakPos(i,1),peakPos(i,2),'Marker','+','Color',[.88 .48 0],'MarkerSize',30)
        hold off;
    end
end