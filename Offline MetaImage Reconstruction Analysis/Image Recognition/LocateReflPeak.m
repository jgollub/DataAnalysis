function peakPos=LocateReflPeak(x_range, y_range,recon_array,recon_fig)


%set function parameters and find max

thres = (max([min(max(recon_array,[],1))  min(max(recon_array,[],2))]));
 %filt = (fspecial('gaussian', 3,1));
filt = (fspecial('average', 3));
edg =3;
res=2;
%use centroid fitting to find max position (between pixels)
%%find peak position

%scale function for plotting

[ysize xsize]=size(recon_array);

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
peakPos=FastPeakFind(recon_array, thres, filt,edg, res);
figure(recon_fig)
for i=1:length(peakPos)/2
    if nargin>3
        figure(recon_fig);
        hold on;
        plot(x_pos(peakPos(i)),y_pos(peakPos(i+1)),'Marker','+','Color',[.88 .48 0],'MarkerSize',30)
        hold off;
    end
end