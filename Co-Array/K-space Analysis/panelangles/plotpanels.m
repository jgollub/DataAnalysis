% function plotpanels = plotpanels(panellist, sampconv)
%
% calculate accessible angles from rectangular regions for
% a list of points
%
% panellist = [ x1 y1 x2 y2 ] rectangles (on z=0 plane)
% sampconv = sampling rate for convolution (spatial units)

function [panelpat,sx,sy] = panelangles(panellist, sampconv)

[r1,c1]=size(panellist);

mincoor = min(panellist,[],1);
maxcoor = max(panellist,[],1);

minpt = [ mincoor(1) mincoor(2) ]; 
maxpt = [ maxcoor(3) maxcoor(4) ];

ctr = 0.5*(minpt+maxpt);
spandif = maxpt-minpt;

minpt = ctr-0.75*spandif;
maxpt = ctr+0.75*spandif;

sampx = ceil((maxpt(1)-minpt(1))/sampconv);
sampy = ceil((maxpt(2)-minpt(2))/sampconv);

panelpat = zeros(sampx,sampy);

vecdifmin = minpt;
vecdifspan = maxpt-minpt;
minspanval = [ 1 1 ];
maxspanval = [ sampx sampy ];

for n1=1:r1
    for x1=panellist(n1,1):sampconv:panellist(n1,3)
        for y1=panellist(n1,2):sampconv:panellist(n1,4)
            vecdif = [ x1 y1 ];
            vecdif = floor(1.5 + maxspanval.*((vecdif - vecdifmin) ./ vecdifspan));
            vecdif = min(max(vecdif,minspanval),maxspanval);
            panelpat(vecdif(1),vecdif(2)) = panelpat(vecdif(1),vecdif(2)) + 1;
        end
    end
end

sx=minpt(1):sampconv:maxpt(1)-sampconv;
sy=minpt(2):sampconv:maxpt(2)-sampconv;

imagesc(sx,sy,panelpat.');colormap(1-gray);colorbar;
title('Pattern');
xlabel('X axis');
ylabel('Y axis');

return;
end



