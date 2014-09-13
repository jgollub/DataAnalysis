% function angles = panelangles(panellist, probelist, pointlist, samp, sampconv)
%
% calculate accessible angles from rectangular regions for
% a list of points
%
% assumes isotropic radiators
% brute force but simple
%
% panellist = [ x1 y1 x2 y2 ] rectangles (on z=0 plane)
% probelist = [ x1 y1 x2 y2 ] rectangles (on z=0 plane) 
% point = [ x y z ] point to evaluate angles to
% samp = sampling rate (in tangent units)
% sampconv = sampling rate for convolution (spatial units)
% freqscal = frequency scaling

function [directions,sx,sy,sz] = panelangles(panellist, probelist, point, samp, sampconv, freqscal)

[r1,c1]=size(panellist);
[r2,c2]=size(probelist);

mincoor = min([panellist; probelist],[],1);
maxcoor = max([panellist; probelist],[],1);

minpt = [ mincoor(1)-point(1) mincoor(2)-point(2) -point(3) ];
minpt = minpt/sqrt(sum(abs(minpt).^2));
maxpt = [ maxcoor(3)-point(1) maxcoor(4)-point(2) point(3) ] ;
maxpt = maxpt/sqrt(sum(abs(maxpt).^2));

ctr = 0.5*(minpt+maxpt);
spandif = maxpt-minpt;

minpt = ctr-2*spandif;
maxpt = ctr+2*spandif;

sampx = ceil((maxpt(1)-minpt(1))/samp);
sampy = ceil((maxpt(2)-minpt(2))/samp);
sampz = ceil((maxpt(3)-minpt(3))/samp);

directions = zeros(sampx,sampy,sampz);

vecdifmin = minpt;
vecdifspan = maxpt-minpt;
minspanval = [ 1 1 1 ];
maxspanval = [ sampx sampy sampz ];

for curfreq=freqscal
    for n1=1:r1
        for n2=1:r2
            for x1=panellist(n1,1):sampconv:panellist(n1,3)
                for y1=panellist(n1,2):sampconv:panellist(n1,4)
                    for x2=probelist(n2,1):sampconv:probelist(n2,3)
                        for y2=probelist(n2,2):sampconv:probelist(n2,4)
                            vec1 = [ x1-point(1) y1-point(2) -point(3) ];
                            vec2 = [ point(1)-x2 point(2)-y2 point(3) ];
                            vec1 = vec1 * (curfreq / sqrt(sum(abs(vec1).^2)));
                            vec2 = vec2 * (curfreq / sqrt(sum(abs(vec2).^2)));
                            vecdif = vec1 - vec2;
                            vecdif = floor(1.5 + maxspanval.*((vecdif - vecdifmin) ./ vecdifspan));
                            vecdif = min(max(vecdif,minspanval),maxspanval);
                            directions(vecdif(1),vecdif(2),vecdif(3)) = directions(vecdif(1),vecdif(2),vecdif(3)) + 1;
                        end
                    end
                end
            end
        end
    end
end

dirx = squeeze(sum(directions,1));
diry = squeeze(sum(directions,2));
dirz = squeeze(sum(directions,3));

sx=minpt(1):samp:maxpt(1)-samp;
sy=minpt(2):samp:maxpt(2)-samp;
sz=minpt(3):samp:maxpt(3)-samp;

subplot(1,3,1);
imagesc(sy,sz,abs(dirx).');colormap(1-gray);colorbar;
title(sprintf('X cross section\nFourier space'));
xlabel('Y axis');
ylabel('Z axis');
subplot(1,3,2);
imagesc(sx,sz,abs(diry).');colormap(1-gray);colorbar;
title(sprintf('Y cross section\nFourier space'));
xlabel('X axis');
ylabel('Z axis');
subplot(1,3,3);
imagesc(sx,sy,abs(dirz).');colormap(1-gray);colorbar;
title(sprintf('Z cross section\nFourier space'));
xlabel('X axis');
ylabel('Y axis');

return;
end
