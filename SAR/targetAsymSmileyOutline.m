function target = targetAsymSmileyOutline(xRes,yRes,xgrid,ygrid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

target=zeros(size(xgrid));

xstep=abs(mean(diff(xgrid(1,:))));
ystep=abs(mean(diff(ygrid(:,1))));

i = ones(round(yRes/ystep),round(xRes/xstep));
o = zeros(size(i));

[m, n] = size(i);


target(floor(end/2-3*n/2):floor(end/2-3*n/2)+(3*n-1),floor(end/2-4*m/2):floor(end/2-4*m/2)+(4*m-1)) = ...
                [ i o i o; ...
                  o o o o; ...
                  i i i i];
end

