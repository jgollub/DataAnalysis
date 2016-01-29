%find maximum sub point
%func = function to search
%p = upsampling
function [xp,yp]=findsubmax(func,p)

n=10;
pfunc=abs(upsamp(func,p));
[ry,rx]=size(pfunc);
% rx
% ry
[yvals,ypos]=max(pfunc);
[xvals,xpos]=max(yvals);
xc=xpos-1;
yc=ypos(xpos)-1;
[xv,yv]=meshgrid(min(max(xc-n:xc+n,1),rx),min(max(yc-n:yc+n,1),ry));
ind=yv*ry+xv+1;
[xx,yy]=meshgrid((-n:n),(-n:n));
rr=sqrt(xx.^2+yy.^2);
rr=(rr<=n).*(pi*rr/(2*(n+1)));
windw=cos(rr).*pfunc(ind);
total=sum(sum(windw));
xp=(sum(sum(windw.*xv))/(total*p)+1);
yp=(sum(sum(windw.*yv))/(total*p)+1);
return;


% Upsample a function using FFT
% funct = function to upsample
% p = number to times to upsample
% nf = upsampled function
function nf = upsamp(funct,p)
[n,m,l]=size(funct);
funct=circshift(fft(fft(funct,[],1),[],2),[floor(n/2) floor(m/2)]);
nf = zeros(n*p,m*p,l);
nf(floor((p-1)*n/2)+1:floor((p+1)*n/2),floor((p-1)*m/2)+1:floor((p+1)*m/2),:)=funct;
nf=(p*p)*ifft(ifft(circshift(nf,[-floor(n*p/2) -floor(m*p/2)]),[],1),[],2);
return;