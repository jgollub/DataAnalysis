function [y,rad] = fitcyl(cr)

cr=cr(:,1:3);
% cr=cr(cr(:,3)>800,:);
%x0 = [0 0 0 0 1 1 0];
x0 =randn(1,7);
optfunc = @(x) cylpenalty(x,cr);
opt = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxIter',10000);
x = lsqnonlin(optfunc,x0,[],[],opt);
transfit = [ norm(optfunc(x0)) norm(optfunc(x)) ]
ncr = dotrans(x,cr);
ctrx = sum(ncr(:,1))/length(ncr(:,1));
ctry = sum(ncr(:,2))/length(ncr(:,2));
rad = sqrt(sum((ncr(:,1)-ctrx).^2+(ncr(:,2)-ctry).^2)/length(ncr(:,1)));

xcir0 = 0.01*rand(1,3);
optfunc = @(x) circlepenalty(x,ncr);
opt = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Jacobian','on','MaxIter',10000);
xcir = lsqnonlin(optfunc,xcir0,[],[],opt);

cirfit = [ norm(optfunc(xcir0)) norm(optfunc(xcir)) ] 
rad=sqrt(xcir(1));

v=[ xcir(2) xcir(3) 0];
xrot = rotvec(v,x(4:6),-x(7));

x(1:3)=x(1:3)-xrot;

ncr2 = dotrans(x,cr);
figure(1);clf;
subplot(1,3,1);
scatter3(cr(:,1),cr(:,2),cr(:,3));
subplot(1,3,2);
scatter3(ncr(:,1),ncr(:,2),ncr(:,3));
subplot(1,3,3);
scatter3(ncr2(:,1),ncr2(:,2),ncr2(:,3));

y=x;
return;

function pts2 = transvec(pts,vec)

[rs,cs]=size(pts);
pts2 = pts + repmat(vec,[rs 1]);
return;

function pts2 = rotvec(pts,axs,ang)

[rs,cs]=size(pts);
axs = axs/norm(axs);
axs = repmat(axs,[rs 1]);
ang = repmat(ang,[rs 3]);
pts2 = pts.*cos(ang) + cross(pts,axs).*sin(ang) + axs.*(repmat(sum(axs.*pts,2),[1 3]).*(1-cos(ang)));
return;

function scr = dotrans(x,cr)
ang = x(7);
axs = x(4:6);
trans = x(1:3);
cr = transvec(cr,trans);
scr = rotvec(cr,axs,ang);
return;

function scr = cylpenalty(x,cr)
cr = dotrans(x,cr);
axs = x(4:6);
scr = [cr(:,1);cr(:,2);1e-3*cr(:,3)];
scr = [scr;100*(sum(abs(axs).^2)-1)];
return;

function [scr,j] = circlepenalty(xcir,cr)
cr(:,1)=cr(:,1)-xcir(2);
cr(:,2)=cr(:,2)-xcir(3);
scr = cr(:,1).^2+cr(:,2).^2-xcir(1);
j = [ repmat(-1,[length(cr(:,1)) 1]) -2*cr(:,1) -2*cr(:,2) ];
return;
