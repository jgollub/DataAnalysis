

h = findobj(gca,'Type','line')
x=get(h,'Xdata');
y=get(h,'Ydata');
figure()
 xmm=x{1};

% xsar=x{2};
xsar=x;
 ymm=y{1};

% ysar=y{2};
ysar=y;
%%
xmm_itp=linspace(min(xmm),max(xmm),100*length(xmm));
xsar_itp=linspace(min(xsar),max(xsar),100*length(xsar));

ymm_itp=interp1(xmm,db2mag(ymm),xmm_itp,'spline');
ymm_itp=ymm_itp/max(ymm_itp);
ysar_itp=interp1(xsar,db2mag(ysar),xsar_itp,'spline');
ysar_itp=ysar_itp/max(ysar_itp);

% ymm_itp=interp1(xmm,ymm,xmm_itp,'spline');
% ymm_itp=ymm_itp./max(ymm_itp);
% ysar_itp=interp1(xsar,ysar,xsar_itp,'spline');
% ysar_itp=ysar_itp./max(ysar_itp);

plot(xsar_itp,db(ysar_itp),'k','lineWidth',2)
hold on
plot(xmm_itp,db(ymm_itp),'r', 'lineWidth',2)
xlim([-.03,0.03]);
ylim([-40,0]);

figure()
plot(xsar_itp,db(ysar_itp),'k','lineWidth',2)
hold on
plot(xmm_itp,db(ymm_itp),'r', 'lineWidth',2)

xlim([-.8,0.8])
ylim([-80,0]);