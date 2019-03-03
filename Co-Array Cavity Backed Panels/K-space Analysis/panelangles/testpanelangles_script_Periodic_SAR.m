sx_tot=[];
sy_tot=[];
sz_tot=[];
dirx_tot=[];
diry_tot=[];
dirz_tot=[];


panellist = [];
for n1=-4:4
    for n2=-4:4
        panellist = [];
%        if (rand(1,1) < 0.2)
            panellist = [ panellist; n1*2 n2*2 n1*2+1 n2*2+1 ];
%          end
%     end
% end

probelist = [];
% for n1=-4:4
%     for n2=-4:4
%          if (rand(1,1) < 0.2)
            probelist = [ probelist; n1*2 n2*2 n1*2+2 n2*2+2 ];
%          end


point = [0 0 5];

% figure(1);clf;
% subplot(1,2,1);
% panelpat = plotpanels(panellist, 0.15);
% title('Panel pattern');
% subplot(1,2,2);
% probepat = plotpanels(probelist, 0.15);
% title('Probe pattern');

figure(2);clf
[directions,sx, sy, sz]= panelangles(panellist, probelist, point, 0.05, 0.15, (30:1:30)/30);

dirx = squeeze(sum(directions,1));
diry = squeeze(sum(directions,2));
dirz = squeeze(sum(directions,3));

    dirx_max=max(dirx_max,dirx);
     diry_max=max(diry_max,diry);
      dirz_max=max(dirz_max,dirz);

    if size(sx_tot)~=size(sx)
    sx_tot=zeros(size(sx));
    sy_tot=zeros(size(sy));
    sz_tot=zeros(size(sz));
end
dirx_tot = dirx+dirx_tot;
diry_tot = diry+diry_tot;
dirz_tot = dirz+dirz_tot;

sx_tot=sx_tot+sx;
sy_tot=sy_tot+sy;
sz_tot=sz_tot+sz;
% figure(3);clf
% [directions,sx, sy, sz] = panelangles(panellist, probelist, point, 0.05, 0.15, (22:1:30)/30);
    end
end

subplot(1,3,1);
hold on;
imagesc(sy_tot,sz_tot,abs(dirx_tot).');colormap(1-gray);colorbar;
title(sprintf('X cross section\nFourier space'));
xlabel('Y axis');
ylabel('Z axis');
drawnow;
subplot(1,3,2);
hold on;
imagesc(sx_tot,sz_tot,abs(diry_tot).');colormap(1-gray);colorbar;
title(sprintf('Y cross section\nFourier space'));
xlabel('X axis');
ylabel('Z axis');
drawnow;
subplot(1,3,3);
hold on;
imagesc(sx_tot,sy_tot,abs(dirz_tot).');colormap(1-gray);colorbar;
title(sprintf('Z cross section\nFourier space'));
xlabel('X axis');
ylabel('Y axis');
drawnow;
