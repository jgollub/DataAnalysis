[y,rad]=fitcyl(data.');

figure(305); clf; hold on;
scatter3(data(1,:)',data(2,:).',data(3,:).',30,'filled','k');
axis equal; xlabel('x'); ylabel('y'); zlabel('z'); title('Cyl. Pts in Panel Coord. Sys')
drawnow;

% translate 
[rs,cs]=size(data(1:3,:).');
dataprime = data(1:3,:).' + repmat(y(1:3),[rs 1]);

%rotate
[rs,cs]=size(dataprime);
axs = y(4:6)/norm(y(4:6));
axs = repmat(axs,[rs 1]);
ang = repmat(y(7),[rs 3]);
dataprime = dataprime.*cos(ang) + cross(dataprime,axs).*sin(ang) + axs.*(repmat(sum(axs.*dataprime,2),[1 3]).*(1-cos(ang)));
% 
figure(305);
scatter3(dataprime(:,1),dataprime(:,2),dataprime(:,3).',30,'filled','k');
axis equal; xlabel('x'); ylabel('y'); zlabel('z'); title('Cyl. Pts in Panel Coord. Sys')
drawnow;
CalPitch=.025;

dtheta=CalPitch/rad;


normCyl=[-1 0 0];
axs = y(4:6)/norm(y(4:6));
normCylPrime= normCyl.*cos(-y(7)) + cross(normCyl,axs).*sin(-y(7)) + axs.*(repmat(sum(axs.*normCyl,2),[1 3]).*(1-cos(-y(7))));
angleshift=-asin(norm(cross(normCylPrime,normCyl))/(norm(normCylPrime)*norm(normCyl)));

theta=angleshift+pi/2:dtheta:3*pi/2+angleshift;
xfit=rad*cos(theta);
yfit=rad*sin(theta);
zfit=min(dataprime(:,3)):CalPitch:max(dataprime(:,3));
zfitpad=repmat(zfit,numel(xfit),1);
xypad=repmat([xfit.',yfit.'],numel(zfit),1);
gridprime=[xypad,zfitpad(:)];

figure(305); hold on;
scatter3(gridprime(:,1),gridprime(:,2),gridprime(:,3),10,'g')
axis equal

%rotate
[rs,cs]=size(gridprime);
axs = y(4:6)/norm(y(4:6));
axs = repmat(axs,[rs 1]);
ang = repmat(-y(7),[rs 3]);
gridprime = gridprime.*cos(ang) + cross(gridprime,axs).*sin(ang) + axs.*(repmat(sum(axs.*gridprime,2),[1 3]).*(1-cos(ang)));
% 
figure(305);
scatter3(gridprime(:,1),gridprime(:,2),gridprime(:,3).',3,'filled','b');
axis equal; xlabel('x'); ylabel('y'); zlabel('z'); title('Cyl. Pts in Panel Coord. Sys')
drawnow;

% translate
[rs,cs]=size(gridprime);
grid = gridprime + repmat(-y(1:3),[rs 1]);

figure(305); hold on;
scatter3(grid(:,1),grid(:,2),grid(:,3).',30,'filled','b');
axis equal; xlabel('x'); ylabel('y'); zlabel('z'); title('Cyl. Pts in Panel Coord. Sys')
drawnow;