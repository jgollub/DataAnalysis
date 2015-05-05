
panellist = [];
for n1=-3:3
    for n2=-3:3
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

figure(1);clf;
subplot(1,2,1);
panelpat = plotpanels(panellist, 0.15);
title('Panel pattern');
subplot(1,2,2);
probepat = plotpanels(probelist, 0.15);
title('Probe pattern');

figure(2);clf
directions = panelangles(panellist, probelist, point, 0.05, 0.15, (30:1:30)/30);

figure(3);clf
directions = panelangles(panellist, probelist, point, 0.05, 0.15, (22:1:30)/30);
    end
end