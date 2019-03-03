
panellist = [];
for n1=-1:2:1
    for n2=-1:2:1
        if (rand(1,1) < 1)
            panellist = [ panellist; n1*2 n2*2 n1*2+1 n2*2+1 ];
        end
    end
end

probelist = [];
for n1=-1:1
    for n2=-1:1
        if (rand(1,1) < 1)
            probelist = [ probelist; n1 n2 n1+.2 n2+.2 ];
        end
    end
end

point = [0 0 5];

figure(1);clf;
subplot(1,2,1);
panelpat = plotpanels(panellist, 0.15);
title('Panel pattern');
subplot(1,2,2);
probepat = plotpanels(probelist, 0.15);
title('Probe pattern');

figure(2);clf
directions = panelangles(panellist, probelist, point, 0.05, 0.05, (17.5:0.0891:26.5)/30);

figure(3); clf
directions = panelangles_freq(panellist, probelist, point, 0.05, 0.001, (17.5:0.0891:26.5)/30);
