% distance=-0.06:-.0001:-0.07;
distance=-0.005:+.0001:+.0005;


 figure(101); clf; clear image FM;
for loop=1:numel(distance)
image{loop}=bp(measurements,X,Y,distance(loop));
FM{loop} = fmeasure(image{loop}, 'LAPE',[]);

 figure(101);  hold on;
%         subplot(1,numel(distance),loop); 
        hold on; cla;
        scatter3(xx(:),yy(:),zz(:),6,20*log10(abs(image{loop}(:))),'filled')
        axis image; colormap('hot');set(gcf,'color','w');
       title([num2str(distance(loop))]);
        xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
        view(90,0); hold on;
        drawnow;
end
figure(102); clf;
merit=[FM{:}];
plot(1:numel(distance),merit)
[~, min_merit]=min(merit);
[~, max_merit]=max(merit);
display(['minimum: ',num2str(distance(min_merit))])
display(['maximum: ',num2str(distance(max_merit))])