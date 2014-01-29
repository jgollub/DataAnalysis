% Set up the movie.
writerObj = VideoWriter('/Users/jgollub/Dropbox/MetaImager Project/MetaImager Data/MetaImager Scenes/test'); % Name it.
writerObj.FrameRate = 1; % How many frames per second.
open(writerObj); 
for frm=1:8
    [CurrentObj3D,Az,El,Z,F,reg]=VoxPlot_Realtime_image_reconstruction_revA3_1of2panel(frm,obj_saved);
    figure(frm+8);
    hold all;
    voxelimage_projections(CurrentObj3D,Az,El,Z)
    
   %make movie
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end
close(writerObj);