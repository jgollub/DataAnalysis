function draw_Kinect_object_scene(xyz,rgb,objs,object_nums)

rgb = double(rgb)./255;

for on=object_nums
        notobj = repmat(~objs(:,:,on),[1,1,3]);
        xyz_obj = xyz;
        rgb_obj = rgb;
        xyz_obj(notobj) = NaN;
        rgb_obj(notobj) = NaN;
        surface(-xyz_obj(:,:,1),xyz_obj(:,:,3),xyz_obj(:,:,2),rgb_obj,'edgecolor','none');

end

axis equal;
axis tight;
axis([-2.5 2.5, 0 5, -1.5 2, 0 5])
%view(-45,+40)