function plot_saved_reconstruction_slices(scene_data,dataset_num,upsample,range_indexes,grid)

if nargin>=4
    ranges = range_indexes;
else
    ranges = 1:size(scene_data.obj_saved(dataset_num).reconstructed,3);
end

if nargin==2
    Z = scene_data.Z;
    Az = scene_data.Az;
    El = scene_data.El;
    obj3D = abs(scene_data.obj_saved(dataset_num).reconstructed);
    plot_range_slices(obj3D,Az,El,Z)
else
    Z = scene_data.Z;
    Az = scene_data.Az;
    El = scene_data.El;
    obj3D = abs(scene_data.obj_saved(dataset_num).reconstructed);
    plot_range_slices(obj3D,Az,El,Z,upsample,ranges,grid)
end

