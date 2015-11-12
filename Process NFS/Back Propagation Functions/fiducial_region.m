function [ search_region ] = fiducial_region(fields,field_pos, fiducial_pos,fiducial_range)

fid_dist_lt =  fiducial_range>sqrt((field_pos.y-fiducial_pos.lt.y).^2+(field_pos.z-fiducial_pos.lt.z).^2);
fid_dist_rt =  fiducial_range>sqrt((field_pos.y-fiducial_pos.rt.y).^2+(field_pos.z-fiducial_pos.rt.z).^2);
fid_dist_lb =  fiducial_range>sqrt((field_pos.y-fiducial_pos.lb.y).^2+(field_pos.z-fiducial_pos.lb.z).^2);
fid_dist_rb =  fiducial_range>sqrt((field_pos.y-fiducial_pos.rb.y).^2+(field_pos.z-fiducial_pos.rb.z).^2);

search_region.lt.fields=fields(fid_dist_lt);
search_region.rt.fields=fields(fid_dist_rt);
search_region.lb.fields=fields(fid_dist_lb);
search_region.rb.fields=fields(fid_dist_rb);

search_region.lt.x=field_pos.x(fid_dist_lt);
search_region.lt.y=field_pos.y(fid_dist_lt);
search_region.lt.z=field_pos.z(fid_dist_lt);

search_region.rt.x=field_pos.x(fid_dist_rt);
search_region.rt.y=field_pos.y(fid_dist_rt);
search_region.rt.z=field_pos.z(fid_dist_rt);

search_region.lb.x=field_pos.x(fid_dist_lb);
search_region.lb.y=field_pos.y(fid_dist_lb);
search_region.lb.z=field_pos.z(fid_dist_lb);

search_region.rb.x=field_pos.x(fid_dist_rb);
search_region.rb.y=field_pos.y(fid_dist_rb);
search_region.rb.z=field_pos.z(fid_dist_rb);

