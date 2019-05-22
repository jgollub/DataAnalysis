function [y_s2,z_s2,im] = gridrefine(y_s,z_s,im_orig,interp_factor)

[Ys,Zs] = meshgrid(y_s,z_s);
y_s2 = linspace(y_s(1),y_s(end),interp_factor*numel(y_s));
z_s2 = linspace(z_s(1),z_s(end),interp_factor*numel(z_s));
[Ys1,Zs1] = meshgrid(y_s2,z_s2);

im = interp2(Ys,Zs,im_orig,Ys1,Zs1,'spline');

end