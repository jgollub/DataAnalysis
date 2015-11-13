% panel_transform transforms and translates panel from initial position at
% the origin using a pre-determined tranformation, T. T is a matrix
% homogenous coordinates
function p1 = panel_transform(p1, T)

if iscell(p1)
    for in=1:length(p1)        
        p1{in} = panel_transform(p1{in}, alpha, beta, gamma);
    end
    return;
end

pos_size = size(p1.x);
l = length(p1.y(:));

%transform the x,y,z
pos_mat = T*[p1.x(:) p1.y(:) p1.z(:) ones(l,1)].';

p1.x = reshape(pos_mat(1,:), pos_size);
p1.y = reshape(pos_mat(2,:), pos_size);
p1.z = reshape(pos_mat(3,:), pos_size);

%transform the dipoles
T_r=[T(1:4,1:3), [0 0 0 1].'];

d1 = p1.dipoles;
d_size = size(d1.x);
l = length(d1.y(:));

m_mat = T_r*[d1.x(:) d1.y(:) d1.z(:) ones(l,1)].';

d1.x = reshape(m_mat(1,:), d_size);
d1.y = reshape(m_mat(2,:), d_size);
d1.z = reshape(m_mat(3,:), d_size);
p1.dipoles = d1;

%transform the feed locations
f_mat = T*[p1.feedLocs ones(size(p1.feedLocs(:,1)))].';
p1.feedLocs = f_mat(1:3,:).';