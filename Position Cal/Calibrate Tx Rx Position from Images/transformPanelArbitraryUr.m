%Transform using arbitrary rotation matrix and translation vector 

function [ p1 ] = transformPanelArbitraryUr(p1,U,r)

if iscell(p1)
    for in=1:length(p1)
            p1{in} = transformPanelArbitraryUr(p1{in}, U, r);
    end
    return;
end

Tu=cat(2,[U(:,:); 0 0 0],[0 0 0 1].');
Tr=cat(2, [eye(3); [0 0 0]],[r(1) r(2) r(3) 1].');
T=Tr*Tu;
%!!!!!!
pos_size = size(p1.x);
l = length(p1.y(:));

%transform alpha beta gamma x,y,z 
pos_mat = T*[p1.x(:) p1.y(:) p1.z(:) ones(l,1)].';

p1.x = reshape(pos_mat(1,:), pos_size);
p1.y = reshape(pos_mat(2,:), pos_size);
p1.z = reshape(pos_mat(3,:), pos_size);

%transform the dipoles
d1 = p1.dipoles;
d_size = size(d1.x);
l = length(d1.y(:));

m_mat = Tu*[d1.x(:) d1.y(:) d1.z(:) ones(l,1)].';

d1.x = reshape(m_mat(1,:), d_size);
d1.y = reshape(m_mat(2,:), d_size);
d1.z = reshape(m_mat(3,:), d_size);
p1.dipoles = d1;

%transform the feed locations
f_mat = T*[p1.feedLocs ones(size(p1.feedLocs(:,1)))].';
p1.feedLocs = f_mat(1:3,:).';

% try to transform the individual antenna orientations, if they exist
try
    p1.u = Tu(1:3,1:3)*p1.u;
    p1.v = Tu(1:3,1:3)*p1.v;
catch
end

