
function [E_values,yyy2, zzz2]=reSquare(values, yy, zz)
yyy=unique(yy);
s_y=length(yyy);
zzz=unique(zz); 
s_z=length(zzz);
[yyy2, zzz2]=meshgrid(yyy,zzz);
sqr_mtx_vec=cat(2,yyy2(:), zzz2(:));
org_vec=cat(2,yy,zz);
[m_log, m_idx]=ismember(org_vec,sqr_mtx_vec,'rows');
m_idx=m_idx(m_idx~=0);
[iz, iy]=ind2sub([s_z,s_y],m_idx);
E_values=zeros(size(yyy2));
for ii=1:length(values(:))
E_values(iz(ii),iy(ii))=values(ii);
end