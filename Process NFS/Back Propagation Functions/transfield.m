
% transvec: 3 X 3   translation vector (in wavelengths)
% sampx:            sampling x spatial frequency (samples/wavelength)
% sampy:            sampling y spatial frequency (samples/wavelength)
% newfld:           translated field
function newfld=transfield(fld,transvec,sampx,sampy)

[n,m,l]=size(fld);
fld=circshift(fld,[floor(n/2) floor(m/2) 0]);
fld=fft(fft(fld,[],1),[],2);
fld=circshift(fld,[floor(n/2) floor(m/2) 0]);
[fx,fy]=ndgrid((-n/2:n/2-1)*(sampx/n),(-m/2:m/2-1)*(sampy/m));
kz2=(1-fx.*fx-fy.*fy);
zr=kz2>0;
fz=sqrt(kz2.*zr);
phs=zr.*exp((-2*pi*1j)*(transvec(1)*fx+transvec(2)*fy+transvec(3)*fz));
newfld=multfld(fld,phs);
newfld=circshift(newfld,[-floor(n/2) -floor(m/2) 0]);
newfld=ifft(ifft(newfld,[],1),[],2);
newfld=circshift(newfld,[-floor(n/2) -floor(m/2) 0]);
return;

% Multiply N X M X L vector field by a scalar field
% fld = field of vectors (N X M X L)
% scalar = scalar field
function newfld=multfld(fld,scalar)
[n,m,l]=size(fld);
newfld=zeros(n,m,l);
for d=1:l
    newfld(:,:,d)=scalar.*fld(:,:,d);
end
return;
