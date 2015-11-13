%
% Rotate Field
%
% fld: N X M X 3    vector field to be rotated
% rotmat: 3 X 3     matrix from old to new coordinates (orthonormal)
% sampx:            sampling x spatial frequency (samples/wavelength)
% sampy:            sampling y spatial frequency (samples/wavelength)
% newfld:           rotated field
function newfld=rotfield(fld,rotmat,sampx,sampy)

[n,m,l]=size(fld);
fld=circshift(fld,[floor(n/2) floor(m/2) 0]);
fld=fft(fft(fld,[],1),[],2);
fld=circshift(fld,[floor(n/2) floor(m/2) 0]);
[fx,fy]=ndgrid((-n/2:n/2-1)*(sampx/n),(-m/2:m/2-1)*(sampy/m));
kz2=(1-fx.*fx-fy.*fy);
zr=kz2>0;
fxyz=zeros(n,m,3);
fxyz(:,:,1)=fx;
fxyz(:,:,2)=fy;
fxyz(:,:,3)=sqrt(kz2.*zr);
fnewxyz=multfld(fxyz,rotmat);
nfx=fnewxyz(:,:,1)*(n/sampx)+(n/2+1);
nfy=fnewxyz(:,:,2)*(m/sampy)+(m/2+1);
newfld=interpupsample(fld,2,nfx,nfy);
newfld=newfld.*repmat(zr,[1 1 l]);
newfld=multfld(newfld,rotmat);
newfld=circshift(newfld,[-floor(n/2) -floor(m/2) 0]);
newfld=ifft(ifft(newfld,[],1),[],2);
newfld=circshift(newfld,[-floor(n/2) -floor(m/2) 0]);
return;

% Upsample a function using FFT
% funct = function to upsample
% p = number to times to upsample
% nf = upsampled function
function nf = upsamp(funct,p)
[n,m,l]=size(funct);
funct=circshift(fft(fft(funct,[],1),[],2),[floor(n/2) floor(m/2)]);
nf = zeros(n*p,m*p,l);
nf(floor((p-1)*n/2)+1:floor((p+1)*n/2),floor((p-1)*m/2)+1:floor((p+1)*m/2),:)=funct;
nf=(p*p)*ifft(ifft(circshift(nf,[-floor(n*p/2) -floor(m*p/2)]),[],1),[],2);
return;

% interp indices upsample
% func = function to interpolate (N X M X L)
% p = upsampling rate to use
% xcoors = x coordinates to upsample
% ycoors = y coordinates to upsample
function ifunc=interpupsample(func,p,xcoors,ycoors)

upfunc=upsamp(func,p);    % upsample function
[n,m,l]=size(upfunc);
xcoors=min(max((xcoors-1)*p,0),n-1);
ycoors=min(max((ycoors-1)*p,0),m-1);
fracxcoors = xcoors-floor(xcoors);
fracycoors = ycoors-floor(ycoors);
[np,mp]=size(xcoors);
ifunc=zeros(np,mp,l);
spac=n*m;
for d=1:l
    ifunc(:,:,d)=(1-fracxcoors).*(1-fracycoors).*upfunc(floor(xcoors)+floor(ycoors)*n+(spac*(d-1)+1)) + ...
                 (1-fracxcoors).*fracycoors.*upfunc(floor(xcoors)+ceil(ycoors)*n+(spac*(d-1)+1)) + ...
                 fracxcoors.*(1-fracycoors).*upfunc(ceil(xcoors)+floor(ycoors)*n+(spac*(d-1)+1)) + ...
                 fracxcoors.*fracycoors.*upfunc(ceil(xcoors)+ceil(ycoors)*n+(spac*(d-1)+1));
end
return;

% Multiply N X M X L vector field by a single matrix
% fld = field of vectors (N X M X L)
% mat = L X L matrix to multiply each vector by
% newfld = new N X M X L matrix
function newfld=multfld(fld,mat)
[n,m,l]=size(fld);
newfld=zeros(n,m,l);
for d=1:l
    newfld(:,:,d)=sum(repmat(reshape(mat(d,:),[1 1 l]),[n m 1]).*fld,3);
end
return;
