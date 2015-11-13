function f = TVnorm_RI(x,nx,ny,nz)

x = reshape(x,nx,ny,nz);

f = 0;
parfor zi = 1:nz
    f = f+TVnorm(real(squeeze(x(:,:,zi))))+TVnorm(imag(squeeze(x(:,:,zi))));
end