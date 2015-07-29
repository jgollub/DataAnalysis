function f = mycalltoTVnewMatrix_RI(x,th,piter,nx,ny,nz)

x = reshape(x,nx,ny,nz);

x_r = real(x);
x_i = imag(x);

f_r = mycalltoTVnewMatrix(x_r,th,piter);
f_i = mycalltoTVnewMatrix(x_i,th,piter);

f = f_r + 1i*f_i;
f = f(:);