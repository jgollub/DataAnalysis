
function p=testcoarray(array1,array2,x,y,n)

[r,c]=size(array1);
[rr,cc]=meshgrid(-r/2:r/2-1,-c/2:c/2-1);
expn=exp((i*2*pi)*(rr*(y/r)+cc*(x/c)));
ftotal=expn*0;
for it=1:n
    spary1=array1.*randn(size(array1));
    spary2=array2.*randn(size(array2));
    ft=ifft2(fft2(spary1).*fft2(spary2));
    val=sum(sum(ft.*expn));
    ftotal=ftotal+conj(ft).*val;
end

figure(1);
ftotali=fft2(ftotal);
imagesc(fftshift(abs(ftotali)));
colormap(gray);
colormap;
