% generate coarray by ad-hoc method
% by DLM

% size1: size of array 1 [rows cols]
% size2: size of array 2 [rows cols]
% dup: array 2 duplicate of array 1 (1=yes 0=no)
% flip: (0=convolution, 1=correlation)
function [best1,best2,bestcvl]=gencoarray(size1,size2,dup,flip)

itr=1000;

numelements = @(x) (size(x,1)+size(x,2));

bestscr = 1e38;
if dup>0
    size2=size1;
end
n=4;
ep=1;
for it=0:itr-1
    zrary=zeros(size1(1)+size2(1),size1(2)+size2(2));
    array1=zrary;
    array2=zrary;
    array1(1:size1(1),1:size1(2))=topelements(rand(size1),numelements(zrary));
    if dup>0
        if flip>0
            array2(1:size1(1),1:size1(2))=fliplr(flipud(array1(1:size1(1),1:size1(2))));
        else
            array2(1:size1(1),1:size1(2))=array1(1:size1(1),1:size1(2));
        end
    else
        array2(1:size2(1),1:size2(2))=topelements(rand(size2),numelements(zrary));
    end
    inneritr=5*(sum(size1)+sum(size2));
    for innert=0:inneritr-1
        cvfft1=fft2(array1);
        cvfft2=fft2(array2);
        cvl=real(ifft2(cvfft1.*cvfft2));
        cvl=cvl.^n;
        cvl=(2^n)*cvl./(1+(2^n)*cvl);
        %figure(1);clf;subplot(1,3,1);imagesc(array1);subplot(1,3,2);imagesc(array2);subplot(1,3,3);imagesc(cvl);drawnow;pause(0.5);
        cvfftl=fft2(cvl);
        if dup>0
            arn = real(ifft2(cvfftl.*cvfft2./(abs(cvfft2).^2+ep)));
            arn = arn(1:size1(1),1:size1(2));
            array1(1:size1(1),1:size1(2)) = topelements(arn,numelements(arn));
            if flip>0
                array2(1:size1(1),1:size1(2))=fliplr(flipud(array1(1:size1(1),1:size1(2))));
            else
                array2(1:size1(1),1:size1(2))=array1(1:size1(1),1:size1(2));
            end
        else            
            if mod(innert,2)==0
                arn = real(ifft2(cvfftl.*cvfft2./(abs(cvfft2).^2+ep)));
                arn = arn(1:size1(1),1:size1(2));
                array1(1:size1(1),1:size1(2)) = topelements(arn,numelements(arn));
            else
                arn = real(ifft2(cvfftl.*cvfft1./(abs(cvfft1).^2+ep)));
                arn = arn(1:size2(1),1:size2(2));
                array2(1:size2(1),1:size2(2)) = topelements(arn,numelements(arn));
            end
        end
    end
    cvfft1=fft2(array1);
    cvfft2=fft2(array2);
    cvfftl=cvfft1.*cvfft2;
    cvl=real(ifft2(cvfftl));
    %score=sum(sum(abs(cvl-1)));
    score=sum(sum(abs(cvl-circshift(cvfftl,[0 1])).^2))+sum(sum(abs(cvl-circshift(cvfftl,[1 0])).^2));
    if score<bestscr
        bestscr=score;
        best1=array1;
        best2=array2;
        bestcvl=cvl;
        figure(2);clf;
        subplot(1,4,1);
        imagesc(best1);colormap(gray);
        subplot(1,4,2);
        imagesc(best2);colormap(gray);
        subplot(1,4,3);
        imagesc(bestcvl);colormap(gray);colorbar;
        subplot(1,4,4);
        imagesc(fftshift(abs(cvfftl)));colormap(gray);colorbar;
        drawnow;pause(0.02);
    end
end
return;

function tp=topelements(x,n)

[r,c]=size(x);
y=reshape(x,[1 r*c]);
y=sort(y);
tp=x > y(1,r*c-n+1);
return;
