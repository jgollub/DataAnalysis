% generate coarray by ad-hoc method
% by DLM

% size1: size of array 1 [rows cols]
% size2: size of array 2 [rows cols]
% dup: array 2 duplicate of array 1 (1=yes 0=no)
% flip: (0=convolution, 1=correlation)
function [best1,best2,bestcvl]=gencoarray(size1,size2,dup,flip)

itr=10000;

%numelements = @(x,f) (size(x,1)+size(x,2)-1);
numelements = @(x,f) floor((2-f)*(size(x,1)+size(x,2)));

bestscr = 1e38;
if dup>0
    size2=size1;
end
n=2;
ep=1;
sizet=size1+size2;
[rr,cc]=ndgrid(-sizet(1)/2:sizet(1)/2-1,-sizet(2)/2:sizet(2)/2-1);
filt = cos((pi/sqrt(sizet(1)^2+sizet(2)^2))*(rr.^2+cc.^2));
%filt=1;
zrary=zeros(sizet(1),sizet(2));
for it=0:itr-1
    array1=zrary;
    array2=zrary;
    array1(1:size1(1),1:size1(2))=topelements(rand(size1),numelements(size1,0));
    if dup>0
        if flip>0
            array2(1:size1(1),1:size1(2))=fliplr(flipud(array1(1:size1(1),1:size1(2))));
        else
            array2(1:size1(1),1:size1(2))=array1(1:size1(1),1:size1(2));
        end
    else
        array2(1:size2(1),1:size2(2))=topelements(rand(size2),numelements(size2,0));
    end
    inneritr=3*sum(sizet);
    for innert=0:inneritr-1
        cvfft1=fft2(array1);
        cvfft2=fft2(array2);
        cvl=real(ifft2(cvfft1.*cvfft2));
        cvl=cvl.^n;
        cvl=(2^n)*cvl./(1+(2^n)*cvl);
        cvl=cvl.*filt;
        cvfftl=fft2(cvl);
        if dup>0
            arn = real(ifft2(cvfftl.*cvfft2./(abs(cvfft2).^2+ep)));
            arn = arn(1:size1(1),1:size1(2));
            array1(1:size1(1),1:size1(2)) = topelements(arn,numelements(arn,innert/inneritr));
            if flip>0
                array2(1:size1(1),1:size1(2))=fliplr(flipud(array1(1:size1(1),1:size1(2))));
            else
                array2(1:size1(1),1:size1(2))=array1(1:size1(1),1:size1(2));
            end
        else            
            if mod(innert,2)==0
                arn = real(ifft2(cvfftl.*cvfft2./(abs(cvfft2).^2+ep)));
                arn = arn(1:size1(1),1:size1(2));
                array1(1:size1(1),1:size1(2)) = topelements(arn,numelements(arn,innert/inneritr));
            else
                arn = real(ifft2(cvfftl.*cvfft1./(abs(cvfft1).^2+ep)));
                arn = arn(1:size2(1),1:size2(2));
                array2(1:size2(1),1:size2(2)) = topelements(arn,numelements(arn,innert/inneritr));
            end
        end
    end
    cvfft1=fft2(array1);
    cvfft2=fft2(array2);
    cvfftl=cvfft1.*cvfft2;
    cvl=real(ifft2(cvfftl));
    %ex=2;score=(sum(sum(abs(cvl-circshift(cvl,[0 1])).^ex))+sum(sum(abs(cvl-circshift(cvl,[1 0])).^ex))); %/(sum(sum(abs(cvl).^ex)));
    ex=2;score=sum(sum(abs(cvl).^ex));
    if score<bestscr
        bestscr=score;
        score
        best1=array1;
        best2=array2;
        bestcvl=cvl;
        figure(2);clf;
        subplot(2,3,1);
        imagesc(best1);colormap(gray);
        title('Pattern 1');
        subplot(2,3,2);
        imagesc(best2);colormap(gray);
        title('Pattern 2');
        subplot(2,3,3);
        bst2=bestcvl(1:size(bestcvl,1)-1,1:size(bestcvl,2)-1);
        bns=0:max(max(bst2));
        hist(reshape(bst2,[1 numel(bst2)]),bns);
        title('Redundancy Histogram');
        subplot(2,3,4);
        imagesc(bestcvl);colormap(gray);colorbar;
        title('Fourier Coverage');
        subplot(2,3,5);
        imagesc(fftshift(abs(cvfftl)));colormap(gray);colorbar;
        title('PSF Raw');
        subplot(2,3,6);
        imagesc(fftshift(abs(fft2(cvl > 0.5))));colormap(gray);colorbar;
        title('PSF Equalized');
        drawnow;pause(0.02);
    end
end
return;

function tp=topelements(x,n)

[r,c]=size(x);
y=reshape(x,[1 r*c]);
y=sort(y);
tp=x >= y(1,max(r*c-n+1,1));
return;
