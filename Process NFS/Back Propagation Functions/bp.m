function [magret]=bp(measurements,X,Y,range)

%load('d4_fid_d5_4mm_samp_tilted_7_25_deg.mat');
sampdist=(abs(X(1,1)-X(1,2)))/1000; %NFS sampling period
%sampdist=0.004;
ztrans=range;
    [ry,rx,nf,p]=size(measurements);
    by=2^(ceil(log(ry)/log(2)));
    bx=2^(ceil(log(rx)/log(2)));
    fmax=26.5e9;
    fmin=17.5e9;
    speedoflight=3e8;
    magaddsq=zeros(by,bx);
    for frno=1:nf
        frcur=fmin+(frno-1)*(fmax-fmin)/(nf-1);
        m=zeros(by,bx,3);
        m(floor((by-ry)/2)+1:floor((by-ry)/2)+ry,floor((bx-rx)/2)+1:floor((bx-rx)/2)+rx,1:p) = ...
               squeeze(measurements(:,:,frno,:));   
           
        wavl=speedoflight/frcur;
        %fld=transfield(m,[0 0 ztrans/wavl],wavl/sampdist,wavl/sampdist);
        fld=transfield(m,[0 0 ztrans/wavl],wavl/sampdist,wavl/sampdist);
        
        switch p
            case 1
        magaddsq=magaddsq+sqrt(abs(fld(:,:,1)).^2+abs(fld(:,:,2)).^2);
            case 2
        magaddsq=magaddsq+sqrt(abs(fld(:,:,1)).^2);
        end
    end

magret=magaddsq(floor((by-ry)/2)+1:floor((by-ry)/2)+ry,floor((bx-rx)/2)+1:floor((bx-rx)/2)+rx);
    
%     figure(1);
%     ndisp=0.3;
%     imagesc(unique(fx),unique(fy),abs(magaddsq).^ndisp);colormap(hot);colorbar;
%     xlabel(sprintf('ztrans=%g',ztrans));%axis xy;
%     drawnow;
    %pause(0.2);
return;