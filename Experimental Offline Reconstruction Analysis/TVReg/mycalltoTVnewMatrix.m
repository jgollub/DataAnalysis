function [img_estimated]=mycalltoTVnewMatrix(mycube,th,piter)


% x = zeros(squeeze(size(mycube(:,:,1))));
img_estimated=zeros(size(mycube));

[uy ux] = size(mycube(:,:,1));

dh = @(x) mycircconv1(x);
dv = @(x) mycircconv3(x);
dht = @(x) mycircconv2(x);
dvt = @(x) mycircconv4(x);
vect = @(x) x(:);
opQ = @(x) [vect(dh(x)) vect(dv(x))] ;
opQt = @(x) dht(reshape(x(:,1),uy,ux))+dvt(reshape(x(:,2),uy,ux));


for i = 1:size(mycube,3)
    x = mycube(:,:,i);
    
    
    img_estimated(:,:,i) = x - projk(x,th/2,opQ,opQt,piter);
end

end

function y = mycircconv1(x)
x = [x zeros(size(x,1),1)] - [zeros(size(x,1),1) x];
x(:,end) = x(:,end)+x(:,1);
y = x(:,2:end);
end


function y = mycircconv2(x)
x = [zeros(size(x,1),1) x] - [x zeros(size(x,1),1)];
x(:,1) = x(:,end)+x(:,1);
y = x(:,1:end-1);
end


function y = mycircconv3(x)
x = [x; zeros(1,size(x,2))]- [zeros(1,size(x,2)); x];
x(end,:) = x(end,:)+x(1,:);
y = x(2:end,:);
end

function y = mycircconv4(x)
x = [zeros(1,size(x,2)); x]-[x; zeros(1,size(x,2))];
x(1,:) = x(end,:)+x(1,:);
y = x(1:end-1,:);
end