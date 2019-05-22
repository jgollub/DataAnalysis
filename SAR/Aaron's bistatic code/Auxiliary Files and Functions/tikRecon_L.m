function im = tikRecon_L(h,g)

addpath(genpath('C:\Program Files\MATLAB\Regtools'))

[U,s,V] = csvd(double(h));
[reg_corner,~,~,~] = l_curve(U,s,g,'Tikh');

D = zeros(size(diag(s)));
for ss = 1:min(size(diag(s)))
    D(ss,ss) = s(ss)/(s(ss)^2 + reg_corner^2);
end

im = V*D'*U'*g;

end
















