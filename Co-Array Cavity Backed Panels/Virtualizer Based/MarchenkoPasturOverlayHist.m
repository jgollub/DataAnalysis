% Marchenko Pastur Distribution

% In Random Matrix Theory, MP law gives the probability density function
% of singular values of large rectangular random matrices;
% when the dimensions of matrix tend to infinity.

% This contribution illustrates the PDF of matrix Y(N,N)=(T^-1)X*X^T, 
% where X is random matrix whose entries X_i,j are independent 
% and identically distributed random variables with zero mean
% and variance s^2. The program is applicable for both uniform and random
% distributions.

% Ref :
% Marchenko,V. A., Pastur, L. A. (1967) "Distribution of eigenvalues for some sets of
% random matrices", Mat. Sb. (N.S.), 72(114):4, 507536

% (c) Youssef KHMOU, Applied Mathematics ,30 January,2015.


function y = MarchenkoPasturOverlayHist(hc,bins,c)

ex=1;
ovshoot=1.0;

s=1;
% Boundaries 
a=(s^2)*(1-sqrt(c))^2;
b=(s^2)*(1+sqrt(c))^2;
b2=ovshoot*(b-a)+a;
hc=ovshoot*(hc-min(hc))*((b-a)/(max(hc)-min(hc)))+a;

[f,lambdalin]=hist(hc,linspace(a,b2,bins));
lambda=linspace(a.^ex,b.^ex,bins);
% Normalization
f=f/sum(f);
% Theoretical pdf
ft=@(lambda,a,b,c) ((1/ex).*(lambda.^(1/ex-1))).*(1./(2*pi*(lambda.^(1/ex))*c*s^(2))).*sqrt((b-(lambda.^(1/ex))).*((lambda.^(1/ex))-a));
F=ft(lambda,a,b,c);
% Processing numerical pdf
F(isnan(F))=0;
F=F/sum(F);
%if (c>1)
%    F=F/c;
%end

% Results
h=bar(lambdalin,f);
set(h,'FaceColor',[.75 .75 .8]);
set(h,'LineWidth',0.25);
xlabel('Singular Value');
ylabel(sprintf('Singular Value Count (%%)'));
title(sprintf('Marchenko Pastur distribution \\lambda=%g',c));
%lmin=min(l);
%lmax=max(l);

hold on;
plot(lambda.^(1/ex),real(F),'g','LineWidth',2);
hold off;
