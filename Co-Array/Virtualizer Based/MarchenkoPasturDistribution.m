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

%clear;
%close all;
N=1000;   %~measurements
T=700;  %~number of voxels

% Ratio of matrix dimensions
c=N/T;

% Sample
x=randn(N,T); % Normal distribution
%x=rand(N,T);  % Uniform distribution

s=std(x(:));

% spectral matrix
r=x*x'/T;

ex=1.0;

%eigenvalues
l=eig(r).^ex;
% Probability Density Function 
% number of points for measurement.
n=50;
% Boundaries 
a=(s^2)*(1-sqrt(c))^2;
b=(s^2)*(1+sqrt(c))^2;

[f,lambda]=hist(l,linspace(a.^ex,b.^ex,n));
% Normalization
f=f/sum(f);
% Theoretical pdf
ft=@(lambda,a,b,c) ((1/ex).*(lambda.^(1/ex-1))).*(1./(2*pi*(lambda.^(1/ex))*c*s^(2))).*sqrt((b-(lambda.^(1/ex))).*((lambda.^(1/ex))-a));
F=ft(lambda,a,b,c);
% Processing numerical pdf
F(isnan(F))=0;
F=F/sum(F);
if (N>T)
    F=F*(T/N);
end

% Results
figure;
h=bar(lambda,f);
set(h,'FaceColor',[.75 .75 .8]);
set(h,'LineWidth',0.25);
xlabel('Eigenvalue \lambda');
ylabel(' Probability Density Function f(\lambda)');
title(' MarchenkoPastur distribution');
lmin=min(l);
lmax=max(l);
%axis([-1 2*lmax 0 max(f)+max(f)/4]);

hold on;
plot(lambda,F,'g','LineWidth',2);
hold off;
