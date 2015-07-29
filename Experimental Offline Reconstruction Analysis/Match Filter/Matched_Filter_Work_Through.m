
data=1:.1:4*pi;
act_sig=sin(data);
meas=sin(data)'+2*randn(length(data),1);
try figure(meas_plot)
catch
    clf;
    meas_plot=figure()
end

plot(meas)

%match filter
R=cov(2*randn(length(data),1));
K=R^(-1)./sqrt(meas'*R^(-1)*meas);
h_hat=K*meas;
h=fliplr(h_hat');
y=filter(h,1,meas);
hold on
plot(y,'-.r*')
plot(act_sig,'g')
