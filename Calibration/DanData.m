data
figure(201)
subplot(5,1,1)

%signal 1
ft_dif12=zeros(1,size(data.ctf,2));

for nd=16:16;
sig1=1;
meas1=data.ctf(sig1,:,nd)
plot(data.f_vec,abs(meas1))
%time signal
num_f=length(data.f_vec);
B=data.f_vec(end)-data.f_vec(1);

df=B/(num_f-1)
ft_meas1=ifftshift(fft(fftshift(meas1)))
subplot(5,1,2)
t=0:1/B:2*1/(2*df);
plot(t,abs(ft_meas1))

%signal 2
sig2=2;
meas2=data.ctf(sig2,:,nd)
subplot(5,1,3)
plot(data.f_vec,abs(meas2))
%time signal

ft_meas2=ifftshift(fft(fftshift(meas2)))
subplot(5,1,4)
plot(t,abs(ft_meas2))

%subtract the two signals
subplot(5,1,5)
ft_dif12=(ft_meas2-ft_meas1)/24+ft_dif12;
plot(t,abs(ft_dif12))
end
%step size