%% ===========================================================================
%generate a random phasor and apply envelop that represents resonant
%behavoir



%bandwidth is 9 Ghz across the the K-band (17.5-26.5 Ghz), at nyquist 1/2B
dt=1/(2*9e9);
df=9e9/101;
T=0:dt:1/df;


%amplitude
Hn=1;
%signal sampling rate
Hn=Hn*1/sqrt(2)*(randn(1,length(T))+1j*randn(1,length(T)));
figure(1); cla; 
subplot(1,2,1); 
    scatter(real(Hn),imag(Hn)); axis equal; axis tight;
subplot(1,2,2)
    plot(T,abs(Hn));
    mean(abs(Hn))
% The Q of the cavity is related to the decay constant by e^-t/tau
Q=1000;
tau=Q/(2*pi*22e9*2);
Hn=Hn.*exp(-T./tau);
hold on; plot(T,abs(Hn),'r');

HFft=fft(Hn);
figure(3); hold on; plot(abs(HFft),'r');