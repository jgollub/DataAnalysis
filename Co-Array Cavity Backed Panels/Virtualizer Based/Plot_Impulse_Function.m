         B=(max(panel.f)- min(panel.f));
          dt=1/(2*B);
          df=B/numel(panel.f);
          T=0:1/(2*B):1/(2*df)-1/(2*B); 
          
          h_t=1/sqrt(2)*(randn(numel(theta),numel(T))...
                         +1j*randn(numel(theta),numel(T)));
            Q_cavity_1=10;
             Q_cavity_2=100;
              Q_cavity_3=1000;
                     
             tau1=2*Q_cavity_1/(2*pi*mean(panel.f));
             tau2=2*Q_cavity_2/(2*pi*mean(panel.f)); 
             tau3=2*Q_cavity_3/(2*pi*mean(panel.f)); 
             
hold on;     figure(); plot(T,abs(h_t(1,:)),'-'); 
             h_t1=bsxfun(@times, h_t, exp(-T./tau1));
              h_t2=bsxfun(@times, h_t, exp(-T./tau2));
               h_t3=bsxfun(@times, h_t, exp(-T./tau3));
               
figure(6); cla
t=colormap(lines);
% subplot(1,2,1);
hold on; plot(T/1e-9,abs(h_t3(1,:)).^2,'-','color',t(3,:),'linewidth',2); 
         plot(T/1e-9,abs(h_t2(1,:)).^2,'-','color',t(2,:),'linewidth',2);
         plot(T/1e-9,abs(h_t1(1,:)).^2,'-','color',t(1,:),'linewidth',2);
         axis tight;
           xlabel('Time (ns)')
          ylabel('Magnitude (a.u.)')
          legend('Q=1000','Q=100','Q=10')
          
          H_f1=fft(h_t1,[],2)*dt;
               H_01=mean(mean(abs(H_f1)));
          H_f1=H_f1./H_01;
          H_f2=fft(h_t2,[],2)*dt;
               H_02=mean(mean(abs(H_f2)));
          H_f2=H_f2./H_02;
          H_f3=fft(h_t3,[],2)*dt;
           H_03=mean(mean(abs(H_f3)));
          H_f3=H_f3./H_03;
          
         figure(7);cla; hold on; plot(panel.f/1e9,abs(H_f3(1,:)).^2,'-','color',t(3,:),'linewidth',1.5);
        plot(panel.f/1e9,abs(H_f2(1,:)).^2,'-','color',t(2,:),'linewidth',1.5);
        plot(panel.f/1e9,abs(H_f1(1,:)).^2,'-','color',t(1,:),'linewidth',1.5);
         axis tight;
           xlabel('Frequency (GHz)')
          ylabel('Magnitude (a.u.)')      
        legend('Q=1000','Q=100','Q=10')
          