   
% g_up=cell();
% g_down=cell();
% g_sig=[];

cycles=20
g_sig=[];
for ii=1:cycles
ft245bitbang(1)
pause(.5)
g=MeasureVna(vna,txSelect,rxSelect,gAverage);
g_up{ii}=g{1,3};
ft245bitbang(0)
pause(.5)
g=MeasureVna(vna,txSelect,rxSelect,gAverage);
g_down{ii}=g{1,3};
g_sig{ii}=g_up{ii}-g_down{ii};
end
g_ave1=zeros(101,1);
ave=numel(g_sig)
for jj=1:cycles
g_ave1=g_sig{jj}+g_ave/ave
end 

g_c1=cell(1,1);
g_c1{1}=g_ave1;

% figure(302); hold on; plot(linspace(0,1/(f(2)-f(1)),101)*3e8/2,abs(fft(g_ave)))

g_sig=[];
for ii=1:cycles
ft245bitbang(1)
pause(.5)
g=MeasureVna(vna,txSelect,rxSelect,gAverage);
g_up{ii}=g{1,3};
ft245bitbang(0)
pause(.5)
g=MeasureVna(vna,txSelect,rxSelect,gAverage);
g_down{ii}=g{1,3};
g_sig{ii}=g_up{ii}-g_down{ii};
end
g_ave2=zeros(101,1);
ave=numel(g_sig)
for jj=1:cycles
g_ave2=g_sig{jj}+g_ave/ave
end 

g_c2=cell(1,1);
g_c2{1}=g_ave2

for i=1:size(g_c1,1)
    for j=1:size(g_c1,2)
        for k=1:size(g_c1,3)
            % Cross correlate
            xc{i,j,k}=abs(ifftshift(ifft(fftshift(g_c1{i,j,k}.*conj(g_c2{i,j,k})))));
            % Score
            merit(i,j,k)=max(xc{i,j,k}(45:55))/mean(xc{i,j,k});
        end
    end
end

figure(); plot(10*log10(abs(xc{1})))

%variance

