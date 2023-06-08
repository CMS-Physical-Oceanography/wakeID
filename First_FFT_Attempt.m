% Pull Data for Analysis
clear all

load ('processed_hatchery_june2022.mat',"time", "P")

% Plot a chart of pressue against time

D1= 'June 08, 2022 15:34:00';     % Time frame for data to be analyzed 
D2= 'June 08, 2022 15:35:00';

t1= datenum(D1);
t2= datenum(D2);

T=[t1,t2];

figure,
plot (time,P,'-')
xlabel("time")
ylabel("Pressure")
title('Pressure vs Time for 1 wake')
set(gca,'xlim',T, 'xtick',t1:1/24/60/10:t2 )
datetick(gca,'x','mm/dd HH:MM:SS', 'keeplimits','keepticks')

% Find Frequency Ranges

ind = find(time>=t1&time<=t2+1/16/86400);
time0 = time(ind);
P0 = P(ind);

fs= 16.66;      
N = 962;
fr=0:fs/N:3;

P1= detrend(P0);
sigp1 = var(P1);
w = hann(N);
P2 = w'.*P1;
sigp2 = var(P2);
P3 = P2*sqrt(sigp1/sigp2);

[Ys,f,t] = spectrogram(P3,hann(N,"symmetric"),.5*N,fr,fs,"power", "xaxis");

Yp = fft(P1)/N; 
Yw = (fft(P3))/N;
dfp= (fs/N);                      % fundamental frequency is 1/(n*dt)
fnp= 1/2;                         % nyquist freq is 1/(2*dt)
if mod(N,2)==0
fp = [(0:(N)/2),-(N/2-1):-1]*(fs/N);         % frequencies 
else 
    fp=[(0:(N-1)/2),-(N-1)/2:-1]*(fs/N);
end

Yp = detrend(Yp) ;
Yw = detrend (Yw) ;
Ys = detrend (Ys) ;

% generate the figure and axes
Spp = Yp.*conj(Yp)/dfp;
Sww = Yw.*conj(Yw)/dfp;
Sss = Ys.*conj(Ys)/(N^2*dfp);

Spp = medfilt1(Spp,10);
Sww = medfilt1(Sww,10);
Sss = medfilt1(Sss,10);

figure,
semilogy(abs(fp),Spp,'-k',abs(fp),Sww,'-r');
xlabel('Frequency (Hz)')
ylabel('Magnitude of Spectra')
title ('One minute of time FFT')

hold on
  semilogy(abs(f),Sss,'--g')
  legend ('Non-Windowed','Manually Windowed','Spectrograms Output')
hold off


%%
%-----------------Redo Process with 5 times the sample size------------

D1= 'June 08, 2022 13:17:30';     % Time frame for data to be analyzed 
D2= 'June 08, 2022 13:22:30';

t1= datenum(D1);
t2= datenum(D2);

T=[t1,t2];

figure,
plot (time,P,'-')
xlabel("time")
ylabel("Pressure")
title('5 minute Pressure vs Time')
set(gca,'xlim',T, 'xtick',t1:1/24/60/10:t2 )
datetick(gca,'x','mm/dd HH:MM:SS', 'keeplimits','keepticks')

ind = find(time>=t1-4/16/86400 & time<=t2+7/16/86400);
time0 = time(ind);
P0 = P(ind);

fs= 16.66;  
N = length(P0);
L = 5;
Olap = 1;
fr = 0:fs/N:3;

P1= detrend(P0);
sigp1 = var(P1);
w = hann(L);
P2 = buffer (P1,L);
numW = ((length(P1))/L);
Pw = [];

for Wx = 1:numW     
   win{Wx} = P2(:,Wx).*w;
   Pw = [ Pw; win{Wx}];
end
 
sigp2 = var(Pw);
P3 = Pw'*sqrt(sigp1/sigp2);

[Ys,f,t] = spectrogram(P3,L,Olap,fr,fs,"power", "xaxis");

Yp = fft(P1)/N; 
Yw = (fft(Pw))./N;
dfp= (fs/N);                      % fundamental frequency is 1/(n*dt)
fnp= 1/2;                         % nyquist freq is 1/(2*dt)
if mod(N,2)==0
    fp = [(0:(N)/2),-(N/2-1):-1]*(fs/N);         % frequencies 
else 
    fp=[(0:(N-1)/2),-(N-1)/2:-1]*(fs/N);
end

Yp = detrend(Yp) ;
Yw = detrend (Yw) ;
Ys = detrend (Ys) ;

% generate the figure and axes
Spp = Yp.*conj(Yp)/dfp;
Sww = Yw.*conj(Yw)/dfp;
Sss = Ys.*conj(Ys)/(N^2*dfp);

Spp = medfilt1(Spp,50);     % Noise reduction by applying a median filter 
Sww = medfilt1(Sww,50);
AvgSss = mean(Sss,2);

figure,
semilogy(abs(fp),Spp,'-k',abs(fp),Sww,'-r');
xlabel('Frequency (Hz)')
ylabel('Magnitude of Spectra')
title ('5 minute FFT')

hold on
  semilogy(abs(f),AvgSss,'--b')
  legend ('Non-Windowed','Manually Windowed','Spectrograms Averaged Output')
hold off

figure,
spectrogram(P3,L,Olap,fr,fs,"power", "xaxis");
title  'Spectrogram for 5 minutes of data'

%%
%--------------------------------------------------------------------
%                     Very Large Set of Data

D1= 'June 08, 2022 15:00:00' ;     % Time frame for data to be analyzed 
D2= 'June 08, 2022 15:30:00';

t1= datenum(D1);
t2= datenum(D2);

T=[t1,t2];

figure,
plot (time,P,'-')
xlabel("time")
ylabel("Pressure")
title('2 hours Pressure vs Time')
set(gca,'xlim',T, 'xtick',t1:1/24/4:t2 )
datetick(gca,'x','mm/dd HH:MM:SS', 'keeplimits','keepticks')


ind = find(time>=t1-1/16/86400 & time<=t2+8/16/86400);
time0 = time(ind);
P0 = P(ind);


fs= 16.66;  
N = length(P0);
L = 43;
Olap = 41;
fr = 0:fs/N:3;

P1= detrend(P0);
sigp1 = var(P1);
w = hann(L);
P2 = buffer (P1,L);
numW = ((length(P1))/L);
Pw = [];

for Wx = 1:numW ;
    
   win{Wx} = P2(:,Wx).*w;
   Pw = [ Pw; win{Wx}];

end

sigp2 = var(Pw);
P3 = Pw'*sqrt(sigp1/sigp2);

[Ys,f,t] = spectrogram(P3,hann(L,"symmetric"),Olap,fr,fs,"power", "xaxis");

Yp = fft(P1)/N; 
Yw = (fft(Pw))./N;
dfp= (fs/N);                      % fundamental frequency is 1/(n*dt)
fnp= 1/2*16.66;                         % nyquist freq is 1/(2*dt)
if mod(N,2)==0
fp = [(0:(N)/2),-(N/2-1):-1]*(fs/N);         % frequencies 
else 
    fp=[(0:(N-1)/2),-(N-1)/2:-1]*(fs/N);
end

Yp = detrend(Yp) ;
Yw = detrend (Yw) ;
Ys = detrend (Ys) ;

% generate the figure and axes
Spp = Yp.*conj(Yp)/dfp;
Sww = Yw.*conj(Yw)/dfp;
Sss = Ys.*conj(Ys)/(N^2*dfp);

Spp = medfilt1(Spp,100);     % Noise reduction by applying a median filter 
Sww = medfilt1(Sww,100);

sumSss = sum(Sss);
[MaxSss, MaxSssI] = max(sumSss);
MaxSssI = Sss(:,MaxSssI);
[MinSss, MinSssI] = min(sumSss);
MinSssI = Sss(:,MinSssI);

AvgSss = mean(Sss,2);        


%Since Spectrogram outputs every single window as a row of the matrix Sss
%then we average these windows together. 

figure,
semilogy(abs(fp),(Spp),'-k',abs(fp),(Sww),'-r',abs(f),(AvgSss),'-b',abs(f),(MaxSssI),'-g',abs(f),MinSssI,'-g');
xlabel('Frequency (Hz)')
ylabel('Magnitude of Spectra')
title ('2 hours FFT')
legend ('Non-Windowed','Manually Windowed','Spectrograms Averaged Output','Min and Max Windows')


figure,
spectrogram(P3,hann(L,"symmetric"),Olap,fr,fs,"power", "xaxis")
title  'Spectrogram for two hours of data'

fig9 = figure ;
imagesc(t,f,log10(Sss))
colormap(bone)
