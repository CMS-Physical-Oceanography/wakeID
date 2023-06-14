% This set of code will give you a set of different resolutions based on
% different window sizes and overlap lengths based on time and sample rate.

close all;

load ('processed_hatchery_june2022.mat',"time", "P")

D1= 'June 08, 2022 15:00:00';     % Time frame for data to be analyzed 
D2= 'June 08, 2022 15:10:00';
    t1= datenum(D1);
    t2= datenum(D2);

ind = find(time>=t1&time<=t2+1/16/86400); % Framing data to match time frame
time0 = time(ind);
seconds = (time0-time0(1))*86400;
P0 = P(ind);
P1= detrend(P0);

fs= 16.66;      % Sample Rate 
N= length(ind);

Ts = [6,8,12,15,20,30] ; % in seconds
dTs = 0.5 ; %in seconds
for ii= 1:length(Ts)
    T= Ts(ii);
    L= round(T*fs);
        Olap= L-round(dTs*fs);
        [s, f, tfs]=spectrogram(P1,L,Olap,0:fs/N:3,fs,"yaxis"); 
        Ns = length(f);
        Fig = figure;
        Fig.Name = sprintf ('Spectorgram L%d Olap%d',L,Olap);
        imagesc(tfs,f,(s.*conj(s)/Ns*fs));
        set(gca,'ylim',[0 1.5]);
        saveas(gcf,sprintf('2ndFigure%d.fig',ii))
end
