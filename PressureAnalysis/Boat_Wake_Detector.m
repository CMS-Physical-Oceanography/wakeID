clear all;
%load ("D:\BoatTracking_Fall2022\ICW_survey_20220922\rbr_pressure_CMSsouth_2022.mat")
% this data runs from '22-Sep-2022 15:24:03' to'11-Oct-2022 10:05:31'.
load("D:\BoatTracking_Fall2022\ICW_survey_20220922\rbr_pressure_CMSsouth_2022.mat")

P=data(2).pres;
time=data(2).time;
clear data

D1= 'September 23, 2022 09:48:32';     % Time frame for data to be analyzed 
D2= 'September 23, 2022 09:58:32';
    t1= datenum(D1);
    t2= datenum(D2);
fs= 16.66;      % Sample Rate 


ind = find(time>=t1&time<=t2+1/16/86400); % Framing data to match time frame
time0 = time(ind);
seconds = (time0-time0(1))*86400;
N= length(ind);

Pmean = mean(P);
F= hanning(round(32*fs)); 
F=F./sum(F);
Pconv= conv(P-Pmean,F,"same");
PerP= P-Pconv-Pmean;   % perturbation pressure

P0 = PerP(ind);

Ts = 32 ; % length of window in seconds
dTs = 0.5 ; % overlap in seconds
L= round(Ts*fs);
Olap= L-round(dTs*fs);
dfp= fs/L;
[s, f, tfs]=spectrogram(P0,hamming(L),Olap,0:fs/L:1.5,fs,"yaxis"); 
S = (s.*conj(s)/(L*fs));

Rmin = min(S(:));
Rmax = max(S(:));
C = 255 - round(((S-0)/(.00001-0)).*255);
C= uint8(C);

Fig = figure;
Fig.Name = sprintf ('Spectorgram L%d Olap%d',L,Olap);
imagesc(C); hold on; axis off
colormap(bone);
caxis([0 255])



%% This will attempt to detect boats based off their cyclical frequency

imwrite(C,"C:\Users\lwlav\OneDrive\Documents\MATLAB\Pressure\TrainingImages\pressure_Glitchy2boats_Sept2323.tif","tiff")


