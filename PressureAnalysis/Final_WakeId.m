%% Final Product of everything
clear all; close all;
load("D:\BoatTracking_Fall2022\ICW_survey_20220922\rbr_pressure_CMSsouth_2022.mat") %Load Data
load("D:\WakeID_Summer2023\imgNet_v4.mat")      %Load Network
load("D:\WakeID_Summer2023\Met_Data\Winds.mat") % Load Met data

P=data(2).pres;             % Take data from surface sensor
T=data(2).time;
clear data

fs = 16.66;
Pmean = mean(P);
F= hanning(round(32*fs)); 
F=F./sum(F);
Pconv= conv(P-Pmean,F,"same");
P0= P-Pconv-Pmean;   % perturbation pressure

%% Transform the data into frequency data to be analyzed and then into image data
% then feed it through the network to ID the boats. 
DataStart = datetime(datestr(T(1)));
DataStartVec = datevec(DataStart);

vect = datetime(datevec(T));

LP= length(P);                  
seg = round((LP/fs)/600);
N = round(LP/seg);

Ts = 32 ; % length of window in seconds
dTs = 0.5 ; % overlap in seconds
L= round(Ts*fs);
Olap= L-round(dTs*fs);
dfp= fs/L;


%% Loop to tally number of boat wakes
 %Date of when the data begins in order to log into spreadsheet
dy= 1;
timeu = unique(dateshift(vect,'start','days'));
Dailynumb = zeros(1,length(timeu));

wkeEind = [];% energy of individual wakes in 10 min window
wkeEseg = [];% energy of all wakes in current 10 min window
wkeEday = zeros(1,length(timeu));

wndEday = zeros(1,length(timeu));
wndEseg = [];
vidNum  = [];% which video are these wakes from

wkeIind = [];% individual wake indices in spectrogram of current video
wkeTind = [];% same with time in seconds relative to dataStart
wkeIseg = 1:N:N*seg;% 10 min segment indices in spectrogram of current video
wkeTseg = T(wkeIseg);% same with time in seconds relative to dataStart

Fig4= figure;
Fig1 = figure;
Fig3=figure;
Fig6 = figure;
Fig5 = figure;
for ii = 1:seg
    t1 = T(1+N*(ii-1):N*ii);
    P1 = P0(1+N*(ii-1):N*ii);
    [s,f,tfs]=spectrogram(P1,hamming(L),Olap,0:fs/L:1.5,fs,"yaxis",'psd');
    Nf = length(f);
    Nt = length(tfs);
    S = (s.*conj(s)/(L*fs));
    C = 255 - round(((S-0)/(.00001-0)).*255);
    C= uint8(C);
    act = activations(imgNet,C,'softmax');
    Bw = act(1:Nf,1:Nt,2)>=.75;
    CleanC= bwareaopen(Bw,100);
    [labim,number] = bwlabel(CleanC);
    labimS = max(labim,[],1);
    startTs = [];
    endTs   = [];
    for ww = 1:max(labimS(:))
        startTs(ww) = find(labimS==ww,1,'first');
        endTs(ww)   = find(labimS==ww,1,'last');
    end
    % create a binary mask for what is a wake vs wind
    whatwake = labimS>0;
    wakeMask = repmat(whatwake,Nf,1);
    % % find start/end indices of wakes
    % iswake = diff(whatwake);
    % firstwkeind=find(iswake~=0,1,'first');
    % lastwkeind=find(iswake~=0,1,'last');
    % if iswake(firstwkeind)==-1
    %     iswake(1)=1;
    % end
    % if iswake(lastwkeind)==1
    %     iswake(end)=-1;
    % end
    % startTs= find(iswake>0)+1;
    % endTs  = find(iswake<0);
    % stats1 = regionprops(CleanC,'centroid');
    % stats2 = regionprops(CleanC,'MajorAxisLength');
    % stats3 = regionprops(CleanC,'MinorAxisLength');
    % centers = cat(1,stats1.Centroid);
    % Fwidth = cat(1,stats3.MinorAxisLength);
    % Twidth = cat(1,stats2.MajorAxisLength);
    vidNum  = cat(1,vidNum, 0*startTs'+ii);
    wkeIind = cat(1,wkeIind,[startTs', endTs']);
    wkeTind = cat(1,wkeTind,[tfs(startTs)', tfs(endTs)']+N*(ii-1)/fs);
    wkeEseg = cat(1,wkeEseg, sum(sum(S.*wakeMask))*dfp*dTs);
    wndEseg = cat(1,wndEseg, sum(sum(S.*~wakeMask))*dfp*dTs);
    E = sum(sum(S))*dfp*dTs;
    for jj = 1:number
        startTind = startTs(jj);
        endTind   = endTs(jj);
        startTsec = tfs(startTind); 
        endTsec   = tfs(endTind); 
        Ef   = sum(S(:,:),1).*dfp;
        Ewke = sum(Ef(startTind:endTind))*dTs;
        indP = P1(round(startTind*dTs*fs):round(endTind*dTs*fs));
        %data.DD(dy).seg(ii).WkeT(jj) = DataStart+startTsec/24/60/60;
        %data.DD(dy).seg(ii).WkeE(jj) = Ewke;
        wkeEind = cat(1,wkeEind,Ewke);
        Dailynumb(dy) = Dailynumb(dy) + number;
        wkeEday(dy) = wkeEday(dy) + Ewke;
        E = E - Ewke;
    end
    wndEday(dy) = wndEday(dy) +E;
    dyswitch = datetime(datevec(T(1+N*(ii-1))));
    dyswitch2 = datetime(datevec(T(1+N*(ii))));
    Dys = dyswitch.Day;
    Dys2 = dyswitch2.Day;
    if Dys~=Dys2
        dy=dy+1;
    end
    
    if ii == seg/2
    DISP('50%')
    end

% Pressure plots
FontSizeTitle=48;
FontsizeTick=28;
FontsizeLabel=28;


% figure(Fig1)
% axs1 = subplot(2,1,1);
% plot (T,(P-Pmean),'-')
% xlabel("$t$ [mm/dd]", 'interpreter','latex','FontSize',FontsizeLabel)
% ylabel("$P$ [m]", 'interpreter','latex','FontSize',FontsizeLabel)
% title("Pressure vs Time", 'September 22 2022 to October 11 2022','FontSize',FontSizeTitle)
% set(axs1,'xlim',[T(1) T(end)],'FontSize',FontsizeTick)
% datetick(axs1,'x','mm/dd', 'keeplimits','keepticks')
% 
% axs2 = subplot(2,1,2);
% plot (T,P0,'-')
% xlabel("$t$ [mm/dd]", 'interpreter','latex','FontSize',FontsizeLabel)
% ylabel("$P'$ [m]", 'interpreter','latex','FontSize',FontsizeLabel)
% title("Perterbation Pressure vs Time",'FontSize',FontSizeTitle)
% set(axs2,'xlim',[T(1) T(end)],'FontSize',FontsizeTick)
% datetick(axs2,'x','mm/dd', 'keeplimits','keepticks')
% xline(datenum(datestr(DataStart))+wndTday(:,1)/86400,'r-','linewidth',.5)

% 
% 
%  figure(Fig6)
%       A = labeloverlay(C,labim);
%       imagesc(A)
%       xline(wkeIind(vidNum==ii,1),'g-','linewidth',1.5)
%       xline(wkeIind(vidNum==ii,2),'r-','linewidth',1.5)
%     xticks([])
% yticks([])

% 
% figure(Fig3)
% plot (indP,'-')
% hold on;
% axs3 = gca;
% xlabel("$t$ [HH:MM:SS]", 'interpreter','latex','FontSize',FontsizeLabel)
% ylabel("$P'$ [m]", 'interpreter','latex','FontSize',FontsizeLabel)
% title("Perterbation Pressure vs Time", 'September 22 2022 16:02-16:04','FontSize',FontSizeTitle)
% set(axs3,'TickLabelInterpreter','latex','TickDir','in','xlim',[1 1458],'FontSize',FontsizeTick,'xtick',1:182.1:1458)
% xticklabels({'16:02:00' '16:02:15' '16:02:30' '16:02:45' '16:03:00' '16:03:15' '16:03:30' '16:03:45' '16:04:00'})
% 
% % geofig =figure;
% geoscatter(34.139610, -77.862707,'red','filled'); hold on;
% geoscatter(34.139956, -77.862401,'yellow','filled')
% axsg = gca;
% title('Camera and Bouy Location','FontSize',FontSizeTitle)
% set(axsg,'Basemap','satellite','ZoomLevel',18,'FontSize',FontsizeTick)
% hold off

% figure(Fig4)
% imagesc(tfs,f,S)
% axs4=gca;
% cbar=colorbar;
% title("Spectrogram from 09/22/23 15:54-16:04", 'Interpreter','latex','FontSize',FontSizeTitle)
% xlabel(axs4,"$t$ [HH:MM]", 'Interpreter','latex','FontSize',FontsizeLabel)
% ylabel(axs4,"$f$ [Hz]", 'Interpreter','latex','FontSize',FontsizeLabel)
% set(axs4,'FontSize',FontsizeTick,'xtick',15.9:113.51:583.58)
% set(cbar,"FontSize",FontsizeTick)
% xticklabels({'15:54' '15:56' '15:58' '16:00' '16:02' '16:04'})
% colormap(flipud(bone))
% cbar.Label.String= "$S_{ss}$ [m$^2$/Hz]";
% cbar.Label.Interpreter = "latex";
% cbar.Label.Rotation = 270;
% cbar.Label.VerticalAlignment = "bottom";
% 
% figure(Fig5)
% axs5 = subplot(1,1,1);
% plot (P1,'-')
% xlabel("$t$ [HH:MM]", 'interpreter','latex','FontSize',FontsizeLabel)
% ylabel("$P'$ [m]", 'interpreter','latex','FontSize',FontsizeLabel)
% title("Perterbation Pressure vs Time",'FontSize',FontSizeTitle)
% set(axs5,'FontSize',FontsizeTick,'xtick',1:1999:9996)
% xticklabels({'15:54' '15:56' '15:58' '16:00' '16:02' '16:04'})
%   xline(wkeInds(vidNum==ii,1)*dTs*fs,'g-','linewidth',1.5)
%       xline(wkeInds(vidNum==ii,2)*dTs*fs,'r-','linewidth',1.5)  

% Fig6=figure;,
% plot (T,P,'-')
% xlabel("time",'FontSize',FontsizeLabel)
% ylabel("Pressure",'FontSize',FontsizeLabel)
% title('Pressure vs Time for 1 wake','FontSize',FontSizeTitle)
% set(gca,'xlim',T, 'xtick',t1:1/24/60/10:t2,'FontSize',FontsizeTick)
% datetick(gca,'x','mm/dd HH:MM:SS', 'keeplimits','keepticks')
% pause(0.1)
end


% save('D:\WakeID_Summer2023\Final_Wake_Data\WveWndData.mat','WveH')WveH = 100*4*sqrt(wkeEseg/Nt/dTs);
WndH = 100*4*sqrt(wndEseg/Nt/dTs);
ActE = 1024*9.8*(wkeEseg/Nt/dTs);

Nfilt = round(3600/(Nt*dTs));
filter= hamming(Nfilt); filter = filter./sum(filter);
WveHs = conv(WveH,filter,'valid');
WndHs = conv(WndH,filter,'valid');
wkeTs = wkeTseg(Nfilt/2:end-Nfilt/2);

dTwind = (MSNBt(2)-MSNBt(1))*86400;
Nfilt = round(3600/dTwind);
filter= hamming(Nfilt); filter = filter./sum(filter);
MSNBwindS = conv(MSNBwindSpeedmps,filter,'valid');
MSNBwindD = conv(MSNBwindDir,filter,'valid');
MSNBts = MSNBt(Nfilt/2:end-Nfilt/2);

figure7 = figure('units','normalized','outerposition',[0 0 1 1]);
axs8 = axes;
plot(MSNBts,abs(cosd(MSNBwindD).*MSNBwindS),'-k','linewidth',2)
axs7=axes;
p2 = plot([0 0],[-100 -101],'-k','linewidth',2);hold on
p0 = plot(wkeTs,WveHs,'-b','linewidth',2);
p1 = plot(wkeTs,WndHs,'-r','linewidth',2);

legend([p0 p1 p2],'Averaged Wake Height','Averaged Wind Wave Height','MSNB Wind Speed')
title("Wind and Wake Wave Heights Against Time",'Interpreter','latex','FontSize',FontSizeTitle)
xlabel("$t$ [mm/dd]", 'Interpreter','latex','FontSize',FontsizeLabel)
ylabel("$H$ [cm]", 'Interpreter','latex','FontSize',FontsizeLabel)
set(axs7, 'xlim',[T(1) T(end)],'tickdir','out','fontsize',20,'box','off','color','none')
datetick(axs7,'x','mm/dd','keeplimits','keepticks')

set(axs8,'color','none','xtick',[],'yaxislocation','right','tickdir','out','fontsize',20,'box','off','position',get(axs7,'position'),'xlim',[T(1) T(end)])
ylabel(axs8,'$U_\mathrm{wind}$ [m/s]','interpreter','latex','Fontsize',FontsizeLabel)
%datetick(axs7,'x','mm/dd','keepli
