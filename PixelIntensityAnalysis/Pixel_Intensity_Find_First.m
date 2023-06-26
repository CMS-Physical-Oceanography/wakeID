%This code will attemp to get a video of the web and find pixel intensity.
%WebPath ="https://stage-ams.srv.axds.co/archive/mp4/uncw/cms_dock_south/2022/06/18/cms_dock_south-2022-06-18-105207Z.mp4";
%options=weboptions; options.CertificateFilename=(''); options.RequestMethod=("Post");
%VidData = websave(CMS_PIXEl_VIDEO.html,WebPath);
clear all;

VidPath = "C:\Users\lwlav\Downloads\cms_dock_south-2022-09-23-134837Z.mp4";
VidData = VideoReader(VidPath);

D1= 'September 23, 2022 09:48:32';     % Time frame for data to be analyzed 
D2= 'September 23, 2022 09:58:32';

t1= datenum(D1);
t2= datenum(D2);

t1= datetime(t1,'ConvertFrom','datenum');
t2= datetime(t2,'ConvertFrom','datenum');

T=t1:seconds(0.2):t2;  %time in frames
frames = 1:VidData.NumFrames ; %frames to be analyzed 

int = zeros(1,length(frames)) ;
for i = 1:length(frames)
    img = read(VidData,frames(i));
    img = rgb2gray(img(700:710,500:520,:));
    int(i) = mean(img(:));
end

%% This code finds and plots the RGB values over the course of the video 
%figure,
%for gg = 1:length(frames)
%    img = read(VidData,frames(gg));
%    A = img(1285:1591,510:1348,:);
%    Rmean(gg) = mean(mean(A(:,:,1)));
%    Gmean(gg) = mean(mean(A(:,:,2)));
%    Bmean(gg) = mean(mean(A(:,:,3)));
%    plot (T(1:VidData.NumFrames),Rmean,'r-')
%    hold on
%    plot (T(1:VidData.NumFrames),Gmean,'g-')
%    plot (T(1:VidData.NumFrames),Bmean,'b-')
%    hold off
%end
% figure; hold on; grid on; box on;
% set(gcf,'position',[279 288 735 295])
% 
% plot(T(1:VidData.NumFrames),intensity,'k');
% title('Pixel Intensity');
% 
% xlabel('Time (UTC)'); ylabel('Pixel Intensity');
% 
% Gmean = zeros(length(frames),1);
% for gg = 5:VidData.NumFrames
%    img = read(VidData,frames(gg));
%     A = img(677:724,477:497,:);
%     Gmean(gg) = mean(mean(A(:,:,2)));
%     if (Gmean(gg)<= Gmean(gg-1)*1.1) && (Gmean(gg)>= Gmean(gg-1)*1.04)
%        xline(T(gg),'b-',{'Boat Wake'})
%     elseif Gmean(gg)>=Gmean(gg-1)*1.1
%            xline(T(gg),'r-',{'Bad Data'})
%     end
%     if (Gmean(gg)>= Gmean(gg-1)*.9) && (Gmean(gg)<= Gmean(gg-1)*.96)
%        xline(T(gg),'b-',{'Boat Wake'})
%     elseif Gmean(gg)<=Gmean(gg-1)*.9
%            xline(T(gg),'r-',{'Bad Data'})
%     end
% end   

%% STFT of a pixel window over a time frame 
I0= detrend(int);
T= 25;
dTs= .5;
fs=5;
L= round(T*fs);
Olap= L-round(dTs*fs);
fr=0:fs/L:1.5;

[sp ,fp, tfp] = spectrogram(I0,hanning(L),Olap,fr,fs);

Sp = (sp.*conj(sp)/(L*fs));

Rmin = min(Sp(:));
Rmax = max(Sp(:));

C = 255-round(((Sp-Rmin)/(Rmax-Rmin)).*255);
C = uint8(C);

figure,
imagesc(C); hold on; axis off;
colormap(bone);
caxis([0 255])
hold off

%%

imwrite(C,"C:\Users\lwlav\OneDrive\Documents\MATLAB\Pixel\TrainingImages\pixel_Glitchy2boats_Sept2323.tif","tiff")
