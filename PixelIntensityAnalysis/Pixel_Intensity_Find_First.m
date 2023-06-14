%This code will attemp to get a video of the web and find pixel intensity.
%WebPath ="https://stage-ams.srv.axds.co/archive/mp4/uncw/cms_dock_south/2022/06/18/cms_dock_south-2022-06-18-105207Z.mp4";
%options=weboptions; options.CertificateFilename=(''); options.RequestMethod=("Post");
%VidData = websave(CMS_PIXEl_VIDEO.html,WebPath);

VidPath ="C:\Users\lwlav\Downloads\cms_dock_south-2022-06-18-173443Z.mp4";
VidData = VideoReader(VidPath);

D1= 'June 18, 2022 13:34:40';      
D2= 'June 18, 2022 13:39:09';

t1= datenum(D1);
t2= datenum(D2);

t1= datetime(t1,'ConvertFrom','datenum');
t2= datetime(t2,'ConvertFrom','datenum');

T=t1:seconds(0.2):t2; %time in frames
frames = 1:VidData.NumFrames;  %frames to be analyzed 
% cd C:/Users/user/my_project

intensity = zeros(1,length(frames)) ;
for i = 1:length(frames)
    img = read(VidData,frames(i));
    img = rgb2gray(img(1240:1452,979:1417,:));
    intensity(i) = mean(img(:));
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

%% Plot Intensity and flag data values by putting a vertical line on 
% sharp jumps that couldn't be explained by boat wakes 

figure; hold on; grid on; box on;
set(gcf,'position',[279 288 735 295])

plot(T(1:VidData.NumFrames),intensity,'k');
title('Pixel Intensity');

xlabel('Time (UTC)'); ylabel('Pixel Intensity');

Gmean = zeros(length(frames),1);
for gg = 2:VidData.NumFrames
   img = read(VidData,frames(gg));
    A = img(1285:1591,510:1348,:);
    Gmean(gg) = mean(mean(A(:,:,2)));
    if Gmean(gg)>= Gmean(gg-1)*1.05
       xline(T(gg),'b-',{'Bad Data'})
    end
    if Gmean(gg)<= Gmean(gg-1)*.95
       xline(T(gg),'b-',{'Bad Data'})
    end
end
hold off
