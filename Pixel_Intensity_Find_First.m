%This code will attemp to get a video of the web and find pixel intensity.

VidPath ="C:\Users\lwlav\Downloads\cms_dock_south-2022-06-18-105207Z.mp4";
VidData = VideoReader(VidPath);

D1= 'June 18, 2022 06:55:06';      
D2= 'June 18, 2022 06:55:21';

t1= datenum(D1);
t2= datenum(D2);

t1= datetime(t1,'ConvertFrom','datenum');
t2= datetime(t2,'ConvertFrom','datenum');

T=t1:seconds(0.2):t2; %time in frames
frames = 925:1000; %frames to be analyzed 
% cd C:/Users/user/my_project

for i = 1:length(frames)-1
    img = read(VidData,frames(i));
    img = rgb2gray(img(630:650,1330:1350,:));
    intensity(i) = mean(img(:));
end

%% Plot Intensity
T=t1:seconds(0.2):t2

figure; hold on; grid on; box on;
set(gcf,'position',[279 288 735 295])

plot(T(1:end-1),intensity,'k','linewidth',2);
title('Pixel Intensity');

xlabel('Time (UTC)'); ylabel('Pixel Intensity');

% plot max intensity
[m, idx] = max(intensity);
scatter(t(idx),intensity(idx),'red','filled');
