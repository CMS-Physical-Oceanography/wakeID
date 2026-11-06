clear all
close all

% netFile = '/Users/derekgrimes/OneDriveUNCW/BoatWakes/Summer2023/imgNet_v4.mat';
netFile = './imgNet_v4.mat';
load(netFile)

rootDir = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/Gannon, Katherine M''s files - Raw data mat files/';

imgDir  = fullfile(rootDir,'OutputFolder','TrainingImages');

files = dir([imgDir,'*.tif']);

figure
for ii=1:length(files)
fin = fullfile(files(ii).folder,files(ii).name);
C = imread(fin);
[Nf,Nt] = size(C);
% $$$ imshow(C)
% $$$ title(num2str(ii))
% $$$ pause(0.5)
% $$$ end

%% get the boat-wake probability using neural network
act = activations(imgNet,C,'softmax');
%% set wake threshold to 75% probability
Bw = act(1:Nf,1:Nt,2)>=.75;
%% isolate connected regions (individual wakes), and trim small regions
CleanC= bwareaopen(Bw,100); % Remove connected regions under 100 pixels
[labim,number] = bwlabel(CleanC); % Label separate detected regions
labimS = max(labim,[],1);   % Collapse labels to one value per time column

A = labeloverlay(C,labim);

imagesc(A)
pause(0.25)
xlabel('Time')
       
end

