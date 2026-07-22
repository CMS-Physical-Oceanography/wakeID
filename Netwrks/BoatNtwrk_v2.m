%% BoakNetwrk
%
% get a structure with training data
imageDir      = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData\TrainingImages";
labelDir      = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData\TrainingLabels";

% The TrainingImages folder contains every spectrogram.
% The TrainingLabels folder only contains images that have been manually
% labeled. To train correctly, use only image files that have a matching
% label file with the same name.
imageFiles = dir(fullfile(imageDir,"*.tif"));
labelFiles = dir(fullfile(labelDir,"*.tif"));

if isempty(imageFiles)
    error("No training images were found in: %s",imageDir);
end
if isempty(labelFiles)
    error("No training labels were found in: %s",labelDir);
end

imageNames = string({imageFiles.name});
labelNames = string({labelFiles.name});

matchedImageFiles = strings(0,1);
matchedLabelFiles = strings(0,1);

for labelNumber = 1:numel(labelNames)
    matchingImageNumber = find(imageNames == labelNames(labelNumber),1);

    if ~isempty(matchingImageNumber)
        matchedImageFiles(end+1,1) = fullfile(imageFiles(matchingImageNumber).folder, ...
            imageFiles(matchingImageNumber).name);
        matchedLabelFiles(end+1,1) = fullfile(labelFiles(labelNumber).folder, ...
            labelFiles(labelNumber).name);
    end
end

if isempty(matchedImageFiles)
    error("No matching image/label filename pairs were found.");
end

fprintf('Found %d labeled image(s) for training.\n',numel(matchedImageFiles));

%             background, nonFront, front
classNames = ["bckgrnd" ,"wake" ];
classLabel = [0         ,255    ];          % Labels from NtwrkTrainer_v2 use 0 for background and 255 for wake.
matchedLabelFilesCell = cellstr(matchedLabelFiles);
matchedImageFilesCell = cellstr(matchedImageFiles);

pxds       = pixelLabelDatastore(matchedLabelFilesCell,classNames,classLabel);
imds       = imageDatastore(matchedImageFilesCell);
%
%
%
% Add one background row and one background column while each pair is read.
% The source TIFF files remain unchanged. Padding gives the network even
% spatial dimensions, so its output is the same size as its training mask.
TrainingData  = transform(combine(imds,pxds),@padTrainingPair);
%
numFiles = length(imds.Files);
Ntrain   = floor(0.75*numFiles);
perm     = randperm(numFiles);
imdsValid  = subset(TrainingData,perm(Ntrain+1:numFiles));
imdsTrain  = subset(TrainingData,perm(1:Ntrain));
%
%
% get input layer size
sampleNumber = min(46,numFiles);
img = readimage(imds,perm(sampleNumber));
lbl = readimage(pxds,perm(sampleNumber));
[originalNy,originalNx,Nt] = size(img);
Ny = originalNy + 1;
Nx = originalNx + 1;
%
% determine the ratios of pixel classes
frac   = [0 0];
for ii = 1:numFiles
    lbl   = readimage(pxds,perm(ii));
    [~,class] = ismember(lbl,{'bckgrnd','wake'});
    Nbg   = numel(find(class==1));        
    Nn    = numel(find(class==2));    
    frac  = frac + [Nbg Nn]./(originalNy*originalNx);
end
frac= frac./numFiles;
unitGain = 1./frac(1);
weights = 1./frac./unitGain;
%
%
numFilters = 32;
filterSize = 3;
numClasses = 2;
layers = [
    imageInputLayer([Ny Nx Nt])
    convolution2dLayer(filterSize,numFilters,'Padding','same')
    batchNormalizationLayer
    reluLayer()
    maxPooling2dLayer(3,'Stride',2,'Padding',[1 1])
    convolution2dLayer(filterSize,numFilters,'Padding','same')
    batchNormalizationLayer    
    reluLayer()
    dropoutLayer(.4)
    transposedConv2dLayer(4,numFilters,'Stride',2,'Cropping','same');
    convolution2dLayer(1,numClasses,'Padding','same');
    batchNormalizationLayer
    softmaxLayer()
    pixelClassificationLayer('Classes',classNames,'ClassWeights',[1 25])
    ];
%
%
%
options = trainingOptions('sgdm',...
                          'InitialLearnRate',0.01,...
                          'MaxEpochs', 12,...
                          'Shuffle','every-epoch',...
                          'ValidationFrequency',3,...
                          'Verbose',1,...
                          'ValidationData',imdsValid,...                          
                          'Plots','training-progress');
%
imgNet = trainNetwork(imdsTrain,layers,options);
%%
% vifsually evaluate skill
% get input layer size
imIn = imds.Files{perm(sampleNumber)}
[~,imNum]= ismember(imIn,imds.Files);
img = readimage(imds,imNum);
lbl = readimage(pxds,imNum);
imgPadded = padarray(img,[1 1],0,'post');
lbl0Padded = semanticseg(imgPadded,imgNet);
actPadded  = activations(imgNet,imgPadded,'softmax');

% Remove the temporary bottom row and right column from the prediction.
lbl0 = lbl0Padded(1:size(img,1),1:size(img,2));
act  = actPadded(1:size(img,1),1:size(img,2),:);
A = labeloverlay(img,lbl);
B = labeloverlay(img,lbl0);
figure,
ax1 =subplot(2,1,1);
imshow(A)
title('Pixel Network')
ax2=subplot(2,1,2);
imshow(B)
hold on, contour(act(:,:,2),[0.9 0.9],'r')
title('Image Network')
linkaxes([ax1,ax2])
set(ax1,'plotboxaspectratio',[3 1 1])
%
%

%%
% plot box sizes
xm = 2.5;
ym = 2;
ag = 1;
pw = 10;
ph = 5;
ppos1 = [xm ym pw ph];
ppos2 = [xm ym+ph+ag pw ph];
ps    = [xm+pw+xm ym+ph+ag+ph+ym];
fig   = figure('units','centimeters');
pos   = get(fig,'position');
pos(3:4) = ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
ax1 = axes('units','centimeters','position',ppos1);
imagesc(B)
hold on, contour(act(:,:,2),[0.85 0.85],'r')
title('Image Network')
xlabel(ax1,'$t$ [seconds]','interpreter','latex')
set(ax1,'ticklabelinterpreter','latex','tickdir','out')

save(fullfile(fileparts(labelDir),"imgNet_v5.mat"),'imgNet')

function paddedPair = padTrainingPair(trainingPair)
% Pad the image with zero-valued pixels along the bottom and right edges.
paddedPair = trainingPair;
paddedPair{1} = padarray(trainingPair{1},[1 1],0,'post');

% Pad the corresponding categorical mask with the background class.
label = trainingPair{2};
backgroundLabel = categorical({'bckgrnd'},categories(label));
paddedLabel = repmat(backgroundLabel,size(label,1)+1,size(label,2)+1);
paddedLabel(1:size(label,1),1:size(label,2)) = label;
paddedPair{2} = paddedLabel;
end
