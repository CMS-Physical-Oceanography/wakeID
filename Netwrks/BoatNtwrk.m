%% BoakNetwrk
%
% get a structure with training data
imageDir      = "C:\Users\lwlav\OneDrive\Documents\MATLAB\Pressure\TrainingImages";
labelDir      = "C:\Users\lwlav\OneDrive\Documents\MATLAB\Pressure\TrainingLabels";

%             background, nonFront, front
classNames = ["bckgrnd" ,"wake" ];
classLabel = [0         ,1      ];          %%%These seem incorrect? Training labels are logical arrays...%%%
pxds       = pixelLabelDatastore(labelDir,classNames,classLabel);
imds       = imageDatastore(imageDir);
%
%
%
TrainingData  = combine(imds,pxds);
%
numFiles = length(imds.Files);
Ntrain   = floor(0.75*numFiles);
perm     = randperm(numFiles);
imdsValid  = subset(TrainingData,perm(Ntrain+1:numFiles));
imdsTrain  = subset(TrainingData,perm(1:Ntrain));
%
%
% get input layer size
img = readimage(imds,perm(46));
lbl = readimage(pxds,perm(46));
[Ny,Nx,Nt] = size(img);
%
% determine the ratios of pixel classes
frac   = [0 0];
for ii = 1:numFiles
    lbl   = readimage(pxds,perm(ii));  %%% Every single point is in the background since range is [0 1] due to logical array.
    [~,class] = ismember(lbl,{'bckgrnd','wake'});
    Nbg   = numel(find(class==1));        
    Nn    = numel(find(class==2));    
    frac  = frac + [Nbg Nn]./(Ny*Nx);
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
imIn = "C:\Users\lwlav\OneDrive\Documents\MATLAB\Pressure\TrainingImages\pressure_Choppy5boats1slow_Sept2223.tif"
[~,imNum]= ismember(imIn,imds.Files);
img = readimage(imds,imNum);
lbl = readimage(pxds,imNum);
lbl0 = semanticseg(img,imgNet);
act  = activations(imgNet,img,'softmax');
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

%save("D:\WakeID_Summer2023\imgNet_v4.mat",'imgNet')
