%% TEST IMGNET_V5 FILTER FUNCTION ON ONE RANDOM TRAINING SPECTROGRAM
%
% This script is now only a tester/display script.
%
% The actual wake-detection filter settings and mask-creation logic live in:
%
%   imgNetfilters_v2.m
%
% This script:
%   1. Picks one random training spectrogram.
%   2. Loads imgNet_v5.
%   3. Calls imgNetfilters_v2.m to create the filtered wake mask.
%   4. Displays the mask over the original spectrogram.
%   5. Displays the row-mean-subtracted probability map.

clear; close all; clc;

%% 1. LOCATE THE NETWORK, TRAINING IMAGES, AND FILTER FUNCTION
scriptFolder = fileparts(mfilename("fullpath"));
addpath(scriptFolder);

trainingDataFolder = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData";
trainingImageFolder = fullfile(trainingDataFolder,"TrainingImages");
displayPreviewFolder = fullfile(trainingDataFolder,"displaypreview");
networkFile = fullfile(trainingDataFolder,"imgNet_v5.mat");

if ~isfolder(trainingImageFolder)
    error("Training image folder not found: %s",trainingImageFolder);
end
if ~isfolder(displayPreviewFolder)
    warning("displaypreview folder not found: %s",displayPreviewFolder);
end
if ~isfile(networkFile)
    error("imgNet_v5 file not found: %s",networkFile);
end
if ~isfile(fullfile(scriptFolder,"imgNetfilters_v2.m"))
    error("imgNetfilters_v2.m was not found in: %s",scriptFolder);
end

loadedNetwork = load(networkFile,"imgNet");
if ~isfield(loadedNetwork,"imgNet")
    error("The network file does not contain a variable named imgNet.");
end
imgNet_v5 = loadedNetwork.imgNet;

% These settings are only used to recreate the display axes. The actual
% filter settings are inside imgNetfilters_v2.m.
samplingFrequency = 16;
windowLengthSeconds = 32;
spectrogramStepSeconds = 0.5;

%% 2. RANDOMLY SELECT ONE TRAINING IMAGE
trainingFiles = dir(fullfile(trainingImageFolder,"*.tif"));
if isempty(trainingFiles)
    error("No TIFF training images were found in: %s",trainingImageFolder);
end

rng("shuffle");
randomImageNumber = randi(numel(trainingFiles));
trainingFileName = string(trainingFiles(randomImageNumber).name);
trainingFile = fullfile(trainingFiles(randomImageNumber).folder, ...
    trainingFiles(randomImageNumber).name);

[~,baseName,~] = fileparts(trainingFileName);

fprintf("Random training image: %s\n",trainingFileName);

trainingImage = imread(trainingFile);
originalRows = size(trainingImage,1);
originalColumns = size(trainingImage,2);

% Load the matching fixed-PSD/less-filtered displaypreview image before
% running imgNetfilters_v2. imgNet still runs on trainingImage, but
% imgNetfilters_v2 uses this displaypreview image for TIFF-value
% post-filtering when it is available.
displayPreviewFile = fullfile(displayPreviewFolder,baseName+".tif");
displayPreviewImage = [];

if isfile(displayPreviewFile)
    displayPreviewImage = imread(displayPreviewFile);
else
    fprintf(["Matching axis-free displaypreview TIFF not found. " ...
        "TIFF-value post-filtering will fall back to the filtered " ...
        "training image.\n"]);
end

%% 3. APPLY IMGNETFILTERS
[~,~,labeledWakes,numberOfWakes,filterInfo,filteredWakeProbability, ...
    seedGrowWakePixels] = ...
    imgNetfilters_v2(trainingImage,imgNet_v5,displayPreviewImage);

% This is the intermediate mask immediately after imgNet seed/grow. It is
% before the later TIFF-value cutoff, bridge removal, row/frequency-span
% checks, and final area filtering.
[seedGrowLabels,numberOfSeedGrowWakes] = bwlabel(seedGrowWakePixels);

fprintf("\n============================================================\n");
fprintf("TOTAL ORIGINAL TIFF IMAGE-VALUE SUM: %.0f\n", ...
    filterInfo.totalImageValueSum);
fprintf("DETECTION MODE SELECTED: %s\n",filterInfo.detectionMode);
fprintf("RAW WAKE PROBABILITY THRESHOLD: %.2f\n", ...
    filterInfo.rawWakeProbabilityThreshold);
fprintf("ROW-MEAN-SUBTRACTED PROBABILITY THRESHOLD: %.2f\n", ...
    filterInfo.filteredProbabilityThreshold);
fprintf("GROWTH PROBABILITY THRESHOLD: %.2f\n", ...
    filterInfo.growthProbabilityThreshold);
fprintf("MAXIMUM ALLOWED TIFF VALUE IN WAKE: %d\n", ...
    filterInfo.maximumWeakTiffValue);
fprintf("TIFF-VALUE FILTERING IMAGE: %s\n", ...
    filterInfo.tiffFilteringImageSource);
fprintf("SECOND DISPLAY IMAGE: %s\n",filterInfo.postFilterImageSource);
fprintf("BRIDGE REMOVAL RADIUS: %d pixels\n", ...
    filterInfo.bridgeRemovalRadiusPixels);
fprintf("MINIMUM REACHED FREQUENCY: %.2f Hz\n", ...
    filterInfo.minimumReachedFrequency);
fprintf("MINIMUM FREQUENCY BIN RANGE: %d bins\n", ...
    filterInfo.minimumFrequencyBinRange);
fprintf("MINIMUM WAKE AREA: %d pixels\n",filterInfo.minimumWakeArea);
fprintf("FILTERED TRAINING IMAGE NON-WHITE PIXELS: %d\n", ...
    filterInfo.trainingNonWhitePixelCount);
fprintf("POST-FILTER IMAGE NON-WHITE PIXELS: %d\n", ...
    filterInfo.postFilterNonWhitePixelCount);
fprintf("SEED/GROW WAKE REGIONS BEFORE FINAL FILTERS: %d\n", ...
    numberOfSeedGrowWakes);
fprintf("DETECTED WAKE REGIONS: %d\n",numberOfWakes);
fprintf("============================================================\n\n");

%% 4. CREATE THE SEED/GROW OVERLAY ON THE FILTERED TRAINING IMAGE
% This top display shows the first mask stage: imgNet seed/grow on the
% filtered image that the network actually sees.
if numberOfSeedGrowWakes > 0
    seedGrowColors = lines(numberOfSeedGrowWakes);
    seedGrowHighlightedImage = labeloverlay(trainingImage,seedGrowLabels, ...
        "Colormap",seedGrowColors,"Transparency",0.45);
else
    seedGrowHighlightedImage = repmat(trainingImage,1,1,3);
end

%% 5. CREATE THE FINAL OVERLAY COLORS
if numberOfWakes > 0
    wakeColors = lines(numberOfWakes);
    highlightedImage = labeloverlay(trainingImage,labeledWakes, ...
        "Colormap",wakeColors,"Transparency",0.45);
else
    % Convert the grayscale image to RGB so it displays the same way as an
    % image produced by labeloverlay, even when no wakes were detected.
    highlightedImage = repmat(trainingImage,1,1,3);
end

%% 6. RECREATE THE CORRESPONDING TIME/FREQUENCY AXES
% The filename begins with the exact segment start date and time, for
% example: spectrogram_2022-09-22_07-00-00_to_2022-09-22_07-09-59.tif
nameParts = regexp(baseName, ...
    '^spectrogram_(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})_to_', ...
    'tokens','once');
if isempty(nameParts)
    error("Could not read the segment start time from filename: %s", ...
        trainingFileName);
end

segmentStart = datetime(nameParts{1}, ...
    "InputFormat","yyyy-MM-dd_HH-mm-ss");
spectrogramSeconds = windowLengthSeconds/2 + ...
    (0:originalColumns-1)*spectrogramStepSeconds;
absoluteSpectrogramTime = segmentStart + seconds(spectrogramSeconds);
absoluteTimeNumber = datenum(absoluteSpectrogramTime);

frequencyStep = samplingFrequency/(windowLengthSeconds*samplingFrequency);
frequency = (0:originalRows-1)*frequencyStep;

%% 7. LOAD THE AXIS-FREE DISPLAYPREVIEW IMAGE
% The displaypreview PNG has labels/margins baked into the raster image, so
% resizing the mask over that full PNG shifts the overlay. Instead, this
% uses the axis-free displaypreview TIFF and recreates the same time/frequency
% axes here. That keeps the mask and spectrogram matrix aligned.
displayPreviewOverlay = [];

if ~isempty(displayPreviewImage)
    if numberOfWakes > 0
        displayPreviewOverlay = labeloverlay(displayPreviewImage, ...
            labeledWakes,"Colormap",wakeColors,"Transparency",0.45);
    else
        displayPreviewOverlay = repmat(displayPreviewImage,1,1,3);
    end
end

%% 8. DISPLAY BOTH MASK OVERLAYS AS VERTICAL SUBPLOTS
overlayFigure = figure("Name","imgNetfilters_v2 overlay comparison: "+baseName, ...
    "Color","white","Units","normalized","OuterPosition",[0 0 1 1]);
overlayLayout = tiledlayout(overlayFigure,2,1, ...
    "TileSpacing","compact","Padding","compact");
title(overlayLayout, ...
    ["imgNetfilters_v2 detections: "+string(numberOfWakes)+" connected region(s)"; ...
    "Seed/grow regions before final filters: "+string(numberOfSeedGrowWakes); ...
    "Mode: "+filterInfo.detectionMode+ ...
    " | Filtered threshold: "+string(filterInfo.filteredProbabilityThreshold)+ ...
    " | Growth threshold: "+string(filterInfo.growthProbabilityThreshold); ...
    trainingFileName],"Interpreter","none");

trainingAxes = nexttile(overlayLayout,1);
imagesc(trainingAxes,absoluteTimeNumber,frequency,seedGrowHighlightedImage);
axis(trainingAxes,"ij");
xlabel(trainingAxes,"Time");
ylabel(trainingAxes,"Frequency [Hz]");
title(trainingAxes, ...
    "Seed/grow mask over filtered training image | non-white pixels: "+ ...
    string(filterInfo.trainingNonWhitePixelCount));
datetick(trainingAxes,"x","HH:MM:SS","keeplimits");

displayAxes = nexttile(overlayLayout,2);
if ~isempty(displayPreviewOverlay)
    imagesc(displayAxes,absoluteTimeNumber,frequency,displayPreviewOverlay);
    axis(displayAxes,"ij");
    title(displayAxes, ...
        "Fully filtered final mask over fixed-PSD displaypreview | non-white pixels: "+ ...
        string(filterInfo.postFilterNonWhitePixelCount));
else
    imagesc(displayAxes,absoluteTimeNumber,frequency,highlightedImage);
    axis(displayAxes,"ij");
    title(displayAxes, ...
        "displaypreview TIFF missing; final mask shown on training image | non-white pixels: "+ ...
        string(filterInfo.postFilterNonWhitePixelCount));
end
xlabel(displayAxes,"Time");
ylabel(displayAxes,"Frequency [Hz]");
datetick(displayAxes,"x","HH:MM:SS","keeplimits");

%% 9. DISPLAY THE MEAN-SUBTRACTED PROBABILITY MAP
figure("Name","mean-subtracted imgNet probability: "+baseName, ...
    "Color","white","Units","normalized","OuterPosition",[0 0 1 1]);
imshow(filteredWakeProbability,[]);
colormap(gca,"parula");
colorbar;
title(["Row-mean-subtracted imgNet wake probability"; ...
    "Created inside imgNetfilters_v2.m"; ...
    trainingFileName],"Interpreter","none");

fprintf("Displayed imgNetfilters_v2 output for: %s\n",trainingFileName);
