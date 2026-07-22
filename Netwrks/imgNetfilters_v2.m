function [tifMask,wakePixels,labeledWakes,numberOfWakes,filterInfo, ...
    filteredWakeProbability,filteredOnlyWakePixels] = imgNetfilters_v2( ...
    trainingImage,imgNet_v5,postFilterImage)
%IMGNETFILTERS_V2 Apply imgNet_v5 with mode-dependent filtering.
%
% trainingImage is the filtered spectrogram image that imgNet sees.
% postFilterImage is optional. When supplied, it should be the matching
% fixed-PSD/less-filtered axis-free displaypreview TIFF used by tester
% scripts for display/counting. The mask itself is seeded, grown, and
% filtered using trainingImage.
%
% Strict mode:
%   Uses the full imgNetfilters-style post-processing pipeline.
%
% Nonstrict mode:
%   Starts wake seeds where row-mean-subtracted probability is high enough,
%   then grows those seeds through connected weaker probability pixels.
%   Then it removes small detections and applies a light radius-based bridge
%   removal. No TIFF cutoff or frequency-bin range filtering is applied in
%   nonstrict mode.
%
% Main output:
%   tifMask  - uint8 mask image. Background is 0 and wake is 255.

%% 1. FILTER SETTINGS
totalImageValueSum = sum(double(trainingImage(:)));
trainingNonWhitePixelCount = nnz(trainingImage ~= 255);
tiffValueModeThreshold = 12000000;

strictModeProbabilityThreshold = 0.97;
nonstrictModeProbabilityThreshold = 0.95;
nonstrictFilteredProbabilityThreshold = 0.2;
nonstrictGrowthProbabilityThreshold = 0.80;
nonstrictMinimumWakeArea = 300;
nonstrictBridgeRemovalRadiusPixels = 2;
nonstrictMinimumReachedRow = 17;
nonstrictMinimumFrequencyBinRange = 5;
nonstrictMaximumWeakTiffValue = 200;

% Strict-mode-only settings.
minimumReachedFrequency = 0.5;
strictMinimumWakeArea = 400;
minimumFrequencyBinRange = 5;
strictFilteredProbabilityThreshold = 0.25;
strictGrowthProbabilityThreshold = 0.01;
strictBridgeRemovalRadiusPixels = 3;
strictMaximumWeakTiffValue = 180;

% These match the spectrogram settings used when the training images were
% created. They are needed to convert image rows into frequency values.
samplingFrequency = 16;
windowLengthSeconds = 32;

%% 2. ACCEPT EITHER A NETWORK VARIABLE OR A .MAT FILE PATH
if isstring(imgNet_v5) || ischar(imgNet_v5)
    loadedNetwork = load(imgNet_v5,"imgNet");
    if ~isfield(loadedNetwork,"imgNet")
        error("The network file does not contain a variable named imgNet.");
    end
    imgNet = loadedNetwork.imgNet;
else
    imgNet = imgNet_v5;
end

%% 2B. ACCEPT AN OPTIONAL DISPLAY IMAGE FOR TESTER/DEBUGGING
if nargin < 3 || isempty(postFilterImage)
    % Fallback keeps older calls working. The optional image is not used for
    % filtering anymore; all filtering is done on trainingImage.
    postFilterImage = trainingImage;
    postFilterImageSource = "filtered training image display";
else
    postFilterImageSource = "fixed-PSD displaypreview image display";
end

if ~isequal(size(postFilterImage,1),size(trainingImage,1)) || ...
        ~isequal(size(postFilterImage,2),size(trainingImage,2))
    error(["postFilterImage must have the same row/column size as " ...
        "trainingImage."]);
end

postFilterNonWhitePixelCount = nnz(postFilterImage ~= 255);

%% 3. CHOOSE STRICT OR NONSTRICT MODE
if totalImageValueSum < tiffValueModeThreshold
    detectionMode = "strict";
    rawWakeProbabilityThreshold = strictModeProbabilityThreshold;
else
    detectionMode = "nonstrict";
    rawWakeProbabilityThreshold = nonstrictModeProbabilityThreshold;
end

%% 4. APPLY IMGNET_V5
originalRows = size(trainingImage,1);
originalColumns = size(trainingImage,2);

% BoatNtwrk_v2 trained the network using one extra row and column so the
% image size would work cleanly. Do the same temporary padding here, then
% remove the extra row/column from the output.
paddedImage = padarray(trainingImage,[1 1],0,"post");
activation = activations(imgNet,paddedImage,"softmax");
activation = activation(1:originalRows,1:originalColumns,:);

% Channel 2 is the wake class probability.
wakeProbability = activation(:,:,2);

%% 5. ALWAYS CALCULATE THE ROW-MEAN-SUBTRACTED PROBABILITY
% In strict mode, this is used for seed/grow. In nonstrict mode, this is
% used as a second required threshold so tiny row means do not create false
% detections.
rowMeanWakeProbability = mean(wakeProbability,2,"omitnan");
filteredWakeProbability = wakeProbability-rowMeanWakeProbability;
filteredWakeProbability(filteredWakeProbability < 0) = 0;

%% 6. CREATE THE MASK
if detectionMode == "nonstrict"
    %% NONSTRICT MODE: ROW-MEAN-SUBTRACTED SEED/GROW
    % A pixel can start a wake seed only if it passes BOTH tests:
    %
    %   1. Raw imgNet wake probability is high enough.
    %   2. Row-mean-subtracted probability is high enough.
    %
    % Once a seed exists, the wake can grow through connected weaker pixels
    % down to nonstrictGrowthProbabilityThreshold.
    filteredProbabilityThreshold = nonstrictFilteredProbabilityThreshold;
    growthProbabilityThreshold = nonstrictGrowthProbabilityThreshold;

    strongWakeSeeds = wakeProbability >= rawWakeProbabilityThreshold & ...
        filteredWakeProbability >= filteredProbabilityThreshold;
    looseWakeGrowthArea = wakeProbability >= growthProbabilityThreshold;
    wakePixels = imreconstruct(strongWakeSeeds,looseWakeGrowthArea);
    filteredOnlyWakePixels = wakePixels;

    minimumWakeArea = nonstrictMinimumWakeArea;
    maximumWeakTiffValue = nonstrictMaximumWeakTiffValue;
    bridgeRemovalRadiusPixels = nonstrictBridgeRemovalRadiusPixels;

    % Nonstrict cleanup order:
    %   1. Apply the maximum weak TIFF value cutoff.
    %   2. Shrink/regrow with a radius filter to remove narrow bridges and
    %      random branches.
    %   3. Keep only blobs that reach the minimum requested image row.
    %   4. Keep only blobs that span enough frequency bins.
    %   5. Remove connected detections below the minimum area.
    wakePixels(trainingImage > maximumWeakTiffValue) = false;
    bridgeRemovalShape = strel("disk",bridgeRemovalRadiusPixels);
    wakePixels = imopen(wakePixels,bridgeRemovalShape);

    % Each remaining nonstrict blob must reach at least the requested image
    % row. Row numbers increase downward in the spectrogram image, so this
    % acts like a minimum downward frequency-bin reach requirement.
    [candidateLabels,numberOfCandidates] = bwlabel(wakePixels);
    rowQualifiedWakePixels = false(size(wakePixels));

    for candidateNumber = 1:numberOfCandidates
        candidatePixels = candidateLabels == candidateNumber;
        rowsUsedByCandidate = find(any(candidatePixels,2));

        if any(rowsUsedByCandidate >= nonstrictMinimumReachedRow)
            rowQualifiedWakePixels = rowQualifiedWakePixels | candidatePixels;
        end
    end

    wakePixels = rowQualifiedWakePixels;

    % Each remaining nonstrict blob must span at least the requested number
    % of frequency bins/rows.
    [candidateLabels,numberOfCandidates] = bwlabel(wakePixels);
    frequencyRangeQualifiedWakePixels = false(size(wakePixels));

    for candidateNumber = 1:numberOfCandidates
        candidatePixels = candidateLabels == candidateNumber;
        candidateRows = find(any(candidatePixels,2));

        if isempty(candidateRows)
            continue
        end

        candidateFrequencyBinRange = max(candidateRows)-min(candidateRows)+1;

        if candidateFrequencyBinRange >= nonstrictMinimumFrequencyBinRange
            frequencyRangeQualifiedWakePixels = frequencyRangeQualifiedWakePixels | ...
                candidatePixels;
        end
    end

    wakePixels = frequencyRangeQualifiedWakePixels;

    % Apply minimum area at the end, after all other nonstrict cleanup.
    wakePixels = bwareaopen(wakePixels,minimumWakeArea);
else
    %% STRICT MODE: FULL FILTER PIPELINE
    filteredProbabilityThreshold = strictFilteredProbabilityThreshold;
    growthProbabilityThreshold = strictGrowthProbabilityThreshold;
    minimumWakeArea = strictMinimumWakeArea;
    maximumWeakTiffValue = strictMaximumWeakTiffValue;
    bridgeRemovalRadiusPixels = strictBridgeRemovalRadiusPixels;

    % Strong pixels start each wake. To become a seed, a pixel must pass
    % both the row-mean-subtracted threshold and the raw imgNet probability
    % threshold.
    strongWakeSeeds = filteredWakeProbability >= filteredProbabilityThreshold & ...
        wakeProbability >= rawWakeProbabilityThreshold;
    looseWakeGrowthArea = filteredWakeProbability >= growthProbabilityThreshold;
    wakePixels = imreconstruct(strongWakeSeeds,looseWakeGrowthArea);
    filteredOnlyWakePixels = wakePixels;

    % Remove weak pixels using the same filtered spectrogram that imgNet
    % analyzed. The displaypreview image is not used for filtering.
    wakePixels(trainingImage > maximumWeakTiffValue) = false;

    % Remove narrow bridges/random branches.
    bridgeRemovalShape = strel("disk",bridgeRemovalRadiusPixels);
    wakePixels = imopen(wakePixels,bridgeRemovalShape);

    % Require each blob to reach at least 0.5 Hz.
    frequencyStep = samplingFrequency/(windowLengthSeconds*samplingFrequency);
    frequencyByRow = (0:originalRows-1)*frequencyStep;

    [candidateLabels,numberOfCandidates] = bwlabel(wakePixels);
    frequencyQualifiedWakePixels = false(size(wakePixels));

    for candidateNumber = 1:numberOfCandidates
        candidatePixels = candidateLabels == candidateNumber;
        rowsUsedByCandidate = any(candidatePixels,2);

        if any(frequencyByRow(rowsUsedByCandidate) >= minimumReachedFrequency)
            frequencyQualifiedWakePixels = frequencyQualifiedWakePixels | ...
                candidatePixels;
        end
    end

    wakePixels = frequencyQualifiedWakePixels;

    % Apply the same filtered-spectrogram TIFF cutoff again after shape
    % cleanup.
    wakePixels(trainingImage > maximumWeakTiffValue) = false;

    % Require a minimum vertical frequency-bin range.
    [candidateLabels,numberOfCandidates] = bwlabel(wakePixels);
    frequencyRangeQualifiedWakePixels = false(size(wakePixels));

    for candidateNumber = 1:numberOfCandidates
        candidatePixels = candidateLabels == candidateNumber;
        candidateRows = find(any(candidatePixels,2));

        if isempty(candidateRows)
            continue
        end

        candidateFrequencyBinRange = max(candidateRows)-min(candidateRows)+1;

        if candidateFrequencyBinRange >= minimumFrequencyBinRange
            frequencyRangeQualifiedWakePixels = frequencyRangeQualifiedWakePixels | ...
                candidatePixels;
        end
    end

    wakePixels = frequencyRangeQualifiedWakePixels;

    % Remove small detections after all other strict filters.
    wakePixels = bwareaopen(wakePixels,minimumWakeArea);
end

%% 7. CREATE LABELS AND TIFF-STYLE MASK
[labeledWakes,numberOfWakes] = bwlabel(wakePixels);
tifMask = uint8(wakePixels)*255;

%% 8. RETURN SETTINGS USED FOR DEBUGGING / RECORD KEEPING
filterInfo = struct();
filterInfo.totalImageValueSum = totalImageValueSum;
filterInfo.trainingNonWhitePixelCount = trainingNonWhitePixelCount;
filterInfo.postFilterNonWhitePixelCount = postFilterNonWhitePixelCount;
filterInfo.detectionMode = detectionMode;
filterInfo.tiffValueModeThreshold = tiffValueModeThreshold;
filterInfo.rawWakeProbabilityThreshold = rawWakeProbabilityThreshold;
filterInfo.filteredProbabilityThreshold = filteredProbabilityThreshold;
filterInfo.growthProbabilityThreshold = growthProbabilityThreshold;
filterInfo.maximumWeakTiffValue = maximumWeakTiffValue;
filterInfo.postFilterImageSource = postFilterImageSource;
filterInfo.tiffFilteringImageSource = "filtered training image";
filterInfo.minimumReachedFrequency = minimumReachedFrequency;
filterInfo.minimumWakeArea = minimumWakeArea;
filterInfo.minimumFrequencyBinRange = minimumFrequencyBinRange;
filterInfo.bridgeRemovalRadiusPixels = bridgeRemovalRadiusPixels;
filterInfo.nonstrictMinimumReachedRow = nonstrictMinimumReachedRow;
filterInfo.nonstrictMinimumFrequencyBinRange = ...
    nonstrictMinimumFrequencyBinRange;
filterInfo.numberOfWakes = numberOfWakes;
filterInfo.nonstrictUsesSeedGrow = detectionMode == "nonstrict";

end
