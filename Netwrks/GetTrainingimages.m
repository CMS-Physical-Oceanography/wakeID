%% CREATE SPECTROGRAM TRAINING IMAGES FOR IMGNET
% This script repeats the pressure preprocessing and spectrogram conversion
% used by Final_WakeId.m. It does not run or train imgNet.
%
% Two matching versions of every daytime spectrogram are saved:
%   TrainingImages    - the exact axis-free uint8 image supplied to imgNet.
%   TimeAxisPreviews  - a reference copy with date/time and frequency axes.
% Both copies use the same filename stem, including their start and end
% date-times, so a training image and its preview are easy to pair.
%
% Use TrainingImages with Ntwrk_Trainer.m. Do not train imgNet with the
% preview figures because their axes, labels, and titles add unrelated pixels.

clear; close all; clc;

%% 1. INPUT AND OUTPUT LOCATIONS
pressureFile = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2022\ICW_survey_20220922\rbr_pressure_CMSsouth_2022.mat";

scriptFolder = fileparts(mfilename("fullpath"));
outputFolder = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData";
trainingImageFolder = fullfile(outputFolder,"TrainingImages");
previewFolder = fullfile(outputFolder,"TimeAxisPreviews");

if ~isfolder(trainingImageFolder)
    mkdir(trainingImageFolder);
end
if ~isfolder(previewFolder)
    mkdir(previewFolder);
end

%% 2. SETTINGS COPIED FROM FINAL_WAKEID
fs = 16;                         % Pressure sampling frequency [Hz]
segmentLengthSeconds = 600;      % One image for each complete 10-minute segment
windowLengthSeconds = 32;        % Length of each sliding Hamming window [s]
spectrogramStepSeconds = 0.5;    % Time between spectrogram columns [s]
maximumFrequency = 1.5;          % Highest frequency displayed [Hz]

firstAllowedTime = hours(6)+minutes(59);  % 06:59
lastAllowedTime = hours(19)+minutes(11);  % 19:11

% Final_WakeId maps PSD values from 0 to 0.00001 into an 8-bit image.
maximumImagePSD = 0.00001;

% Denoise each spectrogram image by subtracting the mean value of each
% frequency row. This removes persistent horizontal background bands, which
% are often caused by wind/noise, while keeping localized wake structure.
applyMeanRowDenoising = true;

%% 3. LOAD AND PREPROCESS THE PRESSURE RECORD
loadedPressure = load(pressureFile,"data");

pressureDbar = loadedPressure.data(2).pres(:); % Surface-sensor pressure [dbar]
rawTime = loadedPressure.data(2).time(:);
clear loadedPressure

% Convert dbar to approximately meters of water, as in Final_WakeId.
pressureMeters = pressureDbar*0.99;

if isdatetime(rawTime)
    pressureTime = rawTime;
else
    pressureTime = datetime(rawTime,"ConvertFrom","datenum");
end

% Remove the record mean and its slowly varying 32-second component.
pressureMean = mean(pressureMeters,"omitnan");
smoothingWindow = hanning(round(32*fs));
smoothingWindow = smoothingWindow/sum(smoothingWindow);
slowPressure = conv(pressureMeters-pressureMean,smoothingWindow,"same");
perturbationPressure = pressureMeters-slowPressure-pressureMean;

%% 4. DEFINE THE 10-MINUTE SEGMENTS AND SPECTROGRAM WINDOW
samplesPerSegment = round(segmentLengthSeconds*fs);
windowLengthSamples = round(windowLengthSeconds*fs);
overlapSamples = windowLengthSamples- ...
    round(spectrogramStepSeconds*fs);
frequencyVector = 0:fs/windowLengthSamples:maximumFrequency;

% The final incomplete segment is skipped so every saved image has the same
% dimensions. Equal dimensions are important when training the network.
segmentStarts = 1:samplesPerSegment: ...
    (numel(perturbationPressure)-samplesPerSegment+1);

%% 5. CREATE AND SAVE DAYTIME SPECTROGRAMS
savedFileNames = strings(0,1);
savedPreviewFileNames = strings(0,1);
savedStartTimes = datetime.empty(0,1);
savedEndTimes = datetime.empty(0,1);
savedRows = zeros(0,1);
savedColumns = zeros(0,1);

previewFigure = figure("Visible","off","Color","white", ...
    "Position",[100 100 1400 500]);

fprintf("Checking %d complete 10-minute segments.\n",numel(segmentStarts));

for segmentNumber = 1:numel(segmentStarts)
    firstSample = segmentStarts(segmentNumber);
    sampleIndex = firstSample:(firstSample+samplesPerSegment-1);

    segmentStart = pressureTime(sampleIndex(1));
    segmentEnd = pressureTime(sampleIndex(end));

    % Only save spectrograms lying completely inside 06:59-19:11. This
    % prevents any before-hours or after-hours pressure from entering an image.
    if timeofday(segmentStart) < firstAllowedTime || ...
            timeofday(segmentEnd) > lastAllowedTime
        continue
    end

    segmentPressure = perturbationPressure(sampleIndex);

    % Create the same PSD spectrogram used as the imgNet input source.
    [~,frequency,spectrogramSeconds,S] = spectrogram( ...
        segmentPressure,hamming(windowLengthSamples),overlapSamples, ...
        frequencyVector,fs,"yaxis","psd");

    % Convert PSD into grayscale image

    noisenumber = sum(sum(S,[],'omitnan'))
    maxS=max(sum(S,[],'omitnan'));
    avgS=mean(S,[2],'omitnan');
    threshold = 0.1;
    if noisenumber <= threshold
        maximumImagePSD = 0.00001;
        %maximumImagePSD = 0.9*maxS;
    else
        maximumImagePSD = 0.00006*sqrt(noisenumber/threshold);
        rowSmoothingBins = 5;
        maximumRowPSDMultiplier = 50;
        normalRowPSD = median(avgS,"omitnan");

        if normalRowPSD <= 0 || isnan(normalRowPSD)
            rowNoiseFactor = ones(size(avgS));
        else
            rowNoiseFactor = avgS./normalRowPSD;
        end

        rowNoiseFactor(rowNoiseFactor < 1) = 1;
        rowNoiseFactor(rowNoiseFactor > maximumRowPSDMultiplier) = maximumRowPSDMultiplier;

        maximumImagePSD = maximumImagePSD.*sqrt(rowNoiseFactor);
        %maximumImagePSD = max(0.85*max(S,[],'all'),1e-4);
    end

    % Row-adaptive PSD cutoff:
    %   maximumImagePSD above is still the baseline cutoff for the whole
    %   image. Here, each frequency row gets its own version of that cutoff.
    %   Rows with higher average PSD are treated as noisier, so their cutoff
    %   is raised. Quiet rows keep the original maximumImagePSD value.

    C = 255-round((S./maximumImagePSD)*255);

    % Subtract the mean from the top 3 frequency bins before median filtering:
    %   Row 1 is the top frequency bin in the saved image display. For each
    %   of the top 3 rows, subtract that row's own average from every pixel
    %   in that same row. This reduces persistent horizontal energy in those
    %   specific rows without changing the rest of the spectrogram.
    topFrequencyRowsToSubtract = 1:min(3,size(C,1));
    topFrequencyRowMeans = mean(C(topFrequencyRowsToSubtract,:),2,"omitnan");
    C(topFrequencyRowsToSubtract,:) = C(topFrequencyRowsToSubtract,:) + ...
        0.15*topFrequencyRowMeans;

    % Median-neighborhood cleanup before contrast enhancement:
    %   Each pixel is replaced by the median value in its local 4x4
    %   neighborhood. Random isolated specks get pulled toward their
    %   surrounding pixels before the contrast stretch makes them stronger.
    medianFilterSize = [3 3];
    C = medfilt2(C,medianFilterSize);

    % if applyMeanRowDenoising
    %     % Baby version:
    %     %   1. Look across time in one frequency row.
    %     %   2. Find the mean/background brightness for that row.
    %     %   3. Subtract that row's background from every pixel in the row.
    %     %
    %     % This removes horizontal bands that stay present across most of the
    %     % image. After subtraction, values are shifted and stretched back
    %     % into the normal 0-255 TIFF image range.
    %     C = double(C);
    %     rowBackground = mean(C,2,"omitnan");
    %     C = C-rowBackground;
    %     C = C-min(C(:));
    % 
    %     if max(C(:)) > 0
    %         C = 255*C/max(C(:));
    %     end
    % end

    % Re-enhance the wake contrast after noise filtering:
    %   1. Ignore the most extreme dark/bright pixels using percentiles.
    %   2. Stretch the remaining values back across the full 0-255 range.
    %   3. Apply mild gamma correction to make weaker wake features easier
    %      to see without letting a few extreme pixels control the image.
    C = double(C);
    lowContrastPercentile = 1;
    highContrastPercentile = 60;
    gammaValue = 2.5;

    lowContrastValue = prctile(C(:),lowContrastPercentile);
    highContrastValue = prctile(C(:),highContrastPercentile);

    if highContrastValue > lowContrastValue
        C = (C-lowContrastValue)/(highContrastValue-lowContrastValue);
        C(C < 0) = 0;
        C(C > 1) = 1;
        C = C.^gammaValue;
        C = 255*C;
    end

    % Values outside [0,255] are clipped by uint8.
    C = uint8(C);

    % Put the exact start and end date-times in both filenames. Colons are
    % replaced by hyphens because Windows does not allow colons in filenames.
    startText = string(segmentStart,"yyyy-MM-dd_HH-mm-ss");
    endText = string(segmentEnd,"yyyy-MM-dd_HH-mm-ss");
    baseName = "spectrogram_"+startText+"_to_"+endText;
    trainingFile = fullfile(trainingImageFolder,baseName+".tif");
    previewFile = fullfile(previewFolder,baseName+".png");

    %% 5A. SAVE THE EXACT IMAGE THAT WILL BE GIVEN TO IMGNET
    imwrite(C,trainingFile,"tif");

    %% 5B. SAVE A SEPARATE HUMAN-READABLE COPY WITH A TIME AXIS
    figure(previewFigure);
    clf(previewFigure);

    absoluteSpectrogramTime = segmentStart+seconds(spectrogramSeconds);
    absoluteTimeNumber = datenum(absoluteSpectrogramTime);

    imagesc(absoluteTimeNumber,frequency,C);
    % Match the matrix orientation seen by imgNet: row 1 (0 Hz) is at the
    % top, and frequency increases downward toward 1.5 Hz.
    axis ij
    colormap(gray(256));
    clim([0 255]);
    xlabel("Time");
    ylabel("Frequency [Hz]");
    title("Pressure spectrogram: "+ ...
        string(segmentStart,"MM/dd/yyyy HH:mm:ss")+" to "+ ...
        string(segmentEnd,"HH:mm:ss")+"noise number =" + noisenumber);
    datetick("x","HH:MM:SS","keeplimits");
    colorbar;

    exportgraphics(previewFigure,previewFile,"Resolution",150);

    %% 5C. RECORD WHICH TIME BELONGS TO EACH TRAINING IMAGE
    savedFileNames(end+1,1) = baseName+".tif";
    savedPreviewFileNames(end+1,1) = baseName+".png";
    savedStartTimes(end+1,1) = segmentStart;
    savedEndTimes(end+1,1) = segmentEnd;
    savedRows(end+1,1) = size(C,1);
    savedColumns(end+1,1) = size(C,2);

    if mod(segmentNumber,100) == 0
        fprintf("Checked %d of %d segments; saved %d images.\n", ...
            segmentNumber,numel(segmentStarts),numel(savedFileNames));
    end
end

%close(previewFigure);

%% 6. SAVE AN INDEX CONNECTING EACH IMAGE TO ITS REAL TIME
trainingImageIndex = table(savedFileNames,savedPreviewFileNames, ...
    savedStartTimes,savedEndTimes,savedRows,savedColumns,"VariableNames", ...
    ["TrainingFileName","PreviewFileName","StartTime","EndTime", ...
    "ImageRows","ImageColumns"]);

writetable(trainingImageIndex,fullfile(outputFolder,"TrainingImageIndex.csv"));
save(fullfile(outputFolder,"TrainingImageIndex.mat"),"trainingImageIndex", ...
    "frequencyVector","fs","windowLengthSeconds", ...
    "spectrogramStepSeconds","maximumImagePSD","applyMeanRowDenoising");

fprintf("Finished. Saved %d training images to:\n%s\n", ...
    height(trainingImageIndex),trainingImageFolder);
fprintf("Timestamped preview figures were saved to:\n%s\n",previewFolder);
