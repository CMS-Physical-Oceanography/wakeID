%% CREATE CONSISTENT DISPLAY PREVIEWS FOR IMGNET MASK REVIEW
% This script uses the same pressure preprocessing, 10-minute segments, and
% spectrogram settings as GetTrainingimages.m.
%
% Difference from GetTrainingimages.m:
%   - No adaptive PSD cutoff.
%   - No median filtering.
%   - No contrast/gamma enhancement.
%   - No special row filtering.
%   - Uses one fixed PSD cutoff for every image:
%
%       maximumImagePSD = 0.00002
%
% These display previews are meant to give you a consistent axis/background
% image for visually reviewing masks. They are not the filtered TIFF images
% that imgNet actually sees.

clear; close all; clc;

%% 1. INPUT AND OUTPUT LOCATIONS
pressureFile = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2022\ICW_survey_20220922\rbr_pressure_CMSsouth_2022.mat";

outputFolder = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData";
displayPreviewFolder = fullfile(outputFolder,"displaypreview");

if ~isfolder(displayPreviewFolder)
    mkdir(displayPreviewFolder);
end

%% 2. SETTINGS COPIED FROM GETTRAININGIMAGES
fs = 16;                         % Pressure sampling frequency [Hz]
segmentLengthSeconds = 600;      % One image for each complete 10-minute segment
windowLengthSeconds = 32;        % Length of each sliding Hamming window [s]
spectrogramStepSeconds = 0.5;    % Time between spectrogram columns [s]
maximumFrequency = 1.5;          % Highest frequency displayed [Hz]

firstAllowedTime = hours(6)+minutes(59);  % 06:59
lastAllowedTime = hours(19)+minutes(11);  % 19:11

% Fixed display cutoff. This is intentionally constant for every preview.
maximumImagePSD = 0.00002;

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

% The final incomplete segment is skipped so every preview has the same
% time span as the matching training image.
segmentStarts = 1:samplesPerSegment: ...
    (numel(perturbationPressure)-samplesPerSegment+1);

%% 5. CREATE AND SAVE DAYTIME DISPLAY PREVIEWS
savedDisplayImageFileNames = strings(0,1);
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

    % Only save previews lying completely inside 06:59-19:11, matching
    % GetTrainingimages.m.
    if timeofday(segmentStart) < firstAllowedTime || ...
            timeofday(segmentEnd) > lastAllowedTime
        continue
    end

    segmentPressure = perturbationPressure(sampleIndex);

    % Create the same PSD spectrogram used as the imgNet input source.
    [~,frequency,spectrogramSeconds,S] = spectrogram( ...
        segmentPressure,hamming(windowLengthSamples),overlapSamples, ...
        frequencyVector,fs,"yaxis","psd");

    % Simple fixed-cutoff grayscale image for consistent preview display.
    C = 255-round((S/maximumImagePSD)*255);
    C = uint8(C);

    % Put the exact start and end date-times in the filename. This matches
    % the stem used by GetTrainingimages.m so files are easy to pair.
    startText = string(segmentStart,"yyyy-MM-dd_HH-mm-ss");
    endText = string(segmentEnd,"yyyy-MM-dd_HH-mm-ss");
    baseName = "spectrogram_"+startText+"_to_"+endText;
    displayImageFile = fullfile(displayPreviewFolder,baseName+".tif");
    previewFile = fullfile(displayPreviewFolder,baseName+".png");

    %% SAVE THE AXIS-FREE FIXED-PSD IMAGE USED FOR EXACT MASK OVERLAYS
    imwrite(C,displayImageFile,"tif");

    %% SAVE A HUMAN-READABLE COPY WITH A CONSISTENT TIME/FREQUENCY AXIS
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
    title("Display preview spectrogram: "+ ...
        string(segmentStart,"MM/dd/yyyy HH:mm:ss")+" to "+ ...
        string(segmentEnd,"HH:mm:ss")+" | fixed max PSD = "+ ...
        string(maximumImagePSD));
    datetick("x","HH:MM:SS","keeplimits");
    colorbar;

    exportgraphics(previewFigure,previewFile,"Resolution",150);

    %% RECORD WHICH TIME BELONGS TO EACH DISPLAY PREVIEW
    savedDisplayImageFileNames(end+1,1) = baseName+".tif";
    savedPreviewFileNames(end+1,1) = baseName+".png";
    savedStartTimes(end+1,1) = segmentStart;
    savedEndTimes(end+1,1) = segmentEnd;
    savedRows(end+1,1) = size(C,1);
    savedColumns(end+1,1) = size(C,2);

    if mod(segmentNumber,100) == 0
        fprintf("Checked %d of %d segments; saved %d display previews.\n", ...
            segmentNumber,numel(segmentStarts),numel(savedPreviewFileNames));
    end
end

% close(previewFigure);

%% 6. SAVE AN INDEX CONNECTING EACH DISPLAY PREVIEW TO ITS REAL TIME
displayPreviewIndex = table(savedDisplayImageFileNames,savedPreviewFileNames, ...
    savedStartTimes,savedEndTimes,savedRows,savedColumns,"VariableNames", ...
    ["DisplayImageFileName","PreviewFileName","StartTime","EndTime", ...
    "ImageRows","ImageColumns"]);

writetable(displayPreviewIndex, ...
    fullfile(displayPreviewFolder,"DisplayPreviewIndex.csv"));
save(fullfile(displayPreviewFolder,"DisplayPreviewIndex.mat"), ...
    "displayPreviewIndex","frequencyVector","fs","windowLengthSeconds", ...
    "spectrogramStepSeconds","maximumImagePSD");

fprintf("Finished. Saved %d display previews to:\n%s\n", ...
    height(displayPreviewIndex),displayPreviewFolder);
