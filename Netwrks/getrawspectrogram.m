%% SAVE RAW PSD SPECTROGRAM MATRICES FOR IMGNET MASK ANALYSIS
% This script uses the same pressure preprocessing, 10-minute segments, and
% spectrogram settings as gettrainingimages_v2.m.
%
% Difference from gettrainingimages_v2.m:
%   - No grayscale image conversion.
%   - No PSD cutoff.
%   - No uint8 scaling.
%   - Saves the raw spectrogram PSD matrix S for each segment.
%
% Each saved .mat file uses the same base filename as the matching
% TrainingImages/displaypreview spectrogram, so imgNet masks can be applied
% directly to S later.

clear; close all; clc;

%% 1. INPUT AND OUTPUT LOCATIONS
pressureFile = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2022\ICW_survey_20220922\rbr_pressure_CMSsouth_2022.mat";

outputFolder = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData";
rawSpectrogramFolder = fullfile(outputFolder,"rawspectrogramdata");

if ~isfolder(rawSpectrogramFolder)
    mkdir(rawSpectrogramFolder);
end

%% 2. SETTINGS COPIED FROM GETTRAININGIMAGES_V2
fs = 16;                         % Pressure sampling frequency [Hz]
segmentLengthSeconds = 600;      % One spectrogram for each complete 10-minute segment
windowLengthSeconds = 32;        % Length of each sliding Hamming window [s]
spectrogramStepSeconds = 0.5;    % Time between spectrogram columns [s]
maximumFrequency = 1.5;          % Highest frequency saved [Hz]

firstAllowedTime = hours(6)+minutes(59);  % 06:59
lastAllowedTime = hours(19)+minutes(11);  % 19:11

%% 3. LOAD AND PREPROCESS THE PRESSURE RECORD
loadedPressure = load(pressureFile,"data");

pressureDbar = loadedPressure.data(2).pres(:); % Surface-sensor pressure [dbar]
rawTime = loadedPressure.data(2).time(:);
clear loadedPressure

% Convert dbar to approximately meters of water, as in Final_WakeId and
% gettrainingimages_v2.
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
overlapSamples = windowLengthSamples-round(spectrogramStepSeconds*fs);
frequencyVector = 0:fs/windowLengthSamples:maximumFrequency;

% The final incomplete segment is skipped so every raw S matrix has the
% same dimensions as the matching TrainingImages/displaypreview spectrogram.
segmentStarts = 1:samplesPerSegment: ...
    (numel(perturbationPressure)-samplesPerSegment+1);

%% 5. CREATE AND SAVE RAW DAYTIME SPECTROGRAM MATRICES
savedRawSpectrogramFileNames = strings(0,1);
savedStartTimes = datetime.empty(0,1);
savedEndTimes = datetime.empty(0,1);
savedRows = zeros(0,1);
savedColumns = zeros(0,1);

fprintf("Checking %d complete 10-minute segments.\n",numel(segmentStarts));

for segmentNumber = 1:numel(segmentStarts)
    firstSample = segmentStarts(segmentNumber);
    sampleIndex = firstSample:(firstSample+samplesPerSegment-1);

    segmentStart = pressureTime(sampleIndex(1));
    segmentEnd = pressureTime(sampleIndex(end));

    % Only save spectrograms lying completely inside 06:59-19:11, matching
    % gettrainingimages_v2.m.
    if timeofday(segmentStart) < firstAllowedTime || ...
            timeofday(segmentEnd) > lastAllowedTime
        continue
    end

    segmentPressure = perturbationPressure(sampleIndex);

    % Create the exact same raw PSD spectrogram matrix used before any
    % display conversion into a .tif image.
    [~,frequency,spectrogramSeconds,S] = spectrogram( ...
        segmentPressure,hamming(windowLengthSamples),overlapSamples, ...
        frequencyVector,fs,"yaxis","psd");

    % Put the exact start and end date-times in the filename. This matches
    % the stem used by TrainingImages and displaypreview.
    startText = string(segmentStart,"yyyy-MM-dd_HH-mm-ss");
    endText = string(segmentEnd,"yyyy-MM-dd_HH-mm-ss");
    baseName = "spectrogram_"+startText+"_to_"+endText;
    rawSpectrogramFileName = baseName+".mat";
    rawSpectrogramFile = fullfile(rawSpectrogramFolder, ...
        rawSpectrogramFileName);

    absoluteSpectrogramTime = segmentStart+seconds(spectrogramSeconds);

    % S is the important raw PSD matrix. The other variables keep enough
    % context to align the matrix with the saved images and masks.
    save(rawSpectrogramFile,"S","frequency","spectrogramSeconds", ...
        "absoluteSpectrogramTime","segmentStart","segmentEnd", ...
        "fs","windowLengthSeconds","spectrogramStepSeconds", ...
        "frequencyVector","pressureFile");

    savedRawSpectrogramFileNames(end+1,1) = rawSpectrogramFileName;
    savedStartTimes(end+1,1) = segmentStart;
    savedEndTimes(end+1,1) = segmentEnd;
    savedRows(end+1,1) = size(S,1);
    savedColumns(end+1,1) = size(S,2);

    if mod(segmentNumber,100) == 0
        fprintf("Checked %d of %d segments; saved %d raw spectrograms.\n", ...
            segmentNumber,numel(segmentStarts), ...
            numel(savedRawSpectrogramFileNames));
    end
end

%% 6. SAVE AN INDEX CONNECTING EACH RAW SPECTROGRAM TO ITS REAL TIME
rawSpectrogramIndex = table(savedRawSpectrogramFileNames,savedStartTimes, ...
    savedEndTimes,savedRows,savedColumns,"VariableNames", ...
    ["RawSpectrogramFileName","StartTime","EndTime","Rows","Columns"]);

writetable(rawSpectrogramIndex, ...
    fullfile(rawSpectrogramFolder,"RawSpectrogramIndex.csv"));
save(fullfile(rawSpectrogramFolder,"RawSpectrogramIndex.mat"), ...
    "rawSpectrogramIndex","frequencyVector","fs","windowLengthSeconds", ...
    "spectrogramStepSeconds");

fprintf("Finished. Saved %d raw spectrogram matrices to:\n%s\n", ...
    height(rawSpectrogramIndex),rawSpectrogramFolder);
