%% GATHER FINAL BOAT PARAMETERS FROM PRESSURE WAKES AND CAMERA CROSSINGS
%
% Goal:
%   Build the final matched boat/wake dataset.
%
% Current version:
%   - Loads manually approved pressure wake times from Boatwaketimes.xlsx.
%   - Loops through filtered TrainingImages spectrograms.
%   - Applies imgNetfilters_v2 to recreate the final wake mask.
%   - Loads the matching raw spectrogram S matrix.
%   - Matches each approved pressure wake to boat_crossing_events.csv using
%     a +/- 15 second search window.
%   - If multiple camera boats are inside the window, picks the closest
%     boat before the pressure wake. If none are before the wake, picks the
%     closest boat after the wake.
%   - Writes a MATLAB structure to final_boat_parameters.mat.
%
% Later:
%   Fill in the pressure-variable placeholder section with slope, total
%   energy, peak frequency, duration, mask area, etc.

clear; close all; clc;

%% 1. PATHS AND SETTINGS
scriptFolder = fileparts(mfilename("fullpath"));
repoFolder = fileparts(scriptFolder);
networkFolder = fullfile(repoFolder,"Netwrks");
addpath(networkFolder);

trainingDataFolder = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData";
trainingImageFolder = fullfile(trainingDataFolder,"TrainingImages");
displayPreviewFolder = fullfile(trainingDataFolder,"displaypreview");
rawSpectrogramFolder = fullfile(trainingDataFolder,"rawspectrogramdata");
networkFile = fullfile(trainingDataFolder,"imgNet_v5.mat");

approvedWakeFile = fullfile(trainingDataFolder,"Boatwaketimes.xlsx");
boatCrossingFile = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\boat_crossing_events.csv";
finalOutputFile = fullfile(trainingDataFolder,"final_boat_parameters.mat");

cameraMatchWindowSeconds = 15;

if ~isfolder(trainingImageFolder)
    error("Training image folder not found: %s",trainingImageFolder);
end
if ~isfolder(displayPreviewFolder)
    error("displaypreview folder not found: %s",displayPreviewFolder);
end
if ~isfolder(rawSpectrogramFolder)
    error("rawspectrogramdata folder not found: %s",rawSpectrogramFolder);
end
if ~isfile(networkFile)
    error("imgNet_v5 file not found: %s",networkFile);
end
if ~isfile(approvedWakeFile)
    error("Approved wake file not found: %s",approvedWakeFile);
end
if ~isfile(boatCrossingFile)
    error("Boat crossing CSV not found: %s",boatCrossingFile);
end
if ~isfile(fullfile(networkFolder,"imgNetfilters_v2.m"))
    error("imgNetfilters_v2.m was not found in: %s",networkFolder);
end

%% 2. LOAD INPUT TABLES AND NETWORK
approvedWakeTimes = readtable(approvedWakeFile,"TextType","string");
boatCrossings = readtable(boatCrossingFile,"TextType","string");

if width(boatCrossings) < 19
    error(["boat_crossing_events.csv must have at least 19 columns. " ...
        "Required columns: time=3, direction=5, length=14, " ...
        "distance=17, speed=19."]);
end

boatCrossingTimeVariable = boatCrossings.Properties.VariableNames{3};
boatCrossingTimes = parseDatetimeColumn( ...
    boatCrossings.(boatCrossingTimeVariable));

loadedNetwork = load(networkFile,"imgNet");
if ~isfield(loadedNetwork,"imgNet")
    error("The network file does not contain a variable named imgNet.");
end
imgNet_v5 = loadedNetwork.imgNet;

requiredApprovedColumns = ["TrainingFileName","BaseName","WakeNumber", ...
    "WakeStartTime"];
missingApprovedColumns = setdiff(requiredApprovedColumns, ...
    string(approvedWakeTimes.Properties.VariableNames));

if ~isempty(missingApprovedColumns)
    error("Boatwaketimes.xlsx is missing required column(s): %s", ...
        strjoin(missingApprovedColumns,", "));
end

approvedWakeTimes.WakeStartTime = parseDatetimeColumn( ...
    approvedWakeTimes.WakeStartTime);

%% 3. LOAD AND SORT TRAINING IMAGES
trainingFiles = dir(fullfile(trainingImageFolder,"*.tif"));
if isempty(trainingFiles)
    error("No TIFF training images were found in: %s",trainingImageFolder);
end

[~,sortOrder] = sort(string({trainingFiles.name}));
trainingFiles = trainingFiles(sortOrder);

finalBoatParameters = initializeFinalBoatParameterStruct();
finalBoatIndex = 0;

fprintf("Processing %d spectrogram(s).\n",numel(trainingFiles));
fprintf("Approved wake rows available: %d\n",height(approvedWakeTimes));
fprintf("Camera crossing rows available: %d\n",height(boatCrossings));

%% 4. LOOP THROUGH EACH SPECTROGRAM
for imageNumber = 1:numel(trainingFiles)
    trainingFileName = string(trainingFiles(imageNumber).name);
    trainingFile = fullfile(trainingFiles(imageNumber).folder, ...
        trainingFiles(imageNumber).name);
    [~,baseName,~] = fileparts(trainingFileName);
    baseName = string(baseName);

    approvedRowsForSpectrogram = approvedWakeTimes( ...
        approvedWakeTimes.TrainingFileName == trainingFileName | ...
        approvedWakeTimes.BaseName == baseName,:);

    if isempty(approvedRowsForSpectrogram)
        continue
    end

    fprintf("\nSpectrogram %d/%d: %s\n",imageNumber,numel(trainingFiles), ...
        trainingFileName);
    fprintf("Approved wake(s) in this spectrogram: %d\n", ...
        height(approvedRowsForSpectrogram));

    %% 4A. LOAD FILTERED IMAGE, DISPLAYPREVIEW IMAGE, AND RAW S DATA
    trainingImage = imread(trainingFile);

    displayPreviewFile = fullfile(displayPreviewFolder,baseName+".tif");
    if isfile(displayPreviewFile)
        displayPreviewImage = imread(displayPreviewFile);
    else
        warning("Missing displaypreview TIFF for %s. Using training image.", ...
            trainingFileName);
        displayPreviewImage = trainingImage;
    end

    rawSpectrogramFile = fullfile(rawSpectrogramFolder,baseName+".mat");
    if ~isfile(rawSpectrogramFile)
        warning("Missing raw spectrogram data for %s. Skipping.", ...
            trainingFileName);
        continue
    end

    rawSpectrogramData = load(rawSpectrogramFile,"S","frequency", ...
        "spectrogramSeconds","absoluteSpectrogramTime");
    S = rawSpectrogramData.S;
    frequency = rawSpectrogramData.frequency;
    spectrogramSeconds = rawSpectrogramData.spectrogramSeconds;
    absoluteSpectrogramTime = rawSpectrogramData.absoluteSpectrogramTime;

    if ~isequal(size(S),size(trainingImage))
        warning(["Raw S matrix size does not match training image for %s. " ...
            "Skipping."],trainingFileName);
        continue
    end

    %% 4B. APPLY IMGNETFILTERS_V2 TO RECREATE FINAL MASK
    [~,~,labeledWakes,numberOfWakes,filterInfo] = ...
        imgNetfilters_v2(trainingImage,imgNet_v5,displayPreviewImage);

    %% 4C. PROCESS EACH APPROVED WAKE IN THIS SPECTROGRAM
    for approvedWakeRowNumber = 1:height(approvedRowsForSpectrogram)
        approvedWakeRow = approvedRowsForSpectrogram(approvedWakeRowNumber,:);
        pressureWakeStartTime = approvedWakeRow.WakeStartTime;
        wakeNumber = approvedWakeRow.WakeNumber;

        if isnan(wakeNumber) || wakeNumber < 1 || wakeNumber > numberOfWakes
            warning(["Approved wake number %g is not valid for %s, which " ...
                "has %d detected wake(s). Skipping this wake."], ...
                wakeNumber,trainingFileName,numberOfWakes);
            continue
        end

        singleWakeMask = labeledWakes == wakeNumber;

        if ~any(singleWakeMask,"all")
            warning("Wake %d mask is empty for %s. Skipping.", ...
                wakeNumber,trainingFileName);
            continue
        end

        %% PLACEHOLDER: CALCULATE PRESSURE-SENSOR WAKE VARIABLES LATER
        % This is where the next step will go. The inputs needed for those
        % calculations are already available here:
        %
        %   singleWakeMask
        %   S
        %   frequency
        %   spectrogramSeconds
        %   absoluteSpectrogramTime
        %
        % Future examples:
        %   pressureVariables.WakeSlope_dt_df_s_per_Hz = ...
        %   pressureVariables.TotalEnergy = ...
        %   pressureVariables.Frequency = ...
        %   pressureVariables.DurationSeconds = ...
        %   pressureVariables.MaskAreaPixels = ...
        pressureVariables = struct();
        pressureVariables.WakeSlope_dt_df_s_per_Hz = ...
            calculateWakeSlopeDtDf(singleWakeMask,spectrogramSeconds, ...
            frequency);
        pressureVariables.RMSHeight = calculateRMSHeight( ...
            singleWakeMask,S,frequency,spectrogramSeconds);
        pressureVariables.TotalEnergy = calculateTotalEnergy( ...
            singleWakeMask,S,frequency,spectrogramSeconds);
        pressureVariables.Frequency = calculateFrequencyProfile( ...
            singleWakeMask,S,frequency);
        pressureVariables.DurationSeconds = NaN;
        pressureVariables.MaskAreaPixels = nnz(singleWakeMask);

        %% MATCH PRESSURE WAKE START TIME TO CAMERA BOAT CROSSING
        [matchedBoatRow,matchedCameraTime,matchTimeOffsetSeconds, ...
            foundCameraMatch] = findMatchingBoatCrossing(boatCrossings, ...
            boatCrossingTimes,pressureWakeStartTime, ...
            cameraMatchWindowSeconds);

        if ~foundCameraMatch
            fprintf("No camera boat crossing within +/- %.1f seconds of " +...
                "%s. Skipping wake %d.\n",cameraMatchWindowSeconds, ...
                string(pressureWakeStartTime,"yyyy-MM-dd HH:mm:ss.SSS"), ...
                wakeNumber);
            continue
        end

        %% SAVE CAMERA, MATCH, FILE, AND PRESSURE DATA INTO STRUCTURE
        finalBoatIndex = finalBoatIndex+1;
        finalBoatParameters(finalBoatIndex,1).Files.TrainingImage = ...
            trainingFile;
        finalBoatParameters(finalBoatIndex,1).Files.DisplayPreview = ...
            displayPreviewFile;
        finalBoatParameters(finalBoatIndex,1).Files.RawSpectrogram = ...
            rawSpectrogramFile;
        finalBoatParameters(finalBoatIndex,1).Files.TrainingFileName = ...
            trainingFileName;
        finalBoatParameters(finalBoatIndex,1).Files.BaseName = baseName;

        finalBoatParameters(finalBoatIndex,1).Match.PressureWakeStartTime = ...
            pressureWakeStartTime;
        finalBoatParameters(finalBoatIndex,1).Match.CameraTime = ...
            matchedCameraTime;
        finalBoatParameters(finalBoatIndex,1).Match.CameraMinusPressureSeconds = ...
            matchTimeOffsetSeconds;
        finalBoatParameters(finalBoatIndex,1).Match.CameraMatchWindowSeconds = ...
            cameraMatchWindowSeconds;

        finalBoatParameters(finalBoatIndex,1).Camera.Direction = ...
            getTableValue(matchedBoatRow,5);
        finalBoatParameters(finalBoatIndex,1).Camera.BoatLengthFt = ...
            getTableValue(matchedBoatRow,14);
        finalBoatParameters(finalBoatIndex,1).Camera.DistanceFromPressureFt = ...
            getTableValue(matchedBoatRow,17);
        finalBoatParameters(finalBoatIndex,1).Camera.SpeedFtPerSec = ...
            getTableValue(matchedBoatRow,19);

        finalBoatParameters(finalBoatIndex,1).Pressure.Mask = singleWakeMask;
        finalBoatParameters(finalBoatIndex,1).Pressure.WakeSlope_dt_df_s_per_Hz = ...
            pressureVariables.WakeSlope_dt_df_s_per_Hz;
        finalBoatParameters(finalBoatIndex,1).Pressure.RMSHeight = ...
            pressureVariables.RMSHeight;
        finalBoatParameters(finalBoatIndex,1).Pressure.TotalEnergy = ...
            pressureVariables.TotalEnergy;
        finalBoatParameters(finalBoatIndex,1).Pressure.Frequency = ...
            pressureVariables.Frequency;
        finalBoatParameters(finalBoatIndex,1).Pressure.DurationSeconds = ...
            pressureVariables.DurationSeconds;
        finalBoatParameters(finalBoatIndex,1).Pressure.MaskAreaPixels = ...
            pressureVariables.MaskAreaPixels;

        fprintf("Matched wake %d at %s to camera boat at %s " +...
            "(camera-pressure = %.3f s).\n",wakeNumber, ...
            string(pressureWakeStartTime,"HH:mm:ss.SSS"), ...
            string(matchedCameraTime,"HH:mm:ss.SSS"), ...
            matchTimeOffsetSeconds);
    end
end

%% 5. SAVE FINAL OUTPUT STRUCTURE
if finalBoatIndex == 0
    warning("No matched boat/wake rows were created.");
else
    save(finalOutputFile,"finalBoatParameters");
    fprintf("\nSaved %d final matched boat parameter structure(s) to:\n%s\n", ...
        numel(finalBoatParameters),finalOutputFile);
end

%% LOCAL FUNCTIONS
function finalBoatParameters = initializeFinalBoatParameterStruct()
    finalBoatParameters = struct();
    finalBoatParameters.Files = struct( ...
        "TrainingImage","", ...
        "DisplayPreview","", ...
        "RawSpectrogram","", ...
        "TrainingFileName","", ...
        "BaseName","");
    finalBoatParameters.Match = struct( ...
        "PressureWakeStartTime",NaT, ...
        "CameraTime",NaT, ...
        "CameraMinusPressureSeconds",NaN, ...
        "CameraMatchWindowSeconds",NaN);
    finalBoatParameters.Camera = struct( ...
        "Direction",missing, ...
        "BoatLengthFt",NaN, ...
        "DistanceFromPressureFt",NaN, ...
        "SpeedFtPerSec",NaN);
    finalBoatParameters.Pressure = struct( ...
        "Mask",[], ...
        "WakeSlope_dt_df_s_per_Hz",NaN, ...
        "RMSHeight",struct("Max_cm",NaN,"Mean_cm",NaN, ...
        "ByTime_cm",[],"TimeStepSeconds",NaN), ...
        "TotalEnergy",struct("cm2s",NaN, ...
        "JPerM2",struct("FullTimeVector",[],"Max",NaN,"Mean",NaN), ...
        "JSecondsPerM2",NaN), ...
        "Frequency",struct("HzVals",[],"m2perHz",[], ...
        "MaxEnergyHz",NaN,"EnergyWeightedMeanHz",NaN, ...
        "BandwidthHz",NaN), ...
        "DurationSeconds",NaN, ...
        "MaskAreaPixels",NaN);
    finalBoatParameters(1) = [];
end

function parsedTimes = parseDatetimeColumn(timeColumn)
    if isdatetime(timeColumn)
        parsedTimes = timeColumn;
        return
    end

    if isnumeric(timeColumn)
        parsedTimes = datetime(timeColumn,"ConvertFrom","datenum");
        return
    end

    timeText = string(timeColumn);
    parsedTimes = NaT(size(timeText));

    inputFormats = [ ...
        "yyyy-MM-dd HH:mm:ss.SSS", ...
        "yyyy-MM-dd HH:mm:ss", ...
        "yyyy-MM-dd'T'HH:mm:ss.SSS", ...
        "yyyy-MM-dd'T'HH:mm:ss", ...
        "MM/dd/yyyy HH:mm:ss.SSS", ...
        "MM/dd/yyyy HH:mm:ss", ...
        "M/d/yyyy HH:mm:ss.SSS", ...
        "M/d/yyyy HH:mm:ss"];

    for formatNumber = 1:numel(inputFormats)
        stillMissing = isnat(parsedTimes) & strlength(timeText) > 0;

        if ~any(stillMissing,"all")
            break
        end

        try
            parsedTimes(stillMissing) = datetime(timeText(stillMissing), ...
                "InputFormat",inputFormats(formatNumber));
        catch
            % Try the next format.
        end
    end

    stillMissing = isnat(parsedTimes) & strlength(timeText) > 0;

    if any(stillMissing,"all")
        try
            parsedTimes(stillMissing) = datetime(timeText(stillMissing));
        catch
            warning("Some datetime values could not be parsed.");
        end
    end
end

function value = getTableValue(tableRow,columnNumber)
    rawValue = tableRow{1,columnNumber};

    if iscell(rawValue)
        value = rawValue{1};
    elseif isstring(rawValue)
        value = rawValue(1);
    elseif isnumeric(rawValue) || islogical(rawValue)
        value = rawValue(1);
    elseif iscategorical(rawValue)
        value = string(rawValue(1));
    else
        value = rawValue;
    end
end

function wakeSlope_dt_df_s_per_Hz = calculateWakeSlopeDtDf( ...
    singleWakeMask,spectrogramSeconds,frequency)

    [wakeRows,wakeColumns] = find(singleWakeMask);
    uniqueWakeRows = unique(wakeRows);

    % Need at least two frequency rows to fit a least-squares centerline.
    if numel(uniqueWakeRows) < 2
        wakeSlope_dt_df_s_per_Hz = NaN;
        return
    end

    medianColumnByRow = zeros(size(uniqueWakeRows));

    for rowNumber = 1:numel(uniqueWakeRows)
        currentRow = uniqueWakeRows(rowNumber);
        medianColumnByRow(rowNumber) = median( ...
            wakeColumns(wakeRows == currentRow));
    end

    % Ordinary least-squares line fit to the wake centerline:
    %   timeColumn = slopePixelsPerFrequencyBin*frequencyRow + intercept
    slopeFit = polyfit(double(uniqueWakeRows),medianColumnByRow,1);
    slopePixelsPerFrequencyBin = slopeFit(1);

    secondsPerColumn = mean(diff(spectrogramSeconds),"omitnan");
    hzPerFrequencyBin = mean(diff(frequency),"omitnan");

    wakeSlope_dt_df_s_per_Hz = ...
        slopePixelsPerFrequencyBin*secondsPerColumn/hzPerFrequencyBin;
end

function rmsHeight = calculateRMSHeight(singleWakeMask,S,frequency, ...
    spectrogramSeconds)

    frequencyBinWidth = mean(diff(frequency),"omitnan");
    timeStepSeconds = mean(diff(spectrogramSeconds),"omitnan");

    wakeColumns = find(any(singleWakeMask,1));

    if isempty(wakeColumns)
        rmsHeight = struct("Max_cm",NaN,"Mean_cm",NaN, ...
            "ByTime_cm",[],"TimeStepSeconds",timeStepSeconds);
        return
    end

    firstWakeColumn = min(wakeColumns);
    lastWakeColumn = max(wakeColumns);
    wakeColumnRange = firstWakeColumn:lastWakeColumn;
    rmsHeightByTime_cm = nan(1,numel(wakeColumnRange));

    for outputColumn = 1:numel(wakeColumnRange)
        timeColumn = wakeColumnRange(outputColumn);
        columnMask = singleWakeMask(:,timeColumn);

        if any(columnMask)
            variance_m2 = sum(S(columnMask,timeColumn),"omitnan")* ...
                frequencyBinWidth;
            rmsHeightByTime_cm(outputColumn) = sqrt(variance_m2)*100;
        end
    end

    rmsHeight = struct();
    rmsHeight.Max_cm = max(rmsHeightByTime_cm,[],"omitnan");
    rmsHeight.Mean_cm = mean(rmsHeightByTime_cm,"omitnan");
    rmsHeight.ByTime_cm = rmsHeightByTime_cm;
    rmsHeight.TimeStepSeconds = timeStepSeconds;
end

function totalEnergy = calculateTotalEnergy(singleWakeMask,S, ...
    frequency,spectrogramSeconds)

    frequencyBinWidth = mean(diff(frequency),"omitnan");
    timeStepSeconds = mean(diff(spectrogramSeconds),"omitnan");
    metersSquaredToCentimetersSquared = 100^2;
    seawaterDensity_kg_per_m3 = 1025;
    gravity_m_per_s2 = 9.81;

    wakeColumns = find(any(singleWakeMask,1));

    if isempty(wakeColumns)
        totalEnergy = struct("cm2s",NaN, ...
            "JPerM2",struct("FullTimeVector",[],"Max",NaN,"Mean",NaN), ...
            "JSecondsPerM2",NaN);
        return
    end

    firstWakeColumn = min(wakeColumns);
    lastWakeColumn = max(wakeColumns);
    wakeColumnRange = firstWakeColumn:lastWakeColumn;
    varianceByTime_m2 = nan(1,numel(wakeColumnRange));
    varianceByTime_cm2 = nan(1,numel(wakeColumnRange));

    for outputColumn = 1:numel(wakeColumnRange)
        timeColumn = wakeColumnRange(outputColumn);
        columnMask = singleWakeMask(:,timeColumn);

        if any(columnMask)
            variance_m2 = sum(S(columnMask,timeColumn),"omitnan")* ...
                frequencyBinWidth;
            varianceByTime_m2(outputColumn) = variance_m2;
            varianceByTime_cm2(outputColumn) = ...
                variance_m2*metersSquaredToCentimetersSquared;
        end
    end

    jPerM2ByTime = seawaterDensity_kg_per_m3*gravity_m_per_s2* ...
        varianceByTime_m2;

    totalEnergy = struct();
    totalEnergy.cm2s = sum(varianceByTime_cm2,"omitnan")* ...
        timeStepSeconds;
    totalEnergy.JPerM2 = struct();
    totalEnergy.JPerM2.FullTimeVector = jPerM2ByTime;
    totalEnergy.JPerM2.Max = max(jPerM2ByTime,[],"omitnan");
    totalEnergy.JPerM2.Mean = mean(jPerM2ByTime,"omitnan");
    totalEnergy.JSecondsPerM2 = sum(jPerM2ByTime,"omitnan")* ...
        timeStepSeconds;
end

function frequencyProfile = calculateFrequencyProfile(singleWakeMask,S, ...
    frequency)

    meanPSDByFrequency_m2PerHz = nan(size(frequency));

    for frequencyRow = 1:size(S,1)
        rowMask = singleWakeMask(frequencyRow,:);

        if any(rowMask)
            meanPSDByFrequency_m2PerHz(frequencyRow) = mean( ...
                S(frequencyRow,rowMask),"omitnan");
        end
    end

    [~,maximumIndex] = max(meanPSDByFrequency_m2PerHz,[],"omitnan");

    frequencyProfile = struct();
    frequencyProfile.HzVals = frequency;
    frequencyProfile.m2perHz = meanPSDByFrequency_m2PerHz;

    if isempty(maximumIndex) || isnan(maximumIndex)
        frequencyProfile.MaxEnergyHz = NaN;
    else
        frequencyProfile.MaxEnergyHz = frequency(maximumIndex);
    end

    validFrequencyBins = ~isnan(meanPSDByFrequency_m2PerHz);
    totalMeanPSD = sum(meanPSDByFrequency_m2PerHz(validFrequencyBins), ...
        "omitnan");

    if totalMeanPSD > 0
        frequencyProfile.EnergyWeightedMeanHz = sum( ...
            frequency(validFrequencyBins).* ...
            meanPSDByFrequency_m2PerHz(validFrequencyBins), ...
            "omitnan")/totalMeanPSD;
        frequencyProfile.BandwidthHz = sqrt(sum( ...
            meanPSDByFrequency_m2PerHz(validFrequencyBins).* ...
            (frequency(validFrequencyBins)- ...
            frequencyProfile.EnergyWeightedMeanHz).^2, ...
            "omitnan")/totalMeanPSD);
    else
        frequencyProfile.EnergyWeightedMeanHz = NaN;
        frequencyProfile.BandwidthHz = NaN;
    end
end

function [matchedBoatRow,matchedCameraTime,matchTimeOffsetSeconds, ...
    foundCameraMatch] = findMatchingBoatCrossing(boatCrossings, ...
    boatCrossingTimes,pressureWakeStartTime,cameraMatchWindowSeconds)

    timeDifferenceSeconds = seconds(boatCrossingTimes-pressureWakeStartTime);
    candidateRows = abs(timeDifferenceSeconds) <= cameraMatchWindowSeconds;
    candidateRows = candidateRows & ~isnat(boatCrossingTimes);

    if ~any(candidateRows)
        matchedBoatRow = table();
        matchedCameraTime = NaT;
        matchTimeOffsetSeconds = NaN;
        foundCameraMatch = false;
        return
    end

    candidateIndices = find(candidateRows);
    candidateOffsets = timeDifferenceSeconds(candidateIndices);

    % Priority 1: closest camera boat before the pressure wake.
    beforeCandidateIndices = candidateIndices(candidateOffsets <= 0);
    beforeCandidateOffsets = timeDifferenceSeconds(beforeCandidateIndices);

    if ~isempty(beforeCandidateIndices)
        [~,bestBeforeNumber] = min(abs(beforeCandidateOffsets));
        matchedIndex = beforeCandidateIndices(bestBeforeNumber);
    else
        % Priority 2: if none are before, use the closest camera boat after.
        afterCandidateIndices = candidateIndices(candidateOffsets > 0);
        afterCandidateOffsets = timeDifferenceSeconds(afterCandidateIndices);
        [~,bestAfterNumber] = min(afterCandidateOffsets);
        matchedIndex = afterCandidateIndices(bestAfterNumber);
    end

    matchedBoatRow = boatCrossings(matchedIndex,:);
    matchedCameraTime = boatCrossingTimes(matchedIndex);
    matchTimeOffsetSeconds = timeDifferenceSeconds(matchedIndex);
    foundCameraMatch = true;
end
