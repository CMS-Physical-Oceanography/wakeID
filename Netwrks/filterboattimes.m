%% FILTER USABLE BOAT WAKE TIMES FROM IMGNET-DETECTED SPECTROGRAMS
%
% Workflow:
%   1. Start with the first TrainingImages spectrogram and move forward.
%   2. Apply imgNetfilters_v2 exactly like imgNet_tester_v2.
%   3. Display the seed/grow mask on the filtered image and the final mask
%      on the fixed-PSD displaypreview image.
%   4. Click the final detected wakes that you want to keep.
%   5. Press the right arrow key to save those wake start times and advance.
%
% Output:
%   Boatwaketimes.xlsx

clear; close all; clc;

%% 1. LOCATE THE NETWORK, TRAINING IMAGES, DISPLAY IMAGES, AND OUTPUT FILE
scriptFolder = fileparts(mfilename("fullpath"));
addpath(scriptFolder);

trainingDataFolder = "C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData";
trainingImageFolder = fullfile(trainingDataFolder,"TrainingImages");
displayPreviewFolder = fullfile(trainingDataFolder,"displaypreview");
networkFile = fullfile(trainingDataFolder,"imgNet_v5.mat");
outputExcelFile = fullfile(trainingDataFolder,"Boatwaketimes.xlsx");

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

%% 2. LOAD AND SORT ALL TRAINING IMAGES
trainingFiles = dir(fullfile(trainingImageFolder,"*.tif"));
if isempty(trainingFiles)
    error("No TIFF training images were found in: %s",trainingImageFolder);
end

[~,sortOrder] = sort(string({trainingFiles.name}));
trainingFiles = trainingFiles(sortOrder);

fprintf("Reviewing %d training spectrograms.\n",numel(trainingFiles));
fprintf("Click usable final wake masks. Press right arrow to save and advance.\n");
fprintf("Press R to clear current selections. Close the figure to stop.\n");

%% 3. REVIEW EACH SPECTROGRAM IN ORDER
for imageNumber = 1:numel(trainingFiles)
    trainingFileName = string(trainingFiles(imageNumber).name);
    trainingFile = fullfile(trainingFiles(imageNumber).folder, ...
        trainingFiles(imageNumber).name);
    [~,baseName,~] = fileparts(trainingFileName);

    trainingImage = imread(trainingFile);
    originalRows = size(trainingImage,1);
    originalColumns = size(trainingImage,2);

    displayPreviewFile = fullfile(displayPreviewFolder,baseName+".tif");
    displayPreviewImage = [];

    if isfile(displayPreviewFile)
        displayPreviewImage = imread(displayPreviewFile);
    else
        % FIX: Changed the brackets [...] to a plus sign +
        fprintf("Matching axis-free displaypreview TIFF not found for %s. " + ...
            "The final mask will be shown on the filtered training image.\n", ...
            trainingFileName);
    end

    %% APPLY IMGNETFILTERS_V2
    [~,~,labeledWakes,numberOfWakes,filterInfo,~,seedGrowWakePixels] = ...
        imgNetfilters_v2(trainingImage,imgNet_v5,displayPreviewImage);
    [seedGrowLabels,numberOfSeedGrowWakes] = bwlabel(seedGrowWakePixels);

    %% RECREATE THE CORRESPONDING TIME/FREQUENCY AXES
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

    %% CREATE OVERLAY IMAGES
    if numberOfSeedGrowWakes > 0
        seedGrowColors = lines(numberOfSeedGrowWakes);
        seedGrowHighlightedImage = labeloverlay(trainingImage, ...
            seedGrowLabels,"Colormap",seedGrowColors,"Transparency",0.45);
    else
        seedGrowHighlightedImage = repmat(trainingImage,1,1,3);
    end

    if numberOfWakes > 0
        wakeColors = lines(numberOfWakes);
        filteredFinalOverlay = labeloverlay(trainingImage,labeledWakes, ...
            "Colormap",wakeColors,"Transparency",0.45);

        if ~isempty(displayPreviewImage)
            displayPreviewOverlay = labeloverlay(displayPreviewImage, ...
                labeledWakes,"Colormap",wakeColors,"Transparency",0.45);
        else
            displayPreviewOverlay = filteredFinalOverlay;
        end
    else
        filteredFinalOverlay = repmat(trainingImage,1,1,3);

        if ~isempty(displayPreviewImage)
            displayPreviewOverlay = repmat(displayPreviewImage,1,1,3);
        else
            displayPreviewOverlay = filteredFinalOverlay;
        end
    end

    %% DISPLAY BOTH MASK OVERLAYS AS VERTICAL SUBPLOTS
    selectedWakeNumbers = [];

    reviewFigure = figure("Name","Boat wake time filter: "+baseName, ...
        "Color","white","Units","normalized","OuterPosition",[0 0 1 1]);
    reviewLayout = tiledlayout(reviewFigure,2,1, ...
        "TileSpacing","compact","Padding","compact");
    title(reviewLayout, ...
        ["Boat wake review "+string(imageNumber)+" of "+ ...
        string(numel(trainingFiles)); ...
        "Click usable wakes. Right arrow = save/next. R = clear selections."; ...
        "Detected final wakes: "+string(numberOfWakes)+ ...
        " | Mode: "+filterInfo.detectionMode; ...
        trainingFileName],"Interpreter","none");

    seedAxes = nexttile(reviewLayout,1);
    imagesc(seedAxes,absoluteTimeNumber,frequency,seedGrowHighlightedImage);
    axis(seedAxes,"ij");
    xlabel(seedAxes,"Time");
    ylabel(seedAxes,"Frequency [Hz]");
    title(seedAxes, ...
        "Seed/grow mask over filtered training image");
    datetick(seedAxes,"x","HH:MM:SS","keeplimits");

    finalAxes = nexttile(reviewLayout,2);
    imagesc(finalAxes,absoluteTimeNumber,frequency,displayPreviewOverlay);
    axis(finalAxes,"ij");
    xlabel(finalAxes,"Time");
    ylabel(finalAxes,"Frequency [Hz]");
    title(finalAxes, ...
        "Final mask over fixed-PSD displaypreview. Click wakes to save.");
    datetick(finalAxes,"x","HH:MM:SS","keeplimits");
    hold(finalAxes,"on");

    fprintf("\nImage %d of %d: %s\n",imageNumber,numel(trainingFiles), ...
        trainingFileName);
    fprintf("Detected %d final wake region(s).\n",numberOfWakes);

    %% CLICK WAKES UNTIL RIGHT ARROW IS PRESSED
    advanceToNextImage = false;

    while ishandle(reviewFigure) && ~advanceToNextImage
        wasKeyPress = waitforbuttonpress;

        if ~ishandle(reviewFigure)
            fprintf("Figure closed. Stopping review.\n");
            return
        end

        if wasKeyPress
            pressedKey = get(reviewFigure,"CurrentKey");

            if strcmp(pressedKey,"rightarrow")
                advanceToNextImage = true;
            elseif strcmpi(pressedKey,"r")
                selectedWakeNumbers = [];
                cla(finalAxes);
                imagesc(finalAxes,absoluteTimeNumber,frequency, ...
                    displayPreviewOverlay);
                axis(finalAxes,"ij");
                xlabel(finalAxes,"Time");
                ylabel(finalAxes,"Frequency [Hz]");
                title(finalAxes, ...
                    "Final mask over fixed-PSD displaypreview. Click wakes to save.");
                datetick(finalAxes,"x","HH:MM:SS","keeplimits");
                hold(finalAxes,"on");
                fprintf("Cleared selections for this spectrogram.\n");
            end
        else
            currentAxes = gca;

            if currentAxes ~= seedAxes && currentAxes ~= finalAxes
                continue
            end

            clickPoint = get(currentAxes,"CurrentPoint");
            clickedTimeNumber = clickPoint(1,1);
            clickedFrequency = clickPoint(1,2);

            [~,clickedColumn] = min(abs(absoluteTimeNumber-clickedTimeNumber));
            [~,clickedRow] = min(abs(frequency-clickedFrequency));

            clickedWakeNumber = labeledWakes(clickedRow,clickedColumn);

            if clickedWakeNumber <= 0
                fprintf("Clicked background; no wake selected.\n");
                continue
            end

            if ismember(clickedWakeNumber,selectedWakeNumbers)
                fprintf("Wake %d was already selected.\n",clickedWakeNumber);
                continue
            end

            selectedWakeNumbers(end+1) = clickedWakeNumber;
            selectedWakeNumbers = unique(selectedWakeNumbers,"stable");

            highlightSelectedWake(finalAxes,labeledWakes,clickedWakeNumber, ...
                absoluteTimeNumber,frequency);

            selectedWakeStartTime = getWakeStartTime(labeledWakes, ...
                clickedWakeNumber,absoluteSpectrogramTime);
            fprintf("Selected wake %d, start time %s.\n",clickedWakeNumber, ...
                string(selectedWakeStartTime,"yyyy-MM-dd HH:mm:ss"));
        end
    end

    if ~ishandle(reviewFigure)
        fprintf("Figure closed. Stopping review.\n");
        return
    end

    %% SAVE SELECTED WAKE START TIMES BEFORE ADVANCING
    numberOfSavedSelections = saveSelectedWakeTimes(outputExcelFile, ...
        trainingFileName,baseName, ...
        selectedWakeNumbers,labeledWakes,absoluteSpectrogramTime, ...
        filterInfo.detectionMode,numberOfWakes);

    if numberOfSavedSelections > 0
        fprintf("Saved %d selected wake time(s) to:\n%s\n", ...
            numberOfSavedSelections,outputExcelFile);
    else
        % FIX: Changed the brackets [...] to a plus sign + to concatenate strings
        fprintf("No wakes selected. Skipped this spectrogram without " + ...
            "changing the spreadsheet.\n");
    end

    close(reviewFigure);
end

fprintf("\nFinished reviewing all spectrograms.\n");
fprintf("Boat wake times are saved in:\n%s\n",outputExcelFile);

%% LOCAL FUNCTIONS
function wakeStartTime = getWakeStartTime(labeledWakes,wakeNumber, ...
    absoluteSpectrogramTime)
    wakeColumns = find(any(labeledWakes == wakeNumber,1));
    wakeStartColumn = min(wakeColumns);
    wakeStartTime = absoluteSpectrogramTime(wakeStartColumn);
end

function highlightSelectedWake(targetAxes,labeledWakes,wakeNumber, ...
    absoluteTimeNumber,frequency)
    selectedMask = labeledWakes == wakeNumber;
    selectedBoundaries = bwboundaries(selectedMask);
    for boundaryNumber = 1:numel(selectedBoundaries)
        boundary = selectedBoundaries{boundaryNumber};
        boundaryRows = boundary(:,1);
        boundaryColumns = boundary(:,2);
        boundaryTimes = interp1(1:numel(absoluteTimeNumber), ...
            absoluteTimeNumber,boundaryColumns,"linear","extrap");
        boundaryFrequencies = interp1(1:numel(frequency),frequency, ...
            boundaryRows,"linear","extrap");
        plot(targetAxes,boundaryTimes,boundaryFrequencies, ...
            "k-","LineWidth",2.5);
    end
end

function numberOfSavedSelections = saveSelectedWakeTimes(outputExcelFile, ...
    trainingFileName,baseName, ...
    selectedWakeNumbers,labeledWakes,absoluteSpectrogramTime,detectionMode, ...
    numberOfDetectedWakes)

    numberOfSavedSelections = 0;
    numberOfSelectedWakes = numel(selectedWakeNumbers);

    % If nothing was selected, treat the right-arrow as a skip. Do not
    % remove or overwrite any previously saved wakes for this spectrogram.
    if numberOfSelectedWakes == 0
        return
    end

    if isfile(outputExcelFile)
        existingBoatWakeTimes = readtable(outputExcelFile, ...
            'TextType','string');
        if any(strcmp(existingBoatWakeTimes.Properties.VariableNames, ...
                "WakeStartTime"))
            existingBoatWakeTimes.WakeStartTime = string( ...
                existingBoatWakeTimes.WakeStartTime);
        end
    else
        existingBoatWakeTimes = table( ...
            strings(0,1),strings(0,1),zeros(0,1),strings(0,1), ...
            zeros(0,1),strings(0,1),zeros(0,1), ...
            'VariableNames',["TrainingFileName","BaseName", ...
            "WakeNumber","WakeStartTime","WakeStartColumn", ...
            "DetectionMode","NumberOfDetectedWakes"]);
    end
    
    % If this spectrogram was already reviewed, replace its old rows with
    % the current selections. This prevents duplicate wake times if you rerun
    % the reviewer from the first image.
    if ~isempty(existingBoatWakeTimes) && ...
            any(strcmp(existingBoatWakeTimes.Properties.VariableNames, ...
            "TrainingFileName"))
        existingBoatWakeTimes( ...
            existingBoatWakeTimes.TrainingFileName == trainingFileName,:) = [];
    end
    
    wakeStartTimes = NaT(numberOfSelectedWakes,1);
    wakeStartTimeText = strings(numberOfSelectedWakes,1);
    wakeStartColumns = zeros(numberOfSelectedWakes,1);
    
    for selectedNumber = 1:numberOfSelectedWakes
        wakeNumber = selectedWakeNumbers(selectedNumber);
        wakeColumns = find(any(labeledWakes == wakeNumber,1));
        wakeStartColumns(selectedNumber) = min(wakeColumns);
        wakeStartTimes(selectedNumber,1) = ...
            absoluteSpectrogramTime(wakeStartColumns(selectedNumber));
        wakeStartTimeText(selectedNumber,1) = string( ...
            wakeStartTimes(selectedNumber), ...
            "yyyy-MM-dd HH:mm:ss.SSS");
    end
    
    newBoatWakeTimes = table( ...
        repmat(trainingFileName,numberOfSelectedWakes,1), ...
        repmat(string(baseName),numberOfSelectedWakes,1), ...
        selectedWakeNumbers(:), ...
        wakeStartTimeText, ...
        wakeStartColumns, ...
        repmat(string(detectionMode),numberOfSelectedWakes,1), ...
        repmat(numberOfDetectedWakes,numberOfSelectedWakes,1), ...
        'VariableNames',["TrainingFileName","BaseName","WakeNumber", ...
        "WakeStartTime","WakeStartColumn","DetectionMode", ...
        "NumberOfDetectedWakes"]);
        
    updatedBoatWakeTimes = [existingBoatWakeTimes; newBoatWakeTimes];
    writetable(updatedBoatWakeTimes,outputExcelFile);
    numberOfSavedSelections = numberOfSelectedWakes;
end
