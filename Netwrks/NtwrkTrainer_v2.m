%% MANUALLY CREATE WAKE LABELS FROM THE GENERATED TRAINING IMAGES
% This is an updated version of Ntwrk_Trainer.m for the images produced by
% imgnettrainer.m.
%
% Workflow:
%   1. Randomly choose 300 images that do not already have saved labels.
%   2. Draw wake regions directly on one enlarged training image.
%   3. Review the combined mask as a transparent red overlay.
%   4. Save, redraw, or skip that image.
%
% The original training image is never moved or changed. The saved label has
% the same filename as its training image, which keeps image/label pairing easy.


clear; close all; clc;


%% 1. LOCATE THE IMAGES CREATED BY IMGNETTRAINER
scriptFolder = fileparts(mfilename("fullpath"));
trainingDataFolder = fullfile(scriptFolder,"SpectrogramTrainingData");
trainingImageFolder = fullfile(trainingDataFolder,"TrainingImages");
trainingLabelFolder = fullfile(trainingDataFolder,"TrainingLabels");


if ~isfolder(trainingImageFolder)
    error("Training image folder not found: %s",trainingImageFolder);
end
if ~isfolder(trainingLabelFolder)
    mkdir(trainingLabelFolder);
end


%% 2. LOAD AND SORT THE AVAILABLE TRAINING IMAGES
trainingImages = imageDatastore(trainingImageFolder, ...
    "FileExtensions",".tif");


% Sorting makes the list chronological because the filenames begin with dates.
[~,sortOrder] = sort(string(trainingImages.Files));
trainingImages.Files = trainingImages.Files(sortOrder);


if isempty(trainingImages.Files)
    error("No TIFF training images were found in: %s",trainingImageFolder);
end


%% 3. RANDOMLY CHOOSE 300 IMAGES THAT HAVE NOT ALREADY BEEN LABELED
numberOfRandomImages = 300;
alreadyLabeled = false(numel(trainingImages.Files),1);


for imageNumber = 1:numel(trainingImages.Files)
    [~,baseName,~] = fileparts(trainingImages.Files{imageNumber});
    expectedLabel = fullfile(trainingLabelFolder,baseName+".tif");
    alreadyLabeled(imageNumber) = isfile(expectedLabel);
end


availableImageNumbers = find(~alreadyLabeled);
if isempty(availableImageNumbers)
    disp("Every training image already has a saved label.");
    return
end


numberToChoose = min(numberOfRandomImages,numel(availableImageNumbers));
if numberToChoose < numberOfRandomImages
    warning(["Only %d unlabeled images remain, so all of them will be " ...
        "used instead of 300."],numberToChoose);
end


% Seed from the current time so each new run produces a different sample.
rng("shuffle");
randomOrder = randperm(numel(availableImageNumbers),numberToChoose);
selectedImageNumbers = availableImageNumbers(randomOrder);


fprintf(['Randomly selected %d unlabeled image(s). Existing masks were ' ...
    'excluded.\n'],numberToChoose);


%% 4. DISPLAY EACH IMAGE AND DRAW ITS WAKE REGION(S)
for selectionNumber = 1:numel(selectedImageNumbers)
    imageNumber = selectedImageNumbers(selectionNumber);
    trainingImageFile = string(trainingImages.Files{imageNumber});
    [~,baseName,~] = fileparts(trainingImageFile);


    labelFile = fullfile(trainingLabelFolder,baseName+".tif");


    trainingImage = imread(trainingImageFile);


    fprintf('\nImage %d of %d: %s\n',selectionNumber, ...
        numel(selectedImageNumbers),baseName);


    imageFinished = false;
    while ~imageFinished
        %% SHOW ONE LARGE IMAGE; STRETCH ITS SHORT FREQUENCY AXIS FOR DRAWING
        labelFigure = figure("Name","Wake label: "+baseName, ...
            "Color","white","Units","normalized", ...
            "OuterPosition",[0 0 1 1]);
        drawingAxes = axes(labelFigure,"Position",[0.04 0.10 0.94 0.82]);


        imageHandle = imagesc(drawingAxes,trainingImage);
        colormap(drawingAxes,gray(256));
        clim(drawingAxes,[0 255]);
        set(drawingAxes,"YDir","reverse"); % Row 1/0 Hz stays at the top
        xlim(drawingAxes,[0.5 size(trainingImage,2)+0.5]);
        ylim(drawingAxes,[0.5 size(trainingImage,1)+0.5]);
        xlabel(drawingAxes,"Spectrogram time column");
        ylabel(drawingAxes,"Frequency row (0 Hz is at the top)");
        title(drawingAxes,["Draw directly over the wake(s)"; baseName], ...
            "Interpreter","none");
        hold(drawingAxes,"on");


        % The red layer is exactly the same number of rows and columns as
        % the training image. Its transparency is updated after every ROI.
        redImage = zeros(size(trainingImage,1),size(trainingImage,2),3);
        redImage(:,:,1) = 1;
        maskOverlay = image(drawingAxes,redImage, ...
            "XData",[1 size(trainingImage,2)], ...
            "YData",[1 size(trainingImage,1)]);
        maskOverlay.AlphaData = zeros(size(trainingImage,1), ...
            size(trainingImage,2));
        maskOverlay.HitTest = "off";
        maskOverlay.PickableParts = "none";


        % Keep the original image interactive beneath the transparent layer.
        imageHandle.HitTest = "off";


        regionResponse = strtrim(input([ ...
            'How many separate wake regions will you draw? ' ...
            'Enter 0 to save a no-wake mask, or K to skip: '],'s'));


        if strcmpi(regionResponse,"k")
            close(labelFigure);
            fprintf('Skipped %s without saving.\n',baseName);
            imageFinished = true;
            continue
        end


        numberOfRegions = str2double(regionResponse);


        if numberOfRegions == 0
            % A completely black mask is a valid no-wake/background label.
            wakeMask = false(size(trainingImage,1),size(trainingImage,2));
            imwrite(uint8(wakeMask)*255,labelFile,"tif");
            close(labelFigure);
            fprintf('Saved no-wake background mask: %s\n',labelFile);
            imageFinished = true;
            continue
        end


        if isnan(numberOfRegions) || ~isscalar(numberOfRegions) || ...
                numberOfRegions < 0 || ...
                numberOfRegions ~= floor(numberOfRegions)
            close(labelFigure);
            error("Enter a nonnegative whole number or K to skip.");
        end


        %% 5. DRAW FREEHAND REGIONS AND PAINT THEIR MASK DIRECTLY OVER THE IMAGE
        wakeMask = false(size(trainingImage,1),size(trainingImage,2));
        for regionNumber = 1:numberOfRegions
            fprintf(['Click and drag around wake region %d of %d. ' ...
                'Release the mouse to finish.\n'],regionNumber,numberOfRegions);


            wakeRegion = drawfreehand(drawingAxes, ...
                "Color","yellow","LineWidth",1.5);


            % Convert the ROI coordinates into a mask with exactly the same
            % dimensions as the underlying spectrogram matrix.
            points = wakeRegion.Position;
            regionMask = poly2mask(points(:,1),points(:,2), ...
                size(trainingImage,1),size(trainingImage,2));
            wakeMask = wakeMask | regionMask;


            % Paint accepted mask pixels red while leaving the spectrogram
            % visible underneath, making alignment easy to inspect one-to-one.
            maskOverlay.AlphaData = 0.40*double(wakeMask);
            drawnow;
        end


        %% 6. REVIEW THE OVERLAY BEFORE WRITING ANYTHING TO DISK
        reviewChoice = lower(strtrim(input( ...
            "Type S to save, R to redraw, or K to skip: ","s")));


        if reviewChoice == "s"
            % Background is 0 and wake is 255 for pixelLabelDatastore.
            imwrite(uint8(wakeMask)*255,labelFile,"tif");
            fprintf('Saved label: %s\n',labelFile);
            close(labelFigure);
            imageFinished = true;
        elseif reviewChoice == "r"
            close(labelFigure);
            fprintf('Redrawing %s.\n',baseName);
        elseif reviewChoice == "k"
            close(labelFigure);
            fprintf('Skipped %s without saving.\n',baseName);
            imageFinished = true;
        else
            close(labelFigure);
            error("Unknown choice. Run the script again and use S, R, or K.");
        end
    end
end


fprintf('\nFinished. Labels are located in:\n%s\n',trainingLabelFolder);


