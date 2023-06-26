%% This Code will allow the user to select parts of a spectrogram and convert
% the spectrogram into a binary matrix with 1's where the user selected a polygon 
% 0's where the polygon hasn't been selected. This can be used to create training 
% weights for machine learning. 

images = imageDatastore("C:\Users\lwlav\OneDrive\Documents\MATLAB\Pixel\TrainingImages","IncludeSubfolders",true);

i=5;

img = imread(images.Files{i});

numROI = 1;

figure,
imagesc(img)

figure,
imshow(img)
for k = 1:numROI

    h(k)= drawassisted(gca);
end

input('Press enter when done')

mask= false(height(img(1)),length(img(2)));
for k = 1:numROI
    mask = mask | createMask(h(k));
end

imwrite(mask,"C:\Users\lwlav\OneDrive\Documents\MATLAB\Pixel\TrainingLabels\Lab_pixelGlass_2boats_Oct0223.tif","tiff")
