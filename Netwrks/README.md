## imgNet wake-detection workflow

This folder contains the scripts used to create, train, test, and apply the `imgNet_v5` wake detector that is eventually used by `PressureAnalysis/gatherwakepressuredata.m`.

### 1. Create the spectrogram images used by imgNet

Run `GetTrainingimages.m`.

This script starts with the raw pressure record and creates 10-minute spectrogram images. These are saved into the `SpectrogramTrainingData` folder:

- `TrainingImages`: the axis-free `.tif` spectrograms that are actually given to imgNet.
- `TimeAxisPreviews`: human-readable preview figures with time and frequency axes.

The important output for network training is `TrainingImages`. The preview images are only for checking what the spectrograms look like.

### 2. Create fixed-display spectrogram previews

Run `gettrainingimages_v2.m`.

This creates the matching fixed-PSD display images and saves them in:

- `displaypreview`

These are not used to train imgNet. They are used later so the final masks can be viewed on a consistent, less-filtered spectrogram image.

### 3. Save the raw PSD spectrogram matrices

Run `getrawspectrogram.m`.

This recreates the same spectrogram segmentation, but instead of saving an image, it saves the raw PSD matrix `S` for each spectrogram in:

- `rawspectrogramdata`

These `.mat` files are used later by `gatherwakepressuredata.m` to calculate pressure-derived wake variables such as RMS height, total energy, frequency metrics, and slope.

### 4. Manually label wake masks for training

Run `NtwrkTrainer_v2.m`.

This script opens random unlabeled images from `TrainingImages` and lets the user manually draw wake masks. The masks are saved in:

- `TrainingLabels`

Each label file has the same filename as its matching training image. The label image uses:

- `0` for background
- `255` for wake

These paired `TrainingImages` and `TrainingLabels` files are the supervised training set.

### 5. Train imgNet

Run `BoatNtwrk_v2.m`.

This script loads only image/label pairs where the filename exists in both `TrainingImages` and `TrainingLabels`. It creates the pixel-label datastore, pads each image and label by one row and one column, defines the convolutional network, trains the network, and saves the result as:

- `imgNet_v5.mat`

The saved variable inside the file is named:

- `imgNet`

This is the network later used for wake detection.

### 6. Test and tune imgNet filtering

Use `imgNet_tester_v2.m`.

This script loads a random training spectrogram, applies `imgNet_v5`, and then applies `imgNetfilters_v2.m`. It displays:

- the seed/grow detection on the filtered training image
- the final filtered mask on the fixed-PSD displaypreview image

The actual filtering parameters are kept inside `imgNetfilters_v2.m`, not inside the tester.

### 7. Final imgNet filtering function

The main reusable detection function is `imgNetfilters_v2.m`.

This function takes:

- the filtered training spectrogram image
- `imgNet_v5` or the loaded `imgNet` network
- optionally, the matching displaypreview image

It applies `imgNet`, creates wake probabilities, performs strict or nonstrict mode filtering, and outputs the final wake mask.

Important outputs include:

- `tifMask`: uint8 mask, where background is `0` and wake is `255`
- `wakePixels`: logical final wake mask
- `labeledWakes`: connected-component labels for detected wakes
- `numberOfWakes`: number of detected wake regions
- `filterInfo`: settings and summary information from the filtering run
- `filteredWakeProbability`: row-mean-subtracted wake probability image
- `filteredOnlyWakePixels`: earlier seed/grow mask before final filtering

This function is the bridge between the trained imgNet and the later pressure-analysis scripts.

### 8. Select usable wake times

Run `filterboattimes.m`.

This script applies `imgNetfilters_v2.m` across the training spectrograms and displays the detected wakes. The user clicks the wakes that are usable. The selected wake start times are saved to:

- `Boatwaketimes.xlsx`

This file tells the later pressure-analysis code which detected wakes should be included.

### 9. Use imgNet output in final pressure/camera analysis

Run `PressureAnalysis/gatherwakepressuredata.m`.

This script loads:

- `TrainingImages`
- `displaypreview`
- `rawspectrogramdata`
- `imgNet_v5.mat`
- `Boatwaketimes.xlsx`
- `boat_crossing_events.csv`

For each approved wake, it applies `imgNetfilters_v2.m` again to recreate the wake mask, applies that mask to the matching raw PSD spectrogram `S`, calculates pressure variables, matches the wake to the camera boat crossing data, and saves the final structure to:

- `final_boat_parameters.mat`

So the core chain is:

```text
GetTrainingimages.m
    -> TrainingImages

NtwrkTrainer_v2.m
    -> TrainingLabels

BoatNtwrk_v2.m
    -> imgNet_v5.mat

imgNetfilters_v2.m
    -> final wake masks

filterboattimes.m
    -> Boatwaketimes.xlsx

gatherwakepressuredata.m
    -> final_boat_parameters.mat
```
