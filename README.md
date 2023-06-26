# WakeID
This repository contains MATLAB code for identifying and analyzing boat wakes using pressure sesnsors and cameras.
## Purpose
By identfying boat wakes through pressure sensors and cameras. We can work towards creating algorithms that can predict erosion on shorelines due to wake energy with just a simple sensor.
The code in this repository is outlined to work with raw data in the form of
 - Pressure 
   - Gathered from a sensor deployed in the ICW just off the CMS pier
     - Code FFTs data and then uses machine learning in order to identify boat wakes automatically. 
 - Pixel Intensity 
   - Gathered from a low cost camera attached to the CMS pier that overlooks the waterway. 
