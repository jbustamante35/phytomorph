--------------------
Introduction to TRiP
--------------------

TRiP is a matlab-based program for estimating circadian period from whole plant image data. TRiP includes a grid-based cropping function that takes each image stack as input and crops the images using grid coordinates to output image stacks for each plant. A motion estimation algorithm is applied to the image stacks and outputs a motion vector for each image over time. The motion vectors are used to estimate circadian period using a single frequency FFT-NLLS analysis. The current code has been tested on Matlab 2014b for Mac.

-------------
TRiP Contents 
-------------

Within the TRiP folder are 3 directories: code/ input/ output/

The code directory contains the TRiP functions:
cropAll.m
errorFunc.m
estimateAll.m
estimateMotion.m
evaluateModel.m
modelFitAll.m
space_time_deriv.m

Also within code/ is the testdata.txt file required for the cropAll function described below. 

The input directory contains 379 images of 9 plants that were imaged every 20 min. 

-----------------------------
Description of TRiP Functions
-----------------------------

1. Make a .txt file containing the cropping coordinates for each image stack including the name of the subdirectory for each image (see testdata.txt for an example). Put the image stack files into the input directory. 

2. run CropAll to generate ../input/<subdir>, each <subdir> contains an image stack for a single plant with the prefix 'crop_' and the name given in the 'txt' file.

NOTE: If you are using the ImageJ record macros function to obtain the coordinates for the crop you will need to include the following code to line 62 of the cropAll.m function in order to convert the ImageJ coordinates to the matlab format where the upper left corner is (1,1). 

% transpose imageJ files
y3 = y2;
y4 = y2 +x2;
x3 = y1;
x4 = y1+ x1;

In addition, replace line 64 with:
imC = im(y3:y4,x3:x4, :); % crop

3. run estimateAll. This will generate a .csv and a .png file for each <subdir> in output. These .csv files contain a single column of the vertical motion as a function of time and the .png files shows a plot of the vertical motion for each plant.

4. run modelFitAll. This takes all the .csv files found in the output directory. This will create a .txt and .png file with the results of model fitting (the frequency of the estimated motion). 

---------------------------
Running TRiP on Sample Data
---------------------------

1. open matlab

2. change directories to TRiP/code

3. at the matlab prompt: 

cropAll(’testdata.txt’);

4. when cropAll is done, you should see 9 folders in TRiP/input. Each of these folders contains 379 frames of a single plant.

5. at the matlab prompt: 

estimateAll;

6. when estimateAll is done, you should see 9 .csv files and 9 .png files in TRiP/ouput. Each .csv file contains the motion trace over time, and each .png file contains a plot of this motion trace over time.

7. at the matlab prompt: 

modelfitAll;

8. when modelFitAll is done, you should see Plant1_model.txt and Plant1_model.png for all 9 plants in TRiP/ouput. These contain the results of fitting a single harmonic to the motion trace.

NOTE: The output of modelFitAll is frequency, not period. We have not included the calculation for period because it depends on the time-series resolution. You can include this calculation in the modelFitAll code or do it separately in a spreadsheet application. To calculate period from frequency, use the following equation:


To calculate period from frequency:

T = No. of frames (based on estimateAll output .csv file)/frequency
Period= T/# of images per hour

Example:

The testdata plants were imaged every 20 min. The period for Plant1:

T = 360/4.730735
P = 76.09/3
Period = 25.4

---------------
TRiP for Octave
---------------

NOTE: Running TRiP on Octave has not been fully optimized. We have only tested the code on Octave-GUI 3.8 (binary installation) run on a Mac OS X Maverick. Currently we have not found an Octave equivalent function for our ModelFit optimization step that outputs equivalent frequency values. We recommend Octave users to take the motion vectors from the EstimateAll '.csv' output files and use other available circadian period estimation platforms (ie. Biodare). We will continue to work on a modelFit code for Octave. Please refer to the GitHub site for updates: http://github.com/KTgreenham/TRiP

* The following Octave packages and their dependencies should be installed and loaded prior to running TRiP:

image
optim

* The estimateAll step is much slower in Octave but does produce the same motion vector output.
* the matlab function 'getframe' is not available in Octave so there will be no '.png' files output with the plots. We are working on implementing an Octave equivalent function but at this time the user will have to generate the plot manually if needed.

----------------------
Running TRiP on Octave
----------------------

1. Open Octave (Octave-Gui 3.8)
2. In the File Browser window, change directories to TRiP/code
3. In the Octave command window:

pkg load image
pkg load optim

4. In the Octave command window: 

cropAll(’testdata.txt’);

5. when cropAll is done, you should see 9 folders in TRiP/input. Each of these folders contains 379 frames of a single plant.
6. In the Octave command window:

estimateAll;

7. When estimateAll is complete, you should see 9 .csv files in TRiP/output. Each .csv file contains the motion vectors over time.
8.The motion vectors are the input to the circadian period estimation method of choice. 