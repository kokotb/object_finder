# object_finder
Automated code for object finding and recording based on image analysis-based decision-making.

## Hardware requirements
The Abberior Stimulated Emission Depletion (STED) system using Imspector software. The code was tested up to Imspector 16.3.16130-w2224 installed on Windows 10 Enterprise 2021 LTSC running Expert (Infinity) line system. It should work the same with the Facility line system (not tested). 

## Software requirements
To run the code python 3.7.10 version was used but the code should work with versions up to 3.10.4 (Imspector specpy library support). 

## Short description
The code enables the selection of the custom number of sites (wells) where custom-sized panoramas are recorded (for example 3 by 5) with set overlap. At each field-of-view (one site, one panorama tile) an overview window with lesser quality (lower number of pixels per range) is recorded. We then calculate nuclear masks, segment the cells, and determine the centres of detected nuclei (objects). A pre-defined number of objects per FOV are selected from the detected objects. At these positions a finer xz is recorded to more accurately determine the centre of the object and z position is adjusted accordingly. Next, we record a finer xy STED image with a smaller FOV and save the recorded data to a specified directory. When a maximum number per frame is reached, we record the next image of the panorama. The procedure is repeated until the maximum number of objects per well is reached. Once this condition is satisfied, the next well is recorded. The code continues until the panorama is finished recording or a maximum number of objects per well are recorded.

