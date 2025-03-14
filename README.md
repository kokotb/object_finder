# object_finder
Automated code for object finding and recording based on image analysis based decision making.

## Short description
The provided coded is modeled for Abberior STED system operated thorugh their proprietary Imspector software, that enables python based microscope control. Beside the simple functionalities such as moving the stage, we implmented custom calculation based autofocus correction that requries xz recordings. 

The code is flexible and enables running any number of windows with different settings that can be used for decision making which additional windows will be recorded. The provided code users to choose any number of sites (wells) in which panoramas of custom size are recorded. At each recording a decision step is implemented where certain type of object in microscope channel is detected. Upon detection, the selected objects are chosen as targets for smaller more detailed field-of-view recording, preceeded by a more accurate z-position determination. A maximum number of analyzed objects per frame and per site are selected.

The code was developed for sparse samples where a faster less acurate field-of-view is recorded to determine the presence and the position of the desired objects.
