	0. GENERAL INFORMATION

This folder contains two scripts:

	- 0_Processing_Camera1_example.py
	- 1_TrackingDataset_3Dreconstruction_example.py

The first (0_Processing_Camera1_example.py) processes an example experiment video (minimal due to space constraints) from an individual camera 
(located at ./VideoExampletrackingV3/Cam 1). The script performs a series of color and shape filters to estimate the location and orientation 
of a plastic element in a video sequence. Continuous calibration is performed to estimate the position of the camera with respect to the world 
reference system. This records a processed video example (Example_output_detection_processed.avi) and saves the raw coordinates in camera 
reference (Example_output_detection.csv) in the following structure:

id, 
FrameID (number of frame), 
centre_x (camera frame coordinate x for the detected object), 
centre_y (camera frame coordinate y for the detected object), 
mainAx_x (component_x of the main axis orientation),
mainAx_y (component_y of the main axis orientation),
rot1 (3 components of the camera rotation vector),
rot2,
rot3,
trans1 (3 components of the camera rotation vector),
trans2,
trans3,
Surfaced (binary value for 1 = plastic connected to the water surface, and 0 = suspended plastic)
 
The second script (1_TrackingDataset_3Dreconstruction_example.py) accepts an stereo-camera tracking output (2 cameras), ./input_test_3D_reconstruction
contains an experiment example. Particle 3D positions are extracted through 3D ray retracing and intersection between the stereo-couple. 

Both scripts require intrinsic and extrinsic camera parameters that are formated accordingly and can be found at in 
./ExtrinsicCalibration
./IntrinsicCalibration

Additionally, the 3D reconstruction requires a synchronization shift between two cameras, pre-computed and found in ./ExtrinsicCalibration/DataFrameSync_Phase2_01_Transport_test _Cups(+ve).xls
