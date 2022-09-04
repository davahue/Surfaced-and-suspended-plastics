
	0. GENERAL INFORMATION

This folder contains the 3D trajectories reconstruction obtained from the image analysis for all experiments.

	1. RAW DATA

The source videos used for the 3D trajectories take several Gb of space and ara available upon reasonable request.

Trajectories obtained from the automatic detection of plastics can be found in the following folders:
- "Cups_trajectories"
- "Films_trajectories"
- "Masks_trajectories"

Each experiment has a file (for instance, "Mask_V3_d5" corresponds to the masks experiment with flow ID V3).
The following information is available within each file:

	FrameId	x	y	z	Area	Surfaced
0	215	43.29393663	-10.71244616	69.98418304	6512.5	1

- FrameId: number of the frame in which a plastic element is detected.
- x: x coordinate, relative to ArUcos position.
- y: y coordinate, relative to ArUcos position.
- z: z coordinate, relative to ArUcos position.
- Area: area of the silhouette of the detected plastic element (mm2).
- Surfaced: flag indicating if the plastic element was in contact with the free surface.

	2. PROCESSED DATA

The following codes are available:
- "Automated_Concentration_profiles_extraction.py": routines to estimate count profiles from trajectories data.
- "Concentration_profiles_faux.py": auxiliary functions supporting the above code. It includes theoretical profile routines.

When running the code "Automated_Concentration_profiles_extraction.py":
- The folder and data file needs to be identified as: 
	- folder, for instance "Cups_trajectories".
	- filename, for instance "Cup_PP_98_def_V5_d5.csv".
- Information from the trajectories is used to compute the C.o.G. count, based on crossing the mid of the observation window.
- Plots are generated, including the KS test and the count profiles (referenced to free surface = 0-position).
- txt files with the following data are produced:
	- tp (time),
	- xp, yp and zp (coordinates relative to ArUco markers, streamwise, bed-normal and traverse),
	- vpx, vpy, vpz (velocities of the plastic elements),
	- surfp (flag for surfaced plastic, = 1),
	- st, pv, Nsamples (information from the KS test, with increasing number of samples, Nsamples),
	- scalars (relevant outcomes of the experiment, including "plane, beta, binsize, Ca, a_surfaced_video, a_surfaced_KS")

These data generated are also contained within the folder "Plastic concentrations".
