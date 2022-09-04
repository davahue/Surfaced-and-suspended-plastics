
	0. GENERAL INFORMATION

This folder contains the data extracted from the plastic elements trajectories and a Monte Carlo analysis comparing theoretical suspended transport profiles.

	1. RAW DATA

Dictionary of the data originating from palstic elements trajectories:
	- tp (time),
	- xp, yp and zp (coordinates relative to ArUco markers, streamwise, bed-normal and traverse),
	- vpx, vpy, vpz (velocities of the plastic elements),
	- surfp (flag for surfaced plastic, = 1),
	- st, pv, Nsamples (information from the KS test, with increasing number of samples, Nsamples),
	- scalars (relevant outcomes of the experiment, including "plane, beta, binsize, Ca, a_surfaced_video, a_surfaced_KS")


The other data contained in this folder correspond to a Monte Carlo analysis comparing two Rouse profiles via a Kolmogorov-Smirnov (KS) test.
In the KS, one profile includes a limitted number of samples.

Dictionary of the data obtained in the Monte Carlo analysis:
	- kmax_list11330 (KS stat): number of samples included in the limitted-samples Rouse profile.
	- st_list11330 (KS stat): KS stat.
	- pv_list11330 (KS stat): p-value corresponding to a given KS test.

	2. PROCESSING THE DATA

The code "Plot KS experimental data.py" loads all data and shows how the empirical concentrations differ from a theoretical suspended concentration profile (Rouse profile), via a KS test and the associated p-value.

The code "Plastic monitoring strategies.py" shows the performance of different riverine monitoring strategies on the estimation of the total plastic budget, based on the experimental observations.

The code Concentration_profiles_faux.py" provides supporting functions to other codes.
